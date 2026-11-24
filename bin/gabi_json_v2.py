#!/usr/bin/env python3
import argparse
import csv
from datetime import datetime
import glob
import json
import logging
import os
import re
from typing import Any, Callable

logger = logging.getLogger(__name__)

REMOVE_KEYS = {"content_curves", "quality_curves", "kmer_count"}
Handler = Callable[[str, re.Match, dict[str, Any]], None]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate bioinformatics pipeline outputs into a single JSON file."
    )
    parser.add_argument("--output", "-o", required=True, help="Output JSON filepath")
    parser.add_argument("--taxon", "-t", required=True, help="Target taxon name")
    parser.add_argument("--yaml", "-y", required=True, help="Software versions YAML file")
    parser.add_argument("--settings", "-x", required=True, help="Pipeline settings JSON file")
    parser.add_argument("--sample", "-s", required=True, help="Sample identifier")
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable detailed debug logging"
    )
    return parser.parse_args()


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def cast_value(val: str, strip_paths: bool = False) -> Any:
    """Safely convert numeric strings to int/float and optionally strip directory paths."""
    if not val:
        return val
    if val.isdigit():
        return int(val)
    try:
        return float(val)
    except ValueError:
        if strip_paths and val.startswith("/") and "/" in val[1:]:
            return os.path.basename(val)
        return val


# =====================================================================
# File Parsers
# =====================================================================


def abricate_to_serotype(entries: list[dict[str, Any]]) -> str:
    h_genes = set()
    o_genes = set()
    for entry in entries:
        gene = entry["GENE"].rsplit("-", 1)[-1]
        if gene.startswith("H"):
            h_genes.add(gene)
        elif gene.startswith("O"):
            o_genes.add(gene)
    return f"{'/'.join(sorted(o_genes))}:{'/'.join(sorted(h_genes))}"


def parse_bracken(filepath: str) -> list[dict[str, str]]:
    data = []
    with open(filepath, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if float(row["fraction_total_reads"]) > 0.01:
                data.append(dict(row))
    return data


def parse_checkm(filepath: str) -> dict[str, Any]:
    with open(filepath, "r") as f:
        header = f.readline().rstrip("\r\n").split("\t")
        values = f.readline().rstrip("\r\n").split("\t")
    return {h: cast_value(v) for h, v in zip(header, values)}


def parse_kraken(filepath: str) -> list[dict[str, Any]]:
    data = []
    with open(filepath, "r") as f:
        for line in f:
            elements = line.strip().split("\t")
            if len(elements) >= 6 and elements[3] == "S":
                percentage = float(elements[0])
                if percentage >= 1.0:
                    data.append({
                        "taxon": " ".join(elements[5:-1]),
                        "percentage": percentage,
                    })
    return data


def parse_json_file(filepath: str) -> Any:
    with open(filepath, "r") as f:
        return json.load(f)


def parse_quast(filepath: str) -> dict[str, Any]:
    data = {}
    with open(filepath, "r") as f:
        for line in f:
            key, value = line.rstrip("\r\n").split("\t", 1)
            data[key] = "0" if value == "0" else cast_value(value)
    return data


def parse_csv(filepath: str) -> list[dict[str, Any]]:
    data = []
    with open(filepath, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({k: cast_value(v) for k, v in row.items()})
    return data


def parse_tabular(filepath: str, max_rows: int | None = None) -> list[dict[str, Any]]:
    data = []
    with open(filepath, "r") as f:
        header_line = f.readline()
        if not header_line:
            return data
        header = header_line.rstrip("\r\n").split("\t")

        for idx, line in enumerate(f):
            if max_rows is not None and idx >= max_rows:
                break
            elements = line.rstrip("\r\n").split("\t")
            row = {
                h: cast_value(elements[i], strip_paths=True)
                for i, h in enumerate(header)
                if i < len(elements)
            }
            data.append(row)
    return data


def parse_mosdepth_global(filepath: str) -> dict[str, float]:
    data = {}
    with open(filepath, "r") as f:
        for line in f:
            chrom, cov, perc = line.rstrip("\r\n").split("\t")
            if chrom == "total" and int(cov) <= 100:
                data[cov] = round(float(perc) * 100, 2)
    return data


def parse_genbank(filepath: str) -> dict[str, str]:
    data = {}
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("LOCUS"):
                data["locus"] = line.split()[1]
            elif line.startswith("DEFINITION"):
                elements = line.split()
                data["definition"] = " ".join(elements[1:-1]).split(", ")[0]
            elif "Assembly:" in line:
                data["assembly"] = line.split("Assembly:", 1)[-1].strip()
    return data


def parse_samtools_stats(filepath: str) -> dict[str, Any]:
    data: dict[str, Any] = {}
    inserts = []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("SN"):
                trimline = line.split(" #", 1)[0]
                elements = trimline.split("\t")
                data[elements[1].rstrip(":")] = elements[2].strip()
            elif line.startswith("IS"):
                elements = line.split("\t")
                pos = int(elements[1])
                if 0 < pos < 1000:
                    inserts.append(int(elements[2]))
    data["insert_sizes"] = inserts
    return data


def parse_nanostat(filepath: str) -> dict[str, Any]:
    data: dict[str, Any] = {}
    with open(filepath, "r") as f:
        for line in f:
            if "\t" in line:
                key, rest = line.split("\t", 1)
                if key.startswith(">Q"):
                    clean_key = key.replace(">", "").replace(":", "")
                    data[clean_key] = int(rest.split()[0])
            elif "Mean read length" in line:
                data["mean_read_length"] = float(line.split()[-1].replace(",", ""))
            elif "Read length N50" in line:
                data["read_length_n50"] = float(line.split()[-1].replace(",", ""))
    return data


def parse_taxonkit(filepath: str) -> dict[str, Any]:
    species_data: dict[str, int] = {}
    genus_data: dict[str, int] = {}
    total_bp = 0

    with open(filepath, "r") as f:
        for line in f:
            elements = line.rstrip("\r\n").split("\t")
            if len(elements) < 4:
                continue

            basepairs = int(elements[2])
            total_bp += basepairs

            for item in elements[3].split("|"):
                if "__" not in item:
                    continue
                level, name = item.split("__", 1)
                if level == "s":
                    species_data[name] = species_data.get(name, 0) + basepairs
                elif level == "g":
                    genus_data[name] = genus_data.get(name, 0) + basepairs

    denom = total_bp or 1
    return {
        "genus": [
            {"genus": g, "basepairs": bp, "fraction": round(bp / denom, 2)}
            for g, bp in genus_data.items()
        ],
        "species": [
            {"species": s, "basepairs": bp, "fraction": round(bp / denom, 2)}
            for s, bp in species_data.items()
        ],
        "length": total_bp,
    }


def parse_yaml(filepath: str) -> dict[str, dict[str, str]]:
    data: dict[str, dict[str, str]] = {}
    current_key = ""
    with open(filepath, "r") as f:
        for raw_line in f:
            line = raw_line.rstrip("\r\n").replace(":", "")
            if not line.strip():
                continue
            if line[0].isspace():
                tool, version = line.strip().split(None, 1)
                data[current_key][tool] = version
            else:
                current_key = line.strip()
                data[current_key] = {}
    return data


def parse_bcftools(filepath: str) -> dict[str, int]:
    data = {}
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("SN"):
                elements = line.rstrip("\r\n").split("\t")
                value = int(elements[-1])
                key = elements[-2].split(" ")[-1][:-1]
                data[key] = value
    return data


def delete_keys_from_dict(target: dict[str, Any], keys_to_remove: set[str]) -> dict[str, Any]:
    for k in keys_to_remove:
        target.pop(k, None)
    for v in target.values():
        if isinstance(v, dict):
            delete_keys_from_dict(v, keys_to_remove)
    return target


def parse_md5(filepath: str) -> dict[str, str]:
    data = {}
    with open(filepath, "r") as f:
        for line in f:
            elements = line.split()
            if len(elements) == 2:
                data[elements[1]] = elements[0]
    return data


def first_line_has(filepath: str, substring: str) -> bool:
    try:
        with open(filepath, "r") as f:
            return substring in f.readline()
    except UnicodeDecodeError:
        return False


# =====================================================================
# Dispatch Handlers & Routing Table
# =====================================================================


def handle_single_row_serotype(tool: str) -> Handler:
    """Factory for single-row tabular serotype tools (ectyper, seqsero2, etc.)."""
    def _handler(filepath: str, match: re.Match, matrix: dict[str, Any]) -> None:
        rows = parse_tabular(filepath, max_rows=1)
        if not rows:
            logger.warning("Empty output in %s file: %s", tool, filepath)
            return
        matrix["serotype"][tool] = rows[0]
    return _handler


def handle_mosdepth_summary(filepath: str, match: re.Match, matrix: dict[str, Any]) -> None:
    plat = (match.group("plat") or "total").lower()
    total_rows = [d for d in parse_tabular(filepath) if d.get("chrom") == "total"]
    if not total_rows:
        logger.warning("No 'total' row found in mosdepth summary: %s", filepath)
        return
    matrix["mosdepth"][plat] = total_rows[0]


def handle_mosdepth_global(filepath: str, match: re.Match, matrix: dict[str, Any]) -> None:
    plat = (match.group("plat") or "total").lower()
    matrix["mosdepth_global"][plat] = parse_mosdepth_global(filepath)


def handle_mlst(filepath: str, match: re.Match, matrix: dict[str, Any]) -> None:
    mlst = parse_json_file(filepath)
    if not mlst:
        logger.warning("Empty MLST JSON file: %s", filepath)
        return
    matrix["mlst"].append(mlst[0])


def handle_busco(filepath: str, match: re.Match, matrix: dict[str, Any]) -> None:
    busco = parse_json_file(filepath)
    busco["dataset"] = os.path.basename(busco["dataset"])
    busco["input_file"] = os.path.basename(busco["input_file"])
    matrix["busco"] = busco


def handle_kaptive(filepath: str, match: re.Match, matrix: dict[str, Any]) -> None:
    kaptive = parse_json_file(filepath)
    kaptive.pop("expected_genes_outside_locus", None)
    matrix["serotype"]["kaptive"] = kaptive


def handle_abricate(filepath: str, match: re.Match, matrix: dict[str, Any]) -> None:
    results = parse_tabular(filepath)
    if not results:
        logger.warning("Empty abricate results (no hits) in: %s", filepath)
        return
    if "ecoh.abricate" in os.path.basename(filepath):
        matrix["serotype"]["ecoh"] = {"serotype": abricate_to_serotype(results)}
    else:
        db = results[0]["DATABASE"]
        matrix["amr"].setdefault("abricate", {})[db] = results


DISPATCH_RULES: list[tuple[re.Pattern, Handler]] = [
    # Bracken & Kraken
    (re.compile(r"\.assembly\.bracken\.tsv$"),
     lambda f, m, mx: mx.update(assembly=parse_bracken(f))),
    (re.compile(r"\.(?P<plat>ILLUMINA|NANOPORE|PACBIO)\.bracken\.tsv$"),
     lambda f, m, mx: mx["bracken"].update({m.group("plat"): parse_bracken(f)})),
    (re.compile(r"(?P<plat>ILLUMINA|NANOPORE).*report_bracken"),
     lambda f, m, mx: mx["kraken"].update({m.group("plat"): parse_kraken(f)})),

    # Mosdepth
    (re.compile(r"(?:(?P<plat>ILLUMINA|NANOPORE|PACBIO)\.)?mosdepth\.summary\.txt$"),
     handle_mosdepth_summary),
    (re.compile(r"(?:(?P<plat>ILLUMINA|NANOPORE|PACBIO)\.)?mosdepth\.global\.dist\.txt$"),
     handle_mosdepth_global),

    # Confindr
    (re.compile(r"(?P<plat>ILLUMINA|NANOPORE).*confindr"),
     lambda f, m, mx: mx["confindr"][m.group("plat").lower()].append(parse_csv(f))),

    # QC & Assembly
    (re.compile(r"\.NanoStats\.txt$"),
     lambda f, m, mx: mx.update(nanostat=parse_nanostat(f))),
    (re.compile(r"\.fastp\.json$"),
     lambda f, m, mx: mx.update(fastp=delete_keys_from_dict(parse_json_file(f), REMOVE_KEYS))),
    (re.compile(r"^report\.tsv$"),
     lambda f, m, mx: mx.update(quast=parse_quast(f))),
    (re.compile(r"\.gbff$"),
     lambda f, m, mx: mx.update(reference=parse_genbank(f))),
    (re.compile(r"\.stats$"),
     lambda f, m, mx: mx.update(samtools=parse_samtools_stats(f))),
    (re.compile(r"bcftools_stats\.txt$"),
     lambda f, m, mx: mx.update(variants=parse_bcftools(f))),
    (re.compile(r"checkm2_report\.tsv$"),
     lambda f, m, mx: mx.update(checkm=parse_checkm(f))),
    (re.compile(r"taxonkit\.txt$"),
     lambda f, m, mx: mx.update(taxonkit=parse_taxonkit(f))),
    (re.compile(r"short_summary.*\.json$"),
     handle_busco),

    # Typing, AMR & Plasmids
    (re.compile(r"mlst\.json$"), handle_mlst),
    (re.compile(r"ectyper\.tsv$"), handle_single_row_serotype("ectyper")),
    (re.compile(r"seqsero2\.tsv$"), handle_single_row_serotype("seqsero2")),
    (re.compile(r"lissero\.tsv$"), handle_single_row_serotype("lissero")),
    (re.compile(r"\.stecfinder\.tsv$"), handle_single_row_serotype("stecfinder")),
    (re.compile(r"\.sistr\.tab$"), handle_single_row_serotype("sistr")),
    (re.compile(r"btyper3\.tsv$"), handle_single_row_serotype("btyper3")),
    (re.compile(r"sccmec"), handle_single_row_serotype("sccmec")),
    (re.compile(r"\.kaptive\.json$"), handle_kaptive),
    (re.compile(r"mobtyper_results\.txt$"),
     lambda f, m, mx: mx.update(plasmids=parse_tabular(f))),
    (re.compile(r"abricate"), handle_abricate),
    (re.compile(r"\.md5$"),
     lambda f, m, mx: mx.update(md5sums=parse_md5(f))),
]


def route_file(filepath: str, matrix: dict[str, Any]) -> bool:
    """Match a filepath against DISPATCH_RULES and parse it into matrix. Returns True if handled."""
    fname = os.path.basename(filepath)
    for pattern, handler in DISPATCH_RULES:
        match = pattern.search(fname)
        if match:
            logger.debug("Matched '%s' with pattern '%s'", fname, pattern.pattern)
            handler(filepath, match, matrix)
            return True

    if first_line_has(filepath, "Protein identifier"):
        logger.debug("Matched '%s' via AMRFinder header fallback", fname)
        matrix["amr"]["amrfinder"] = parse_tabular(filepath)
        return True

    logger.debug("Skipping unrecognized file: %s", filepath)
    return False


# =====================================================================
# Main Entrypoint
# =====================================================================


def main(sample: str, taxon: str, yaml_file: str, settings_json: str, output: str) -> None:
    files = sorted(os.path.abspath(f) for f in glob.glob("*/*") if os.path.isfile(f))
    logger.info("Found %d candidate files for sample '%s'", len(files), sample)

    matrix: dict[str, Any] = {
        "date": datetime.today().strftime("%Y-%m-%d"),
        "sample": sample,
        "taxon": taxon,
        "quast": {},
        "mlst": [],
        "confindr": {"illumina": [], "nanopore": []},
        "confindr_nanopore": [],
        "serotype": {},
        "mosdepth": {},
        "reference": {},
        "mosdepth_global": {},
        "kraken": {},
        "bracken": {},
        "amr": {},
        "assembly": [],
        "checkm": {},
        "software": parse_yaml(yaml_file),
        "taxonkit": {},
        "pipeline_settings": parse_json_file(settings_json),
    }

    parsed_count = 0
    for filepath in files:
        try:
            if route_file(filepath, matrix):
                parsed_count += 1
        except Exception:
            logger.exception("Failed to parse file: %s", filepath)

    with open(output, "w") as fo:
        json.dump(matrix, fo, indent=4, sort_keys=True)

    logger.info(
        "Successfully parsed %d/%d files and wrote output to %s",
        parsed_count,
        len(files),
        output,
    )


if __name__ == "__main__":
    args = parse_args()
    configure_logging(args.verbose)
    main(args.sample, args.taxon, args.yaml, args.settings, args.output)
