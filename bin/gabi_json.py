#!/usr/bin/env python3
import argparse
import csv
from datetime import datetime
import glob
import json
import os
import re
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate bioinformatics pipeline outputs into JSON")
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--taxon", "-t", required=True)
    parser.add_argument("--yaml", "-y", required=True)
    parser.add_argument("--settings", "-x", required=True)
    parser.add_argument("--sample", "-s", required=True)
    return parser.parse_args()


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


def abricate_to_serotype(entries: list[dict[str, Any]]) -> str:
    h_genes = set()
    o_genes = set()
    for entry in entries:
        gene = entry["GENE"].rsplit("-", 1)[-1]
        if gene.startswith("H"):
            h_genes.add(gene)
        elif gene.startswith("O"):
            o_genes.add(gene)
    return f"{'/'.join(o_genes)}:{'/'.join(h_genes)}"


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


def main(sample: str, taxon: str, yaml_file: str, settings_json: str, output: str) -> None:
    files = sorted(os.path.abspath(f) for f in glob.glob("*/*") if os.path.isfile(f))
    remove_keys = {"content_curves", "quality_curves", "kmer_count"}

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

    platforms = ("ILLUMINA", "NANOPORE", "PACBIO")

    for filepath in files:
        fname = os.path.basename(filepath)

        # 1. Bracken & Kraken
        if fname.endswith(".assembly.bracken.tsv"):
            matrix["assembly"] = parse_bracken(filepath)
        elif any(fname.endswith(f".{p}.bracken.tsv") for p in platforms):
            plat = next(p for p in platforms if fname.endswith(f".{p}.bracken.tsv"))
            matrix["bracken"][plat] = parse_bracken(filepath)
        elif "report_bracken" in fname:
            for plat in ("ILLUMINA", "NANOPORE"):
                if plat in fname:
                    matrix["kraken"][plat] = parse_kraken(filepath)
                    break

        # 2. Mosdepth (summary & global dist)
        elif fname.endswith("mosdepth.summary.txt"):
            plat = next((p.lower() for p in platforms if f"{p}.mosdepth.summary.txt" in fname), "total")
            total_rows = [d for d in parse_tabular(filepath) if d.get("chrom") == "total"]
            if total_rows:
                matrix["mosdepth"][plat] = total_rows[0]
        elif fname.endswith("mosdepth.global.dist.txt"):
            plat = next((p.lower() for p in platforms if f"{p}.mosdepth.global.dist.txt" in fname), "total")
            matrix["mosdepth_global"][plat] = parse_mosdepth_global(filepath)

        # 3. Confindr
        elif "confindr" in fname:
            for plat in ("ILLUMINA", "NANOPORE"):
                if plat in fname:
                    matrix["confindr"][plat.lower()].append(parse_csv(filepath))
                    break

        # 4. QC & Assembly Stats
        elif fname.endswith(".NanoStats.txt"):
            matrix["nanostat"] = parse_nanostat(filepath)
        elif fname.endswith(".fastp.json"):
            matrix["fastp"] = delete_keys_from_dict(parse_json_file(filepath), remove_keys)
        elif fname == "report.tsv":
            matrix["quast"] = parse_quast(filepath)
        elif fname.endswith(".gbff"):
            matrix["reference"] = parse_genbank(filepath)
        elif fname.endswith(".stats"):
            matrix["samtools"] = parse_samtools_stats(filepath)
        elif fname.endswith("bcftools_stats.txt"):
            matrix["variants"] = parse_bcftools(filepath)
        elif fname.endswith("checkm2_report.tsv"):
            matrix["checkm"] = parse_checkm(filepath)
        elif fname.endswith("taxonkit.txt"):
            matrix["taxonkit"] = parse_taxonkit(filepath)
        elif "short_summary" in fname and fname.endswith(".json"):
            busco = parse_json_file(filepath)
            busco["dataset"] = os.path.basename(busco["dataset"])
            busco["input_file"] = os.path.basename(busco["input_file"])
            matrix["busco"] = busco

        # 5. Typing, AMR & Plasmids
        elif fname.endswith("mlst.json"):
            mlst = parse_json_file(filepath)
            if mlst:
                matrix["mlst"].append(mlst[0])
        elif fname.endswith("ectyper.tsv"):
            matrix["serotype"]["ectyper"] = parse_tabular(filepath, max_rows=1)[0]
        elif fname.endswith("seqsero2.tsv"):
            matrix["serotype"]["seqsero2"] = parse_tabular(filepath, max_rows=1)[0]
        elif fname.endswith("lissero.tsv"):
            matrix["serotype"]["lissero"] = parse_tabular(filepath, max_rows=1)[0]
        elif fname.endswith(".stecfinder.tsv"):
            matrix["serotype"]["stecfinder"] = parse_tabular(filepath, max_rows=1)[0]
        elif fname.endswith(".sistr.tab"):
            matrix["serotype"]["sistr"] = parse_tabular(filepath, max_rows=1)[0]
        elif fname.endswith("btyper3.tsv"):
            matrix["serotype"]["btyper3"] = parse_tabular(filepath, max_rows=1)[0]
        elif "sccmec" in fname:
            matrix["serotype"]["sccmec"] = parse_tabular(filepath, max_rows=1)[0]
        elif fname.endswith(".kaptive.json"):
            kaptive = parse_json_file(filepath)
            kaptive.pop("expected_genes_outside_locus", None)
            matrix["serotype"]["kaptive"] = kaptive
        elif fname.endswith("mobtyper_results.txt"):
            matrix["plasmids"] = parse_tabular(filepath)
        elif "ecoh.abricate" in fname:
            results = parse_tabular(filepath)
            if results:
                matrix["serotype"]["ecoh"] = {"serotype": abricate_to_serotype(results)}
        elif "abricate" in fname:
            results = parse_tabular(filepath)
            if results:
                abricate_db = results[0]["DATABASE"]
                matrix["amr"].setdefault("abricate", {})[abricate_db] = results
        elif fname.endswith(".md5"):
            matrix["md5sums"] = parse_md5(filepath)
        elif first_line_has(filepath, "Protein identifier"):
            matrix["amr"]["amrfinder"] = parse_tabular(filepath)

    with open(output, "w") as fo:
        json.dump(matrix, fo, indent=4, sort_keys=True)


if __name__ == "__main__":
    args = parse_args()
    main(args.sample, args.taxon, args.yaml, args.settings, args.output)
