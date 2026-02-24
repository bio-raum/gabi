#!/usr/bin/env python
from datetime import datetime
import os
import glob
import json
import re
import argparse

parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--output", "-o")
parser.add_argument("--taxon", "-t")
parser.add_argument("--yaml", "-y")
parser.add_argument("--sample", "-s")


args = parser.parse_args()


def parse_bracken(lines):
    header = lines.pop(0).strip().split("\t")
    data = []
    for line in lines:
        this_data = {}
        elements = line.strip().split("\t")
        for idx, h in enumerate(header):
            this_data[h] = elements[idx]

        # only keep if fraction over 1%
        if (float(this_data["fraction_total_reads"]) > 0.01):
            data.append(this_data)

    return data


def parse_checkm(lines):
    data = {}
    header = lines.pop(0).split("\t")
    elements = lines.pop(0).split("\t")
    for h in header:
        value = elements.pop(0)
        if re.match(r"^[0-9]*$", value):
            value = int(value)
        elif re.match(r"^[0-9]*\.[0-9]*$", value):
            value = float(value)

        data[h] = value

    return data


def parse_kraken(lines):
    data = []
    for line in lines:
        elements = line.strip.split("\t")
        level = elements[3]
        percentage = float(elements[0])
        taxon = " ".join(elements[5:-1])

        if level == "S" and percentage >= 1.0:
            entry = {
                "taxon": taxon,
                "percentage": percentage
            }
            data.append(entry)
    return data


def parse_json(lines):
    data = json.loads(" ".join(lines))
    return data


def parse_quast(lines):
    data = {}
    for line in lines:
        key, value = line.split("\t")
        if re.match(r"^[0-9]*$", value):
            value = int(value)
        elif re.match(r"^[0-9]*\.[0-9]*$", value):
            value = float(value)
        data[key] = value
    return data


def parse_csv(lines):
    header = lines.pop(0).strip().split(",")
    data = []
    for line in lines:
        this_data = {}
        elements = line.strip().split(",")
        for idx, h in enumerate(header):
            entry = elements[idx]
            if re.match(r"^[0-9]*$", entry):
                entry = int(entry)
            elif re.match(r"^[0-9]*\.[0-9]*$", entry):
                entry = float(entry)
            this_data[h] = entry
        data.append(this_data)
    return data


def parse_tabular(lines):
    header = lines.pop(0).strip().split("\t")
    data = []
    for line in lines:
        this_data = {}
        elements = line.strip().split("\t")
        for idx, h in enumerate(header):
            if idx < len(elements):
                entry = elements[idx]
                # value is an integer
                if re.match(r"^[0-9]+$", entry):
                    entry = int(entry)
                # value is a float
                elif re.match(r"^[0-9]+\.[0-9]+$", entry):
                    entry = float(entry)
                # value is a file path (messes up md5 fingerprinting)
                elif re.match(r"^\/.*\/.*$", entry):
                    entry = entry.split("/")[-1]
                this_data[h] = entry
        data.append(this_data)

    return data


def parse_mosdepth_global(lines):
    data = {}
    for line in lines:
        chr, cov, perc = line.split("\t")
        if chr == "total" and int(cov) <= 100:
            data[cov] = round(float(perc) * 100, 2)
    return data


def parse_genbank(lines):
    data = {}
    for line in lines:
        elements = line.split()
        if re.match(r"^LOCUS.*", line):
            data['locus'] = elements[1]
        elif re.match(r"^DEFINITION.*", line):
            definition = " ".join(elements[1:-1]).split(", ")[0]
            data["definition"] = definition
        elif re.match(r".*Assembly\:.*", line):
            assembly = line.split("Assembly:")[-1].strip()
            data["assembly"] = assembly
    return data


def parse_samtools_stats(lines):
    data = {}
    inserts = []

    for line in lines:
        if re.match(r"^SN.*", line):
            trimline = line.split(" #")[0]
            elements = trimline.split("\t")
            data[elements[1].replace(":", "")] = elements[2]
        elif re.match(r"^IS.*", line):
            elements = line.split("\t")
            pos = int(elements[1])
            if pos > 0 and pos < 1000:
                inserts.append(int(elements[2]))

    data["insert_sizes"] = inserts
    return data


def parse_nanostat(lines):
    data = {}
    for line in lines:
        elements = line.split("\t")
        if (len(elements) > 1):
            key = elements[0]
            if re.match(r"^>Q", key):
                key = key.replace(">", "").replace(":", "")
                value = int(elements[1].split()[0])
                data[key] = value
        else:
            elements = line.split()
            if re.search("Mean read length", line):
                data["mean_read_length"] = float(elements[-1].replace(",", ""))
            elif re.search("Read length N50", line):
                data["read_length_n50"] = float(elements[-1].replace(",", ""))
    return data


def parse_taxonkit(lines):
    data = {}
    species_data = {}
    genus_data = {}
    sum = 0
    for line in lines:

        elements = line.split("\t")
        basepairs = int(elements[2])

        if len(elements) < 4:
            continue

        contig, this_id, this_length, taxstring = elements

        sum += basepairs

        taxdata = {}
        for e in taxstring.split("|"):
            level, t = e.split("__")
            taxdata[level] = t

        if "s" in taxdata:
            species = taxdata["s"]
            if species in species_data:
                species_data[species] += basepairs
            else:
                species_data[species] = basepairs

        if "g" in taxdata:
            genus = taxdata["g"]
            if genus in genus_data:
                genus_data[genus] += basepairs
            else:
                genus_data[genus] = basepairs

    data["genus"] = []

    for genus, length in genus_data.items():
        data["genus"].append({"genus": genus, "basepairs": length, "fraction": round((length / sum), 2)})

    data["species"] = []
    for species, length in species_data.items():
        data["species"].append({"species": species, "basepairs": length, "fraction": round((length / sum), 2)})

    data["length"] = sum

    return data


def parse_yaml(lines):

    data = {}
    key = ""

    for line in lines:

        line = line.replace(":", "")
        if re.match(r"^\s+.*", line):
            tool, version = line.strip().split()
            data[key][tool] = version
        else:
            key = line.strip()
            data[key] = {}

    return data


def parse_bcftools(lines):
    data = {}

    for line in lines:
        if re.search("^SN", line):
            elements = line.split("\t")
            value = int(elements[-1])
            key = elements[-2].split(" ")[-1][:-1]
            data[key] = value
    return data


def delete_keys_from_dict(dict_del, lst_keys):
    for k in lst_keys:
        try:
            del dict_del[k]
        except KeyError:
            pass
    for v in dict_del.values():
        if isinstance(v, dict):
            delete_keys_from_dict(v, lst_keys)

    return dict_del


def main(sample, taxon, yaml_file, output):

    files = [os.path.abspath(f) for f in glob.glob("*/*")]
    date = datetime.today().strftime('%Y-%m-%d')

    remove_keys = [ 
        "content_curves",
        "quality_curves",
        "kmer_count"
    ]

    with open(yaml_file, "r") as f:
        yaml_lines = [line.rstrip() for line in f]

    versions = parse_yaml(yaml_lines)

    matrix = {
        "date": date,
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
        "software": versions,
        "taxonkit": {}
    }

    for file in files:

        with open(file, "r") as f:
            lines = [line.rstrip() for line in f]

        if re.search(".assembly.bracken.tsv", file):
            matrix["assembly"] = parse_bracken(lines)
        elif re.search(".ILLUMINA.bracken.tsv", file):
            matrix["bracken"]["ILLUMINA"] = parse_bracken(lines)
        elif re.search(".NANOPORE.bracken.tsv", file):
            matrix["bracken"]["NANOPORE"] = parse_bracken(lines)
        elif re.search(r"^.*ILLUMINA.*report_bracken.*", file):
            matrix["kraken"]["ILLUMINA"] = parse_kraken(lines)
        elif re.search(".*NANOPORE*.*report_bracken.*", file):
            matrix["kraken"]["NANOPORE"] = parse_kraken(lines)
        elif re.search(".NanoStats.txt", file):
            matrix["nanostat"] = parse_nanostat(lines)
        elif re.search(".fastp.json", file):
            fastp = parse_json(lines)
            fastp_clean = delete_keys_from_dict(fastp, remove_keys)
            matrix["fastp"] = fastp_clean
        elif re.search(".*mlst.json", file):
            mlst = parse_json(lines)
            matrix["mlst"].append(mlst[0])
        elif re.search(r"^.*ILLUMINA.*confindr.*", file):
            matrix["confindr"]["illumina"].append(parse_csv(lines))
        elif re.search(r"^.*NANOPORE.*confindr.*", file):
            confindr = parse_csv(lines)
            matrix["confindr"]["nanopore"].append(confindr)
        elif re.search(r".*/report.tsv", file):
            matrix["quast"] = parse_quast(lines)
        elif "Protein identifier" in lines[0]:
            matrix["amr"]["amrfinder"] = parse_tabular(lines)
        elif re.search("ectyper.tsv", file):
            ectyper = parse_tabular(lines)
            matrix["serotype"]["ectyper"] = ectyper[0]
        elif re.search(".*seqsero2.tsv", file):
            seqsero = parse_tabular(lines)
            matrix["serotype"]["seqsero2"] = seqsero[0]
        elif re.search(".*lissero.tsv", file):
            lissero = parse_tabular(lines)
            matrix["serotype"]["lissero"] = lissero[0]
        elif re.search(".stecfinder.tsv", file):
            stecfinder = parse_tabular(lines[0:2])
            matrix["serotype"]["stecfinder"] = stecfinder[0]
        elif re.search(".kaptive.json", file):
            kaptive = parse_json(lines)
            kaptive.pop('expected_genes_outside_locus', None)
            matrix["serotype"]["kaptive"] = kaptive
        elif re.search("ILLUMINA.mosdepth.summary.txt", file):
            mosdepth = [d for d in parse_tabular(lines) if d['chrom'] == "total"]
            matrix["mosdepth"]["illumina"] = mosdepth[0]
        elif re.search("NANOPORE.mosdepth.summary.txt", file):
            mosdepth = [d for d in parse_tabular(lines) if d['chrom'] == "total"]
            matrix["mosdepth"]["nanopore"] = mosdepth[0]
        elif re.search("PACBIO.mosdepth.summary.txt", file):
            mosdepth = [d for d in parse_tabular(lines) if d['chrom'] == "total"]
            matrix["mosdepth"]["pacbio"] = mosdepth[0]
        elif re.search(".*mosdepth.summary.txt", file):
            mosdepth = [d for d in parse_tabular(lines) if d['chrom'] == "total"]
            matrix["mosdepth"]["total"] = mosdepth[0]
        elif re.search("ILLUMINA.mosdepth.global.dist.txt", file):
            matrix["mosdepth_global"]["illumina"] = parse_mosdepth_global(lines)
        elif re.search("NANOPORE.mosdepth.global.dist.txt", file):
            matrix["mosdepth_global"]["nanopore"] = parse_mosdepth_global(lines)
        elif re.search("PACBIO.mosdepth.global.dist.txt", file):
            matrix["mosdepth_global"]["pacbio"] = parse_mosdepth_global(lines)
        elif re.search(".*mosdepth.global.dist.txt", file):
            matrix["mosdepth_global"]["total"] = parse_mosdepth_global(lines)
        elif re.search(".sistr.tab", file):
            matrix["serotype"]["sistr"] = parse_tabular(lines)[0]
        elif re.search(".gbff$", file):
            matrix["reference"] = parse_genbank(lines)
        elif re.search(".stats$", file):
            matrix["samtools"] = parse_samtools_stats(lines)
        elif re.search("mobtyper_results.txt", file):
            matrix["plasmids"] = parse_tabular(lines)
        elif re.search("abricate", file):
            results = parse_tabular(lines)
            if not results:
                continue
            abricate_db = results[0]["DATABASE"]
            if ("abricate" not in matrix["amr"]):
                matrix["amr"]["abricate"] = {}
            matrix["amr"]["abricate"][abricate_db] = results
        elif re.search("btyper3.tsv", file):
            matrix["serotype"]["btyper3"] = parse_tabular(lines)[0]
        elif re.search("sccmec", file):
            matrix["serotype"]["sccmec"] = parse_tabular(lines)[0]
        elif re.search("bcftools_stats.txt", file):
            matrix["variants"] = parse_bcftools(lines)
        elif re.search("checkm2_report.tsv", file):
            matrix["checkm"] = parse_checkm(lines)
        elif re.search("taxonkit.txt", file):
            matrix["taxonkit"] = parse_taxonkit(lines)
        elif re.search(r".*short_summary.*json", file):
            busco = parse_json(lines)
            dataset = busco["dataset"]
            ifile = busco["input_file"]
            busco["dataset"] = dataset.split("/")[-1]
            busco["input_file"] = ifile.split("/")[-1]
            matrix["busco"] = busco

    with open(output, "w") as fo:
        json.dump(matrix, fo, indent=4, sort_keys=True)


if __name__ == '__main__':
    main(args.sample, args.taxon, args.yaml, args.output)
