#!/usr/bin/env python
import argparse
import json
import os
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from concurrent.futures import ThreadPoolExecutor


def read_json_files(file_paths):
    data = []
    for file_path in file_paths:
        with open(file_path, 'r') as f:
            data.append(json.load(f))
    return data


def extract_data(data):
    samples = []
    samples_by_technology = {"NANOPORE": [], "ILLUMINA": [], "PACBIO": []}
    bracken_data_all = {}
    serotypes_all = {}
    mlst_all = {}
    insert_sizes_all = {}
    min_insert_size_length = 1000
    busco_data_all = []
    status = {
        "pass": "pass",
        "warn": "warn",
        "fail": "fail",
        "missing": "missing"
    }

    for jdata in data:
        sample = jdata["sample"]
        samples.append(sample)
        qc = jdata["qc"]
        contaminant_status = check_contaminant_status(jdata)
        fastp_q30, fastp_q30_status = check_fastp_q30(jdata)
        nanostat_q15, nanostat_read_n50 = check_nanostat(jdata)
        taxon_count = check_taxon_count(jdata)
        bracken_data_all, samples_by_technology = check_bracken(jdata, bracken_data_all, samples_by_technology)
        checkm_status = check_checkm(jdata)
        taxonkit_genus, taxonkit_genus_status = check_taxonkit_genus(jdata)
        samtools, quast, checkm, samtools_insert_size = check_samtools_and_quast(jdata)
        coverage, coverage_40, coverage_40_illumina, coverage_40_nanopore, coverage_40_pacbio = check_coverage(jdata)
        n50, genome_fraction, contigs, plasmids, assembly, assembly_status = check_assembly(jdata)
        serotype_data = check_serotype(jdata)
        busco_data_all = check_busco(jdata, busco_data_all)
        amrfinder_data, amr_counts = check_amrfinder(jdata)
        mlst_data = check_mlst(jdata)
        data.append({
            "sample": sample,
            "status": check_status(qc),
            "quality_illumina": fastp_q30,
            "quality_illumina_status": fastp_q30_status,
            "quality_nanopore": nanostat_q15,
            "nanopore_n50": nanostat_read_n50,
            "taxon_count_illumina": taxon_count["ILLUMINA"],
            "taxon_count_nanopore": taxon_count["NANOPORE"],
            "taxon_count_pacbio": taxon_count["PACBIO"],
            "coverage": coverage,
            "coverage_40": coverage_40,
            "n50": n50,
            "fraction": genome_fraction,
            "contigs": contigs,
            "plasmids": plasmid_count,
            "assembly": assembly,
            "assembly_status": assembly_status,
            "confindr_illumina_status": contaminant_status["illumina"],
            "confindr_nanopore_status": contaminant_status["nanopore"],
            "quast": quast,
            "checkm": checkm,
            "checkm_status": checkm_status,
            "taxonkit_genus": taxonkit_genus,
            "taxonkit_genus_status": taxonkit_genus_status,
            "mlst": mlst_data,
            "amrfinder": amrfinder_data,
            "amrcount": amr_counts,
            "serotype": serotype_data,
            "busco": busco_data_all
        })
        insert_sizes_all[sample] = samtools_insert_size
    return data, bracken_data_all, samples_by_technology, insert_sizes_all


def check_contaminant_status(jdata):
    contaminant_status = {
        "illumina": "missing",
        "nanopore": "missing"
    }
    if "contamination_illumina" in jdata:
        contaminant_status["illumina"] = jdata["contamination_illumina"]
    if "contamination_nanopore" in jdata:
        contaminant_status["nanopore"] = jdata["contamination_nanopore"]
    return contaminant_status


def check_fastp_q30(jdata):
    fastp_q30, fastp_q30_status = "missing", "missing"
    if "fastp" in jdata:
        fastp_q30, fastp_q30_status = jdata["fastp"]["q30"], "pass" if jdata["fastp"]["q30"] >= 80 else "warn"
    return fastp_q30, fastp_q30_status


def check_nanostat(jdata):
    nanostat_q15, nanostat_read_n50 = "missing", "missing"
    if "nanostat" in jdata:
        nanostat_q15, nanostat_read_n50 = jdata["nanostat"]["q15"], jdata["nanostat"]["n50"]
    return nanostat_q15, nanostat_read_n50


def check_taxon_count(jdata):
    taxon_count = {"ILLUMINA": "missing", "NANOPORE": "missing", "PACBIO": "missing"}
    if "taxon_count_illumina" in jdata:
        taxon_count["ILLUMINA"] = jdata["taxon_count_illumina"]
    if "taxon_count_nanopore" in jdata:
        taxon_count["NANOPORE"] = jdata["taxon_count_nanopore"]
    if "taxon_count_pacbio" in jdata:
        taxon_count["PACBIO"] = jdata["taxon_count_pacbio"]
    return taxon_count


def check_bracken(jdata, bracken_data_all, samples_by_technology):
    if "bracken" in jdata:
        for platform, bracken_data in jdata["bracken"].items():
            if platform not in bracken_data_all:
                bracken_data_all[platform] = []
                samples_by_technology[platform] = []
            for taxon, count in bracken_data.items():
                bracken_data_all[platform].append({"taxon": taxon, "count": count})
                samples_by_technology[platform].append(jdata["sample"])
    return bracken_data_all, samples_by_technology

def check_checkm(jdata):
    return "missing" if "checkm" not in jdata else jdata["checkm"]


def check_taxonkit_genus(jdata):
    taxonkit_genus = "missing"
    if "taxonkit_genus" in jdata:
        taxonkit_genus = jdata["taxonkit_genus"]
    return taxonkit_genus, "missing"


def check_samtools_and_quast(jdata):
    samtools, quast, checkm, samtools_insert_size = "missing", "missing", "missing", []
    if "samtools" in jdata and "quast" in jdata:
        samtools = {
            "q30": jdata["samtools"]["q30"],
            "n50": jdata["quast"]["n50"]
        }
        quast = {
            "completeness": jdata["quast"]["completeness"],
            "duplicated": jdata["quast"]["duplicated"],
            "missing": jdata["quast"]["missing"],
            "fragmented": jdata["quast"]["fragmented"]
        }
        checkm = "missing" if "checkm" not in jdata else jdata["checkm"]
        samtools_insert_size = jdata["samtools"]["insert_size"]
    return samtools, quast, checkm, samtools_insert_size


def check_coverage(jdata):
    coverage, coverage_40, coverage_40_illumina, coverage_40_nanopore, coverage_40_pacbio = "missing", "missing", "missing", "missing", "missing"
    if "coverage" in jdata:
        coverage = jdata["coverage"]
        if "40x" in jdata["coverage"]:
            coverage_40 = jdata["coverage"]["40x"]
            if "illumina" in jdata["coverage"]["40x"]:
                coverage_40_illumina = jdata["coverage"]["40x"]["illumina"]
            if "nanopore" in jdata["coverage"]["40x"]:
                coverage_40_nanopore = jdata["coverage"]["40x"]["nanopore"]
            if "pacbio" in jdata["coverage"]["40x"]:
                coverage_40_pacbio = jdata["coverage"]["40x"]["pacbio"]
    return coverage, coverage_40, coverage_40_illumina, coverage_40_nanopore, coverage_40_pacbio


def check_assembly(jdata):
    n50, genome_fraction, contigs, plasmids, assembly, assembly_status = "missing", "missing", "missing", "missing", "missing", "missing"
    if "assembly" in jdata:
        n50 = jdata["assembly"]["n50"]
        genome_fraction = jdata["assembly"]["fraction"]
        contigs = jdata["assembly"]["contigs"]
        plasmids = jdata["assembly"]["plasmids"]
        assembly = jdata["assembly"]["size"]
        assembly_status = "missing" if "assembly_status" not in jdata else jdata["assembly"]["assembly_status"]
    return n50, genome_fraction, contigs, plasmids, assembly, assembly_status


def check_serotype(jdata):
    serotype_data = "missing"
    if "serotype" in jdata:
        serotype_data = jdata["serotype"]
    return serotype_data


def check_busco(jdata, busco_data_all):
    if "busco" in jdata:
        busco_data_all.append({
            "sample": jdata["sample"],
            "completeness": jdata["busco"]["completeness"],
            "missing": jdata["busco"]["missing"],
            "fragmented": jdata["busco"]["fragmented"],
            "duplicated": jdata["busco"]["duplicated"]
        })
    return busco_data_all


def check_amrfinder(jdata):
    amrfinder_data, amr_counts = [], 0
    if "amrfinder" in jdata:
        for amr_entry in jdata["amrfinder"]:
            amrfinder_data.append({
                "gene_symbol": amr_entry["Gene symbol"],
                "amr_class": amr_entry["Class"],
                "amr_subclass": amr_entry["Subclass"],
                "element_type": amr_entry["Element type"],
                "element_subtype": amr_entry["Element subtype"],
                "sequence_name": amr_entry["Sequence name"]
            })
        amr_counts = len(set(entry["Class"] for entry in amrfinder_data))
    return amrfinder_data, amr_counts


def check_mlst(jdata):
    mlst_data = "missing"
    if "mlst" in jdata:
        mlst_data = jdata["mlst"][0]
    return mlst_data


def check_status(qc):
    return "missing" if "status" not in qc else qc["status"]


def plot_bracken(bracken_data_all, samples_by_technology):
    plots = {}
    for platform, data in bracken_data_all.items():
        df = pd.DataFrame(data, index=samples_by_technology[platform])
        plot = px.bar(df, orientation='h', labels={"index": "Samples", "value": "Percentage"}, height=(len(samples_by_technology[platform]) * 25))
        plots[platform] = plot.to_html(full_html=False)
    return plots


def plot_insert_sizes(insert_sizes_all):
    df = pd.DataFrame(insert_sizes_all)
    plot = px.line(df, labels={"index": "Basepairs", "value": "Count"})
    return plot.to_html(full_html=False)


def plot_busco(busco_data_all):
    df = pd.DataFrame(busco_data_all, index=[item["sample"] for item in busco_data_all])
    plot = px.bar(df, orientation='h', labels={"index": "Samples", "value": "Percentage"}, height=(len(busco_data_all) * 25))
    return plot.to_html(full_html=False)


def parse_yaml(yaml_path):
    software = {}
    with open(yaml_path, "r") as yfile:
        lines = [line.rstrip() for line in yfile]
        current_module = ""
        rmod = re.compile('^[A-Za-z0.*/]')
        for line in lines:
            if (rmod.match(line)):
                current_module = line.split(":")[0]
                software[current_module] = []
            else:
                s, v = line.strip().split()
                software[current_module].append(line.strip())
    return software


def render_template(data, template_path, output_path):
    with open(template_path) as template_file:
        j2_template = Environment().from_string(template_file.read())
        output = j2_template.render(data)
        with open(output_path, "w", encoding="utf-8") as output_file:
            output_file.write(output)


def main(yaml_path, template_path, output_path):
    file_paths = [f for f in os.listdir(".") if f.endswith(".json")]
    with ThreadPoolExecutor() as executor:
        data = list(executor.map(read_json_files, [file_paths]))
    data, bracken_data_all, samples_by_technology, insert_sizes_all = extract_data(data)
    plots = {
        **plot_bracken(bracken_data_all, samples_by_technology),
        **{"Insertsizes": plot_insert_sizes(insert_sizes_all)},
        **{"Busco": plot_busco(data)}
    }
    data.update(plots)
    software = parse_yaml(yaml_path)
    data["packages"] = software
    render_template(data, template_path, output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process JSON files and generate an HTML dashboard.')
    parser.add_argument('--input', help='Path to the YAML file', required=True)
    parser.add_argument('--template', help='Path to the Jinja2 template', required=True)
    parser.add_argument('--output', help='Path to the output HTML file', required=True)
    args = parser.parse_args()
    main(args.input, args.template, args.output)