#!/usr/bin/env python
import plotly.express as px
from jinja2 import Template
import datetime
import pandas as pd
import os
import json
import getpass
import re
import argparse


parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--input", help="An input option")
parser.add_argument("--template", help="A JINJA2 template")
parser.add_argument("--version", help="Pipeline version")
parser.add_argument("--call", help="Command line call")
parser.add_argument("--wd", help="work directory")
parser.add_argument("--output")

args = parser.parse_args()

status = {
    "pass": "pass",
    "warn": "warn",
    "fail": "fail",
    "missing": "missing"
}


def check_status(key, qc):
    if key in qc["fail"]:
        return status["fail"]
    elif key in qc["warn"]:
        return status["warn"]
    elif key in qc["missing"]:
        return status["missing"]
    elif key in qc["pass"]:
        return status["pass"]

    print(f"Unknown key found: {key}")
    return status["missing"]


def main(yaml, template, output, version, call, wd):

    # Read all the JSON files we see in this folder
    json_files = [pos_json for pos_json in os.listdir('.') if pos_json.endswith('.json') and "AQUAMIS" not in pos_json]
    json_files.sort()

    data = {}

    data["user"] = getpass.getuser()
    data["date"] = datetime.datetime.now()
    data["version"] = version
    data["call"] = call
    data["wd"] = wd

    data["summary"] = []

    samples = []

    bracken_data_all = {}
    serotypes_all = {}
    mlst_all = {}
    insert_sizes_all = {}
    min_insert_size_length = 1000
    busco_data_all = []

    for json_file in json_files:

        rtable = {}
        messages = []

        with open(json_file) as f:
            jdata = json.load(f)
            f.close

            taxon = jdata["taxon"]
            sample = jdata["sample"]
            samples.append(sample)

            # Pre-computed QC data
            qc = jdata["qc"]
            this_status = qc["call"]
            messages = qc["messages"]

            #############################################
            # Check for contaminated reads using confindr
            #############################################

            contaminated = {
                "illumina": {"contaminated": "-", "confindr_status": status["missing"]},
                "nanopore": {"contaminated": "-", "confindr_status": status["missing"]}
            }

            for platform, confindr in jdata["confindr"].items():
                contaminated[platform]["contaminated"] = "-"

                contaminated[platform]["confindr_status"] = check_status(f"confindr_{platform.lower()}", qc)

                for set in confindr:

                    for read in set:
                        if ":" in read["Genus"]:
                            contaminated[platform]["contaminated"] = read["Genus"]
                        else:
                            contaminated[platform]["contaminated"] = read["NumContamSNVs"]

            ########################
            # Read quality via FastP
            ########################

            fastp_q30_status = check_status("fastp_q30_rate", qc)
            fastp_q30 = "-"
            if "fastp" in jdata:
                fastp_summary = jdata["fastp"]["summary"]
                fastp_q30 = (round(fastp_summary["after_filtering"]["q30_rate"], 2) * 100)

            ##########################
            # Read stats from NanoStat
            ##########################

            nanostat_q15 = "-"
            nanostat_read_n50 = "-"

            if "nanostat" in jdata:
                nanostat_data = jdata["nanostat"]
                nanostat_q15 = int(nanostat_data["Q15"])
                # nanostat_mean_read_length = nanostat_data["mean_read_length"]
                nanostat_read_n50 = nanostat_data["read_length_n50"]

            ####################
            # Get Bracken results
            ####################

            taxon_count = {
                "ILLUMINA": {
                    "count": "-",
                    "status": "missing"
                },
                "NANOPORE": {
                    "count": "-",
                    "status": "missing"
                }
            }

            taxon_count["ILLUMINA"]["status"] = check_status("read_hit1_genus_fraction_ILLUMINA", qc)
            taxon_count["NANOPORE"]["status"] = status["missing"]

            if "bracken" in jdata:

                for platform, bracken in jdata["bracken"].items():

                    tcount = 0

                    # Rather than defining it at the beginning, we check if we need this platform in the results
                    # This avoids having to filter all platforms that werent in this analysis
                    if platform not in bracken_data_all:
                        bracken_data_all[platform] = []

                    bracken_results = {}
                    for tax in bracken:
                        this_taxon = tax["name"].replace('"', '')
                        tcount += 1
                        # The Bracken results are all in quotes, so we need to clean that up and convert to precentage
                        tperc = round((float(tax["fraction_total_reads"].replace('"', '')) * 100), 2)

                        bracken_results[this_taxon] = tperc

                    taxon_count[platform]["count"] = tcount

                    bracken_data_all[platform].append(bracken_results)

            ####################
            # Get CheckM results
            ####################

            checkm = jdata["checkm"]["Contamination"]
            checkm_status = check_status("checkm_contamination", qc)

            #####################
            # Taxonkit results
            #####################
            taxonkit_genus = sorted(jdata["taxonkit"]["genus"], key=lambda d: d['fraction'], reverse=True)
            taxonkit_genus_fraction = round(float(taxonkit_genus[0]["fraction"]) * 100, 2)
            taxonkit_genus_status = check_status("taxonkit_genus_fraction", qc)

            ####################
            # Get samtools stats
            ####################

            samtools = {"mean_insert_size": "-", }
            if ("samtools" in jdata):
                insert_size = float(jdata["samtools"]["insert size average"])
                insert_stdv = float(jdata["samtools"]["insert size standard deviation"])
                samtools["mean_insert_size"] = f"{insert_size} (+/-{insert_stdv})"
                inserts = [int(item) for item in jdata["samtools"]["insert_sizes"]]
                insert_sizes_all[sample] = inserts
                if (len(inserts) < min_insert_size_length):
                    min_insert_size_length = len(inserts)

            ####################
            # Get assembly stats
            ####################

            assembly = round((int(jdata["quast"]["Total length"]) / 1000000), 2)
            assembly_status = check_status("quast_assembly", qc)

            contigs = jdata["quast"]["# contigs"]
            contigs_status = check_status("quast_contigs", qc)

            n50 = round((int(jdata["quast"]["N50"]) / 1000), 2)
            n50_status = check_status("quast_n50", qc)

            genome_fraction = "-"
            genome_fraction_status = status["missing"]

            quast = {}

            quast["duplication"] = "-"
            quast["misassembled"] = "-"
            quast["mismatches"] = "-"

            if "Genome fraction (%)" in jdata["quast"]:
                genome_fraction = round(float(jdata["quast"]["Genome fraction (%)"]), 2)
                quast["duplication"] = jdata["quast"]["Duplication ratio"]
                quast["misassembled"] = jdata["quast"]["# misassembled contigs"]
                quast["mismatches"] = jdata["quast"]["# mismatches per 100 kbp"]

            # Highlight if a reference coverage is less than 90%
            # This might indicate a problem with our assembly (or the mash database...)
            quast["N50"] = n50
            quast["size"] = jdata["quast"]["Total length (>= 0 bp)"]
            quast["N"] = jdata["quast"]["# N's per 100 kbp"]
            quast["largest_contig"] = round((int(jdata["quast"]["Largest contig"]) / 1000), 2)
            quast["contigs_1k"] = jdata["quast"]["# contigs (>= 1000 bp)"]
            quast["contigs_5k"] = jdata["quast"]["# contigs (>= 5000 bp)"]
            quast["size_1k"] = round(float(int(jdata["quast"]["Total length (>= 1000 bp)"]) / 1000000), 2)
            quast["size_5k"] = round(float(int(jdata["quast"]["Total length (>= 5000 bp)"]) / 1000000), 2)
            quast["gc"] = float(jdata["quast"]["GC (%)"])
            quast["gc_status"] = check_status("quast_gc", qc)
            quast["duplication_ratio"] = round(float(jdata["quast"]["Duplication ratio"]), 4)
            quast["duplication_status"] = check_status("quast_duplication", qc)

            #################
            # Get serotype(s)
            #################

            serotype_data = {}
            if "serotype" in jdata:
                serotypes = jdata["serotype"]
                for stool, sresults in serotypes.items():
                    # we skip these tools
                    if (stool in ["stecfinder", "seqsero2"]):
                        continue
                    pathotype = ""
                    pathogenes = ""
                    comment = ""
                    if (stool == "ectyper"):
                        serotype = sresults["Serotype"]
                        pathogenes = sresults["PathotypeGenes"]
                        pathotype = "" if sresults["Pathotype"] == "ND" else sresults["Pathotype"]
                        comment = sresults["StxSubtypes"]
                    elif (stool == "stecfinder"):
                        serotype = sresults["Serotype"]
                        pathogenes = sresults["stx type"]
                    elif (stool == "seqsero2"):
                        serotype = sresults['serogroup']
                        pathogenes = ""
                    elif (stool == "sistr"):
                        serotype = sresults['serogroup']
                        pathogenes = ""
                        pathotype = sresults["serovar"]
                    elif (stool == "kaptive"):
                        serotype = sresults["best_match"]
                    elif (stool == "lissero"):
                        serotype = sresults["SEROTYPE"]
                        pathogenes = ""
                    elif (stool == "btyper3"):
                        serotype = sresults["Adjusted_panC_Group(predicted_species)"]
                        pathogenes = sresults["Bt(genes)"]
                    elif (stool == "sccmec"):
                        serotype = sresults["subtype"]
                        pathogenes = "" if sresults["mecA"] == "-" else "mecA"
                    stool_name = f"{stool} ({taxon})"
                    pathogenes = [f"<a href=https://www.uniprot.org/uniprotkb?query={gene}+AND+(taxonomy_id%3A2) target=_new>{gene}</a>" for gene in pathogenes.split(",")]
                    serotype_data = {"tool": stool, "serotype": serotype, "genes": pathogenes, "pathotype": pathotype, "comment": comment}
                    if (stool_name in serotypes_all):
                        serotypes_all[stool_name].append({"sample": sample, "serotype": serotype, "genes": pathogenes})
                    else:
                        serotypes_all[stool_name] = [{"sample": sample, "serotype": serotype, "genes": pathogenes}]

            # Reference genome
            reference = jdata["reference"]

            ##############
            # Busco scores
            ##############

            busco = jdata["busco"]
            busco_completeness_status = check_status("busco_completeness", qc)
            busco_duplicates_status = check_status("busco_duplicates", qc)
            busco_status = status["missing"]
            if status["warn"] in [busco_completeness_status, busco_duplicates_status]:
                busco_status = status["warn"]
            if status["fail"] in [busco_completeness_status, busco_duplicates_status]:
                busco_status = status["fail"]
            if busco_status == status["missing"]:
                if status["pass"] in [busco_completeness_status, busco_duplicates_status]:
                    busco_status = status["pass"]

            busco_total = int(busco["dataset_total_buscos"])
            busco_completeness = round(((int(busco["C"])) / int(busco_total)), 2) * 100
            busco_fragmented = round((int(busco["F"]) / busco_total), 2) * 100
            busco_missing = round((int(busco["M"]) / busco_total), 2) * 100
            busco_duplicated = round((int(busco["D"]) / busco_total), 2) * 100
            busco["completeness"] = round(busco_completeness, 2)
            busco["duplicated"] = round(busco_duplicated, 2)
            busco_data_all.append({"Complete": busco_completeness, "Missing": busco_missing, "Fragmented": busco_fragmented, "Duplicated": busco_duplicated})

            ##############
            # Amrfinder
            ##############

            amrfinder_data = []
            amr_classes = []
            if "amrfinder" in jdata["amr"]:
                adata = jdata["amr"]["amrfinder"]
                for amr_entry in adata:
                    amrfinder_data.append(
                        {
                            "gene_symbol": amr_entry["Gene symbol"],
                            "amr_class": amr_entry["Class"],
                            "amr_subclass": amr_entry["Subclass"],
                            "element_type": amr_entry["Element type"],
                            "element_subtype": amr_entry["Element subtype"],
                            "sequence_name": amr_entry["Sequence name"]
                        }
                    )
                    if (amr_entry["Element type"] == "AMR"):
                        if (amr_entry["Class"] not in amr_classes):
                            amr_classes.append(amr_entry["Class"])

                amrfinder_data = sorted(amrfinder_data, key=lambda x: x['gene_symbol'])
            # Count unique AMR classes found
            amr_counts = len(amr_classes)

            ##############
            # MLST types
            ##############

            mlst = jdata["mlst"]
            mlst_data = mlst[0]

            for mentry in mlst:
                sequence_type = mentry["sequence_type"]
                scheme = mentry["scheme"]

                scheme_name = f"{scheme} ({taxon})"
                if (scheme_name in mlst_all):
                    mlst_all[scheme_name].append({"sample": sample, "sequence_type": sequence_type})
                else:
                    mlst_all[scheme_name] = [{"sample": sample, "sequence_type": sequence_type}]

            ##############
            # Get coverage(s)
            ##############

            # Set defaults in case this platform wasnt used
            coverage = "-"
            coverage_status = check_status("coverage_total_mean", qc)

            coverage_illumina = "-"
            coverage_illumina_status = check_status("coverage_illumina_mean", qc)

            coverage_nanopore = "-"
            coverage_nanopore_status = check_status("coverage_nanopore_mean", qc)

            coverage_pacbio = "-"
            coverage_pacbio_status = check_status("coverage_pacbio_mean", qc)

            # mean coverages
            if "total" in jdata["mosdepth"]:
                coverage = float(jdata["mosdepth"]["total"]["mean"])

            if "illumina" in jdata["mosdepth"]:
                coverage_illumina = float(jdata["mosdepth"]["illumina"]["mean"])

            if "nanopore" in jdata["mosdepth"]:
                coverage_nanopore = float(jdata["mosdepth"]["nanopore"]["mean"])

            if "pacbio" in jdata["mosdepth"]:
                coverage_pacbio = float(jdata["mosdepth"]["pacbio"]["mean"])

            # fraction covered at 40X
            coverage_40 = "-"
            coverage_40_status = check_status("coverage_total_40x", qc)
            coverage_40_illumina = "-"
            coverage_40_illumina_status = check_status("coverage_illumina_40x", qc)
            coverage_40_nanopore = "-"
            coverage_40_nanopore_status = check_status("coverage_nanopore_40x", qc)
            coverage_40_pacbio = "-"
            coverage_40_pacbio_status = check_status("coverage_pacbio_40x", qc)

            if "mosdepth_global" in jdata:
                if "illumina" in jdata["mosdepth_global"]:
                    if "40" in jdata["mosdepth_global"]["illumina"]:
                        coverage_40_illumina = jdata["mosdepth_global"]["illumina"]["40"]

                if "nanopore" in jdata["mosdepth_global"]:
                    if "40" in jdata["mosdepth_global"]["nanopore"]:
                        coverage_40_nanopore = jdata["mosdepth_global"]["nanopore"]["40"]

                if "pacbio" in jdata["mosdepth_global"]:
                    if "40" in jdata["mosdepth_global"]["pacbio"]:
                        coverage_40_pacbio = jdata["mosdepth_global"]["pacbio"]["40"]

                if "total" in jdata["mosdepth_global"]:
                    if "40" in jdata["mosdepth_global"]["total"]:
                        coverage_40 = jdata["mosdepth_global"]["total"]["40"]

            #########################
            # sample-level dictionary
            #########################

            if not messages:
                messages.append("No values outside of expected range(s).")

            rtable = {
                "sample": sample,
                "messages": "<br> ".join(messages),
                "reference": reference,
                "status": this_status,
                "samtools": samtools,
                "taxon": taxon,
                "busco": busco,
                "quality_illumina": fastp_q30,
                "quality_illumina_status": fastp_q30_status,
                "quality_nanopore": nanostat_q15,
                "nanopore_n50": nanostat_read_n50,
                "busco_status": busco_status,
                "taxon_count_illumina": taxon_count["ILLUMINA"]["count"],
                "taxon_count_illumina_status": taxon_count["ILLUMINA"]["status"],
                "taxon_count_nanopore": taxon_count["NANOPORE"]["count"],
                "taxon_count_nanopore_status": taxon_count["NANOPORE"]["status"],
                "coverage": coverage,
                "coverage_status": coverage_status,
                "coverage_illumina": coverage_illumina,
                "coverage_illumina_status": coverage_illumina_status,
                "coverage_nanopore": coverage_nanopore,
                "coverage_nanopore_status": coverage_nanopore_status,
                "coverage_pacbio": coverage_pacbio,
                "coverage_pacbio_status": coverage_pacbio_status,
                "coverage_40": coverage_40,
                "coverage_40_status": coverage_40_status,
                "coverage_40_illumina": coverage_40_illumina,
                "coverage_40_illumina_status": coverage_40_illumina_status,
                "coverage_40_nanopore": coverage_40_nanopore,
                "coverage_40_nanopore_status": coverage_40_nanopore_status,
                "coverage_40_pacbio": coverage_40_pacbio,
                "coverage_40_pacbio_status": coverage_40_pacbio_status,
                "n50": n50,
                "n50_status": n50_status,
                "fraction": genome_fraction,
                "fraction_status": genome_fraction_status,
                "contigs": contigs,
                "contigs_status": contigs_status,
                "assembly": assembly,
                "assembly_status": assembly_status,
                "contamination_illumina": contaminated["illumina"]["contaminated"],
                "confindr_illumina_status": contaminated["illumina"]["confindr_status"],
                "contamination_nanopore": contaminated["nanopore"]["contaminated"],
                "confindr_nanopore_status": contaminated["nanopore"]["confindr_status"],
                "quast": quast,
                "checkm": checkm,
                "checkm_status": checkm_status,
                "taxonkit_genus": taxonkit_genus_fraction,
                "taxonkit_genus_status": taxonkit_genus_status,
                "mlst": mlst_data,
                "amrfinder": amrfinder_data,
                "amrcount": amr_counts,
                "serotype": serotype_data
            }

        data["summary"].append(rtable)

    #############
    # Plots
    #############

    # Bracken abundances
    if "bracken" in jdata:
        if "ILLUMINA" in bracken_data_all:
            kdata = pd.DataFrame(data=bracken_data_all["ILLUMINA"], index=samples)
            plot_labels = {"index": "Samples", "value": "Percentage"}
            h = (len(samples) * 25) if len(samples) > 10 else (200 + len(samples) * 50)
            fig = px.bar(kdata, orientation='h', labels=plot_labels, height=h)
            data["Bracken_ILLUMINA"] = fig.to_html(full_html=False)
        if "NANOPORE" in bracken_data_all:
            print("Creating Bracken ONT graph")
            kdata = pd.DataFrame(data=bracken_data_all["NANOPORE"], index=samples)
            plot_labels = {"index": "Samples", "value": "Percentage"}
            h = (len(samples) * 25) if len(samples) > 10 else (200 + len(samples) * 50)
            fig = px.bar(kdata, orientation='h', labels=plot_labels, height=h)
            data["Bracken_NANOPORE"] = fig.to_html(full_html=False)
        if "PACBIO" in bracken_data_all:
            print("Creating Bracken Pacbio graph")
            kdata = pd.DataFrame(data=bracken_data_all["PACBIO"], index=samples)
            plot_labels = {"index": "Samples", "value": "Percentage"}
            h = (len(samples) * 25) if len(samples) > 10 else (200 + len(samples) * 50)
            fig = px.bar(kdata, orientation='h', labels=plot_labels, height=h)
            data["Bracken_PACBIO"] = fig.to_html(full_html=False)

    # Insert size distribution
    if "samtools" in jdata:
        # Crop all insert size histograms to the shortest common length
        insert_sizes_all_cropped = {}
        for s, ins in insert_sizes_all.items():
            list_end = min_insert_size_length - 1
            insert_sizes_all_cropped[s] = ins[:list_end]

        plot_labels = {"index": "Basepairs", "value": "Count"}
        hdata = pd.DataFrame(insert_sizes_all_cropped)
        hfig = px.line(hdata, labels=plot_labels)
        data["Insertsizes"] = hfig.to_html(full_html=False)

    # busco score
    if busco_data_all:
        # Draw the busco stats graph
        bdata = pd.DataFrame(data=busco_data_all, index=samples)
        plot_labels = {"index": "Samples", "value": "Percentage"}
        h = (len(samples) * 25) if len(samples) > 10 else (200 + len(samples) * 50)
        print(h)
        fig = px.bar(bdata, orientation='h', labels=plot_labels, height=h)
        data["Busco"] = fig.to_html(full_html=False)

    data["serotypes"] = serotypes_all

    data["mlst"] = mlst_all

    ##############################
    # Parse the versions YAML file
    ##############################

    software = {}
    current_module = ""
    rmod = re.compile('^[A-Za-z0.*/]')
    with open(yaml, "r") as yfile:
        lines = [line.rstrip() for line in yfile]
        for line in lines:
            if (rmod.match(line)):
                current_module = line.split(":")[0]
                software[current_module] = []
            else:
                s, v = line.strip().split()
                software[current_module].append(line.strip())

    data["packages"] = software

    ########################
    # Render Jinja2 template
    ########################

    with open(output, "w", encoding="utf-8") as output_file:
        with open(template) as template_file:
            j2_template = Template(template_file.read())
            output_file.write(j2_template.render(data))


if __name__ == '__main__':
    main(args.input, args.template, args.output, args.version, args.call, args.wd)
