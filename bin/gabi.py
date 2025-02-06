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
parser.add_argument("--references", help="Reference values for various taxa")
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


def main(yaml, template, output, reference, version, call, wd):

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

    kraken_data_all = {"ILLUMINA": [], "NANOPORE": []}
    serotypes_all = {}
    mlst_all = {}
    insert_sizes_all = {}
    min_insert_size_length = 1000
    busco_data_all = []

    with open(reference) as r:
        ref_data = json.load(r)["thresholds"]
        r.close

    for idx, json_file in enumerate(json_files):

        rtable = {}

        with open(json_file) as f:
            jdata = json.load(f)
            f.close

            # Track the sample status
            this_status = status["pass"]
            taxon = jdata["taxon"]
            genus, species = taxon.split(" ")

            # The reference data has thresholds for genus and species level; but not always
            # We take the species first, genus second (if any) and iterate over this list in
            # the check functions. Whatever hits first, gets returned (species, when in doubt)
            this_refs = []
            if species in ref_data:
                this_refs.append(ref_data[species])
            elif genus in ref_data:
                this_refs.append(ref_data[genus])
            else:
                this_refs = [{}]

            #############################################
            # Check for contaminated reads using confindr
            #############################################
            if "confindr" in jdata:
                contaminated = "-"
                confindr = jdata["confindr"]
                confindr_status = status["missing"]

                for set in confindr:

                    if confindr_status == status["missing"]:
                        confindr_status = status["pass"]

                    for read in set:
                        if read["ContamStatus"] == "True":
                            contaminated = True

                            if "PercentContam" in read:
                                if (read["PercentContam"] == "ND"):
                                    perc = "ND"
                                    if ":" in read["Genus"]:
                                        perc = read["Genus"]
                                    this_status = status["fail"]
                                    confindr_status = status["fail"]
                                    contaminated = perc
                                else:
                                    perc = float(read["PercentContam"])

                                    if (perc > contaminated):
                                        contaminated = perc

                                    if (perc >= 10.0):
                                        confindr_status = status["fail"]
                                    elif (perc > 0.0 and confindr_status == status["pass"]):
                                        confindr_status = status["warn"]
                                    else:
                                        if confindr_status == status["pass"]:
                                            confindr_status = status["warn"]
                                        contaminated = "ND"

                            else:
                                # If no percentages are given, we may still have a inter-species
                                # contamination scenario, which we need to report
                                if ":" in read["Genus"]:
                                    contaminated = read["Genus"]
                                    confindr_status = status["fail"]
                                else:
                                    contaminated = "ND"
                                    confindr_status = status["warn"]

            # All the relevant values and optional status classes
            sample = jdata["sample"]
            samples.append(sample)

            ########################
            # Read quality via FastP
            ########################

            fastp_q30 = "-"
            fastp_q30_status = status["missing"]

            if "fastp" in jdata:
                fastp_q30_status = status["pass"]
                fastp_summary = jdata["fastp"]["summary"]
                fastp_q30 = (round(fastp_summary["after_filtering"]["q30_rate"], 2) * 100)
                if fastp_q30 < 85:
                    fastp_q30_status = status["warn"]

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
            # Get Kraken results
            ####################

            taxon_status = status["missing"]
            taxon_count = "-"
            taxon_count_status = status["missing"]

            if "kraken" in jdata:

                for platform, kraken in jdata["kraken"].items():
                    taxon_perc = float(kraken[0]["percentage"])
                    if taxon_perc >= 90.0:
                        taxon_status = status["pass"]
                    elif taxon_perc >= 70.0:
                        taxon_status = status["warn"]
                    else:
                        taxon_status = status["fail"]

                    taxon_count = 0
                    taxon_count_status = status["pass"]

                    kraken_results = {}
                    for tax in kraken:
                        this_taxon = tax["taxon"]
                        tperc = float(tax["percentage"])

                        kraken_results[this_taxon] = tperc

                        if (tperc > 5.0):
                            taxon_count += 1

                    kraken_data_all[platform].append(kraken_results)

                    if (taxon_count > 3):
                        taxon_count_status = status["fail"]
                    elif (taxon_count > 1):
                        taxon_count_status = status["warn"]

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
            assembly_status = check_assembly(this_refs, int(jdata["quast"]["Total length"]))

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

            contigs = int(jdata["quast"]["# contigs"])
            contigs_status = check_contigs(this_refs, int(jdata["quast"]["# contigs"]))

            n50 = round((int(jdata["quast"]["N50"]) / 1000), 2)
            n50_status = check_n50(this_refs, int(jdata["quast"]["N50"]))

            quast["size"] = jdata["quast"]["Total length (>= 0 bp)"]
            quast["N"] = jdata["quast"]["# N's per 100 kbp"]
            quast["largest_contig"] = round((int(jdata["quast"]["Largest contig"]) / 1000), 2)
            quast["contigs_1k"] = jdata["quast"]["# contigs (>= 1000 bp)"]
            quast["contigs_5k"] = jdata["quast"]["# contigs (>= 5000 bp)"]
            quast["size_1k"] = round(float(int(jdata["quast"]["Total length (>= 1000 bp)"]) / 1000000), 2)
            quast["size_5k"] = round(float(int(jdata["quast"]["Total length (>= 5000 bp)"]) / 1000000), 2)
            quast["gc"] = float(jdata["quast"]["GC (%)"])
            quast["gc_status"] = check_gc(this_refs, float(jdata["quast"]["GC (%)"]))

            #################
            # Get serotype(s)
            #################

            if "serotype" in jdata:
                serotypes = jdata["serotype"]
                for sentry in serotypes:
                    for stool, sresults in sentry.items():
                        if (stool == "ectyper"):
                            serotype = sresults["Serotype"]
                        elif (stool == "Stecfinder"):
                            serotype = sresults["Serotype"]
                        elif (stool == "SeqSero2"):
                            serotype = f"{sresults['Predicted serotype']} ({sresults['Predicted antigenic profile']})"
                        elif (stool == "Sistr"):
                            serotype = f"{sresults['serovar']} ({sresults['serogroup']})"
                        elif (stool == "Lissero"):
                            serotype = sresults["SEROTYPE"]

                    stool_name = f"{stool} ({taxon})"
                    if (stool_name in serotypes_all):
                        serotypes_all[stool_name].append({"sample": sample, "serotype": serotype})
                    else:
                        serotypes_all[stool_name] = [{"sample": sample, "serotype": serotype}]

            # Reference genome
            reference = jdata["reference"]

            ##############
            # Busco scores
            ##############

            busco = jdata["busco"]
            busco_status = status["missing"]
            busco_total = int(busco["dataset_total_buscos"])
            busco_completeness = round(((int(busco["C"])) / int(busco_total)), 2) * 100
            busco_fragmented = round((int(busco["F"]) / busco_total), 2) * 100
            busco_missing = round((int(busco["M"]) / busco_total), 2) * 100
            busco_duplicated = round((int(busco["D"]) / busco_total), 2) * 100
            busco["completeness"] = busco_completeness
            busco_data_all.append({"Complete": busco_completeness, "Missing": busco_missing, "Fragmented": busco_fragmented, "Duplicated": busco_duplicated})

            if (busco_completeness > 90.0):
                busco_status = status["pass"]
            elif (busco_completeness > 80.0):
                busco_status = status["warn"]
            else:
                busco_status = status["fail"]

            # Warn if there are duplications in the gene set and busco wasnt already failed
            if (busco_duplicated > 5.0) & (busco_status != status["fail"]):
                busco_status = status["warn"]

            ##############
            # MLST types
            ##############

            mlst = jdata["mlst"]

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

            coverage = "-"
            coverage_status = status["missing"]

            coverage_illumina = "-"
            coverage_illumina_status = status["missing"]

            coverage_nanopore = "-"
            coverage_nanopore_status = status["missing"]

            coverage_pacbio = "-"
            coverage_pacbio_status = status["missing"]

            # mean coverages

            if "total" in jdata["mosdepth"]:
                coverage = float(jdata["mosdepth"]["total"]["mean"])
                if coverage >= 40.0:
                    coverage_status = status["pass"]
                elif coverage >= 20.0:
                    coverage_status = status["warn"]
                else:
                    coverage_status = status["fail"]

            if "illumina" in jdata["mosdepth"]:
                coverage_illumina = float(jdata["mosdepth"]["illumina"]["mean"])
                if coverage_illumina >= 40.0:
                    coverage_illumina_status = status["pass"]
                elif coverage_illumina >= 20.0:
                    coverage_illumina_status = status["warn"]
                else:
                    coverage_illumina_status = status["fail"]

            if "nanopore" in jdata["mosdepth"]:
                coverage_nanopore = float(jdata["mosdepth"]["nanopore"]["mean"])
                if coverage_nanopore >= 40.0:
                    coverage_nanopore_status = status["pass"]
                elif coverage_nanopore >= 20.0:
                    coverage_nanopore_status = status["warn"]
                else:
                    coverage_nanopore_status = status["fail"]

            if "pacbio" in jdata["mosdepth"]:
                coverage_pacbio = float(jdata["mosdepth"]["pacbio"]["mean"])
                if coverage_pacbio >= 40.0:
                    coverage_pacbio_status = status["pass"]
                elif coverage_pacbio >= 20.0:
                    coverage_pacbio_status = status["warn"]
                else:
                    coverage_pacbio_status = status["fail"]

            # fraction covered at 40X

            coverage_40 = "-"
            coverage_40_status = status["missing"]
            coverage_40_illumina = "-"
            coverage_40_illumina_status = status["missing"]
            coverage_40_nanopore = "-"
            coverage_40_nanopore_status = status["missing"]
            coverage_40_pacbio = "-"
            coverage_40_pacbio_status = status["missing"]

            if "mosdepth_global" in jdata:
                if "illumina" in jdata["mosdepth_global"]:
                    coverage_40_illumina = jdata["mosdepth_global"]["illumina"]["40"]
                    if coverage_40_illumina < 90:
                        coverage_40_illumina_status = status["warn"]
                    else:
                        coverage_40_illumina_status = status["pass"]
                if "nanopore" in jdata["mosdepth_global"]:
                    coverage_40_nanopore = jdata["mosdepth_global"]["nanopore"]["40"]
                    if coverage_40_nanopore < 90:
                        coverage_40_nanopore_status = status["warn"]
                    else:
                        coverage_40_nanopore_status = status["pass"]
                if "pacbio" in jdata["mosdepth_global"]:
                    coverage_40_pacbio = jdata["mosdepth_global"]["pacbio"]["40"]
                    if coverage_40_pacbio < 90:
                        coverage_40_pacbio_status = status["warn"]
                    else:
                        coverage_40_pacbio_status = status["pass"]
                if "total" in jdata["mosdepth_global"]:
                    coverage_40 = jdata["mosdepth_global"]["total"]["40"]
                    if coverage_40 < 90:
                        coverage_40_status = status["warn"]
                    elif coverage_40 < 75:
                        coverage_40_status = status["fail"]
                    else:
                        coverage_40_status = status["pass"]

            ######################################
            # Set the overall status of the sample
            ######################################

            # The metrics that by themselves determine overall status:
            for estatus in [confindr_status, taxon_count_status, assembly_status]:
                # if any one metric failed, the whole sample failed
                if estatus == status["fail"]:
                    this_status = estatus
                # if a metric is dubious, the entire sample is dubious, unless it already failed or warned
                elif (estatus == status["warn"]) & (this_status == status["pass"]):
                    this_status = estatus

            # The other metrics should at most warn, but never fail the sample
            for estatus in [busco_status, contigs_status]:
                if (estatus != status["missing"]) & (this_status != status["fail"]) & (estatus != status["pass"]):
                    this_status = status["warn"]

            #########################
            # sample-level dictionary
            #########################

            rtable = {
                "sample": sample,
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
                "taxon_status": taxon_status,
                "taxon_count": taxon_count,
                "taxon_count_status": taxon_count_status,
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
                "contamination": contaminated,
                "confindr_status": confindr_status,
                "quast": quast,
            }

        data["summary"].append(rtable)

    #############
    # Plots
    #############

    # Kraken abundances
    if "kraken" in jdata:
        if "ILLUMINA" in jdata["kraken"]:
            kdata = pd.DataFrame(data=kraken_data_all["ILLUMINA"], index=samples)
            plot_labels = {"index": "Samples", "value": "Percentage"}
            h = len(samples) * 25 if len(samples) > 10 else 400
            fig = px.bar(kdata, orientation='h', labels=plot_labels, height=h)
            data["Kraken_ILLUMINA"] = fig.to_html(full_html=False)
        if "NANOPORE" in jdata["kraken"]:
            kdata = pd.DataFrame(data=kraken_data_all["NANOPORE"], index=samples)
            plot_labels = {"index": "Samples", "value": "Percentage"}
            h = len(samples) * 25 if len(samples) > 10 else 400
            fig = px.bar(kdata, orientation='h', labels=plot_labels, height=h)
            data["Kraken_NANOPORE"] = fig.to_html(full_html=False)

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
        h = len(samples) * 25 if len(samples) > 10 else 400
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


def check_assembly(refs, query):

    for ref in refs:

        if "Total length" in ref:

            ref_intervals = [int(x) for x in ref["Total length"][0]["interval"]]

            # assembly falls between the allowed sizes
            if (any(x >= query for x in ref_intervals) and any(x <= query for x in ref_intervals)):
                return status["pass"]
            elif (any(x >= (query * 0.8) for x in ref_intervals) and any(x <= (query * 1.2) for x in ref_intervals)):
                return status["warn"]
            else:
                return status["fail"]

    return status["missing"]


def check_contigs(refs, query):

    for ref in refs:

        if "# contigs (>= 0 bp)" in ref:

            ref_intervals = [int(x) for x in ref["# contigs (>= 0 bp)"][0]["interval"]]

            if (any(x >= query for x in ref_intervals)):
                return status["pass"]
            elif (any((x * 1.2) >= query for x in ref_intervals)):
                return status["warn"]
            else:
                return status["fail"]

    return status["missing"]


def check_n50(refs, query):

    for ref in refs:

        if "N50" in ref:

            ref_intervals = [int(x) for x in ref["N50"][0]["interval"]]

            if (any(x <= query for x in ref_intervals)):
                return status["pass"]
            elif (any((x * 0.8) <= query for x in ref_intervals)):
                return status["warn"]
            else:
                return status["fail"]

        return status["missing"]


def check_gc(refs, query):

    for ref in refs:

        ref_intervals = [float(x) for x in ref["GC (%)"][0]["interval"]]

        # check if gc falls within expected range, or range +/- 5% - else fail
        if (any(x >= query for x in ref_intervals) and any(x <= query for x in ref_intervals)):
            return status["pass"]
        elif (any(x >= (query * 0.95) for x in ref_intervals) and any(x <= (query * 1.05) for x in ref_intervals)):
            return status["warn"]
        else:
            return status["fail"]

    return status["missing"]


if __name__ == '__main__':
    main(args.input, args.template, args.output, args.references, args.version, args.call, args.wd)
