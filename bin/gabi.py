#!/usr/bin/env python
import plotly.express as px
from jinja2 import Template
import pandas as pd
import os,json,re
import argparse


parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--input", help="An input option")
parser.add_argument("--template", help="A JINJA2 template")
parser.add_argument("--output")
args = parser.parse_args()

status = {
    "pass": "pass",
    "warn": "warn",
    "fail": "fail",
    "missing": "missing"
}

# curated thresholds for relevant bacterial taxa
references = {
    "genomes" : {
        "Escherichia coli": {
            "size": 4.8,
            "max_contigs": 500,
            "n50": 200
        },
        "Listeria monocytogenes": {
            "size": 2.9,
            "max_contigs": 300,
            "n50": 200
        },
        "Salmonella enterica": {
            "size": 4.5,
            "max_contigs": 300,
            "n50": 200
        },
        "Staphylococcus aureus": {
            "size": 2.8,
            "max_contigs": 150,
            "n50": 200
        },
        "Bacillus cereus": {
            "size": 6.0,
            "max_contigs": 300,
            "n50": 200
        },
        "Campylobacter jejuni": {
            "size": 1.6,
            "max_contigs": 300,
            "n50": 150
        },
        "Campylobacter lari": {
            "size": 1.6,
            "max_contigs": 300,
            "n50": 150
        },
        "Campylobacter coli": {
            "size": 1.6,
            "max_contigs": 300,
            "n50": 150
        },
        "Klebsiella pneumoniae": {
            "size": 5.7,
            "max_contigs": 300,
            "n50": 200
        }
    },
    "taxon": {
        "pass": 80.0,
        "warn": 60.0,
        "fail": 0.0
    }
}

def main(yaml, template, output):

    json_files = [pos_json for pos_json in os.listdir('.') if pos_json.endswith('.json')]

    data = {}
    data["summary"] = [ ]

    samples = []

    kraken_data_all = []

    for idx, json_file in enumerate(json_files):

        rtable = {}

        with open(json_file) as f:
            jdata = json.load(f)
            f.close

            # Track the sample status
            this_status = status["pass"]

            # Check for contaminated reads using confindr
            contaminated = 0
            confindr = jdata["confindr"]
            confindr_status = status["pass"]

            for set in confindr:
                for read in set:
                    if read["ContamStatus"] == "True":
                        contaminated = True
                        perc = float(read["PercentContam"])

                        if (perc > contaminated):
                            contaminated = perc

                        if (perc >= 10.0):
                            confindr_status = status["fail"]
                            this_status = status["fail"]
                        elif (perc > 0.0 and confindr_status == status["pass"]):
                            confindr_status = status["warn"]
                            if (this_status == status["pass"]):
                                this_status = status["warn"]

            # All the relevant values and optional status classes
            sample = jdata["sample"]
            samples.append(sample)

            # Get Kraken results
            taxon = jdata["kraken"][0]["taxon"]
            taxon_perc = float(jdata["kraken"][0]["percentage"])
            if taxon_perc >= 80.0:
                taxon_status = status["pass"]
            elif taxon_perc >= 60.0:
                taxon_status = status["warn"]
                if (this_status == status["pass"]):
                    this_status = status["warn"]
            else:
                taxon_status = status["fail"]
                # probably not an outright fail for the whole analysis
                if (this_status == status["pass"]):
                    this_status = status["warn"]

            taxon_count = 0
            taxon_count_status = status["pass"]

            kraken_results = {}
            for tax in jdata["kraken"]:
                this_taxon = tax["taxon"]
                tperc = float(tax["percentage"])

                kraken_results[this_taxon] = tperc
               
                if (tperc > 10.0 ):
                    taxon_count += 1
            
            kraken_data_all.append(kraken_results)

            if (taxon_count > 3):
                taxon_count_status = status["fail"]
                this_status = status["fail"]
            elif (taxon_count > 1 ):
                taxon_count_status = status["warn"]
                if (this_status == status["pass"]):
                    this_status = status["warn"]

            # Get assembly stats
            assembly = round((int(jdata["quast"]["Total length"])/1000000),2)
            assembly_status = check_assembly(taxon,assembly)

            genome_fraction = round(float(jdata["quast"]["Genome fraction (%)"]),2)

            if (genome_fraction > 90):
                genome_fraction_status = status["pass"]
            elif (genome_fraction > 80):
                genome_fraction_status = status["warn"]
                if (status == status["pass"]):
                    this_status = status["warn"]
            else:
                genome_fraction_status = status["fail"]
                this_status = status["fail"]

            contigs = int(jdata["quast"]["# contigs"])
            contigs_status = check_contigs(taxon,contigs)

            n50 = round((int(jdata["quast"]["N50"])/1000),2)
            n50_status = check_n50(taxon,n50)

            quast = {}
            quast["size"] = jdata["quast"]["Total length (>= 0 bp)"]
            quast["duplication"] = jdata["quast"]["Duplication ratio"]
            quast["N"] = jdata["quast"]["# N's per 100 kbp"]
            quast["mismatches"] = jdata["quast"]["# mismatches per 100 kbp"]
            quast["largest_contig"] = jdata["quast"]["Largest contig"]
            quast["misassembled"] = jdata["quast"]["# misassembled contigs"]

            # Get coverage(s)
            coverage_illumina = "-"
            coverage_illumina_status = status["missing"]

            coverage_nanopore = "-"
            coverage_nanopore_status = status["missing"]

            coverage_pacbio = "-"
            coverage_pacbio_status = status["missing"]

            coverage = float(jdata["mosdepth"]["total"]["mean"])
            if coverage >= 40:
                coverage_status = status["pass"]
            elif coverage >= 20:
                coverage_status = status["warn"]
            else:
                coverage_status = status["pass"]
            
            if ("illumina" in jdata["mosdepth"]):

                coverage_illumina = float(jdata["mosdepth"]["illumina"]["mean"])
                if coverage_illumina >= 40:
                    coverage_illumina_status = status["pass"]
                elif coverage_illumina >= 20:
                    coverage_illumina_status = status["warn"]
                else:
                    coverage_illumina_status = status["fail"]

            if ("nanopore" in jdata["mosdepth"]):
                coverage_nanopore = float(jdata["mosdepth"]["nanopore"]["mean"])
                if coverage_nanopore >= 40:
                    coverage_nanopore_status = status["pass"]
                elif coverage_nanopore >= 20:
                    coverage_nanopore_status = status["warn"]
                else:
                    coverage_nanopore_status = status["fail"]

            if ("pacbio" in jdata["mosdepth"]):
                coverage_pacbio = float(jdata["mosdepth"]["pacbio"]["mean"])
                if coverage_pacbio >= 40:
                    coverage_pacbio_status = status["pass"]
                elif coverage_pacbio >= 20:
                    coverage_pacbio_status = status["warn"]
                else:
                    coverage_pacbio_status = status["fail"]
            
            # sample-level dictionary
            rtable = { 
                "sample": sample,
                "status": this_status,
                "taxon": taxon,
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
                "quast": quast
            }
        
        data["summary"].append(rtable)    

    kdata = pd.DataFrame(data=kraken_data_all,index=samples)
    plot_labels = { "index": "Samples", "value": "Percentage"}
    fig = px.bar(kdata,orientation='h',labels=plot_labels)

    data["Kraken"] = fig.to_html(full_html=False)

    # Parse the versions YAML file
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
                s,v = line.strip().split()
                software[current_module].append(line.strip())
            
    data["packages"] = software

    with open(output, "w", encoding="utf-8") as output_file:
        with open(template) as template_file:
            j2_template = Template(template_file.read())
            output_file.write(j2_template.render(data))

def check_assembly(taxon,size):
    if taxon in references["genomes"]:
        ref = float(references["genomes"][taxon]["size"])
        if (size >= (ref*0.9)) and (size <= (ref*1.1)):
            return status["pass"]
        elif (size >= (ref*0.8) and size <= (ref*1.2)):
            return status["warn"]
        else:
            return status["fail"]
    else:
        return status["missing"]

def check_contigs(taxon,contigs):
    if taxon in references["genomes"]:
        ref = int(references["genomes"][taxon]["max_contigs"])
        if contigs <= ref:
            return status["pass"]
        elif contigs <= (ref*1.1):
            return status["warn"]
        else:
            return status["fail"]
    else:
        return status["missing"]

def check_n50(taxon,n50):
    if taxon in references["genomes"]:
        ref = int(references["genomes"][taxon]["n50"])
        if n50 >= ref:
            return status["pass"]
        elif n50 >= (ref*0.9):
            return status["warn"]
        else:
            return status["fail"]
    else:
        return status["missing"] 

if __name__ == '__main__':
    main(args.input, args.template, args.output)
