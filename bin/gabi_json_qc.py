#!/usr/bin/env python
import json
import argparse

parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--input", "-i")
parser.add_argument("--refs", "-r")
parser.add_argument("--output", "-o")


args = parser.parse_args()


status = {
    "pass": "pass",
    "warn": "warn",
    "fail": "fail",
    "missing": "missing"
}


def bracken_stats(hits):

    genus_stats = {}
    for tax in hits:
        this_taxon = tax["name"].replace('"', '')
        genus = this_taxon.split(" ")[0]
        # The Bracken results are all in quotes, so we need to clean that up and convert to precentage
        tperc = round(float(tax["fraction_total_reads"].replace('"', '')), 2)

        if genus in genus_stats:
            genus_stats[genus] += tperc
        else:
            genus_stats[genus] = tperc

    return genus_stats


def parse_json(j):
    with open(j, "r") as json_file:
        data = json.loads(json_file.read())

    return data


def check(key, refs, query):

    for ref in refs:
        if key in ref:
            thresholds = sorted([float(x) for x in ref[key][0]["interval"]])
            bins = ref[key][0]["binscore"]

            # Must be smaller than this reference value
            if bins == [0, 1]:
                if query < thresholds[0]:
                    return status["pass"]
                elif query < (thresholds[0] * 1.1):
                    return status["warn"]
                else:
                    return status["fail"]
            # Must be larger than this reference value
            elif bins == [1, 0]:
                if query > thresholds[0]:
                    return status["pass"]
                elif (query * 1.1) > thresholds[0]:
                    return status["warn"]
                else:
                    return status["fail"]
            # Must be within this interval
            elif bins == [1, 0, 1]:
                if (any(x >= query for x in thresholds) and any(x <= query for x in thresholds)):
                    return status["pass"]
                elif (any(x >= (query * 0.95) for x in thresholds) and any(x <= (query * 1.05) for x in thresholds)):
                    return status["warn"]
                else:
                    return status["fail"]
            # Must be within these intervals, with outlier limits
            elif bins == [2, 1, 0, 1, 2]:
                low, low_ok, high_ok, high = thresholds

                if ((query >= low_ok) and (query <= high_ok)):
                    return status["pass"]
                elif (query < low):
                    return status["fail"]
                elif (query < low_ok):
                    return status["warn"]
                elif (query > high):
                    return status["fail"]
                elif (query > high_ok):
                    return status["warn"]
            # Must be above these limits, with outlier
            elif bins == [2, 1, 0]:
                if query >= thresholds[-1]:
                    return status["pass"]
                elif query >= thresholds[0]:
                    return status["warn"]
                else:
                    return status["fail"]
            # Must be below these values, with outlier
            elif bins == [0, 1, 2]:
                if query >= thresholds[-1]:
                    return status["fail"]
                else:
                    return status["warn"]
            elif bins == [2, 0]:

                if query >= thresholds[-1]:
                    return status["pass"]
                else:
                    return status["fail"]

    return status["missing"]


def main(input, refs, output):

    qc_calls = {"fail": [], "warn": [], "pass": [], "missing": [], "messages": []}

    data = parse_json(input)

    ref_data = parse_json(refs)["thresholds"]

    taxon = data["taxon"]
    genus = taxon.split(" ")[0]

    this_refs = [{}]
    if taxon in ref_data:
        this_refs.append(ref_data[taxon])

    if genus in ref_data:
        this_refs.append(ref_data[genus])

    this_refs.append(ref_data["all Species"])

    # confindr
    for platform, confindr in data["confindr"].items():
        if len(confindr) > 0:

            platform_contaminated = status["missing"]
            contaminated = status["missing"]

            for set in confindr:
                if platform_contaminated == "missing":
                    platform_contaminated = status["pass"]

                if contaminated == status["missing"]:
                    contaminated = status["pass"]

                for read in set:
                    contam_type = "intra-species"

                    if read["BasesExamined"] == 0:
                        contaminated = status["missing"]
                        m = "No ConfindR databases for this species available, contamination check skipped."
                        if m not in qc_calls["messages"]:
                            qc_calls["messages"].append(m)
                        platform_contaminated = status["missing"]
                    else:
                        if ":" in read["Genus"]:
                            contam_type = "inter-species"

                        if read["ContamStatus"] == "True":
                            contaminated = check("NumContamSNVs", this_refs, int(read["NumContamSNVs"]))

                            if contaminated == status["fail"]:
                                m = f"ConfindR: Contamination ({contam_type}) detected in {platform} reads {read['Sample']}"
                                if m not in qc_calls["messages"]:
                                    qc_calls["messages"].append(m)
                            else:
                                m = f"ConfindR: Low levels of contamination ({contam_type}) detected in {platform} reads {read['Sample']}"
                                if m not in qc_calls["messages"]:
                                    qc_calls["messages"].append(m)
                        else:
                            contaminated = status["pass"]

                        if contaminated == status["fail"]:
                            platform_contaminated = status["fail"]
                        elif contaminated == status["warn"] and platform_contaminated != status["fail"]:
                            platform_contaminated = status["warn"]

            qc_calls[platform_contaminated].append(f"confindr_{platform}")

    # fastp
    if "fastp" in data:
        fastp_q30_rate = data["fastp"]["summary"]["after_filtering"]["q30_rate"]
        fastp_q30_status = check("q30_rate_after", this_refs, fastp_q30_rate)
        qc_calls[fastp_q30_status].append("fastp_q30_rate")
        if fastp_q30_rate < 0.75:
            qc_calls["messages"].append("Fastp: Illumina Q30 fraction below 75%")

    # bracken
    for platform, bracken in data["bracken"].items():

        genus_stats = bracken_stats(bracken)

        # Sort abundances from high to low
        abundances = sorted(genus_stats.items(), key=lambda x: x[1], reverse=True)
        first_hit = abundances[0]
        first_hit_status = check("read_hit1_species_fraction", this_refs, first_hit[1])
        first_hit_status = status["warn"] if first_hit_status == status["fail"] else first_hit_status
        qc_calls[first_hit_status].append(f"read_hit1_genus_fraction_{platform}")
        if first_hit_status == status["warn"]:
            qc_calls["messages"].append("Bracken: Read abundance of dominant genus below threshold - possible contamination issue.")

        # If there is a second genus detected
        if len(abundances) > 1:
            second_hit = abundances[1]
            second_hit_status = check("read_hit2_genus_fraction", this_refs, second_hit[1])
            second_hit_status = status["warn"] if second_hit_status == status["fail"] else second_hit_status
            qc_calls[second_hit_status].append(f"read_hit2_genus_fraction_{platform}")
        else:
            qc_calls[status["missing"]].append(f"read_hit2_genus_fraction_{platform}")

    # Assembly taxonomic composition
    genus_stats = bracken_stats(data["assembly"])

    abundances = sorted(genus_stats.items(), key=lambda x: x[1], reverse=True)

    first_hit = abundances[0]
    first_hit_status = check("contig_hit1_genus_fraction", this_refs, first_hit[1])
    qc_calls[first_hit_status].append("contig_genus_species_fraction")
    if first_hit_status == status["fail"]:
        qc_calls["messages"].append("Bracken: Assembly composition abundance of dominant species below threshold - possible contamination issue.")

    # Assembly checkM
    checkm_stats = data["checkm"]
    checkm_status = check("checkm_contamination", this_refs, checkm_stats["Contamination"])
    qc_calls[checkm_status].append("checkm_contamination")

    if checkm_status == status["fail"]:
        qc_calls["messages"].append(("CheckM: Assembly contains predicted contamination above threshold."))
        
    # If there is a second genus detected
    if len(abundances) > 1:
        second_hit = abundances[1]
        second_hit_status = check("contig_hit1_genus_fraction", this_refs, second_hit[1])
        qc_calls[second_hit_status].append(f"contig_hit1_genus_fraction{platform}")
    else:
        qc_calls[status["missing"]].append("contig_hit1_genus_fraction")

    # Taxonkit assembly QC
    taxonkit_genus = sorted(data["taxonkit"]["genus"], key=lambda d: d['fraction'], reverse=True)
    first_genus = taxonkit_genus[0]
    contig_first_genus_fraction = float(first_genus["fraction"])
    taxonkit_genus_status = check("contig_hit1_genus_fraction", this_refs, contig_first_genus_fraction)
    qc_calls[taxonkit_genus_status].append("taxonkit_genus_fraction")

    if taxonkit_genus_status == status["fail"]:
        qc_calls["messages"].append(("Taxonkit: Majority genus fraction below threshold."))

    # quast
    assembly = int(data["quast"]["Total length"])
    assembly_status = check("Total length", this_refs, assembly)
    qc_calls[assembly_status].append("quast_assembly")
    if assembly_status == status["fail"]:
        qc_calls["messages"].append("Quast: Size of the assembly outside of reference range")
    elif assembly_status == status["warn"]:
        qc_calls["messages"].append("Quast: Size of the assembly slightly outside of reference range")

    contigs = int(data["quast"]["# contigs"])
    contigs_status = check("# contigs (>= 0 bp)", this_refs, contigs)
    qc_calls[contigs_status].append("quast_contigs")
    if contigs_status == status["warn"]:
        qc_calls["messages"].append("Quast: Number of contigs slightly outside of reference range")
    elif contigs_status == status["fail"]:
        qc_calls["messages"].append("Quast: Number of contigs well outside of reference range")

    n50_status = check("N50", this_refs, int(data["quast"]["N50"]))
    qc_calls[n50_status].append("quast_n50")
    if n50_status == status["warn"]:
        qc_calls["messages"].append("Quast: N50 of this assembly slightly outside of reference range")
    elif n50_status == status["fail"]:
        qc_calls["messages"].append("Quast: N50 of this assembly well outside of reference range")

    gc_status = check("GC (%)", this_refs, float(data["quast"]["GC (%)"]))
    qc_calls[gc_status].append("quast_gc")
    if gc_status == status["warn"]:
        qc_calls["messages"].append("Quast: GC value of this assembly slightly outside of reference range")
    elif gc_status == status["fail"]:
        qc_calls["messages"].append("Quast: GC value of this assembly well outside of reference range")

    duplication_status = check("Duplication ratio", this_refs, data["quast"]["Duplication ratio"])
    qc_calls[duplication_status].append("quast_duplication")
    if duplication_status == status["warn"]:
        qc_calls["messages"].append("Quast: Duplication ratio of assembly slightly outside of reference range")
    elif duplication_status == status["fail"]:
        qc_calls["messages"].append("Quast: Duplication ratio of this assembly well outside of reference range")

    # busco
    if "busco" in data:
        busco = data["busco"]
        busco_total = int(busco["dataset_total_buscos"])
        busco_completeness = round(((int(busco["C"])) / int(busco_total)), 2)
        busco_completeness_status = check("busco_single", this_refs, busco_completeness)
        qc_calls[busco_completeness_status].append("busco_completeness")
        if busco_completeness_status == status["warn"]:
            qc_calls["messages"].append("Busco: Assembly may be incomplete")
        elif busco_completeness_status == status["fail"]:
            qc_calls["messages"].append("Busco: Assembly likely incomplete")

        busco_duplicates = round((int(busco["D"]) / busco_total), 2)
        busco_duplicates_status = check("busco_duplicates", this_refs, busco_duplicates)
        qc_calls[busco_duplicates_status].append("busco_duplicates")
        if busco_duplicates_status == status["warn"]:
            qc_calls["messages"].append("Busco: Assembly contains some duplications")
        elif busco_duplicates_status == status["fail"]:
            qc_calls["messages"].append("Busco: Assembly contains many duplications")

    # mosdepth total
    for platform, pdata in data["mosdepth"].items():
        coverage = float(pdata["mean"])
        coverage_status = check("assembly_coverageDepth", this_refs, coverage)
        qc_calls[coverage_status].append(f"coverage_{platform}_mean")

        if coverage_status == status["fail"]:
            qc_calls["messages"].append(f"Mosdepth: {platform} coverage below threshold!")
        elif coverage_status == status["warn"]:
            qc_calls["messages"].append(f"Mosdepth: {platform} coverage below threshold!")

    # mosdepth global
    for platform, pdata in data["mosdepth_global"].items():
        coverage_status = status["missing"]
        if "40" in pdata:
            coverage = pdata["40"]
            if coverage < 90:
                coverage_status = status["warn"]
                qc_calls["messages"].append(f"Mosdepth: Less than 90% of assembly coveraged at 40X by {platform} reads - this may be too low")
            elif coverage < 50:
                coverage_status = status["fail"]
                qc_calls["messages"].append(f"Mosdepth: Less than 50% of assembly coveraged at 40X by {platform} reads - this is likely too low")
            else:
                coverage_status = status["pass"]
        qc_calls[coverage_status].append(f"coverage_{platform}_40x")

    data["qc"] = qc_calls

    data["qc"]["call"] = status["missing"]

    fail_categories = ["confindr_illumina", "confindr_nanopore", "quast_contigs", "coverage_total_mean", "read_hit1_species_fraction_ILLUMINA", "read_hit1_species_fraction_NANOPORE", "taxonkit_genus_fraction"]

    # overall qc ruling
    for category in fail_categories:
        if category in qc_calls["fail"]:
            data["qc"]["call"] = status["fail"]

    if data["qc"]["call"] != status["fail"]:
        if len(data["qc"]["warn"]) > 0:
            data["qc"]["call"] = status["warn"]
        else:
            data["qc"]["call"] = status["pass"]

    with open(output, "w") as fo:
        json.dump(data, fo, indent=4, sort_keys=True)


if __name__ == '__main__':
    main(args.input, args.refs, args.output)
