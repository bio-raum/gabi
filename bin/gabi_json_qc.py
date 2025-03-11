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
                                m = f"Contamination ({contam_type}) detected in {platform} reads {read['Sample']}"
                                if m not in qc_calls["messages"]:
                                    qc_calls["messages"].append(m)
                            else:
                                m = f"Low levels of contamination ({contam_type}) detected in {platform} reads {read['Sample']}"
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
            qc_calls["messages"].append("Illumina Q30 fraction below 75%")

    # bracken
    for platform, bracken in data["bracken"].items():
        taxon_count = 0
        taxon_count_same = 0
        taxon_count_status = status["pass"]

        for tax in bracken:
            this_taxon = tax["name"].replace('"', '')
            genus = this_taxon.split(" ")[0]
            # The Bracken results are all in quotes, so we need to clean that up and convert to precentage
            tperc = round((float(tax["fraction_total_reads"].replace('"', '')) * 100), 2)

            if (tperc > 5.0):
                taxon_count += 1
                if genus in taxon:
                    taxon_count_same += 1
        # Some genera will yield multiple taxa hits, so we only trigger this if the signals come from different genera
        if ((taxon_count - taxon_count_same) > 3):
            taxon_count_status = status["fail"]
            qc_calls["messages"].append(f"More than three taxa detected in {platform} read data!")
        elif ((taxon_count - taxon_count_same) > 1):
            taxon_count_status = status["warn"]
            qc_calls["messages"].append(f"More than one taxon detected in the {platform} read data!")

        qc_calls[taxon_count_status].append(f"taxon_{platform.lower()}_count")

    # quast
    assembly = int(data["quast"]["Total length"])
    assembly_status = check("Total length", this_refs, assembly)
    qc_calls[assembly_status].append("quast_assembly")
    if assembly_status == status["fail"]:
        qc_calls["messages"].append("Size of the assembly outside of reference range")
    elif assembly_status == status["warn"]:
        qc_calls["messages"].append("Size of the assembly slightly outside of reference range")

    contigs = int(data["quast"]["# contigs"])
    contigs_status = check("# contigs (>= 0 bp)", this_refs, contigs)
    qc_calls[contigs_status].append("quast_contigs")
    if contigs_status == status["warn"]:
        qc_calls["messages"].append("Number of contigs slightly outside of reference range")
    elif contigs_status == status["fail"]:
        qc_calls["messages"].append("Number of contigs well outside of reference range")

    n50_status = check("N50", this_refs, int(data["quast"]["N50"]))
    qc_calls[n50_status].append("quast_n50")
    if n50_status == status["warn"]:
        qc_calls["messages"].append("N50 of this assembly slightly outside of reference range")
    elif n50_status == status["fail"]:
        qc_calls["messages"].append("N50 of this assembly well outside of reference range")

    gc_status = check("GC (%)", this_refs, float(data["quast"]["GC (%)"]))
    qc_calls[gc_status].append("quast_gc")
    if gc_status == status["warn"]:
        qc_calls["messages"].append("GC value of this assembly slightly outside of reference range")
    elif gc_status == status["fail"]:
        qc_calls["messages"].append("GC value of this assembly well outside of reference range")

    duplication_status = check("Duplication ratio", this_refs, data["quast"]["Duplication ratio"])
    qc_calls[duplication_status].append("quast_duplication")
    if duplication_status == status["warn"]:
        qc_calls["messages"].append("Duplication ratio of assembly slightly outside of reference range")
    elif duplication_status == status["fail"]:
        qc_calls["messages"].append("Duplication ratio of this assembly well outside of reference range")

    # busco
    if "busco" in data:
        busco = data["busco"]
        busco_total = int(busco["dataset_total_buscos"])
        busco_completeness = round(((int(busco["C"])) / int(busco_total)), 2)
        busco_completeness_status = check("busco_single", this_refs, busco_completeness)
        qc_calls[busco_completeness_status].append("busco_completeness")
        if busco_completeness_status == status["warn"]:
            qc_calls["messages"].append("Assembly may be incomplete")
        elif busco_completeness_status == status["fail"]:
            qc_calls["messages"].append("Assembly likely incomplete")

        busco_duplicates = round((int(busco["D"]) / busco_total), 2)
        busco_duplicates_status = check("busco_duplicates", this_refs, busco_duplicates)
        qc_calls[busco_duplicates_status].append("busco_duplicates")
        if busco_duplicates_status == status["warn"]:
            qc_calls["messages"].append("Assembly contains some duplications")
        elif busco_duplicates_status == status["fail"]:
            qc_calls["messages"].append("Assembly contains many duplications")

    # mosdepth total
    for platform, pdata in data["mosdepth"].items():
        coverage = float(pdata["mean"])
        coverage_status = check("assembly_coverageDepth", this_refs, coverage)
        qc_calls[coverage_status].append(f"coverage_{platform}_mean")

        if coverage_status == status["fail"]:
            qc_calls["messages"].append(f"{platform} coverage below threshold!")
        elif coverage_status == status["warn"]:
            qc_calls["messages"].append(f"{platform} coverage below threshold!")

    # mosdepth global
    for platform, pdata in data["mosdepth_global"].items():
        coverage_status = status["missing"]
        if "40" in pdata:
            coverage = pdata["40"]
            if coverage < 90:
                coverage_status = status["warn"]
                qc_calls["messages"].append(f"Less than 90% of assembly coveraged at 40X by {platform} reads - this may be too low")
            elif coverage < 50:
                coverage_status = status["fail"]
                qc_calls["messages"].append(f"Less than 50% of assembly coveraged at 40X by {platform} reads - this is likely too low")
            else:
                coverage_status = status["pass"]
        qc_calls[coverage_status].append(f"coverage_{platform}_40x")

    data["qc"] = qc_calls

    data["qc"]["call"] = status["missing"]

    fail_categories = ["confindr_illumina", "confindr_nanopore", "quast_contigs", "coverage_total_mean"]

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
