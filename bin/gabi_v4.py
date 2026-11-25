#!/usr/bin/env python
import shutil
import subprocess
import tempfile
from pathlib import Path
from plotly.offline import get_plotlyjs
import argparse
import json
import os
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from concurrent.futures import ThreadPoolExecutor
import re
import plotly.express as px

def read_json_files(file_paths):
    """Read multiple JSON files and return a list of data dictionaries."""
    data = []
    for file_path in file_paths:
        with open(file_path, 'r') as f:
            data.append(json.load(f))
    return data

def get_nested_value(data, keys, default="missing"):
    """Safely get nested dictionary values."""
    for key in keys:
        if isinstance(data, dict) and key in data:
            data = data[key]
        else:
            return default
    return data

def check_contaminant_status(jdata):
    """Check contaminant status for both Illumina and Nanopore."""
    return {
        "illumina": jdata.get("contamination_illumina", "missing"),
        "nanopore": jdata.get("contamination_nanopore", "missing")
    }

def check_fastp_q30(jdata):
    """Check fastp Q30 values."""
    fastp = jdata.get("fastp", {})
    q30 = fastp.get("q30", "missing")
    status = "pass" if q30 >= 80 else "warn" if q30 != "missing" else "missing"
    return q30, status

def check_nanostat(jdata):
    """Check nanostat values."""
    nanostat = jdata.get("nanostat", {})
    return nanostat.get("q15", "missing"), nanostat.get("n50", "missing")

def check_taxon_count(jdata):
    """Check taxon counts for different platforms."""
    return {
        "ILLUMINA": jdata.get("taxon_count_illumina", "missing"),
        "NANOPORE": jdata.get("taxon_count_nanopore", "missing"),
        "PACBIO": jdata.get("taxon_count_pacbio", "missing")
    }

def check_bracken(jdata, bracken_data_all, samples_by_technology):
    """Process Bracken data."""
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
    """Check CheckM status."""
    return jdata.get("checkm", "missing")

def check_taxonkit_genus(jdata):
    """Check taxonkit genus."""
    genus = jdata.get("taxonkit_genus", "missing")
    return genus, "missing"

def check_samtools_and_quast(jdata):
    """Process samtools and quast data."""
    if "samtools" in jdata and "quast" in jdata:
        samtools = {
            "q30": jdata["samtools"].get("q30", "missing"),
            "n50": jdata["quast"].get("n50", "missing")
        }
        quast = {
            "completeness": jdata["quast"].get("completeness", "missing"),
            "duplicated": jdata["quast"].get("duplicated", "missing"),
            "missing": jdata["quast"].get("missing", "missing"),
            "fragmented": jdata["quast"].get("fragmented", "missing")
        }
        checkm = jdata.get("checkm", "missing")
        insert_size = jdata["samtools"].get("insert_size", [])
        return samtools, quast, checkm, insert_size
    
    return "missing", "missing", "missing", []

def check_coverage(jdata):
    """Check coverage data."""
    coverage_data = jdata.get("coverage", {})
    coverage_40 = coverage_data.get("40x", "missing")
    
    if isinstance(coverage_40, dict):
        return (
            "missing",
            coverage_40,
            coverage_40.get("illumina", "missing"),
            coverage_40.get("nanopore", "missing"),
            coverage_40.get("pacbio", "missing")
        )
    
    return "missing", coverage_40, "missing", "missing", "missing"

def check_assembly(jdata):
    """Check assembly data."""
    assembly = jdata.get("assembly", {})
    return (
        assembly.get("n50", "missing"),
        assembly.get("fraction", "missing"),
        assembly.get("contigs", "missing"),
        assembly.get("plasmids", "missing"),
        assembly.get("size", "missing"),
        jdata.get("assembly_status", "missing")
    )

def check_serotype(jdata):
    """Check serotype data."""
    return jdata.get("serotype", "missing")

def check_busco(jdata, busco_data_all):
    """Process BUSCO data."""
    if "busco" in jdata:
        busco_data_all.append({
            "sample": jdata["sample"],
            "completeness": jdata["busco"].get("completeness", "missing"),
            "missing": jdata["busco"].get("missing", "missing"),
            "fragmented": jdata["busco"].get("fragmented", "missing"),
            "duplicated": jdata["busco"].get("duplicated", "missing")
        })
    return busco_data_all

def check_amrfinder(jdata):
    """Process AMRFinder data."""
    amrfinder_data = []
    amr_counts = 0
    
    if "amrfinder" in jdata:
        for amr_entry in jdata["amrfinder"]:
            amrfinder_data.append({
                "gene_symbol": amr_entry.get("Gene symbol", ""),
                "amr_class": amr_entry.get("Class", ""),
                "amr_subclass": amr_entry.get("Subclass", ""),
                "element_type": amr_entry.get("Element type", ""),
                "element_subtype": amr_entry.get("Element subtype", ""),
                "sequence_name": amr_entry.get("Sequence name", "")
            })
        amr_counts = len(set(entry["Class"] for entry in amrfinder_data if entry["Class"]))
    
    return amrfinder_data, amr_counts

def check_mlst(jdata):
    """Check MLST data."""
    mlst = jdata.get("mlst", [])
    return mlst[0] if mlst else "missing"

def check_status(qc):
    """Check QC status."""
    return qc.get("status", "missing")

def plot_bracken(bracken_data_all, samples_by_technology):
    """Create Bracken plots."""
    plots = {}
    for platform, data in bracken_data_all.items():
        if data:  # Only create plot if there's data
            df = pd.DataFrame(data, index=samples_by_technology[platform])
            plot = px.bar(df, orientation='h', labels={"index": "Samples", "value": "Percentage"}, height=(len(samples_by_technology[platform]) * 25))
            plots[platform] = plot.to_html(full_html=False, include_plotlyjs=False)
    return plots

def plot_insert_sizes(insert_sizes_all):
    """Create insert size plot."""
    if insert_sizes_all:
        df = pd.DataFrame(insert_sizes_all)
        plot = px.line(df, labels={"index": "Basepairs", "value": "Count"})
        return plot.to_html(full_html=False, include_plotlyjs=False)
    return ""

def plot_busco(busco_data_all):
    """Create BUSCO plot."""
    if busco_data_all:
        df = pd.DataFrame(busco_data_all, index=[item["sample"] for item in busco_data_all])
        plot = px.bar(df, orientation='h', labels={"index": "Samples", "value": "Percentage"}, height=(len(busco_data_all) * 25))
        return plot.to_html(full_html=False, include_plotlyjs=False)
    return ""

def parse_yaml(yaml_path):
    """Parse YAML file to extract software information."""
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

def compile_tailwind_css(html_path: Path, tailwind_bin: str = "tailwindcss") -> str:
    """Compile only the Tailwind CSS utility classes present in html_path."""
    if not shutil.which(tailwind_bin) and not Path(tailwind_bin).exists():
        raise FileNotFoundError(
            f"Tailwind CLI binary '{tailwind_bin}' not found on PATH. "
            "Install via `conda install -c conda-forge tailwindcss`."
        )

    tmp_dir = html_path.parent
    input_css = tmp_dir / "input.css"
    output_css = tmp_dir / "output.css"

    input_css.write_text(
        f'@import "tailwindcss";\n@source "{html_path.resolve()}";\n',
        encoding="utf-8",
    )

    subprocess.run(
        [tailwind_bin, "-i", str(input_css), "-o", str(output_css), "--minify"],
        check=True,
        capture_output=True,
        text=True,
    )
    return output_css.read_text(encoding="utf-8")


def render_template(data, template_path, output_path, alpine_path=None, tailwind_bin="tailwindcss"):
    """Render HTML template with data, compiled Tailwind CSS, and inlined JS assets."""
    template_file = Path(template_path)
    env = Environment(loader=FileSystemLoader(str(template_file.parent)))
    j2_template = env.get_template(template_file.name)

    alpine_js = Path(alpine_path).read_text(encoding="utf-8") if alpine_path and Path(alpine_path).exists() else ""
    base_context = {
        **data,
        "alpine_js": alpine_js,
        "plotly_js": get_plotlyjs(),
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        raw_html = Path(tmpdir) / "raw_report.html"
        raw_html.write_text(j2_template.render(**base_context, compiled_css=""), encoding="utf-8")
        compiled_css = compile_tailwind_css(raw_html, tailwind_bin=tailwind_bin)

    final_output = j2_template.render(**base_context, compiled_css=compiled_css)
    Path(output_path).write_text(final_output, encoding="utf-8")


def extract_data(data_list):
    """Extract and process all data from JSON files."""
    samples = []
    samples_by_technology = {"NANOPORE": [], "ILLUMINA": [], "PACBIO": []}
    bracken_data_all = {}
    insert_sizes_all = {}
    busco_data_all = []
    
    processed_data = []

    for jdata in data_list:
        sample = jdata["sample"]
        samples.append(sample)
        
        # Extract all data using helper functions
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
        n50, genome_fraction, contigs, plasmid_count, assembly, assembly_status = check_assembly(jdata)
        serotype_data = check_serotype(jdata)
        busco_data_all = check_busco(jdata, busco_data_all)
        amrfinder_data, amr_counts = check_amrfinder(jdata)
        mlst_data = check_mlst(jdata)
        
        # Add processed data to main list
        processed_data.append({
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
    
    return processed_data, bracken_data_all, samples_by_technology, insert_sizes_all

def main(yaml_path, template_path, output_path):
    """Main function to process JSON files and generate HTML dashboard."""
    file_paths = [f for f in os.listdir(".") if f.endswith(".json")]
    
    # Read all JSON files
    with ThreadPoolExecutor() as executor:
        data = list(executor.map(read_json_files, [file_paths]))
    
    # Extract data from JSON files
    data, bracken_data_all, samples_by_technology, insert_sizes_all = extract_data(data[0])
    
    # Create plots
    plots = {
        **plot_bracken(bracken_data_all, samples_by_technology),
        **{"Insertsizes": plot_insert_sizes(insert_sizes_all)},
        **{"Busco": plot_busco(busco_data_all)}
    }
    
    data[0].update(plots)
    software = parse_yaml(yaml_path)
    data[0]["packages"] = software
    render_template(data[0], template_path, output_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process JSON files and generate an HTML dashboard.')
    parser.add_argument('--input', help='Path to the YAML file', required=True)
    parser.add_argument('--template', help='Path to the Jinja2 template', required=True)
    parser.add_argument('--output', help='Path to the output HTML file', required=True)
    args = parser.parse_args()
    main(args.input, args.template, args.output)
