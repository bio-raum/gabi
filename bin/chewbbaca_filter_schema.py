#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import argparse
import pathlib
import os

parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--list", help="List of valid alleles")
parser.add_argument("--schema", help="Folder with chewbbaca schema")
parser.add_argument("--output", help="Name of output folder")

args = parser.parse_args()


def main(list, schema, output):

    pathlib.Path(output).mkdir(parents=True, exist_ok=True)

    alleles = glob.glob(schema + "/*/*.fasta")

    with open(list) as file:
        valid = [line.rstrip() for line in file]

    for allele in alleles:
        allele_name = os.path.basename(allele)
        # There probably is a pythonier way to do this
        if allele_name in valid:
            os.system("cp " + allele + " " + output + "/")


if __name__ == '__main__':
    main(args.list, args.schema, args.output)
