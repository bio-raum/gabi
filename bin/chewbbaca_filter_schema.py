#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import argparse


parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--list", help="List of valid alleles")
parser.add_argument("--schema", help="Folder with chewbbaca schema")
parser.add_argument("--output", help="Name of output folder")

args = parser.parse_args()


def main(list, schema, output):

    alleles = glob.glob("*/*.fasta")
    
    with open(list) as file:
        valid = [line.rstrip() for line in file]
 
    for keep in valid:
        


if __name__ == '__main__':
    main(args.list, args.schema, args.output)
