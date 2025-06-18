#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--output", "-o")
parser.add_argument("--input", "-i")

args = parser.parse_args()

def make_id(count):

    base = "contig"
    if count < 10:
        suffix = "0000"
    elif count < 100:
        suffix = "000"
    elif count < 1000:
        suffix = "00"
    else:
        suffix = "0"
    
    return f"{base}{suffix}{count}"

def main(input, output):

    data = []

    records = SeqIO.index(input, "fasta")

    # Iterate over each dict record and build a list of hashes to sort on
    for k, record in records.items():
        entries = record.description.split(" ")
        entries.pop(0)
        meta = {}
        for entry in entries:
            key,value = entry.split("=")
            meta[key] = value
            
        data.append( { "id": record.id, "name": record.name , "length": int(meta["len"]), "coverage": float(meta["cov"]), "description": " ".join(entries) } )
    
    # sort the entries by both sequence length and coverage (which is what Shovill gets wrong!)
    sorted_keys = sorted(data, key = lambda x: (x["length"], x["coverage"]), reverse=True )

    sequences = []
    counter = 0
    
    # use the properly sorted contig names and build new SeqRecords
    for key in sorted_keys:
        counter += 1
        seq = records[key["id"]].seq
        this_id = make_id(counter)
        r = SeqRecord(
            seq,
            id=this_id,
            description=key["description"]
        )
        sequences.append(r)

    SeqIO.write(sequences, output, "fasta")


if __name__ == '__main__':
    main(args.input, args.output)
