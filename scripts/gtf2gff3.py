#!/usr/bin/env python

import argparse
import gffutils
import sys
import os

def gtf_to_gff3(gtf_file, gff3_file):
    """
    Converts a GTF file to GFF3 format without losing information.

    Args:
        gtf_file (str): Path to the input GTF file.
        gff3_file (str): Path to the output GFF3 file.
    """
    # Create a temporary database
    db_fn = os.path.splitext(gtf_file)[0] + '.db'
    if os.path.exists(db_fn):
        os.remove(db_fn)

    print(f"Creating a database from {gtf_file}...")
    try:
        db = gffutils.create_db(
            gtf_file,
            dbfn=db_fn,
            force=True,
            keep_order=True,
            disable_infer_genes=False,
            disable_infer_transcripts=False,
            merge_strategy='merge',
            sort_attribute_values=True
        )
    except Exception as e:
        print(f"Error creating database: {e}")
        sys.exit(1)

    print(f"Writing GFF3 to {gff3_file}...")
    try:
        with open(gff3_file, 'w') as out_handle:
            for feature in db.all_features(order_by=('seqid', 'start', 'end', 'strand')):
                out_handle.write(str(feature) + '\n')
    except Exception as e:
        print(f"Error writing GFF3 file: {e}")
        sys.exit(1)

    # Remove temporary database file
    if os.path.exists(db_fn):
        os.remove(db_fn)
    if os.path.exists(db_fn + '-shm'):
        os.remove(db_fn + '-shm')
    if os.path.exists(db_fn + '-wal'):
        os.remove(db_fn + '-wal')

    print("Conversion completed successfully.")

def main():
    parser = argparse.ArgumentParser(description="Convert GTF to GFF3 without losing information.")
    parser.add_argument('-i', '--input', required=True, help='Input GTF file')
    parser.add_argument('-o', '--output', required=True, help='Output GFF3 file')
    args = parser.parse_args()

    gtf_to_gff3(args.input, args.output)

if __name__ == '__main__':
    main()