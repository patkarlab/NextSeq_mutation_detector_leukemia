#!/usr/bin/python3

import sqlite3
import sys

def calculate_frequency_for_variants(db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # Unique sample count
    cursor.execute('SELECT SAMPLE_NAME FROM Variants;')
    rows = cursor.fetchall()
    unique_samples = set()
    for row in rows:
        sample_names = row[0].split(',')
        unique_samples.update(sample_names)
    total_unique_samples = len(unique_samples)
 
    cursor.execute('SELECT CHROM, POS, REF, ALT, SAMPLE_NAME FROM Variants;')
    variants = cursor.fetchall()

    for variant in variants:
        chrom, pos, ref, alt, sample_names = variant
        samples_for_variant = set(sample_names.split(','))
        
        # Calculate frequency
        frequency = len(samples_for_variant) / total_unique_samples

        print(f"Variant: Chrom:{chrom}, Pos:{pos}, Ref:{ref}, Alt:{alt}")
        print(f"Number of Samples: {len(samples_for_variant)}")
        print(f"Frequency: {frequency:.4f}\n")

    conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: {} <database_name>".format(sys.argv[0]))
        sys.exit(1)

    db_name = sys.argv[1]
    calculate_frequency_for_variants(db_name)
