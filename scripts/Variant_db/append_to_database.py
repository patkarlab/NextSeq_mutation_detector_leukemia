#!/usr/bin/python3
import sqlite3
import sys
import os

def append_to_database(db_name, vcf_files):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    for vcf_file in vcf_files:
        sample_id = os.path.basename(vcf_file).replace('.somaticseq.vcf', '')

        with open(vcf_file, 'r') as vcf_file:
            for line in vcf_file:
                if not line.startswith('#'):
                    data = line.strip().split('\t')
                    chrom, pos, ref, alt = data[0], int(data[1]), data[3], data[4]
                    gt_field = data[9].split(":")[0]
                    alt_num = gt_field.split('/').count('1')

                    cursor.execute('SELECT SAMPLE_NAME, ALT_NUM FROM Variants WHERE CHROM=? AND POS=? AND REF=? AND ALT=?;', (chrom, pos, ref, alt))
                    existing_data = cursor.fetchone()

                    if existing_data:
                        existing_sample_ids, existing_alt_nums = existing_data

                        if sample_id not in existing_sample_ids:
                            # Appending new sample ID to the existing ones and counting the number of sample IDs
                            new_sample_ids = f"{existing_sample_ids},{sample_id}"
                            new_alt_nums = existing_alt_nums + alt_num

                            cursor.execute('UPDATE Variants SET SAMPLE_NAME = ?, ALT_NUM = ? WHERE CHROM=? AND POS=? AND REF=? AND ALT=?;', (new_sample_ids, new_alt_nums, chrom, pos, ref, alt))
                    else:
                        cursor.execute('INSERT INTO Variants (CHROM, POS, REF, ALT, SAMPLE_NAME, ALT_NUM) VALUES (?, ?, ?, ?, ?, ?);', (chrom, pos, ref, alt, sample_id, alt_num))

#    cursor.execute('SELECT CHROM, POS, REF, ALT, SAMPLE_NAME, ALT_NUM FROM Variants;')
#   variants_data = cursor.fetchall()

#    for variant in variants_data:
#        chrom, pos, ref, alt, sample_id, alt_num = variant

        # Calculating frequency
#        num_samples = len(sample_id.split(','))  # Calculate num_samples for each variant
#        frequency = alt_num / (2 * num_samples)

#        print(f"Chrom: {chrom}, Pos: {pos}, Ref: {ref}, Alt: {alt}, Frequency: {frequency}")

    conn.commit()
    conn.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: {} <database_name> <vcf_file1> <vcf_file2> ...".format(sys.argv[0]))
        sys.exit(1)

    db_name = sys.argv[1]
    vcf_files = sys.argv[2:]

    append_to_database(db_name, vcf_files)
    print("Variants from {} VCF files appended to the 'Variants' table in '{}'".format(len(vcf_files), db_name))
