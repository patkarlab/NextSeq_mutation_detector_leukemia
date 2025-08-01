#!/usr/bin/env python3

import sys
import pandas as pd

args = sys.argv

positions_file = args[1]
details_file = args[2]
out_df = args[3]

positions_df = pd.read_csv(positions_file)
details_df = pd.read_csv(details_file)

# Extract position and remove the right shift
details_df['Start'] = details_df.iloc[:, 14].str.split(':').str[1].astype(int) - 1

# Rename columns in details_df to match positions_df
details_df = details_df.rename(columns={
    'wt_coverage': 'REF_COUNT',
    'occurrence': 'ALT_COUNT'
})

# Merge the two DataFrames based on the three columns
merged_df = pd.merge(
    positions_df,
    details_df,
    on=['Start', 'REF_COUNT', 'ALT_COUNT'],
    how='inner'
)

# Drop unwanted columns
merged_df = merged_df.drop(columns=['VAF', 'reference_pos', 'FILTER', 'SOMATIC_FLAG'])

# Remove rows with '-' in the Alt column
merged_df = merged_df[merged_df['Alt'] != '-']

# Keep rows where 'Size' is greater than or equal to 3
merged_df = merged_df[merged_df['size'] >= 3]

# Save the merged output to a CSV file
merged_df.to_csv(out_df, index=False)

print(f"Merged file saved as: {out_df}")

