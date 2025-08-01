#!/usr/bin/env python3
import pandas as pd
import sys

args = sys.argv

input_csv_file = args[1]
filtered_csv_file = args[2]

df = pd.read_csv(input_csv_file)

# Check if the 'reference_pos' column exists
if 'reference_pos' in df.columns:
    filtered_df = df[df['reference_pos'].notna()]  # Filter rows where 'reference_pos' is not NaN

    filtered_df.to_csv(filtered_csv_file, index=False)
    print(f"Filtered details saved to: {filtered_csv_file}")
else:
    print("The input CSV file does not have a 'reference_pos' column.")

