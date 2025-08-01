#!/usr/bin/env python3

import pandas as pd
import sys

args = sys.argv

filename = args[1]
outfile = args[2]

df = pd.read_csv(filename, dtype='unicode')
x = df['Otherinfo1']
discarded_column = df.columns.get_loc('Otherinfo2')

data = dict()
somatic_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'FILTER', 'SOMATIC_FLAG', 'REF_COUNT', 'ALT_COUNT', 'VAF%', 
                'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Gene_full_name.refGene', 
                'Function_description.refGene', 'Disease_description.refGene', 'cosmic84', 'PopFreqMax']

data.setdefault('VAF%', [])
data.setdefault('REF_COUNT', [])
data.setdefault('ALT_COUNT', [])

for row in x:
    rowitems = row.split('\t')
    info = rowitems[10].split(';')
    
    data['REF_COUNT'].append(info[2].split('=')[1])
    data['ALT_COUNT'].append(info[1].split('=')[1])
    data['VAF%'].append(float(info[4].split('=')[1]))  # Convert to float for multiplication

# Multiply VAF% by 100
data['VAF%'] = [value * 100 for value in data['VAF%']]

df1 = df.iloc[:, :5]
df2 = pd.DataFrame(data, columns=data.keys())
df3 = df.iloc[:, 5:discarded_column]

horizontal_stack = pd.concat([df1, df2, df3], axis=1)
horizontal_stack['cosmic84'] = horizontal_stack['cosmic84'].str.replace(',', ';')
horizontal_stack['AAChange.refGene'] = horizontal_stack['AAChange.refGene'].str.replace(',', ';')
horizontal_stack.replace(to_replace='.', value='-1', inplace=True)
horizontal_stack=horizontal_stack.reindex(columns = somatic_cols)
horizontal_stack.to_csv(outfile, index=False)
