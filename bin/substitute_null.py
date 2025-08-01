#!/usr/bin/env python3

import sys
import pandas as pd
import os

input1 = sys.argv[1]
outputFile = sys.argv[2]

if os.path.getsize(input1) != 0:
	df1 = pd.read_csv(input1, sep='\t')
	
	if not df1.empty:
		df1.fillna(value= -1 , inplace=True)
		df1.to_csv(outputFile, sep = '\t', header=True, index=False)
	else :
		pass
else:
	with open(outputFile, 'w') as file:
		pass
