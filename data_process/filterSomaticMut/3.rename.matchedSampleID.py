"""
  Hua Sun
  4/29/2019

  Rename normal sample of meta3 based on paired table

  python3 rename.matchedSampleID.py --table pairedName.table --meta3 input.tsv -o output
  
  --table <file>
  --meta3 <file>
  -o|--out <file>

"""

import argparse
import pandas as pd
import re
import os


parser = argparse.ArgumentParser()

parser.add_argument('-t', '--table', type=str, required=True, help='pair samples')
parser.add_argument('-m', '--meta3', type=str, required=True, help='meta3 form')
parser.add_argument('-o', '--out', type=str, required=True, help='outfile')

args = parser.parse_args()



df_pairSample = pd.read_csv(args.table, sep = "\t", low_memory=False)

df_pairSample = df_pairSample[['Normal', 'Tumor']]
df_pairSample.columns = ['Matched_Norm_Sample_Barcode', 'Tumor_Sample_Barcode']

df_meta3 = pd.read_csv(args.meta3, sep = "\t", low_memory=False)
df_meta3['Matched_Norm_Sample_Barcode'] = df_meta3['Tumor_Sample_Barcode'].map(df_pairSample.set_index('Tumor_Sample_Barcode')['Matched_Norm_Sample_Barcode'])

df_meta3.to_csv(args.out, sep = "\t", index = False, header = True)


