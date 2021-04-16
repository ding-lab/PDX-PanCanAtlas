"""

    Hua Sun
    2021-03-15; 1/25/2020; 4/28/2020

    Filter unknown mutations from PDX mutations
        1. keep human tumor somatic mutations from same case
        2. keep TCGA & COSMIC somatic mutations

    USAGE
    =======
    # filter PDX sample mut based on ref loci & human tumor from matched case
    python3 run.py -m meta3.tsv --info sample.info -r ref_db.loci -o outdir
    or
    python3 run.py -m meta3.tsv --info sample.info -r ref_db.loci -o outdir --outname filtered_pdxMut.meta3.tsv
  

"""

import argparse
import pandas as pd
import glob
import re
import os
import io


# collect input arguments
parser = argparse.ArgumentParser()

parser.add_argument('--info', type=str, required=True, help='sample info')
parser.add_argument('-m', '--meta3', type=str, required=True, help='meta3 form')
parser.add_argument('-r', '--ref', type=str, default='', help='ref loci')
parser.add_argument('-o', '--outdir', type=str, default='./filteredPDXmut', help='outdir')
parser.add_argument('--outname', type=str, default='filtered_PDX_somaticMut.meta3.tsv', help='out file name')

args = parser.parse_args()



'''
    Main
'''

def main():
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    if args.ref == '':
        print('[Error] Please enter reference loci data ...')
        quit()

    FilterPDX_somaticMut(args.info, args.ref, args.meta3, args.outdir)

    f_filteredMut = f'{args.outdir}/{args.outname}'
    MergeAllFilesToOne(f'{args.outdir}/temp', f_filteredMut)

    os.system(f'rm -rf {args.outdir}/temp')

    outlog = f'{args.outdir}/for_check.removed_PDX_somaticMut.meta3.log'
    ExtractRemovedMut(args.meta3, f_filteredMut, outlog)



'''
    Set Func
'''

def FilterPDX_somaticMut(f_info, f_ref_loci, f_meta3, outdir):
    info = pd.read_csv(args.info, sep='\t', usecols=range(0,3))
    info.columns = ['CaseID', 'ID', 'Group']
    info = info.loc[info['Group'].isin(['Human_Tumor','PDX'])]

    ref_loci = pd.read_csv(f_ref_loci, sep='\t', header=None, low_memory=False)
    ref_loci.columns = ['Chr', 'Pos', 'Ref']

    meta3 = pd.read_csv(f_meta3, sep='\t', low_memory=False)


    temp_dir = f'{outdir}/temp'
    if not os.path.isdir(f'{outdir}/temp'):
        os.mkdir(temp_dir)

    for case_id in info['CaseID'].unique():
        case_info = info.loc[info['CaseID']==case_id]
        
        print(f'[INFO] Filter PDX mutations for patient {case_id} ... ')
        FilterMutationPerCase(case_id, case_info, ref_loci, meta3, temp_dir)




''' Filter PDX mut per patient '''
def FilterMutationPerCase(caseName, case_info, ref_loci, meta3, outdir):
    
    target_loci = ref_loci

    if 'Human_Tumor' in case_info['Group'].unique():
        ht_id = case_info['ID'].loc[case_info['Group']=='Human_Tumor'].tolist()
        ht_mut = meta3.loc[meta3['Tumor_Sample_Barcode'].isin(ht_id)]
        
        ht_mut.to_csv(f'{outdir}/{caseName}.ht_meta3.tsv', sep='\t', index=False)

        ht_loci = ht_mut[['Chromosome', 'Start_position', 'Reference_Allele']].drop_duplicates()
        ht_loci.columns = ['Chr', 'Pos', 'Ref']

        target_loci = target_loci.append(ht_loci)
        target_loci.drop_duplicates(inplace=True)
    
    target_loci['tag'] = target_loci['Chr'].astype(str) + '_' + target_loci['Pos'].astype(str) + '_' + target_loci['Ref'].astype(str)

    pdx_id = case_info['ID'].loc[case_info['Group']=='PDX'].tolist()
    pdx_mut = meta3.loc[meta3['Tumor_Sample_Barcode'].isin(pdx_id)]
    pdx_mut['tag'] = pdx_mut['Chromosome'].astype(str) + '_' + pdx_mut['Start_position'].astype(str) + '_' + pdx_mut['Reference_Allele'].astype(str)
    
    pdx_mut_filtered = pdx_mut.loc[pdx_mut['tag'].isin(target_loci['tag'])]
    pdx_mut_filtered = pdx_mut_filtered.drop(['tag'], axis=1)
    pdx_mut_filtered.to_csv(f'{outdir}/{caseName}.pdx_filtered_meta3.tsv', sep='\t', index=False)




''' Merge all files to one '''
def MergeAllFilesToOne(file_path, outfile):
    pattern = f'{file_path}/*.tsv'
    file_list = glob.glob(pattern)
    
    df = pd.read_csv(file_list[0], sep='\t', low_memory=False)

    for file in file_list[1:len(file_list)]:
        df_csv = pd.read_csv(file, sep='\t', low_memory=False)
        df = df.append(df_csv)
    
    df.to_csv(outfile, sep='\t', index=False)

    status = df.groupby(['Tumor_Sample_Barcode']).size().reset_index(name='Count').sort_values(['Count'])
    status.to_csv(f'{outfile}.status', sep='\t', index=False)




''' Extract removed mutation '''
def ExtractRemovedMut(f_orgMut, f_filteredMut, outfile):
    df = pd.read_csv(f_orgMut, sep='\t', low_memory=False)
    filteredMut = pd.read_csv(f_filteredMut, sep='\t', low_memory=False)
    df = df.append(filteredMut)

    df_log = df[~df.duplicated(['Tumor_Sample_Barcode', 'Chromosome', 'Start_position', 'Reference_Allele', 'Tumor_Seq_Allele2'], keep=False)]

    df_log.to_csv(outfile, sep='\t', index=False)






if __name__ == '__main__':
    main()

