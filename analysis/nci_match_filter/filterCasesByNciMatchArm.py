"""
    
    Hua Sun

    8/22/2020
    
    Filter false positive NCI-MATCH trial cases via exclusion_case table
"""


import os
import argparse
import pandas as pd


# collect input arguments
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='alteration input file')
parser.add_argument('--tarArm', type=str, required=True, help='1.alt.annoArm.table')
parser.add_argument('--remCase', type=str, required=True, help='2.exclusionCase.table')
parser.add_argument('--info', type=str, required=True, help='3.sampleInfo')
parser.add_argument('--type', type=str, required=True, help='set extract item')
parser.add_argument('-o', '--outdir', type=str, default='.', help='outfile')

args = parser.parse_args()



'''
    Main
'''

def main():
    # read
    f_in = args.input
    f_tarArm = args.tarArm
    f_remCase = args.remCase
    f_info = args.info
    item = args.type
    outdir = args.outdir

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # read data
    df_arm = pd.read_csv(f_tarArm, sep="\t")
    df_remCase = pd.read_csv(f_remCase, sep="\t")
    df_remCase['flag'] = df_remCase['StudyArm'] + '_' + df_remCase['Formal_CaseID']
    df_info = pd.read_csv(f_info, sep="\t")
    df_info = df_info[['Analysis_ID', 'Formal_CaseID', 'ModelID_proper_v2', 'Group', 'CancerType']]
    df_info = df_info.rename(columns={'ModelID_proper_v2':'ModelID_proper'})
    # 'Analysis_ID', 'Formal_CaseID', 'ModelID_proper', 'Group', 'CancerType'

    # only PDX info
    df_pdx_info = df_info[df_info['Group']=='PDX']
    
    
    if (item == 'meta3'):
        ExtractTargetMut(f_in, df_arm, df_pdx_info, df_remCase, outdir)

    if (item == 'cnv_amp'):
        ExtractTargetCopyNumber('amp', f_in, df_arm, df_pdx_info, df_remCase, outdir)
    
    if (item == 'cnv_del'):
        ExtractTargetCopyNumber('del', f_in, df_arm, df_pdx_info, df_remCase, outdir)

    if (item == 'fusion'):
        FilterFusion(f_in, df_arm, df_pdx_info, df_remCase, outdir)
    




'''
    Mutation
'''

# filter PDXmodels from Mutations 
def ExtractTargetMut(f_in, df_arm, df_info, df_remCase, outdir):
    df_meta3 = pd.read_csv(f_in, sep="\t", low_memory=False)
    df_meta3 = df_meta3[['Tumor_Sample_Barcode', 'Hugo_Symbol', 'HGVSp_Short']]
    df_meta3.drop_duplicates(inplace=True)
    df_meta3.rename(columns={'Tumor_Sample_Barcode':'Sample', 'Hugo_Symbol':'Gene', 'HGVSp_Short':'aaMut'}, inplace=True)

    ##------------- Part1
    # 1. extract mutation arm
    df_arm_tar = df_arm[df_arm['Event'] == 'Mutation']
    targetGeneList = list(df_arm_tar['Gene'].unique())
    df_tar_mut = df_meta3[df_meta3['Gene'].isin(targetGeneList)]

    # 2-1. format form for mut form
    df_tar_mut = df_tar_mut.merge(df_arm_tar[['Gene', 'StudyArm']], on='Gene', how='left')
    # SampleID  Gene  aaMut  StudyArm

    # 2-2. ---> specific mut. for same gene have multiple arm BRAF-V600 & non-V600
    # the BRAF-non-V600 --> 'R'
    df_tar_mut.loc[(df_tar_mut['Gene']=='BRAF') & (~df_tar_mut['aaMut'].str.contains('p.V600')), 'StudyArm'] = 'Z1L'


    # 2-3. add case/model/group/cancertype
    df_tar_mut = df_tar_mut.merge(df_info, left_on='Sample', right_on='Analysis_ID')
    # Sample  Gene  aaMut  StudyArm Analysis_ID  Formal_CaseID  ModelID_proper  Group  CancerType

    
    # 3. filter
    df_tar_mut['flag'] = df_tar_mut['StudyArm'] + '_' + df_tar_mut['Formal_CaseID']
    df_filtered_mut = df_tar_mut[~df_tar_mut['flag'].isin(df_remCase['flag'])]
    
    # 4. outupt
    df_tar_mut.to_csv(outdir + '/1.somaticMut.with_arm.out', sep='\t', index=False)
    df_filtered_mut.to_csv(outdir + '/1.somaticMut.with_arm.filtered.out', sep='\t', index=False)

    ##------------- Part2
    # 5. summary
    summary_PDXmodel_by_gene('Mut', df_filtered_mut, outdir)
    # recurrent
    cal_recurrentMut(df_filtered_mut, outdir)

    ##------------- Part3
    # 6. make table for plot heatmap
    PDXmodel_by_gene_for_complexHeatmap('Mut', df_filtered_mut, outdir)




##============= Recurrent mutations
def cal_recurrentMut(df, outdir):
    df2 = df[['CancerType', 'ModelID_proper', 'StudyArm', 'Gene', 'aaMut']].drop_duplicates().\
                        groupby(['CancerType', 'StudyArm', 'Gene', 'aaMut']).size().reset_index(name="N_pdxmodel_perCancer")
    
    df3 = df[['CancerType', 'ModelID_proper', 'StudyArm', 'Gene', 'aaMut']].drop_duplicates().\
                        groupby(['StudyArm', 'Gene', 'aaMut']).size().reset_index(name="total_pdxmodel")
    df2 = df2.merge(df3, on=['StudyArm', 'Gene', 'aaMut'], how='left')

    df_recurrent = df2[df2['total_pdxmodel'] > 1]
    df_recurrent['GeneMut'] = df_recurrent['Gene'] + ':' + df_recurrent['aaMut']
    df_recurrent.to_csv(outdir + '/2.somaticMut.filtered.recurrentMut.out', sep='\t', index=False)

    return df_recurrent




'''
    Copy number
'''

## Filter copy number
def ExtractTargetCopyNumber(cn_event, f_in, df_arm, df_info, df_remCase, outdir):
    df_cn = pd.read_csv(f_in, sep="\t", low_memory=False)
    df_cn = df_cn[['sample', 'gene', 'segment_log2', 'segment_cn']]

    ##------------- Part1
   # 1. extract mutation arm
    if cn_event == 'amp':
       df_arm_tar = df_arm[df_arm['Event'] == 'Amplification']
       df_cn_alt = df_cn[df_cn['segment_cn'] > 5]
    
    if cn_event == 'del':
       df_arm_tar = df_arm[df_arm['Event'] == 'Deletion']
       df_cn_alt = df_cn[df_cn['segment_log2'] < -1.3]


    targetGeneList = list(df_arm_tar['Gene'].unique())
    df_tar_cn = df_cn_alt[df_cn_alt['gene'].isin(targetGeneList)]

    # 2-1. format form for cn form
    df_tar_cn = df_tar_cn.merge(df_arm_tar[['Gene', 'StudyArm']], left_on='gene', right_on='Gene', how='left')
    # SampleID  Gene   StudyArm

    # 2-2. add case/model/group/cancertype
    df_tar_cn = df_tar_cn.merge(df_info, left_on='sample', right_on='Analysis_ID')
    # Sample  Gene   StudyArm Analysis_ID  Formal_CaseID  ModelID_proper  Group  CancerType

    # 3. filter
    df_tar_cn['flag'] = df_tar_cn['StudyArm'] + '_' + df_tar_cn['Formal_CaseID']
    df_filtered_cn = df_tar_cn[~df_tar_cn['flag'].isin(df_remCase['flag'])]
    
    # 4. outupt
    if cn_event == 'amp':
        df_tar_cn.to_csv(outdir + '/1.cn_amp.with_arm.out', sep='\t', index=False)
        df_filtered_cn.to_csv(outdir + '/1.cn_amp.filtered.out', sep='\t', index=False)
        ##------------- Part2
        # 5. summary
        summary_PDXmodel_by_gene('Amp', df_filtered_cn, outdir)
        ##------------- Part3
        # 6. make table for plot heatmap
        PDXmodel_by_gene_for_complexHeatmap('Amp', df_filtered_cn, outdir)


    if cn_event == 'del':
        df_tar_cn.to_csv(outdir + '/1.cn_del.with_arm.out', sep='\t', index=False)
        df_filtered_cn.to_csv(outdir + '/1.cn_del.filtered.out', sep='\t', index=False)
        ##------------- Part2
        # 5. summary
        summary_PDXmodel_by_gene('Del', df_filtered_cn, outdir)
        ##------------- Part3
        # 6. make table for plot heatmap
        PDXmodel_by_gene_for_complexHeatmap('Del', df_filtered_cn, outdir)
       

    
    
    
'''
    Fusion
'''

##=== Filter fusion
def FilterFusion(f_in, df_arm, df_info, df_remCase, outdir):

    df_fusion = pd.read_csv(f_in, sep="\t", low_memory=False)

    ##------------- Part1
    # 1. extract fusion by arm
    df_arm_tar = df_arm[df_arm['Event'] == 'Fusion']
    targetGeneList = list(df_arm_tar['Gene'].unique())
    df_tar_fusion = extractTargetFusion(targetGeneList, df_fusion)
 
    # 2-1. format form 
    df_tar_fusion = df_tar_fusion.merge(df_arm_tar[['Gene', 'StudyArm']], on='Gene', how='left')
    # Sample	Gene	FusionName	FusionType	StudyArm

    # 2-2. add case/model/group/cancertype
    df_tar_fusion = df_tar_fusion.merge(df_info, left_on='Sample', right_on='Analysis_ID')
    # Sample  Gene  FusionName	FusionType	StudyArm Analysis_ID  Formal_CaseID  ModelID_proper  Group  CancerType

    # 3-1. filter
    df_tar_fusion['flag'] = df_tar_fusion['StudyArm'] + '_' + df_tar_fusion['Formal_CaseID']
    df_filtered_fusion = df_tar_fusion[~df_tar_fusion['flag'].isin(df_remCase['flag'])]

    # 3-2 second filter for fusion --> FusionType == '.'
    df_filtered_fusion_2 = df_filtered_fusion[df_filtered_fusion['FusionType']!='.']
    
    # 4. outupt
    df_tar_fusion.to_csv(outdir + '/1.fusion.with_arm.out', sep='\t', index=False)
    df_filtered_fusion_2.to_csv(outdir + '/1.fusion.with_arm.filtered.out', sep='\t', index=False)

    
    ##------------- Part2
    # 5. summary
    summary_PDXmodel_by_gene('Fusion', df_filtered_fusion_2, outdir)

    ##------------- Part3
    # 6. make table for plot heatmap
    PDXmodel_by_gene_for_complexHeatmap('Fusion', df_filtered_fusion_2, outdir)



##======================= Extract target Fusion
def extractTargetFusion(geneList, df_starFusion):
    max_row = df_starFusion.shape[0]
    
    # create new dataframe
    d_matchedFusion = pd.DataFrame(columns=['Sample', 'Gene', 'FusionName', 'FusionType'])

    
    df_starFusion.reset_index(drop=True)
    
    for i in range(0, max_row):
        
        fusionName = df_starFusion.loc[i,'FusionName']
        fusionGeneList = fusionName.split('--')
        
        for gene in geneList:
            if gene in fusionGeneList:
                sample = df_starFusion.loc[i, 'Sample']
                fusionType = df_starFusion.loc[i, 'PROT_FUSION_TYPE']
                fusionType = df_starFusion.loc[i, 'PROT_FUSION_TYPE']
                rec = {'Sample':sample, 'Gene':gene, 'FusionName':fusionName, 'FusionType':fusionType}
                d_matchedFusion = d_matchedFusion.append(rec, ignore_index=True)  
            else:
                next
    
    return d_matchedFusion





##-------------------------------------------------------------------------------------------------##



'''
    Summary PDXmodel by gene
'''


def summary_PDXmodel_by_gene(altType, df, outdir):

    df_gene_NumModel_per_cancer = df[['CancerType', 'ModelID_proper', 'StudyArm', 'Gene']].drop_duplicates().\
                                        groupby(['CancerType', 'StudyArm', 'Gene']).size().reset_index(name="N_pdxmodel_perCancer")
        
    df_gene_NumModel_per_cancer['Event'] = altType

    df_gene_NumModel_per_cancer.to_csv(outdir + '/summaryPDXmodel.' + altType + '.mutation.byGene.out', sep="\t", index=False)
    # CancerType    StudyArm    Gene    N_PDXmodel   Event   

    return df_gene_NumModel_per_cancer





'''
    Format form for plot
'''

def PDXmodel_by_gene_for_complexHeatmap(altEvent, df, outdir):

    if altEvent=='Mut':
        df_gene_model_per_cancer = df[['Formal_CaseID', 'ModelID_proper', 'Gene', 'CancerType', 'StudyArm']].drop_duplicates()
        df_gene_model_per_cancer = df_gene_model_per_cancer.rename(columns={'Formal_CaseID':'CaseID', 'ModelID_proper':'PDXModelID', 'Hugo_Symbol':'Gene'})
        df_gene_model_per_cancer['Event'] = 'Mut'

    if (altEvent=='Amp') | (altEvent=='Del'):
        df_gene_model_per_cancer = df[['Formal_CaseID', 'ModelID_proper', 'gene', 'CancerType', 'StudyArm']].drop_duplicates()
        df_gene_model_per_cancer = df_gene_model_per_cancer.rename(columns={'Formal_CaseID':'CaseID', 'ModelID_proper':'PDXModelID', 'gene':'Gene'})
        if altEvent=='Amp':
            df_gene_model_per_cancer['Event'] = 'Amp'
        if altEvent=='Del':
            df_gene_model_per_cancer['Event'] = 'Del'

    if altEvent=='Fusion':
        df_gene_model_per_cancer = df[['Formal_CaseID', 'ModelID_proper', 'Gene', 'CancerType', 'StudyArm']].drop_duplicates()
        df_gene_model_per_cancer = df_gene_model_per_cancer.rename(columns={'Formal_CaseID':'CaseID', 'ModelID_proper':'PDXModelID'})
        df_gene_model_per_cancer['Event'] = 'Fusion'

    
    df_gene_model_per_cancer.to_csv(outdir + '/' + altEvent + '.nciMatchTrial.filtered.forTable.out', sep="\t", index=False)

    df_simple = df_gene_model_per_cancer[['PDXModelID', 'Gene', 'CancerType']]
    df_simple.to_csv(outdir + '/' + altEvent + '.nciMatchTrial.filtered.forComplexheatmapPlot.in', sep="\t", index=False)




if __name__ == '__main__':
    main()

