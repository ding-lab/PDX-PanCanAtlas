'''

    Hua Sun
    8/14/2020


    1) For each target, calculate the % of passages display positive signal for the target.
    2) Draw a curve representing the distribution of the % from 1.

    python3 run.py --var sample.alt.tsv --info data.info --drug nci-match_eay131.arm_drug.index --outdir outdir

'''


import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go


parser = argparse.ArgumentParser()
parser.add_argument('-v','--var', type=str, help='sample alteration data')
parser.add_argument('-i','--info', type=str, help='data info')
parser.add_argument('-d','--drug', type=str, help='drug and study arm info')
parser.add_argument('-o','--outdir', type=str, help='out dir')

args = parser.parse_args()


def main():
    # outdir
    outdir = args.outdir

    # read data
    df_var = pd.read_csv(args.var, sep='\t')
    df_var.columns = ['analysisID', 'studyArm']
    print(df_var.shape)

    df_info = pd.read_csv(args.info, sep='\t')
    df_info.columns = ['caseID', 'pdxmodelID', 'sampleID', 'passage', 'analysisID', 'dataType']
 
    # merge event & info
    df_var = df_var.merge(df_info, on='analysisID')
  
    # read drug data
    df_drug = pd.read_csv(args.drug, sep='\t')
    df_drug.columns = ['studyArm', 'drug']
    df_var = df_var.merge(df_drug, on='studyArm')
    # analysisID  studyArm  caseID  pdxmodelID  sampleID passage dataType drug

    print(df_var.shape)



     df_res = estimate_targetArm_covered_passage(df_var, df_info)

    outfile = outdir + '/estimated.targetArm_covered_passage.out'
    df_res.to_csv(outfile, sep='\t', index=False)
    
    # visulization
    plotly_lineplot(df_res, outdir+'/estimated.targetArm_covered_passage.pdf')
    



def estimate_targetArm_covered_passage(df_var, df_info):
    targetArmList = list(df_var['studyArm'].unique())


    df_rec = pd.DataFrame(columns=['StudyArm', 'TotalUniqPassage', 'DetectedUniqPassage', 'TotalSize_uniqPassage', 'DetectedSize_uniqPassage', 'DetectedPercentage'])
    for targetArm in targetArmList:

        # drugname 
        drugname = list(df_var[df_var['studyArm']==targetArm]['drug'].unique())[0]

        # varianted unique passage
        var_passage = list(df_var[df_var['studyArm']==targetArm]['passage'].unique())
        var_passage.sort()
        var_num_passage = len(var_passage)

        # total passages, which target arm related all passages from raw pdx-models
        pdxmodelList = list(df_var[df_var['studyArm']==targetArm]['pdxmodelID'].unique())
        totalPassage = list(df_info[df_info['pdxmodelID'].isin(pdxmodelList)]['passage'].unique())
        totalPassage.sort()
        totalPassage_size = len(totalPassage)

        # var % for unique passage
        freq = round(var_num_passage/totalPassage_size*100, 2)

        # format output list as str
        totalPassage = ','.join(totalPassage)
        var_passage = ','.join(var_passage)

        # record
        drug_with_arm = drugname + ' (' + targetArm + ')'
        rec = {'StudyArm':drug_with_arm, 'TotalUniqPassage':totalPassage, 'DetectedUniqPassage':var_passage, 'TotalSize_uniqPassage':totalPassage_size, \
             'DetectedSize_uniqPassage':var_num_passage ,  'DetectedPercentage':freq}
        df_rec = df_rec.append(rec, ignore_index=True)
    
    df_rec.sort_values('DetectedPercentage', ascending=False, inplace=True)

    return df_rec



def plotly_lineplot(df, outfile):

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df.StudyArm, y=df.DetectedPercentage,
                         line = dict(color='firebrick', width=4)))

    fig.update_layout(
                   xaxis_title='Drug (target Arm)',
                   yaxis_title='Coverage of PDX passages (%)')
    
    # change figure size                  
    fig.update_layout(
        width=500,
        height=300,
        font_size=8,
        margin=go.layout.Margin(
            l=10, #left margin
            r=10, #right margin
            b=150, #bottom margin
            t=10  #top margin
        )
    )

    fig.write_image(outfile)







if __name__ == '__main__':
    main()
