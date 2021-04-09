'''

    Hua Sun
    8/9/2020

    python3 run.py --var sample.alt.tsv --info data.info --outdir outdir

'''


import pandas as pd
import numpy as np
import os
import argparse
import plotly.graph_objects as go

parser = argparse.ArgumentParser()
parser.add_argument('-v','--var', type=str, help='sample alteration data')
parser.add_argument('-i','--info', type=str, help='data info')
parser.add_argument('-o','--outdir', type=str, help='out dir')

args = parser.parse_args()





def main():
    # outdir
    outdir = args.outdir

    # read data
    df_var = pd.read_csv(args.var, sep='\t')
    df_var.columns = ['analysisID', 'studyArm']

    df_info = pd.read_csv(args.info, sep='\t')
    df_info.columns = ['caseID', 'pdxmodelID', 'sampleID', 'passage', 'analysisID', 'dataType']
    # caseID	pdxmodelID	sampleID	passage	analysisID	dataType


    # merge event & info
    df_var = df_var.merge(df_info, on='analysisID')
    # analysisID  studyArm  caseID  pdxmodelID  sampleID passage dataType
    #df_var.to_csv(args.out, sep='\t', index=False)



    ##=========== S1: Assess PDX-Model level 
    df_res = estimatePDX_model(df_var, df_info)

    outfile = outdir + '/estimated.PDX_model.results.out'
    df_res.to_csv(outfile, sep='\t', index=False)
    
    # visulization
    plotly_pie_chart_stable(df_res, outdir+'/estimated.PDX_model.note.pdf')
    plotly_pie_chart_multipleArm(df_res, outdir+'/estimated.PDX_model.studyarm.pdf')


    ##=========== S2: Assess PDX-sample level 
    # visulization
    plot_withArmPDXs_from_1760pdx(df_var, df_info, outdir+'/estimated.PDX_sample.pdf')





def estimatePDX_model(df_var, df_info):
    pdxmodelList = list(df_var['pdxmodelID'].unique())

    df_rec = pd.DataFrame(columns=['StudyArm', 'ArmEventPerModel', 'CaseID', 'PDXmodelID', 'TotalSizeUniqPassage', 'TotalUniqPassage', 'VarSizeUniqPassage', 'VarUniqPassage', 'FreqVarInPassage', 'Catalog'])
    for model in pdxmodelList:
        varEvent = list(df_var[df_var.pdxmodelID==model]['studyArm'].unique())
        Num_of_varEvent = len(varEvent)
        caseid = ''.join(df_var[df_var.pdxmodelID==model]['caseID'].unique())
        # total unique passage
        total_passage = list(df_info[df_info.pdxmodelID==model]['passage'].unique())
        total_passage.sort()
        total_num_passage = len(total_passage)

        # varianted unique passage
        var_passage = list(df_var[df_var.pdxmodelID==model]['passage'].unique())
        var_passage.sort()
        var_num_passage = len(var_passage)

        # var % for unique passage
        freq = round(var_num_passage/total_num_passage*100, 2)

        # summary
        note = 'Not_defined'
        if freq >= 60 :
            note = 'Stable'
        else:
            note = 'Unstable'
        
        if (total_num_passage == 2) and (freq == 50):
            note = 'Equivocal'
        
        if (total_num_passage > 3 and ['P0','P1'] in total_passage):
            if ['P0','P1'] in var_passage:
                note = 'Early_passage'
            elif ('P0' not in var_passage) and ('P1' not in var_passage) and (var_num_passage > 1):
                note = 'Later_passage'
        

        # format output list as str
        total_passage = ','.join(total_passage)
        var_passage = ','.join(var_passage)
        varEvent = ','.join(varEvent)

        # record
        rec = {'StudyArm':varEvent, 'ArmEventPerModel':Num_of_varEvent, 'CaseID':caseid, 'PDXmodelID':model, \
             'TotalSizeUniqPassage':total_num_passage, 'TotalUniqPassage':total_passage, \
             'VarSizeUniqPassage':var_num_passage , 'VarUniqPassage':var_passage, 'FreqVarInPassage':freq, 'Catalog':note}
        df_rec = df_rec.append(rec, ignore_index=True)
    
    return df_rec





def plot_withArmPDXs_from_1760pdx(df_var, df_info, outname):
    pdxmodelist = list(df_var['pdxmodelID'].unique())
    total_target_n = len(df_info[df_info['pdxmodelID'].isin(pdxmodelist)]['sampleID'].unique())
    print('\n' + 'total_tar_sample: ' + str(total_target_n))

    total_var_sample_size = len(df_var.sampleID.unique())
    wt = total_target_n - total_var_sample_size
    d = {'Catalog':['WT', 'Arm-event'], 'Freq':[wt, total_var_sample_size]}
    df_rec = pd.DataFrame(data=d)
    print('total_var_sample: ' + str(total_var_sample_size) + '\n')

    fig = go.Figure(data=[go.Pie(labels=df_rec.Catalog, values=df_rec.Freq, textinfo='label+percent',
                             insidetextorientation='radial'
                            )])
    fig.update_layout(
        width=400,
        height=400
    )

    fig.write_image(outname)






# must including column-name: Catalog
def plotly_pie_chart_stable(df, outname):
    table = df.groupby('Catalog').size().reset_index(name="Freq")
    print(table)

    fig = go.Figure(data=[go.Pie(labels=table.Catalog, values=table.Freq, textinfo='label+percent',
                             insidetextorientation='radial'
                            )])
    fig.update_layout(
        width=400,
        height=400
    )

    fig.write_image(outname)
    


# must including column-name: StudyArm
def plotly_pie_chart_multipleArm(df, outname):
    table = pd.Series([str(x.count(',')+1)+'-Arm' for x in df.StudyArm]).value_counts().reset_index(name='Freq')
    table.columns = ['arm', 'freq']
    fig = go.Figure(data=[go.Pie(labels=table.arm, values=table.freq, textinfo='label+percent',
                             insidetextorientation='radial'
                            )])
    fig.update_layout(
        width=400,
        height=400
    )
    
    fig.write_image(outname)






if __name__ == '__main__':
    main()

