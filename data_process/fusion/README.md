
STAR-Fusion pipeline
---------------------


Install
---------
conda install -c bioconda star-fusion=1.6.0

Detail document for STAR-Fusion 
https://github.com/STAR-Fusion/STAR-Fusion/wiki


Usage
--------

```
$STAR_FUSION --genome_lib_dir $CTAT_LIB_DIR \
	--left_fq $FQ1 \
	--right_fq $FQ2 \
	--FusionInspector validate \
	--examine_coding_effect \
	--output_dir $starFusionOutDir
```

Run in MGI-server
--------

```
1. set the 'call.starFusion.full.v2.sh' location to the 'run.mgi.call.starFusion.byDir.sh'
2. sh run.mgi.call.starFusion.byDir.sh ./RNA/3.pdx_rna_filtered_fq
```

[Note] Filter non-cancer and normal fusions (from TCGA), please use the file 'https://github.com/ding-lab/PDX-PanCanAtlas/blob/master/data_process/appendix_db/filterFusionSet_noncancer_tcgaNor.v1.tsv'



Contact
-------------
Hua Sun, <hua.sun@wustl.edu>


