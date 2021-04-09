
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



Contact
-------------
Hua Sun, <hua.sun@wustl.edu>


