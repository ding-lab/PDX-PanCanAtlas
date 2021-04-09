MSIsensor for tumor-normal pipeline
-----------------------


Install
--------

conda install -c bioconda msisensor=0.5

Detail document https://github.com/ding-lab/msisensor/blob/master/README_msisensor.md


Usage
--------

`msisensor msi -d ${microsatellite} -n ${normal_bam} -t ${tumor_bam} -l 1 -q 1 -o ${OUT}/${sample}.prefix`



MSIsensor2 tumor-only pipeline
----------------------


Install
--------

https://github.com/niu-lab/msisensor2


Usage
--------

`msisensor2 msi -M $models_hg38 -t $bam -l 1 -q 1 -o $out/$sample.msi`



Contact
-------------
Hua Sun, <hua.sun@wustl.edu>






