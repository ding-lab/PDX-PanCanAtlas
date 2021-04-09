

A mouse reads filtering pipeline from WGS/WES and RNA-seq data


Install
-------
   * conda install -c anaconda openjdk=8.0.152
   * conda install -c bioconda bwa=0.7.17
   * conda install -c bioconda samtools=1.5
   * conda install -c bioconda picard=2.17.11
   * conda install -c bioconda ngs-disambiguate=2018.05.03
   * conda install -c bioconda star=2.6.1a
	
	* Download PICARD from Broad Institute (2.17.11)
	* Java (jre) install 8.0.152

Reference
-------

To easy set references, you can download human and mouse references (References - 2020-A (July 7, 2020). It includes genome and STAR DB.) from https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest


Note: In the PDX-Pilot, we used GENCODE GRCh38.v29 and GRCm38.vM19 genome reference.


Set config.ini
-------

To run script, it must set `config.ini`. Please refer to `config/config.pdx.mgi.v2.ini`



Usage
-----

        Version 1.0
        Usage:  pdx_disam_kit.sh -p <command> [options]
	
	Note: Please set path for "scriptDir=" and "config=" to pdx_disam_kit.sh script before use.

Key commands:

        dnaFull            Do full step from *.fq.gz to create new bam 
        rnaFull            Do full step from *.fq.gz to create new *.fq.gz
        
        humanBam           Separate step - Only create human bam
        mouseBam           Separate step - Only create mouse bam
        disambiguate       Separate step - Do disambiguate

Options:
        
        -p  <command>      Key command
        -n  <name>         Any name for building folder
        -1  <fastq.gz>     Fastq gzip file
        -2  <fastq.gz>     Fastq gzip file
        -o  <directory>    Output directory
        -d  <directory>    Directory. It use to separate step.


Example
-------

1. WGS/WES full step
        
        sh pdx_disam_kit.sh -p dnaFull -n sample -1 fq1.gz -2 fq2.gz -o outdir

2. RNA full step
        
        sh pdx_disam_kit.sh -p rnaFull -n sample -1 fq1.gz -2 fq2.gz -o outdir


3. Separate steps
   
   * Make human to bam
   
     `sh pdx_disam_kit.sh -p humanBam -n sample -1 fq1.gz -2 fq2.gz -o outdir`

   * Make mouse to bam
   
     `sh pdx_disam_kit.sh -p mouseBam -n sample -1 fq1.gz -2 fq2.gz -o outdir`

   * Disambiguate
   
     `sh pdx_disam_kit.sh -p disambiguate -n sample -d pdx_wxs_dir`


Input
-------
* WGS/WES/RNA-seq
  
  ```
  Pair fastq gzip files (e.g. *.1.fq.gz / *.2.fq.gz)
  ```

Output
-------
* WGS/WES data output
        
        *.disam.reAlign.remDup.bam
        *.disam.reAlign.remDup.bam.bai

* RNA-seq data output
        
        *.disam.rna-seq.1.fastq.gz
        *.disam.rna-seq.2.fastq.gz


### Docker file

* docker/Dockerfile



### CWL version

The CWL version developed by Matthew Wyczalkowski
(https://github.com/ding-lab/MouseTrap2)



Contact
-------
Hua Sun, <hua.sun@wustl.edu>
