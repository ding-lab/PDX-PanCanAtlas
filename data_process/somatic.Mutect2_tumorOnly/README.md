
## Mutect2 - call tumor-only somatic mutations

GRCh38 version



Install
-------

* [jre1.8.0_152](https://www.oracle.com/java/technologies/javase/javase8-archive-downloads.html)
* conda install -c bioconda samtools=1.5  (or later version)
* [PICARD v2.7.11](https://github.com/broadinstitute/picard/releases/tag/2.17.11)
* [SNPSIFT v4_3t](https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip/download)
* [GATK v4.1.2.0](https://github.com/broadinstitute/gatk/releases/download/4.1.2.0/gatk-4.1.2.0.zip)
	* Or [Docker image for GATK](https://hub.docker.com/r/broadinstitute/gatk/)
   

Note: In the PDX-pilot project, we used GATK v4.1.2. Please use the latest GATK version based on your project.
<https://github.com/broadinstitute/gatk/releases>




Prepare DB for Mutect2
-------------------

* [GECODE - GRCh38_v29](https://www.gencodegenes.org/human/release_29.html)

	> GRCh38.primary_assembly.genome.fa
	
* [VEP v85](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_download.html#old_vep_tool)
	
	> to be consistant with somaticwrapper 1.5 pipeline

* Download `af-only-gnomad.hg38.vcf.gz` and `af-only-gnomad.hg38.vcf.gz.tbi` from <https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false>

* Download Pool normal VCF - gatk4_mutect2_4136_pon.vcf.tar
	<https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files>

* Download dbSNP (b151) and filter somatic mutation sites using COSMIC_v89 VCF (download both of CosmicCoding and NonCoding VCFs from COSMIC). The filtered file is too big (*.vcf.gz >15GB) to share.

* [Download common_biallelic DB](https://sourceforge.net/projects/mutect2-data/)




Set config.ini
-------

To run script, it must set `config.ini`. Please refer to `config.mutect2.mgi.ini`



Usage
-----


Note: Please set "scriptDir=" and/or "config=" into the  `somaticMut.tumor-only.mutect2.sh` before running pipeline. If it only sets "scriptDir=" then you have to use '-c' parameter to set config.ini file every time.

```
sh somaticMut.tumor-only.mutect2.sh -c <config.ini> -p <programname> -n <name> -b <bam> -o <outdir>

# run-full steps by one (-p s0)
sh somaticMut.tumor-only.mutect2.sh -p full -n sample -b bam -o outdir

# re-run
sh somaticMut.tumor-only.mutect2.sh -p full -r yes -n sample -b bam -o outdir

```




Contact
-------
Hua Sun, <hua.sun@wustl.edu>
