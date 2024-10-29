# EVCpipeline 
## Workflow Introduction

# RNA-seq_variant-calling
This is the workflow for RNA-seq germline variant calling based on [GATK RNAseq short variant discovery workflows](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) and [VAtools](https://vatools.readthedocs.io/en/latest/) for data cleaning.
## 
<embed src="[[https://drive.google.com/viewerng/viewer?embedded=true&url=YOUR_PDF_URL](https://github.com/Tina04021997/EVCpipeline/blob/main/workflow_logo/v0.1.pdf)](https://github.com/Tina04021997/EVCpipeline/blob/main/workflow_logo/v0.1.pdf)" width="500" height="375">

## Input data
- Paired-end fastq files
## Output data
-  tsv files transformed from annotated vcf files

## Reference data 
**mm10 resource bundle** see [Create GATK mm10 resource bundle](https://github.com/igordot/genomics/blob/master/workflows/gatk-mouse-mm10.md).

- If you encounter problems while concatenating dbSNP VCF files, try this:
```
bgzip -c vcf_chr_number.vcf > vcf_chr_number.vcf.gz
tabix  vcf_chr_number.vcf.gz
bcftools concat vcf_chr_number.vcf.gz vcf_chr_number.vcf.gz -Oz -o dbSNP.vcf.gz 
```
