# EVCpipeline 
## Workflow Introduction
<img src="https://github.com/Tina04021997/EVCpipeline/blob/main/workflow_logo/v0.1.png" width="95%" height="95%">

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
