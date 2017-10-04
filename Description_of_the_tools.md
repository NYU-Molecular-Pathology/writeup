# Description of the tools

Various software packages are used during the targeted exome sequencing and variant calling analysis used by the 580 gene panel. The bulk of these tools are scripted in the [`sns`](https://github.com/igordot/sns)(https://github.com/igordot/sns) pipeline, with extra custom steps included in the [`snsxt`](https://github.com/NYU-Molecular-Pathology/snsxt)(https://github.com/NYU-Molecular-Pathology/snsxt) pipeline. Automatic processing of sample data is accomplished with the [`lyz`](https://github.com/NYU-Molecular-Pathology/lyz)(https://github.com/NYU-Molecular-Pathology/lyz) program. 

# Trimmomatic Version 0.33

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/

http://www.usadellab.org/cms/?page=trimmomatic

A flexible read trimming tool for Illumina NGS data. Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.

## Settings Used

- `PE`: Paired End

- `-threads 8`: CPU threads

- `ILLUMINACLIP:/ref/contaminants/trimmomatic.fa:2:30:10:1:true`: Cut adapter and other illumina-specific sequences from the read.

- `TRAILING:5`: TRAILING: Cut bases off the end of a read, if below a threshold quality

- `SLIDINGWINDOW:4:15`: Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15

- `MINLEN:35`: Drop the read if it is below a specified length

## Sample Command

```
java -Xms16G -Xmx16G -jar /local/apps/trimmomatic/0.33/trimmomatic-0.33.jar PE 
-threads 8 SampleID_input_forward.fq.gz SampleID_input_reverse.fq.gz SampleID_output_forward_paired.fq.gz SampleID_output_forward_unpaired.fq.gz SampleID_output_reverse_paired.fq.gz SampleID_output_reverse_unpaired.fq.gzILLUMINACLIP:/ref/contaminants/trimmomatic.fa:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35 
```


# BWA Version: 0.7.13-r1126 

https://github.com/lh3/bwa

http://bio-bwa.sourceforge.net/bwa.shtml

BWA is a software package for mapping DNA sequences against a large reference genome, such as the human genome.

## Settings Used

- `bwa mem`: Align 70bp-1Mbp query sequences with the BWA-MEM algorithm.

- `-M`: Mark shorter split hits as secondary (for Picard compatibility).

- `-v 1`: Control the verbose level of the output. 1 for outputting errors only

- `-t 8`: Number of threads

- `-R ... `: Complete read group header line


## Sample Command

```
bwa mem -M -v 1 -t 8 -R '@RG\tID:SampleID\tSM:SampleID\tLB:SampleID\tPL:ILLUMINA' /ref/hg19/BWA/genome.fa SampleID_R1.trim.fastq.gz SampleID_R2.trim.fastq.gz 
```


# Sambamba v0.6.6

https://github.com/lomereiter/sambamba

http://lomereiter.github.io/sambamba/

http://lomereiter.github.io/sambamba/docs/sambamba-view.html

Tools for working with SAM/BAM data. 

Sambamba is a high performance modern robust and fast tool (and library), written in the D programming language, for working with SAM and BAM files. Current functionality is an important subset of samtools functionality, including view, index, sort, markdup, and depth.

## Settings Used

### View

- `view`: tool for extracting information from SAM/BAM/CRAM files

- `--sam-input`: Specify that the input is SAM file (default is BAM for all operations).

- `--nthreads=8`: Number of threads to use.

- `--filter='mapping_quality>=10'`: Set custom filter for alignments.

- `--format=bam`: Specify output format. FORMAT must be one of sam, bam, cram, or json (in lowercase). Default is SAM.

- `--compression-level=0`: Set compression level for BAM output, a number from 0 to 9.

- ``: 

### Sample Command

```
/software/sambamba/sambamba_v0.6.6 view --sam-input --nthreads=8 --filter='mapping_quality>=10' --format=bam --compression-level=0 /dev/stdin 
```

## Sort

- `sort`: tool for sorting BAM files

- `--nthreads=8`: Number of threads to use.

- `--memory-limit=16GB`: Sets an upper bound for used memory

- `--out`: Output file name

### Sample Command

```
/software/sambamba/sambamba_v0.6.6 sort --nthreads=8 --memory-limit=16GB --out=SampleID.bam /dev/stdin
```

## Markdup

- `markdup`: finding duplicate reads in BAM file

- `--nthreads=8`: Number of threads to use.

- `--remove-duplicates`: remove duplicates instead of just marking them

- `--hash-table-size 525000`: size of hash table for finding read pairs

- `--overflow-list-size 525000`: size of the overflow list where reads, thrown away from the hash table, get a second chance to meet their pairs


### Sample Command

```
/software/sambamba/sambamba_v0.6.6 markdup --remove-duplicates --nthreads 8 --hash-table-size 525000 --overflow-list-size 525000 SampleID.bam SampleID.dd.bam
```

## Flagstat

- `flagstat`: getting flag statistics from BAM file

### Sample Command

```
/software/sambamba/sambamba_v0.6.6 flagstat SampleID.dd.bam > SampleID.dd.bam.flagstat.txt
```

# GATK version 3.6-0-g89b7209

https://software.broadinstitute.org/gatk/

https://software.broadinstitute.org/gatk/documentation/

https://software.broadinstitute.org/gatk/documentation/tooldocs/3.6-0/

https://software.broadinstitute.org/gatk/documentation/tooldocs/3.6-0/org_broadinstitute_gatk_engine_CommandLineGATK.php

Genome Analysis Toolkit

Variant Discovery in High-Throughput Sequencing Data

Developed in the Data Sciences Platform at the Broad Institute, the toolkit offers a wide variety of tools with a primary focus on variant discovery and genotyping.

## Settings Used

### DepthOfCoverage

https://software.broadinstitute.org/gatk/documentation/tooldocs/3.6-0/org_broadinstitute_gatk_tools_walkers_coverage_DepthOfCoverage.php

- `-T DepthOfCoverage`: Assess sequence coverage by a wide array of metrics, partitioned by sample, read group, or library

- `-dt NONE`: Type of read downsampling to employ at a given locus

- `-rf BadCigar`: Filters to apply to reads before analysis

- `-nt 8`: Number of data threads to allocate to this analysis

- `--logging_level ERROR`: Set the minimum level of logging

- `--omitIntervalStatistics`: Do not calculate per-interval statistics

- `--omitLocusTable`: Do not calculate per-sample per-depth counts of loci

- `--omitDepthOutputAtEachBase`: Do not output depth of coverage at each base

- `-ct 10 -ct 100`: Coverage threshold (in percent) for summarizing statistics

- `-mbq 20`: Minimum quality of bases to count towards depth

- `-mmq 20`: Minimum mapping quality of reads to count towards depth

- `--reference_sequence`: Reference sequence file

- `--outputFormat csv`: The format of the output file

- `--input_file`: Input file containing sequence data (BAM or CRAM)

- `--out`: An output file created by the walker. Will overwrite contents if file exists

#### Sample Command

```
java -Xms16G -Xmx16G -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T DepthOfCoverage -dt NONE -rf BadCigar -nt  --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence /ref/hg19/genome.fa --input_file SampleID.dd.bam --outputFormat csv --out SampleID.genome
```

### RealignerTargetCreator

https://software.broadinstitute.org/gatk/documentation/tooldocs/3.6-0/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php

- `-T RealignerTargetCreator`: Define intervals to target for local realignment

- `-dt NONE`: Type of read downsampling to employ at a given locus

- `--logging_level ERROR`: Set the minimum level of logging

- `-nt 8`: Number of data threads to allocate to this analysis

- `--reference_sequence`: Reference sequence file

- `--known`: Input VCF file with known indels

- `--intervals`: One or more genomic intervals over which to operate

- `--interval_padding`: Amount of padding (in bp) to add to each interval

- `--input_file`: Input file containing sequence data (BAM or CRAM)

- `--out`: An output file created by the walker. Will overwrite contents if file exists

#### Sample Command

```
java -Xms16G -Xmx16G -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -dt NONE --logging_level ERROR -nt 8 --reference_sequence /ref/hg19/genome.fa -known /ref/hg19/gatk-bundle/1000G_phase1.indels.hg19.vcf -known /ref/hg19/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf --intervals /ifs/data/molecpathlab/NGS580_targets.bed --interval_padding 10 --input_file SampleID.dd.bam 
--out SampleID.intervals
```

### IndelRealigner

https://software.broadinstitute.org/gatk/documentation/tooldocs/3.6-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php

- `-T IndelRealigner`: Perform local realignment of reads around indels

- `-dt NONE`: Type of read downsampling to employ at a given locus

- `--logging_level ERROR`: Set the minimum level of logging

- `--reference_sequence`: Reference sequence file

- `--maxReadsForRealignment 50000`: Max reads allowed at an interval for realignment

- `--known`: Input VCF file with known indels

- `--targetIntervals`: Intervals file output from RealignerTargetCreator

- `--input_file`: Input file containing sequence data (BAM or CRAM)

- `--out`: An output file created by the walker. Will overwrite contents if file exists

#### Sample Command

```
java -Xms16G -Xmx16G -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T IndelRealigner -dt NONE --logging_level ERROR --reference_sequence /ref/hg19/genome.fa --maxReadsForRealignment 50000 -known /ref/hg19/gatk-bundle/1000G_phase1.indels.hg19.vcf -known /ref/hg19/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf -targetIntervals SampleID.intervals --input_file SampleID.dd.bam --out SampleID.dd.ra.bam
```


### BaseRecalibrator

https://software.broadinstitute.org/gatk/documentation/tooldocs/3.6-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php


- `-T BaseRecalibrator`: Detect systematic errors in base quality scores

- `--logging_level ERROR`: Set the minimum level of logging

- `-nt 8`: Number of data threads to allocate to this analysis

- `-rf BadCigar`: Filters to apply to reads before analysis

- `--reference_sequence`: Reference sequence file

- `--knownSites`: A database of known polymorphic sites

- `--intervals`: One or more genomic intervals over which to operate

- `--interval_padding`: Amount of padding (in bp) to add to each interval

- `--input_file`: Input file containing sequence data (BAM or CRAM)

- `--out`: An output file created by the walker. Will overwrite contents if file exists

- `-BQSR`: Base quality score recalibration

#### Sample Command

```
java -Xms16G -Xmx16G -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T BaseRecalibrator --logging_level ERROR -nct 8 -rf BadCigar --reference_sequence /ref/hg19/genome.fa -knownSites /ref/hg19/gatk-bundle/1000G_phase1.indels.hg19.vcf -knownSites /ref/hg19/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf -knownSites /ref/hg19/gatk-bundle/dbsnp_138.hg19.vcf --intervals /ifs/data/molecpathlab/NGS580_targets.bed --interval_padding 10 --input_file SampleID.dd.ra.bam 
--out SampleID.table1.txt

java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T BaseRecalibrator --logging_level ERROR -nct 8 -rf BadCigar --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa -knownSites /ifs/home/id460/ref/hg19/gatk-bundle/1000G_phase1.indels.hg19.vcf -knownSites /ifs/home/id460/ref/hg19/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf -knownSites /ifs/home/id460/ref/hg19/gatk-bundle/dbsnp_138.hg19.vcf --intervals /ifs/data/molecpathlab/NGS580_targets.bed --interval_padding 10 --input_file SampleID.dd.ra.bam -BQSR SampleID.table1.txt --out SampleID.table2.txt
```

### AnalyzeCovariates

https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php

- `-T AnalyzeCovariates`: Create plots to visualize base recalibration results

- `--logging_level ERROR`: Set the minimum level of logging

- `--reference_sequence`: Reference sequence file

- `-before`: file containing the BQSR first-pass report file

- `-after`: file containing the BQSR second-pass report file

- `-csv`: location of the csv intermediate file

- `-plots`: location of the output report


#### Sample Command

```
java -Xms16G -Xmx16G -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T AnalyzeCovariates 
--logging_level ERROR --reference_sequence /ref/hg19/genome.fa -before SampleID.table1.txt -after SampleID.table2.txt -csv SampleID.csv -plots SampleID.pdf
```

### PrintReads

https://software.broadinstitute.org/gatk/documentation/tooldocs/3.6-0/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php

- `-T PrintReads`: Write out sequence read data (for filtering, merging, subsetting etc)

- `--logging_level ERROR`: Set the minimum level of logging

- `--reference_sequence`: Reference sequence file

- `-nt 8`: Number of data threads to allocate to this analysis

- `-rf BadCigar`: Filters to apply to reads before analysis

- `--input_file`: Input file containing sequence data (BAM or CRAM)

- `--out`: An output file created by the walker. Will overwrite contents if file exists

- `-BQSR`: Base quality score recalibration

#### Sample Command

```
java -Xms16G -Xmx16G -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T PrintReads --logging_level ERROR -nct 8 -rf BadCigar --reference_sequence /ref/hg19/genome.fa -BQSR SampleID.table1.txt --input_file SampleID.dd.ra.bam --out SampleID.dd.ra.rc.bam
```

### HaplotypeCaller

https://software.broadinstitute.org/gatk/documentation/tooldocs/3.6-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

- `-T HaplotypeCaller`: Call germline SNPs and indels via local re-assembly of haplotypes

- `-dt NONE`: Type of read downsampling to employ at a given locus

- `--logging_level ERROR`: Set the minimum level of logging

- `-nct 8`: Number of CPU threads to allocate per data thread

- `--max_alternate_alleles 3`: Maximum number of alternate alleles to genotype

- `--standard_min_confidence_threshold_for_calling 50`: The minimum phred-scaled confidence threshold at which variants should be called

- `--reference_sequence`: Reference sequence file

- `--intervals`: One or more genomic intervals over which to operate

- `--interval_padding`: Amount of padding (in bp) to add to each interval

- `--input_file`: Input file containing sequence data (BAM or CRAM)

- `--out`: An output file created by the walker. Will overwrite contents if file exists

#### Sample Command

```
java -Xms16G -Xmx16G -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -dt NONE --logging_level ERROR -nct 8 --max_alternate_alleles 3 --standard_min_confidence_threshold_for_calling 50 --reference_sequence /ref/hg19/genome.fa --intervals /ifs/data/molecpathlab/NGS580_targets.bed --interval_padding 10 --input_file SampleID.dd.ra.rc.bam --out SampleID.original.vcf

```

### MuTect2

https://software.broadinstitute.org/gatk/documentation/tooldocs/3.6-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php

- `-dt NONE`: Type of read downsampling to employ at a given locus

- `--logging_level WARN`: Set the minimum level of logging

- `--standard_min_confidence_threshold_for_calling 30`: The minimum phred-scaled confidence threshold at which variants should be called

- `--max_alt_alleles_in_normal_count 10`: Threshold for maximum alternate allele counts in normal

- `--max_alt_allele_in_normal_fraction 0.05`: Threshold for maximum alternate allele fraction in normal

- `--max_alt_alleles_in_normal_qscore_sum 40`: Threshold for maximum alternate allele quality score sum in normal

- `--reference_sequence`: Reference sequence file

- `--dbsnp`: dbSNP file. rsIDs from this file are used to populate the ID column of the output.

- `--cosmic`: VCF file of COSMIC sites

- `--intervals`: One or more genomic intervals over which to operate

- `--interval_padding`: Amount of padding (in bp) to add to each interval

- `--input_file`: Input file containing sequence data (BAM or CRAM)

- `--out`: An output file created by the walker. Will overwrite contents if file exists

#### Sample Command

```
java -Xms16G -Xmx16G -jar /software/GenomeAnalysisTK/GenomeAnalysisTK-3.6-0/GenomeAnalysisTK.jar -T MuTect2 -dt NONE --logging_level WARN --standard_min_confidence_threshold_for_calling 30 --max_alt_alleles_in_normal_count 10 --max_alt_allele_in_normal_fraction 0.05 --max_alt_alleles_in_normal_qscore_sum 40 --reference_sequence /ref/hg19/genome.fa --dbsnp /ref/hg19/gatk-bundle/dbsnp_138.hg19.vcf --cosmic /ref/hg19/CosmicCodingMuts_v73.hg19.vcf --intervals /ifs/data/molecpathlab/NGS580_WES/NGS580_targets.bed --interval_padding 10 --input_file:tumor SampleID-T.dd.ra.rc.bam --input_file:normal SampleID-N.dd.ra.rc.bam --out SampleID-T_SampleID-N.original.vcf
```


# LoFreq version 2.1.2

http://csb5.github.io/lofreq/

https://github.com/CSB5/lofreq

https://www.ncbi.nlm.nih.gov/pubmed/23066108

LoFreq is a fast and sensitive variant-caller for inferring SNVs and indels from next-generation sequencing data. It makes full use of base-call qualities and other sources of errors inherent in sequencing (e.g. mapping or base/indel alignment uncertainty), which are usually ignored by other methods or only used for filtering.

## Settings Used

- `call-parallel`: Call variants in parallel

- `--call-indels`: also call indels

- `--pp-threads 8`: the number of threads

- `--ref`: reference sequence

- `--bed`: target regions

- `--out`: output file

### Sample Command

```
lofreq call-parallel --call-indels --pp-threads 8 --ref /ref/hg19/genome.fa --bed /ifs/data/molecpathlab/NGS580_targets.pad10.bed --out SampleID.original.vcf SampleID.dd.ra.rc.bam
```

# bcftools version 1.3

https://samtools.github.io/bcftools/

https://samtools.github.io/bcftools/bcftools.html

http://www.htslib.org/download/

BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed.

## Settings Used

- `index`: index VCF/BCF

- `norm`: normalize indels

- `view`: subset, filter and convert VCF and BCF files

- `--multiallelics`: split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+)

- `-both`: abbreviation of "-c indels  -c snps"

- `-c snps`: any SNP records are compatible, regardless of whether the ALT alleles match or not. For duplicate positions, only the first SNP record will be considered and appear on output.

- `-c indels`: all indel records are compatible, regardless of whether the REF and ALT alleles match or not. For duplicate positions, only the first indel record will be considered and appear on output.

- `--output-type v`: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion.

- `--fasta-ref`: reference sequence in fasta format

- `--exclude 'DP<5' `: exclude sites for which EXPRESSION is true. 

### Sample Command

```
bcftools index SampleID.original.vcf.bgz 

bcftools norm --multiallelics -both --output-type v SampleID.original.vcf.bgz | bcftools norm --fasta-ref /ref/hg19/genome.fa --output-type v - | bcftools view --exclude 'DP<5' --output-type v > SampleID.vcf
```

# ANNOVAR version 2016-02-01 00:11:18 -0800 (Mon,  1 Feb 2016)

http://annovar.openbioinformatics.org/en/latest/

ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes (including human genome hg18, hg19, hg38, as well as mouse, worm, fly, yeast and many others)

## Settings Used

- `--format vcf4old`: input format (default: pileup)

- `--includeinfo`: include supporting information in output

- `--outfile`: specify the output file name.

- `--buildver hg19`: genome build version

- `--protocol refGene,snp138,snp138NonFlagged,exac03,esp6500siv2_all,1000g2015aug_all,cosmic80,cadd13gt10,fathmm`: comma-delimited string specifying database protocol

- `--operation g,f,f,f,f,f,f,f,f`: comma-delimited string specifying type of operation

- `--argument '--splicing_threshold 10',,,,,,,,`: comma-delimited strings as optional argument for each operation

- `--nastring .`: string to display when a score is not available (default: null)

- `--remove`: remove all temporary files

### Sample Command

```
perl /software/annovar/annovar-160201/convert2annovar.pl --format vcf4old --includeinfo SampleID.vcf > SampleID.avinput

perl /software/annovar/annovar-160201/table_annovar.pl SampleID.avinput /ref/annovar/hg19/ --outfile SampleID --buildver hg19 --protocol refGene,snp138,snp138NonFlagged,exac03,esp6500siv2_all,1000g2015aug_all,cosmic80,cadd13gt10,fathmm --operation g,f,f,f,f,f,f,f,f --argument '--splicing_threshold 10',,,,,,,, --nastring . --remove
```

# Rsamtools 

http://bioconductor.org/packages/release/bioc/html/Rsamtools.html

This package provides an interface to the 'samtools', 'bcftools', and 'tabix' utilities (see 'LICENCE') for manipulating SAM (Sequence Alignment / Map), FASTA, binary variant call (BCF) and compressed indexed tab-delimited (tabix) files.
