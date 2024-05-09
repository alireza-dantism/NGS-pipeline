# NGS pipeline
Just a number of steps to analyze data

# Step 1 (Linux Commands)

## 1-1) Preparing Data

Searching project number: PRJNA1106990 in https://www.ncbi.nlm.nih.gov/sra

In the opened page under study section click on all runs

Then we need all SRR numbers

_Project number is in paper_

**NCBI Run Browser:**

https://www.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&display=metadata
https://www.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR000002&display=metadata

## 1-2) Download, install and use of SRAToolkit library

Open SRA Toolkit bin Directory:

```linux
./vdb-config --interactive
```
We run this command because it set a some configuration that do not download unnecessary files during prefetch command

```linux
./prefetch -p SRR—Accession-Number O ~/Desktop
```
prefetch <accession> will download the <accession> run file and its dependencies
* Input => accession number such as SRR0000001
* Output => .sra file of accession number such as SRR0000001.sra

```linux
./fastq-dump SRR—Accession-Number.sra —split-3 --outdir ~/Desktop
```
can be used to convert the prefetched Runs from compressed SRA format to fastq
* Input => .sra file from previous step
* Output => .fastq files
    * SRR0000001_1.fastq => left pair
    * SRR0000001_2.fastq => right pair
    * SRR0000001.fastq => unpaired 
 
 **we can remove unpaired files**

## 1-3) Download, install and use of FastQC library

FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

Website: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

```linux
fastqc -f fastq SRR0000001_1.fastq SRR0000001_2.fastq
```


* Input => _1.fastq & _2.fasta files
* Output => Visualizing quality of samples based on .html files

FastQC Manual: https://mugenomicscore.missouri.edu/PDF/FastQC_Manual.pdf

## 1-4) Download, install and use of Trimmomatic library

Trimmomatic Manual: V0.32 => http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

Commands: http://www.usadellab.org/cms/?page=trimmomatic

Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop Illumina (FASTQ) data as well as to remove adapters. These adapters can pose a real problem depending on the library preparation and downstream application.

```linux
java -jar <path to trimmomatic.jar> PE *_1.fastq *_2.fastq -baseout ../trimmingOutput/*_trimmed_reads.fastq LEADING:25 TRAILING:25 SLIDINGWINDOW:4:25 MINLEN:50
```

* Input => _1.fastq & _2.fast1
* Output => 
    * _trimmed_1.fastq
    * _untrimmed_1.fastq
    * _trimmed_2.fastq
    * _untrimmed_2.fastq

**Surviving both should be more than 80%**

More advanced:

```linux
for i in *_1.fastq; do echo $i; java -jar ../../Trimmomatic-0.39/trimmomatic-0.39.jar PE $i ${i/_1/_2} -baseout ../trimmingOutput/${i/_chrX*/} LEADING:25 TRAILING:25 SLIDINGWINDOW:4:20 MINLEN:50; done
```

---  

# Step 2 (Linux Commands)

## 2-1) 



