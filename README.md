# NGS pipeline
Just a number of steps to analyze data

# Step 1 (Linux Commands)

## 1-1) Preparing Data

Searching project number: PRJNA000001 in https://www.ncbi.nlm.nih.gov/sra

In the opened page under study section click on all runs

Then we need all SRR numbers

_Project number is in paper_

**NCBI Run Browser:**

https://www.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&display=metadata
https://www.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR000002&display=metadata

## 1-2) Download, install and use of SRAToolkit library

The SRA Toolkit and SDK from NCBI is a collection of tools and libraries for using data in the INSDC Sequence Read Archives.

https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

Open SRA Toolkit bin Directory:

```linux
./vdb-config --interactive
```
We run this command because it set a some configuration that do not download unnecessary files during prefetch command

```linux
./prefetch -p SRR—Accession-Number -O ~/Desktop
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

## 2-1) Download, install and use of HISAT2 library

HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) to a population of human genomes as well as to a single reference genome. 
HISAT2 uses a large set of small GFM indexes that collectively cover the whole genome. These small indexes (called local indexes), combined with several alignment strategies, enable rapid and accurate alignment of sequencing reads. This new indexing scheme is called a Hierarchical Graph FM index (HGFM).

Website: https://daehwankimlab.github.io/hisat2/

Download Link: https://daehwankimlab.github.io/hisat2/download/

HISAT2 HowTo: https://daehwankimlab.github.io/hisat2/howto/

```linux
hisat2 -q -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR*_chrX_1.fastq -2 chrX_data/samples/ERR1*_chrX_2.fastq -S ERR*_chrX.sam
```

* Input = > 
    * ‘-x’ the indexed reference genome directories (use prefix file name)
    * ‘-1’ and ‘-2’ are used to denote our fwd and rev samples in a paired-end alignment
    * ‘-S’ is used to denote that we would like our output in SAM format
    * ‘-p 8’ denoting the use of 8 threads, ‘–dta’ is used to generate output SAM files
* Output => .sam file


More advanced:

Here is a bash script for the above HISAT2 command called hisat2.sh that will run all the .fastq.gz files for you simultaneously.

```linux
#!/usr/bin/bash

#bash script for hisat2; align all .fastq.gz files to indexed reference genome to generate .sam files

SAMPLES="ERR188044 ERR188104 ERR188234 ERR188245 ERR188257"

for SAMPLE in $SAMPLES; do
    hisat2 -p 11 --dta -x ~/chrX_data/indexes/chrX_tran -1 ~/chrX_data/samples/${SAMPLE}_chrX_1.fastq.gz -2 ~/chrX_data/samples/${SAMPLE}_chrX_2.fastq.gz -S ${SAMPLE}_chrX.sam
done
```

— OR —

```linux
for i in *_1P; do echo $i; ../../hisat2-2.2.1/hisat2 -q -x ../indexes/chrX_tran -1 $i -2 ${i/_1/_2} --add-chrname -S ../hisat2Output/${i/_1P/.sam}; done
```

If we have not index files use below first:

```linux
hisat2-build [options]* <reference_in> <ht2_base>
```

hisat2-build outputs a set of 6 files with suffixes .1.ht2, .2.ht2, .3.ht2, .4.ht2, .5.ht2, .6.ht2, .7.ht2, and .8.ht2. In the case of a large index these suffixes will have a ht2l termination. These files together constitute the index: they are all that is needed to align reads to that reference. The original sequence FASTA files are no longer used by HISAT2 once the index is built.

Helper: https://open.bioqueue.org/home/knowledge/showKnowledge/sig/hisat2-build

## 2-2) Filtering the .sam file

Samtools website: https://www.htslib.org/

Samtools is a set of programs for interacting with high-throughput sequencing data. It is helpful for converting SAM, BAM and CRAM files. One of the most used commands is the “samtools view,” which takes .BAM/.SAM files as input and converts them to .SAM/.BAM, respectively. The “-S” and “-b” commands are used. The alignment of fastq files occurs in random order with the position in the reference genome. Therefore, in “samtools sort,” the BAM files sorting is performed. The “-o” command indicates the output file.

In the sam file second column we have flags and we have to convert these number ro binary format and check if it consists of 4, remove that read

```linux
samtools view -b -F 4 .sam > mapped.bam
```

* Input => .sam file (it consist of unmapped data)
* Output => .bam file which consists of only mapped data

| Flag        | Chr           | Description  |
| ------------- |:-------------:| :----- |
| 0x0001      | p       | the read is paired in sequencing |
| 0x0002      | P       | the read is mapped in a proper pair |
| 0x0004      | u       | the query sequence itself is unmapped |
| 0x0008      | U       | the mate is unmapped |
| 0x0010      | r       | strand of the query (1 for reverse) |
| 0x0020      | R      | strand of the mate |
| 0x0040      | 1       | the read is the first read in a pair |
| 0x0080      | 2       | the read is the second read in a pair |
| 0x0100      | s       | the alignment is not primary |
| 0x0200      | f       | the read fails platform/vendor quality checks |
| 0x0400      | d       | the read is either a PCR or an optical duplicate |

Helper: https://www.biostars.org/p/56246/

Other information:

IGV (Integrative Genomics Viewer) - BAM (Binary Alignment Map)

A BAM file is the binary version of a SAM file and is used to assemble aligned reads into transcripts using StringTie and is also the preferred file format for viewing in IGV. 

Sorting:
```linux
samtools sort -@ 1 -o map/ERR188044_chrX.bam map/ERR188044_chrX.sam
```

## 2-3) Counting by htseq-count

https://htseq.readthedocs.io/en/master/install.html

Given a file with aligned sequencing reads and a list of genomic features, htseq-count can be used to count how many reads map to each feature.

**Requires .gtf file (Genome Transfer Format)**

```linux
htseq-count --format bam sorted_alignment_file.bam genome_annotation > .count
```

* Input => .bam file , gtf file
* Output => .count file

Advanced:

```linux
for i in *.sam; do echo $i; htseq-count $i ../genes/chrX.gtf > ../countsOutput/${i/sam/count}; done
```

# Step 3 (DESeq2 in R)

## 3-1) Install and load DESeq2 in R

```R
install.packages("BiocManager")
require(DESeq2)
```

## 3-1) Preparing Data

Sample of Count file:

| GeneID        | Count           |
| ------------- |:-------------:|
| NM_000032	| 0 |
| NM_000033	| 86 |
| NM_000044	| 1 |
| NM_000047	| 0 |
| ... | ... |

```R
Counts = read.table('ERR188044.count', sep = '\t')[,2]
```

* Output => retrieve only second column values

| Count           |
|:-------------:|
| 0 |
| 86 |
| 1 |
| 0 |
| ... |

## 3-1-2) Merging all Count files

```R
for (i in dir("/project-result/3-countsOutput/", full.names = T)[-1]) {
  Counts = cbind((Counts), read.table(i)[,2])
}
```

```R
rownames(Counts) = read.table('/3-countsOutput/ERR188044.count', sep = '\t')[,1]
```

```R
Counts = Counts[-grep("__",rownames(Counts)),] 
```

* Output => removing the lines in data which are meta data in starts with ‘__’

```R
CountsMatrix = as.matrix(Counts)
```

* Output => Convert Counts into the Matrix

Tip: as.matrix() mixed logical-integer will give a nteger matrix

```R
Class = read.csv("geuvadis_phenodata.csv")[,3]
```

* Output = >  [1] "YRI" "YRI" "YRI" "GBR" "GBR" "YRI" "GBR" "GBR" "GBR" "GBR" "YRI" "YRI"

```R
Class = data.frame(Population = as.factor(Class))
```

* Output =>
        * Population
    * 1         YRI
    * 2         YRI
    * 3         YRI
    * 4         GBR
    * 5         GBR
    * 6         YRI
    * 7         GBR
    * 8         GBR
    * 9         GBR
    * 10        GBR
    * 11        YRI
    * 12        YRI



