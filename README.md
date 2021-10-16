# L1PD (LINE-1 Pattern Detection)

## Overview
L1PD is a tool for efficiently detecting human LINE-1s in subject genomes.  The underlying algorithm uses fixed pre-determined probes that were generated using the GRCh38 reference as well as the LINE-1 database [L1Base2](http://l1base.charite.de/l1base.php).

L1PD may be executed in one of three modes:
* Genome mode
* BAM/CRAM mode
* FASTQ mode

BAM/CRAM mode automatically invokes Genome mode, and FASTQ mode automatically invokes BAM/CRAM mode, as shown in the following flowchart:

![L1PD flowchart](https://user-images.githubusercontent.com/14218905/137537174-ea06cc63-2e69-4f61-bc0e-bd5bfcc7cce3.jpg)

## Modes
### Genome mode
The main purpose of L1PD is to detect human LINE-1s in subject genomes, which is precisely what Genome mode is for.  You provide a subject genome in FASTA format and L1PD will use its algorithm to detect LINE-1s in that genome, providing its output in GFF3 format.  L1PD also generates a histogram in PNG format of where patterns were detected, comparing the patterns from the subject genome with the patterns present in [GRCh38DH](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/) (the version of GRCh38 used by the 1000 Genomes Project to account for decoy sequences, alternative haplotypes and EBV). 

### BAM/CRAM mode
If you don't have the subject genome, but you have an alignment file in BAM or CRAM format as well as the corresponding reference genome (used to align and generate the BAM/CRAM file), then BAM/CRAM mode can use that alignment file and reference genome to generate the subject genome.  Once the subject genome has been generated, L1PD proceeds to use that genome in Genome mode.

Variant calling is performed using BCFtools; examples may be found at [Samtools - Workflows](https://www.htslib.org/workflow/).

### FASTQ mode
If you don't have the subject genome, or an alignment file, but you have paired-end reads in FASTQ format, then FASTQ mode will take the reads along with a reference genome and generate the subject genome.  This is accomplished by aligning the reads to the reference using BWA and processing the results with Samtools to generate a BAM file.  Once the BAM file has been generated, L1PD proceeds to use that alignment file and the reference in BAM/CRAM mode.

The pipeline used is similar to the one found at [Samtools - Workflows](https://www.htslib.org/workflow/).

## Requirements
### Genome mode
* Python 3
  * Matplotlib
  * Numpy
* [mrFAST](https://github.com/BilkentCompGen/mrfast/)
### BAM/CRAM mode
* All requirements of Genome mode
* [BCFtools](https://www.htslib.org/) 1.11 or newer
### FASTQ mode 
* All requirements of BAM/CRAM mode
* [Samtools](https://www.htslib.org/) 1.11 or newer
* [BWA](http://bio-bwa.sourceforge.net/)

## Usage
### Genome mode
<pre><code>L1PD.sh -g <i>genome.fa</i> [-e <i>edit_distance</i>] [-d <i>diff_threshold</i>] [-m <i>min_amt_probes</i>] [-p <i>prefix</i>] [-h]</code></pre>

### BAM/CRAM mode
<pre><code>L1PD.sh -b <i>alignment.bam/.cram</i> -r <i>reference.fa</i> [-t <i>num_threads</i>] [-e <i>edit_distance</i>] [-d <i>diff_threshold</i>] [-m <i>min_amt_probes</i>] [-p <i>prefix</i>] [-h]</code></pre>

### FASTQ mode
<pre><code>L1PD.sh -f <i>file1.fq</i> <i>file2.fq</i> -r <i>reference.fa</i> [-t <i>num_threads</i>] [-e <i>edit_distance</i>] [-d <i>diff_threshold</i>] [-m <i>min_amt_probes</i>] [-p <i>prefix</i>] [-h]</code></pre>

### Mode options
<dl>
   <dt>-g</dt><dd>Subject genome. [<em>Genome mode only</em>]</dd>
   <dt>-b</dt><dd>Alignment file in BAM/CRAM format. [<em>Genome mode only</em>]</dd>
   <dt>-r</dt><dd>Reference file in FASTA format. [<em>BAM/CRAM or FASTQ mode</em>]</dd>
   <dt>-f</dt><dd>Two paired-end read files in FASTQ format. [<em>FASTQ mode only</em>]</dd>
   <dt>-e</dt><dd>Maximum allowed edit distance when mapping probes to subject genome with mrFAST. [<em>Default: 20</em>]</dd>
   <dt>-d</dt><dd>Difference threshold allowed when comparing expected distance and actual distance between probes. [<em>Default: 700</em>]</dd>
   <dt>-m</dt><dd>Minimum amount of probes needed to establish a pattern. [<em>Default: 9</em>]</dd>
   <dt>-t</dt><dd>Additional number of threads to be used. [<em>Default: 0</em>]<br>Note that Genome mode ignores this argument, if specified.</dd>
</dl>
