# L1PD (LINE-1 Pattern Detection)

L1PD is an implementation of an algorithm designed to efficiently detect human LINE-1s in subject genomes.  The algorithm uses fixed pre-determined probes that were generated using the GRCh38 reference as well as the LINE-1 database L1Base2 (http://l1base.charite.de/l1base.php).

L1PD may be executed in one of three modes:
* Genome mode
* BAM/CRAM mode
* FASTQ mode

BAM/CRAM mode automatically invokes Genome mode, and FASTQ mode automatically invokes BAM/CRAM mode, as shown in the following flowchart:

![L1PD flowchart](https://user-images.githubusercontent.com/14218905/137537174-ea06cc63-2e69-4f61-bc0e-bd5bfcc7cce3.jpg)

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
