# L1PD (LINE-1 Pattern Detection)

L1PD is an implementation of an algorithm designed to efficiently detect human LINE-1s in subject genomes.  The algorithm uses fixed pre-determined probes that were generated using the GRCh38 reference as well as the LINE-1 database L1Base2 (http://l1base.charite.de/l1base.php).

## Requirements

* Python 3
  * Matplotlib
  * Numpy
* [mrFAST](https://github.com/BilkentCompGen/mrfast/)

If starting from BAM/CRAM mode, [BCFtools](https://www.htslib.org/) 1.11 or newer is also required.

If starting from FASTQ mode, all requirements for Genome and BAM/CRAM modes must be met.  Additional requirements:
* [Samtools](https://www.htslib.org/) 1.11 or newer
* [BWA](http://bio-bwa.sourceforge.net/)
