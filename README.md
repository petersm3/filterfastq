# Summary
Filter out sequences in Illumina [HiSeq](https://www.illumina.com/systems/sequencing-platforms/hiseq-3000-4000.html)/[MiSeq](https://www.illumina.com/systems/sequencing-platforms/miseq.html) gzipped [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files (single or paired-end) that match a reference genome, e.g., [PhiX](https://www.illumina.com/content/dam/illumina-support/documents/products/technotes/technote_phixcontrolv3.pdf) (Control Libraries)

# Disclaimer
* For removal of contamination it is recommended that you use other programs, e.g., JGI's [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) "Kmer filtering" feature as part of [Data Preprocessing](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/data-preprocessing/)
* <code>filterfastq.pl</code> was originally written as a proof of concept in 2015 (and used in production on a limited basis) due to Illumina's [elimination](https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote-hiseq-low-diversity.pdf) of the dedicated [PhiX](https://www.illumina.com/content/dam/illumina-support/documents/products/technotes/technote_phixcontrolv3.pdf) control lane (#8) on their [HiSeq](https://www.illumina.com/systems/sequencing-platforms/hiseq-3000-4000.html) 3000/4000 instruments
  * e.g., Sequencing a non-indexed sample (whole genome) with a large PhiX spike-in (e.g., 10%) may have required filtering of the final set of FASTQ reads
* <code>filterfastq.pl</code> uses command-line [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (<code>blastn</code>) to compare sequences in a serial fashion, which can take a long time
  * e.g., Approximately 11.75 days to process 264,724,686 single end reads, which had a recorded PhiX spike-in at 0.6%; 0.521% reads were filtered out by the script
# Dependencies
* Reference genome to filter against (in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format), e.g., [phiX174](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta) ([Phi X 174 bacteriophage](https://en.wikipedia.org/wiki/Phi_X_174))
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata) (<code>makeblastdb</code> and <code>blastn</code>)
* [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) (<code>fastq_to_fasta</code>)
* [GNU core utilities](https://www.gnu.org/software/coreutils/) (e.g., <code>split</code>)
* [GNU sed](https://www.gnu.org/software/sed/)

# Usage
```text
./filterfastq.pl -h

Usage: ./filterfastq.pl -r reference_genome.fa -i fastq_input_dir -o fastq_output_dir [-t threads] [-e evalue]

-r, --reference
       FASTA reference genome for sequences to be filtered against
-i, --input
       Directory containing the gzipped FASTQ file(s) to be filtered
       Sets of paired-end reads need to be in the same directory and are
       identified by their filenames containing either *_R1_* and *_R2_*
-o, --output
       Directory to store the filtered FASTQ and intermediary files
       Script will create output directory but not parents, e.g., mkdir -p
-t, --threads
       Optional number of threads for blastn to use (default is 1)
-e, --evalue
       Optional evalue for blastn to use (default is 1e-15)
-q, --quiet
       Do not write script progress to standard output
-h, --help
       This usage information.
```
