# RecoveryUnclassifiedReads

The package was specifically created to work with Oxford Nanopore Technologies (ONT) 16S sequencing projects multiplexed samples. Prior to running the scripts, basic data analysis such as basecalling and demultiplexing is required. 

The scripts work with the unclassified.fastq file from SQK-16S024 sequencing experiments in order to associate additional reads to each barcoded sample. The barcodes sequences are provided by ONT at https://community.nanoporetech.com/technical_documents/chemistry-technical-document/v/chtd_500_v1_revaq_07jul2016/barcode-sequences. 

All files corresponding to each barcode and the unclassified.fastq file should be placed within a directory. The downloaded package should be placed in a dedicated directory along the .fastq files. Simply enter the package directory and run the ./script_16S_unclassified_reads.sh script. Additionally, SeqKit, ncbi-blast+, NanoPlot and several R libraries should be pre-installed.
