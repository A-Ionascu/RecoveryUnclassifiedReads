#!/bin/bash

CWD="$(pwd)"


echo "Preparing unclassified reads ..."
cd ..
mkdir unclassified
mv unclassified.fastq unclassified
cd unclassified
seqkit fq2fa unclassified.fastq > unclassified.fasta
cd $CWD
cp ONT_16S_all_barcodes.fasta ../unclassified

echo "Performing BLAST ..." 
cd ../unclassified
blastn -subject unclassified.fasta -query ONT_16S_all_barcodes.fasta -word_size 24 -outfmt 6 > blast_all_barcodes.csv

echo "Filtering BLAST results ..."
cd $CWD
cp select_reads_part1.R ..
cd ..
Rscript select_reads_part1.R 2> /dev/null
rm select_reads_part1.R
# mv select_reads_part1.R $CWD/pachet_16S_unclassified_recovery

echo "Selecting unique reads ..."
cd unclassified
sed -i 's/\"//g' unique_all_barcodes.txt
grep -F -A 1 -f unique_all_barcodes.txt unclassified.fasta > reads_with_barcodes_tmp.fasta
grep -v -e '--' reads_with_barcodes_tmp.fasta > reads_with_barcodes.fasta
rm reads_with_barcodes_tmp.fasta
blastn -subject reads_with_barcodes.fasta -query ONT_16S_all_barcodes.fasta -word_size 24 -outfmt 6 > blast_unique_barcodes.csv

echo "Creating tables for each barcode ..."
cd $CWD
cp select_reads_part2.R ..
cd ..
Rscript select_reads_part2.R 2> /dev/null
rm select_reads_part2.R
# mv select_reads_part2.R $CWD/pachet_16S_unclassified_recovery

echo "Selecting reads for each barcode ..."
cd unclassified
sed -i 's/\"//g' *.lst

mkdir ../recovered_reads
mv *.lst ../recovered_reads
cd ../recovered_reads

for i in *.lst
do 
seqtk subseq ../unclassified/unclassified.fastq $i > ${i%.lst}".fastq"
done

echo "Creating concatenated reads ..."
mkdir ../concatenated
cp *.fastq ../concatenated
cd ..
cp *.fastq concatenated
cd concatenated


for i in SQK*.fastq 
do
for j in barcode*.fastq
do

if [[ $i == *$j* ]]
then
cat $i $j > "concatenated_"$j
fi

done
done

rm barcode*.fastq
rm SQK*.fastq

cat *.fastq > concatenated_all_barcodes.fastq

cd ../unclassified
cp unclassified.fastq ..


echo "Making NanoPlot analysis ..."
cd ..
mkdir nanoplot
cp concatenated/*.fastq nanoplot
cd nanoplot

for i in *.fastq
do
mkdir "nanoplot_${i%.fastq}"
cp $i "nanoplot_${i%.fastq}"
cd "nanoplot_${i%.fastq}"
NanoPlot --huge -t 6 --tsv_stats --only-report --fastq $i
awk -F " " '{print $0}' NanoStats.txt > NanoStats.csv
rm $i
cd ..
cp ../../pachet_16S_unclassified_recovery/nanoplot_analysis.R ./
done

Rscript nanoplot_analysis.R
rm nanoplot_analysis.R

echo "Preparing Epi2ME folders ..."
cd .. 
mkdir epi2me
cp concatenated/*.fastq epi2me
cd epi2me

for i in *.fastq
do
mkdir "epi2me_${i%.fastq}"
cp $i "epi2me_${i%.fastq}"
done

cd ..

echo "Finished !"







