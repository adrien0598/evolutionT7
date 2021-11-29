#!/bin/bash

#########################################################
##### script for analysis of raw ngs output ############
#########################################################

# step 1
# quality analysis and auto trimming 
# using fastp

# step 2
# from raw paried reads to bam file
# mapping with bwa (better than minimap 2 for our read lenght)

# step 3
# from bam file to sorted and filtered bam file
# using samtools

# step 4
# variant calling from sorted bam file
# using varscan with --min-coverage 1 --min-reads2 1 --min-avg-qual 35 --min-var-freq 0.0000000001 --p-value 0.99


cd X204SC21093372-Z01-F003/raw_data # F003 is the second batch we send to NGS (7 samples)
for i in D*; do
    name=$i.txt
    cd $i
    for j1 in *1.fq.gz; do
        doc1=$j1
    done
    for j2 in *2.fq.gz; do
        doc2=$j2
    done
    cd ../../..
    doc1=X204SC21093372-Z01-F003/raw_data/$i/$doc1
    doc2=X204SC21093372-Z01-F003/raw_data/$i/$doc2
    
    #fastp
    fastp -i $doc1 -I $doc2 -o X204SC21093372-Z01-F003/fastp_data/$j1 -O X204SC21093372-Z01-F003/fastp_data/$j2 --verbose
    
    #alignment BWA-MEM
    ./bwa/bwa index ref.fa
    ./bwa/bwa mem ref.fa X204SC21093372-Z01-F003/fastp_data/$j1 X204SC21093372-Z01-F003/fastp_data/$j2 | \
    samtools fixmate -O bam,level=0 -m - - | \
    samtools sort -l 1 -@8 -o sort.bam -T /tmp/example_prefix
    samtools markdup -O bam,level=1 sort.bam BAM_BWA/$i.bam
    samtools mpileup -f ref.fa BAM_BWA/$i.bam >myPileup.pileup
    
    #Variant call
    java -jar VarScan.v2.3.9.jar pileup2snp myPileup.pileup --min-coverage 1 --min-reads2 1 --min-avg-qual 35 --min-var-freq 0.0000000001 --p-value 0.99 >variant/$name
    
    cd X204SC21093372-Z01-F003/raw_data
done

cd X204SC21093372-Z01-F001/raw_data # F003 is the first batch we send to NGS (18 samples)
for i in A*; do
    name=$i.txt
    cd $i
    for j1 in *1.fq.gz; do
        doc1=$j1
    done
    for j2 in *2.fq.gz; do
        doc2=$j2
    done
    cd ../../..
    doc1=X204SC21093372-Z01-F001/raw_data/$i/$doc1
    doc2=X204SC21093372-Z01-F001/raw_data/$i/$doc2
    
    #fastp
    fastp -i $doc1 -I $doc2 -o X204SC21093372-Z01-F001/fastp_data/$j1 -O X204SC21093372-Z01-F001/fastp_data/$j2 --verbose
    
    #alignment BWA-MEM
    ./bwa/bwa index ref.fa
    ./bwa/bwa mem ref.fa X204SC21093372-Z01-F001/fastp_data/$j1 X204SC21093372-Z01-F001/fastp_data/$j2 | \
    samtools fixmate -O bam,level=0 -m - - | \
    samtools sort -l 1 -@8 -o sort.bam -T /tmp/example_prefix
    samtools markdup -O bam,level=1 sort.bam BAM_BWA/$i.bam
    samtools mpileup -f ref.fa BAM_BWA/$i.bam >myPileup.pileup
    
    #Variant call
    java -jar VarScan.v2.3.9.jar pileup2snp myPileup.pileup --min-coverage 1 --min-reads2 1 --min-avg-qual 35 --min-var-freq 0.0000000001 --p-value 0.99 >variant/$name
    
    cd X204SC21093372-Z01-F001/raw_data
done
