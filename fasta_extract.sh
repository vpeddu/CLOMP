#!/bin/bash

#HPIV3 reads
cat 5780_S28_L001_R1_001_assignments.txt | grep "\t11216" | cut -f1,3 | tr "M0" ">M0" | tr "\t" "\n" > hpiv3.fasta
#Rhinovirus A reads
cat SC5683_S2_L001_R1_001_assignments.txt | grep "\t147711" | cut -f1,3 | tr "M0" ">M0" | tr "\t" "\n" > RhVA.fasta
#Rhinovirus C reads
cat SC5698_S4_L001_R1_001_assignments.txt | grep "\t463676" | cut -f1,3 | tr "M0" ">M0" | tr "\t" "\n" > RhVC.fasta
#moraxella reads 
cat SC5688_S3_L001_R1_001_assignments.txt| grep "\t480" | cut -f1,3 | tr "M0" ">M0" | tr "\t" "\n" > MC_WA6_UW3.fasta
cat sc5727_S26_L001_R1_001_assignments.txt| grep "\t480" | cut -f1,3 | tr "M0" ">M0" | tr "\t" "\n" > MC_WA9_UW6.fasta
cat sc5746_S25_L001_R1_001_assignments.txt | grep "\t480" | cut -f1,3 | tr "M0" ">M0" | tr "\t" "\n" > MC_WA8_UW5.fasta

