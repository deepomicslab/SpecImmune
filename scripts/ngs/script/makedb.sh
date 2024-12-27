#!/bin/bash
bin=$(cd `dirname $0`; pwd)/../bin
db=$(cd `dirname $0`; pwd)/../db
license=$bin/novoalign.lic
if [ -f "$license" ];then
	$bin/novoindex  -k 14 -s 1 $db/CYP/ref/CYP.select.ndx $db/CYP/ref/CYP.select.fasta
	$bin/novoindex  -k 14 -s 1 $db/HLA/ref/HLA.extend.ndx $db/HLA/ref/HLA.extend.fasta
	$bin/novoindex  -k 14 -s 1 $db/HLA/ref/KIR.extend.select.ndx $db/HLA/ref/KIR.extend.select.fasta
fi
$bin/bowtie2-build $db/CYP/ref/CYP.select.fasta $db/CYP/ref/CYP.select.fasta
$bin/bowtie2-build $db/HLA/ref/HLA.extend.fasta $db/HLA/ref/HLA.extend.fasta
$bin/bowtie2-build $db/HLA/ref/KIR.extend.select.fasta $db/HLA/ref/KIR.extend.select.fasta



