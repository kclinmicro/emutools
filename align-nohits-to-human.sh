#!/bin/bash
sampath=$1
# Change the below to "sci run" after installing SciCommander
# (https://github.com/samuell/scicommander) in order to avoid re-running
# completed steps, and tracking lineage.
executor="bash -c"
#executor="sci run"

if [[ -z ${sampath} ]]; then
    echo "Usage: align-nohits-to-human.sh <.sam-file>";
    exit 1
fi

samfile=$(basename ${sampath})
fqfile=${samfile%.sam}.unalign.fq
reffile=GCF_009914755.1_T2T-CHM13v2.0_genomic.fna

echo "--------------------------------------------------------------------------------";
echo "-> Extracting unaligned sequences from ${sampath} ..."
echo "--------------------------------------------------------------------------------";
${executor} "samtools fastq -f 4 ${sampath} > ${fqfile}"

echo "--------------------------------------------------------------------------------";
echo "-> Download human genome ..."
echo "--------------------------------------------------------------------------------";
${executor} "curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/${reffile}.gz > ${reffile}.gz"

echo "--------------------------------------------------------------------------------";
echo "-> Unpack genome ..."
echo "--------------------------------------------------------------------------------";
${executor} "zcat ${reffile}.gz > ${reffile}"

echo "--------------------------------------------------------------------------------";
echo "-> Index reference genome ..."
echo "--------------------------------------------------------------------------------";
${executor} "samtools faidx ${reffile} --fai-idx ${reffile}.fai"

echo "--------------------------------------------------------------------------------";
echo "-> Aligning nohits sequences in ${fqfile} to human genome ..."
echo "--------------------------------------------------------------------------------";
humalnsam=${fqfile%.fq}.aln_human.sam
${executor} "minimap2 -ax lr:hq ${reffile} ${fqfile} > ${humalnsam}"

echo "--------------------------------------------------------------------------------";
echo "-> Converting to bam ..."
echo "--------------------------------------------------------------------------------";
humalnbam=${fqfile%.fq}.aln_human.bam
${executor} "samtools view -b ${humalnsam} > ${humalnbam}"

echo "--------------------------------------------------------------------------------";
echo "-> Sorting bam ..."
echo "--------------------------------------------------------------------------------";
humalnbamsrt=${fqfile%.fq}.aln_human.sorted.bam
${executor} "samtools sort ${humalnbam} > ${humalnbamsrt}"

echo "--------------------------------------------------------------------------------";
echo "-> Indexing bam ..."
echo "--------------------------------------------------------------------------------";
humalnbamsrt=${fqfile%.fq}.aln_human.sorted.bam
${executor} "samtools index ${humalnbamsrt} -o ${humalnbamsrt}.bai"

echo "--------------------------------------------------------------------------------";
echo "-> Determining most common chromosome location ..."
echo "--------------------------------------------------------------------------------";

# Determine the most commonly occuring chromosome / location
genomeloc=$(samtools view ${humalnbamsrt} | awk -F"\t" '{ print $3 ":" $4 }' | sort | uniq -c | sort -nr | head -n 1 | awk '{ print $2 }')
echo "Most common location: ${genomeloc}"

# Format it as a chromosome region for IGV
genomereg=$(echo ${genomeloc} | awk -F: '{ print $1 ":" $2-100 "-" $2 + 1500 }')

echo "--------------------------------------------------------------------------------";
echo "-> Plotting alignment..."
echo "--------------------------------------------------------------------------------";
igvscript=${humalnbamsrt%.bam}.igv
plotfile=${humalnbamsrt%.bam}.igv.png
cat << END > ${igvscript}
new genome ${reffile}
load ${humalnbamsrt}
goto ${genomereg}
snapshot ${plotfile}
END

${executor} "igv -b ${igvscript} # outfile: ${plotfile}"
