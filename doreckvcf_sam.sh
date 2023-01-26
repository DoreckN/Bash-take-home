#!/bin/bash

# Download the files (sample.vcf.gz and sample.sam)from
#https://drive.google.com/drive/folders/11UD52i99CaCSBEJFNb8Y1afo9p3hL8cL?usp=sharing

# Manipulating vcf files

#gunzip the vcf file
#   gunzip sample.vcf.gz 

# 3. How many samples are in the file
    bcftools query -l sample.vcf | wc -l
    
# 4. How many variants are in the file
    bcftools query -f '%ALT\n' sample.vcf | wc -l

# 5. How would you extract the chromosome, position, QualByDepth and RMSMappingQuality fields? 
#Save the output to a tab-delimited file
    bcftools query -f '%CHROM\t%POS\t%INFO/QD\t%INFO/MQ\n' sample.vcf > fields_output.txt

# 6. Extract data that belongs to chromosomes 2,4 and MT
    awk '$1 ~ /^(2|4|MT)$/ {print $0}' sample.vcf > extracted_data.vcf

# 7. Print out variants that do not belong to chr20:1-30000000
    awk '$1 != "20" || ($1 == "chr20" && ($2 < 1 || $2 > 30000000)) \
        {print $1, $2, $4, $5}' sample.vcf > var20_variants.vcf

#8. Extract variants that belong to SRR13107019
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -s SRR13107019 sample.vcf > output19.txt 

#9. Filter out variants with a QualByDepth above 7
    awk -F '\t' '{if ($6>=7) print $0}' sample.vcf > minQD7.vcf 

#10. How many contigs are referred to in the file. Check the header section
    grep -c "^##contig" sample.vcf # answer is 2211

#12. Extract data on the read depth of called variants for sample SRR13107018
    bcftools query -f '%DP\n' -s SRR13107018 sample.vcf > read_depth18.vcf

#13. Extract data on the allele frequency of alternate alleles. Combine this data with the
#chromosome and position of the alternate allele
    bcftools query -f '%CHROM\t%POS\t%AF\n' sample.vcf > allele_freq.vcf


# Manipulating sam files

#3. How many samples are in the file
    samtools view -H sample.sam | grep -c '@RG' # answer is 249

#4. How many alignments are in the file
    samtools view -c -F 4 sample.sam # answer is 35511

#5. Get summary statistics for the alignments in the file
    samtools flagstat sample.sam > sam_stat.txt

#6. Count the number of fields in the file
    awk '{print NF}' sample.sam

#7. Print all lines in the file that have @SQ and sequence name tag beginning with NT_
    grep '@SQ.*NT_' sample.sam 

#8. Print all lines in the file that have @RG and LB tag beginning with Solexa
    grep "@RG.*LB:Solexa" sample.sam

#9.Extract primarily aligned sequences and save them in another file
    samtools view -f 2 -b sample.sam > primary_aligned.sam

#10. Extract alignments that map to chromosomes 1 and 3. Save the output in BAM
#format 
    awk '$1 !~/^@/ && ($3 == "1" || $3 == "3")' sample.sam | samtools view -Sb sample.sam\
         > align1_3.bam

#11. How would you obtain unmapped reads from the file
    samtools view -f 4 sample.sam > unmapped_reads.sam

#12. How many reads are aligned to chromosome 4
    grep -c "^4\t" sample.sam # there are no aligned reads

#14. Extract all optional fields of the file and save them in “optional_fields.txt”
    awk '{for(i=11;i<=NF;i++) print $i}' sample.sam > optional_fields.txt

