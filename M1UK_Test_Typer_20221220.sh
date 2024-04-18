source /etc/profile.d/modules.sh
bindir="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/share/CIDT/bin/"

PATH=$PATH:$HOME/.local/bin:$HOME/bin:/opt/sge/bin/lx-amd64
export PATH="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/JanOw_Dependencies:$PATH"
export PERL5LIB=/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/JanOw_Dependencies/perl_libs

module purge
module load perl/5.12.3
module load Python/2.7
module load bowtie2/2.1.0
module load samtools/1.8
module load freebayes/1.1.0
module load tabix/0.2.6
module load minimap2/2.4
module load snp-sites/2.3.2
module load prokka/1.8
module load BEDTools/2.27.1

i=$SGE_TASK_ID

SAMPLE_LINE=$(cat $WKDIR"/ID_fqdir.txt" | awk 'NR==var' var=$i | sed 's/"//g')
  # $WKDIR is an env variable set by the qsub call
LABID=$(echo $SAMPLE_LINE | awk '{print $1}')
fq1=$(echo $SAMPLE_LINE | awk '{print $2}')
fq2=$(echo $SAMPLE_LINE | awk '{print $3}')

ref1=$WKDIR"/ref.fa"
ref1_idx=$WKDIR"/ref.mmi"
ref1_plo=$WKDIR"/ref.plo"

WK_DIR=$WKDIR"/$LABID/"
mkdir -p $WK_DIR
cd $WK_DIR

echo "LABID is " $LABID
echo "fastq files are in directory:" $fq1 $fq2
echo "Output directory is " $WK_DIR

## Start doing work

minimap2 -ax sr  $ref1_idx "$fq1" "$fq2"  > test0.sam

samtools view -bu -S test0.sam -T $ref1  -o test0.bam
  samtools sort test0.bam -o test0.sorted.bam
  samtools index test0.sorted.bam


  module load bcftools/1.10.2
  bcftools mpileup -I -f $ref1 -O b test0.sorted.bam > test0.bcf
  bcftools call -cv  --ploidy-file $ref1_plo -O v test0.bcf > test0.var.vcf

  samtools depth -aa test0.sorted.bam > temp1.txt
  cat temp1.txt | awk '{i=i+$3};END {print i/NR}' > Mean_Depth.txt

  cat temp1.txt | awk '$3<=3 {print $2" N"}' > out2.txt
  cat test0.var.vcf |  awk '$6 >=20' | grep "AF1=1;" | \
    sed 's/DP=/DP=\t/' | sed 's/;VDB/\t;VDB/' | awk '$9 >=10' | \
    awk 'length($5)==1' | awk '$1!~"#" {print $2" "$5}' >> out2.txt

  cat $WKDIR"/refgenome2.vcf"  "out2.txt" | \
    sort -nk1,1 | sed 's/\..* / /'  |  awk  '!_[$1]++'  > temp2.txt
  echo ">""$LABID" > "$LABID"".fasta"
  echo $(awk '{print $2}' temp2.txt) | sed 's/ //g' | fold >> "$LABID"".fasta"

  cat test0.var.vcf > "$LABID".var.vcf

  rm -f *.sam *.bam test* temp* out*


f2="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/share/YL_yqh8/M1UK/\
M1UK_SNP27_MGAS5005POS_REF_ALT.txt"

cat $f2 | awk '{print LABID"\t"$1-1"\t"$1}' LABID="$LABID" > SNP27.bed
bedtools getfasta -fi "$LABID"".fasta" -bed SNP27.bed -tab | \
  sed 's/:.*-/\t/' > "$LABID"".SNP27.txt"

n1=$(paste $f2 "$LABID"".SNP27.txt" | awk '$3==$6' | wc -l)
n2=$(paste $f2 "$LABID"".SNP27.txt" | awk '$2==$6' | wc -l)
n3=$(paste $f2 "$LABID"".SNP27.txt" | awk '$1==$5' | wc -l)
echo "$LABID" "$n1" "$n2" "$n3"> "$LABID""_Match_M1UK.txt"

rm -rf  SNP27.bed "$LABID".fasta.fai

module purge
