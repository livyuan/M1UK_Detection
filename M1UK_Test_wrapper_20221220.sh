source /etc/profile.d/modules.sh
bindir="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/share/CIDT/bin/"

## STEP 0 Setting up batch working directory

while getopts :i:o: option
do
    case $option in
        i) SAMPLE_LIST=$OPTARG;;
        o) WKDIR=$OPTARG;;
    esac
done

  #SAMPLE_LIST is a text file with 3 columns seperated by comma; No header;
  # Column 1 : LABID; Column 2: path to fastq 1; Column 3: path to fastq 2
  # No space or comma in LABID or path
  # Each LABID should be unique

  #WKDIR is a directory where output files will be saved

echo "Input file is " $SAMPLE_LIST
echo "Output directory is " $WKDIR

mkdir -p $WKDIR
mkdir -p $WKDIR"/qsub_log/"

## STEP 1 Prepare sample list and ref genome 
cd $WKDIR
dos2unix $SAMPLE_LIST
cat $SAMPLE_LIST | sed 's/,/ /g'> ID_fqdir.txt
dos2unix ID_fqdir.txt

module purge
module load bowtie2/2.1.0
module load samtools/1.8
module load minimap2/2.4

#M1ref_5005="/scicomp/home/yqh8/download/GAS_Complete_Genomes_57/\
#GCA_000011765.2_ASM1176v2/GCA_000011765.2_ASM1176v2_genomic.fna"

M1ref_5005="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/share/CIDT/database/\
GCA_000011765.2_ASM1176v2_genomic.fna"

echo ">MGAS5005" > ref.fa
grep -v ">" "$M1ref_5005" >> ref.fa
ref1="$WKDIR""/ref.fa"
rm -f $ref1.*
bowtie2-build $ref1 $ref1
samtools faidx $ref1

chr=$(cat $ref1 |awk 'NR==1'| sed 's/>//')
cat $ref1 |awk 'NR>1' | sed 's/\(.\)/\1 \n/g' | sed '/^$/d' > temp1.txt
n1=$(cat "temp1.txt" | wc -l)
echo "$(seq 1 $n1)" | awk '{print chr" "$1}' chr=$chr > temp2.txt
paste temp2.txt temp1.txt > temp3.txt
cat temp3.txt | awk '{print $1"\t"$2"\t.\t"$3"\t.\t.\t.\t.\t."}' > temp4.txt
printf "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n" > refgenome.vcf
cat temp4.txt >>  refgenome.vcf
rm -f temp*
grep -v "#" refgenome.vcf | awk '{print $2".5 "$4}' > refgenome2.vcf

ref1_idx=$WKDIR"/ref.mmi"
minimap2 -d $ref1_idx  $ref1

ref1_plo=$WKDIR"/ref.plo"
echo "* * * * 1" > $ref1_plo

##STEP 2: run individual typer 
cd $WKDIR

n1=$(cat ID_fqdir.txt | wc -l)
qsub_file="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/share/CIDT/scripts/\
M1UK_Test_Typer_20221220.sh"

qsub -sync y  -t 1-$n1:1 \
  -l h_vmem=8G,h_rt=71:00:00 -pe smp 2 \
  -e $WKDIR"/qsub_log/" -o $WKDIR"/qsub_log/" \
  -v WKDIR=$WKDIR \
  $qsub_file

## STEP 3: compile output

cd $WKDIR
#out_file=Result_Table_$(date '+%Y%m%d').txt
#echo "LABID,SNP_Match_M1UK,SNP_Not_Match_M1UK,SNP_Total" > $out_file
#cat $(find "./" | grep "_Match_M1UK\.txt") | sed 's/ /,/g'>> $out_file

out_file=M1UK_Result_Table_$(date '+%Y%m%d').csv
echo "LABID,SNP_Match_M1UK,M1UK_Result" > $out_file
cat $(find "./" | grep "_Match_M1UK\.txt") | awk '{Res="Negative";\
if ($2==27) {Res="Positive"}; print $1","$2","Res}' >>$out_file

grep "ERROR" ./qsub_log/* > ERROR_Log.txt
grep "fatal" ./qsub_log/* >> ERROR_Log.txt

