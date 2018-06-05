#from: https://legacy.gitbook.com/book/ycl6/methylation-sequencing-analysis/details
#updated with more recent hg19 files
#NB running on OSX so gsort -V used, else use a sort -V or some chr aware sorting method

##here we try to create EPIC array annotation-style BED format file

##format we want to annotate with is:
# chr, probe_start, probe_end, 'anno'

##'anno' = probe_id;gene_name;gene_ID;transcript_id;anno_location;relation_to_island;rmsk

##anno_location is specified as:
#exon, intron, 3_UTR, 5_UTR, TSS1500, TSS200
#this is based on EPIC array manifest "UCSC_RefGene_Group"

##also included a repetitive element not present in EPIC (rmsk)

##staging input:
SCRIPTDIR=$(dirname $0)
CURDIR=${SCRIPTDIR}/../data
BEDDIR=${SCRIPTDIR}/../bed
if [[ ! -d $CURDIR ]];then
  echo "$CURDIR does not exist, creating it"
  mkdir -p $CURDIR
  if [[ ! -d $BEDDIR ]];then
    mkdir -p $BEDDIR
  fi
fi
cd $CURDIR

####################
## DATA DOWNLOADS ##
####################
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip
wget http://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-v1-0-missing-legacy-cpg-b3-vs-b2-annotations.zip

wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeRegTfbsClustered.txt.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz

##(g)unzip those files
unzip infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip
unzip infinium-methylationepic-v1-0-missing-legacy-cpg-b3-vs-b2-annotations.zip

gunzip gencode.v27.annotation.gtf.gz chromInfo.txt.gz cpgIslandExt.txt.gz wgEncodeRegTfbsClustered.txt.gz

#########################################
##create annotation files in bed format##
#########################################
##Infinium EPIC bed:

##making MethylationEPIC_v-1-0_B4.csv into MethylationEPIC_v-1-0_B4.bed
##want 4 ';' delim fields, along with 1,2,3 of chr, start, end
##anno delim ';' -> probe;gene;transcript;feature_type
##NB if your chips start 2011...something... then you don't want to use the Legacy B2 Missing probes
##mine start 2007...something... so I do use it
cp MethylationEPIC_v-1-0_B4.csv MethylationEPIC_v-1-0_B4.Legacy_B2.csv
tail -n+8 "MethylationEPIC Missing Legacy CpG (v1.0_B3 vs. v1.0_B2) Annotations.csv" >> MethylationEPIC_v-1-0_B4.Legacy_B2.csv
tail -n+8 MethylationEPIC_v-1-0_B4.Legacy_B2.csv | perl -ane '@s=split(/\,/);
  $cn=$s[11]; $cn=~s/chr//;
  $st=$s[12]; $en=$st; $en++;
  $str=$s[14];
  if($str eq "R"){$stro="-"} if($str eq "F"){$stro="+"}
  if(($str ne "F") && ($str ne "R")){$stro="*"}
  $gn=$s[15];
  $tn=$s[16];
  $bd=$s[17];
  $in=$s[19];
  if($in eq ""){$in="Open_Sea";}
  $gn=~s/\;/\,/g; $tn=~s/\;/\,/g; $bd=~s/\;/\,/g; $in=~s/\;/\,/g;
  if($cn ne "CHR"){print "$cn\t$s[12]\t$en\t$stro\t$s[0];$gn;$tn;$bd;$in\n";}' | gsort -V | perl -ane 'if(scalar(@F)==5){print $_;}' > MethylationEPIC_v-1-0_B4.Legacy_B2.bed

##repetitive elements (see: https://bit.ly/2H819it)
##NB that this includes an annotated "ELEMENT-TYPE_<line-number>" for traceback
echo "Making Repetitive Elements BED"
gunzip -c rmsk.txt.gz | perl -ane 'chomp; print "$F[5]\t$F[6]\t$F[7]\t$F[11]|$F[9]|$F[12]|$F[10]|$F[11]_$.\n";' | grep -v "?" | gzip > hg19.rmskRM327.NR.bed.gz

##CpG islands: naming is CpG name and length separated by underscore (see: https://bit.ly/2kDogbQ)
##make shore, shelf annotations:
##shore is 2kb+(south, downstream),-(north, upstream), shelves another 2kb out from those:
echo "Making CpG BED"
perl -ane 'chomp; print "$F[1]\t$F[2]\t$F[3]\tCpG$F[5]_$F[6]\n";' cpgIslandExt.txt | gsort -V | perl -ane 'if($F[0]=~m/\_/){next;} chomp; $sru=$F[1]-2000; $sfu=$sru-2000; $srd=$F[2]+2001; $sfd=$srd+2000; $sfue=$sru-1; $sfde=$srd+1; $srue=$F[1]-1; $srde=$F[2]+1; $nsf="$F[0]\t$sfu\t$sfue\t$F[3]_N_Shelf"; $nsr="$F[0]\t$sru\t$srue\t$F[3]_N_Shore"; $isl="$F[0]\t$F[1]\t$F[2]\t$F[3]_Island"; $ssr="$F[0]\t$srde\t$srd\t$F[3]_S_Shore"; $ssf="$F[0]\t$sfde\t$sfd\t$F[3]_S_Shelf";  print "$nsf\n$nsr\n$isl\n$ssr\n$ssf\n";' | gsort -V | bedtools merge -i - -c 4 -o distinct | gzip > cpgIslandExt.relation.bed.gz

# ##TFBS: naming is name, filter on scores above 500 (max. 1000, see https://bit.ly/2HbglLJ)
# perl -ane 'if($F[5] >= 500){chomp;print "$F[1]\t$F[2]\t$F[3]\t$F[4]\n";}' wgEncodeRegTfbsClustered.txt | gsort -V | gzip > wgEncodeRegTfbsClustered.bed.gz

##GTF reformat (see: https://bit.ly/2LdHmkk, NB we don't use that version, but same schema)
##exon: naming is gene_name, gene_id, gene_type, strand
echo "Making Exon BED"
perl -ane 'if($F[2] eq "exon"){chomp; $o="$F[0]\t$F[3]\t$F[4]\texon|$F[15]|$F[9]|$F[13]|$F[6]\n"; $o=~s/\"//g; $o=~s/\;//g; print $o;}' gencode.v27.annotation.gtf | gsort -V | bedtools merge -i - -c 4 -o distinct | gzip > gencode.v27.annotation.exon_merged.bed.gz

##intron: naming as exon
echo "Making Intron BED"
perl -ane 'if($F[2] eq "gene"){chomp; $o="$F[0]\t$F[3]\t$F[4]\tintron|$F[13]|$F[9]|$F[11]|$F[6]\n"; $o=~s/\"//g; $o=~s/\;//g; print $o;}' gencode.v27.annotation.gtf | gsort -V | bedtools subtract -a - -b gencode.v27.annotation.exon_merged.bed.gz | gsort -V | bedtools merge -i - -c 4 -o distinct | gzip > gencode.v27.annotation.intron_merged.bed.gz

##script from: https://davetang.org/muse/2013/01/18/defining-genomic-regions/ Defining the 5′ and 3′ UTRs section; I put it into define_UTR.pl, and changed output, delimiters to match above
echo "Making UTR BED"
perl ${SCRIPTDIR}/define_UTR.pl gencode.v27.annotation.gtf | gsort -V | gsort -V | bedtools merge -i - -c 4 -o distinct | gzip > gencode.v27.annotation.UTR_merged.bed.gz

# ##not using intergenic, exon, intron UTR enough; included for completeness
# ##chromInfo.txt throws error about no valid entries, so use only first 2 cols:
# cut -f 1,2 chromInfo.txt | grep -v "_" | gsort -V > chromInfo.12.txt
#
# # ##intergenic: naming as exon
# # perl -ane 'if($F[2] eq "gene"){print "$F[0]\t$F[3]\t$F[4]\n";}' gencode.v27.annotation.gtf | gsort -V | bedtools complement -i - -g chromInfo.12.txt | gzip > gencode.v27.annotation.intergenic.bed.gz
echo "Making Transcript BED"
perl -ane 'if($F[2] eq "transcript"){chomp; $o="$F[0]\t$F[3]\t$F[4]\t$F[15]|$F[9]|$F[11]|$F[6]\n"; $o=~s/\"//g; $o=~s/\;//g; print $o;}' gencode.v27.annotation.gtf | gsort -V | bedtools merge -i - -c 4 -o distinct | gzip > gencode.v27.annotation.transcript.bed.gz

echo "Making TSS BEDs"
gunzip -c gencode.v27.annotation.transcript.bed.gz | perl -ane 'if($F[3]=~m/\|+/){
    if($F[3]=~m/\|-/){next;}
    $l=$F[1]-1000; $r=$F[1]+500;
  }
  if($F[3]=~m/\|-/){
    if($F[3]=~m/\|+/){next;}
    $l=$F[1]-500; $r=$F[1]+1000;
  }
  if($l < 0){
    $l=0;
  }
  print "$F[0]\t$l\t$r\t$F[3]\n";' | gsort -V | bedtools merge -i - -c 4 -o distinct | gzip > gencode.v27.annotation.TSS1500_merged.bed.gz

gunzip -c gencode.v27.annotation.transcript.bed.gz | perl -ane 'if($F[3]=~m/\|+/){
    if($F[3]=~m/\|-/){next;}
    $l=$F[1]-100; $r=$F[1]+100;
  }
  if($F[3]=~m/\|-/){
    if($F[3]=~m/\|+/){next;}
    $l=$F[1]-100; $r=$F[1]+100;
  }
  if($l < 0){
    $l=0;
  }
  print "$F[0]\t$l\t$r\t$F[3]\n";' | gsort -V | bedtools merge -i - -c 4 -o distinct | gzip > gencode.v27.annotation.TSS200_merged.bed.gz

##put into single bed:
echo "Moving BEDs to BED dir"
for x in $(ls *bed.gz);do
  nm=$(echo $x | sed 's/\.gz//'); echo $nm; gunzip -c $x > $BEDDIR/$nm;
done
