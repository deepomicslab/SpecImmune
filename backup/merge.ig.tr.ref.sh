hg38=/home/wangmengyao/packages/kourami-0.9.6/resources/hs38DH.fa
chm13=/mnt/disk2_workspace/wangmengyao/chm13v2.0.fa

echo ">chr7_trb" >merge.fa
samtools faidx $hg38 chr7:142298000-142815000 |grep -v ">chr" >>merge.fa
echo ">chr9_trb" >> merge.fa
samtools faidx $hg38 chr9:33615000-33790000 |grep -v ">chr" >>merge.fa
echo ">chr14_tra" >> merge.fa
samtools faidx $hg38 chr14:21620000-22553000 |grep -v ">chr" >> merge.fa
echo ">chr7_trg" >> merge.fa
samtools faidx $hg38 chr7:38248000-38370000 |grep -v ">chr" >> merge.fa
#echo ">chr14_trd" >> merge.fa
#samtools faidx $hg38 chr14:22095000-22470000 |grep -v ">chr" >>merge.fa
echo ">chr7_alt_trb" >> merge.fa
samtools faidx $hg38 chr7_KI270803v1_alt:360000-540000 |grep -v ">chr" >> merge.fa

echo ">chr14_igh" >>merge.fa
samtools faidx $chm13 chr14:100130000-101160000 |grep -v ">chr" >> merge.fa
echo ">chr15_igh" >> merge.fa
samtools faidx $hg38 chr15:19962000-19995000 chr15:21714600-21760000 chr15:22170000-22195200|grep -v ">chr"|sed ':z ;N;s/\n/ / ; t z ;' |awk 'gsub(" ","")'  >> merge.fa
#samtools faidx $chm13 chr15:19065000-20034000 |grep -v ">chr" >> merge.fa
echo ">chr16_igh" >> merge.fa
samtools faidx $chm13 chr16:32338000-32458000 chr16:33278000-33388000 chr16:33958000-34188000 chr16:35388000-35588000 |grep -v ">chr" |sed ':z ;N;s/\n/ / ; t z ;' |awk 'gsub(" ","")' >> merge.fa
#samtools faidx $chm13 chr16:32338000-35588000 |grep -v ">chr" >> merge.fa
echo ">chr2_igk" >> merge.fa
samtools faidx $chm13 chr2:87340000-89340000 chr2:89740000-89940000 chr2:90440000-90840000 |grep -v ">chr"|sed ':z ;N;s/\n/ / ; t z ;' |awk 'gsub(" ","")' >> merge.fa

#echo ">chr22_igk" >> merge.fa
#samtools faidx $chm13 chr22:17580000-17612000 |grep -v ">chr" >> merge.fa
echo ">chr22_igl" >> merge.fa
samtools faidx $chm13 chr22:22434000-23350000 |grep -v ">chr" >> merge.fa

mv merge.fa  merge.IG.TR.ref.fasta
samtools faidx merge.IG.TR.ref.fasta
bwa index merge.IG.TR.ref.fasta