#!/usr/bin/perl -w
use FindBin qw($Bin);
use Getopt::Long;

my ($sample, $gene, $dir, $pop, $help);

GetOptions(
           "s=s"     =>      \$sample,
           "g=s"     =>      \$gene,
           "i=s"     =>      \$dir,
           "p=s"     =>      \$pop,
           "h"       =>      \$help
);
my $usage = <<USE;
Usage:
description: combination selection
usage: perl $0 [options]
        Common options:
        -s       <tr>    sample name
        -g       <tr>    gene name "HLA_A|HLA_B|HLA_C|HLA_DPA1|HLA_DPB1|HLA_DQA1|HLA_DQB1|HLA_DRB1"
        -i       <tr>    the directory of phased sequence
        -p       <tr>    population information "Asian|Black|Caucasian|Unknown|nonuse"
        -help|?           print help information
e.g.:
        perl $0 -s samplename -g HLA_A -i indir -p Unknown 
USE
die $usage unless ($sample && $dir && $pop && $gene) ;
print "parameter:\tsample:$sample\tdir:$dir\tpop:$pop\tgene:$gene\n";
my $k = 2;
my (%hashp, %hashpp, %hashg, %hashc, %hash,%hashdd);
my $db="$Bin/../db/";
my $bin="$Bin/../bin";
my $fadir=$dir;
my $workdir = "$dir/tmp";
`mkdir  -p $workdir`;

my $class = $gene;
my $gene_class = "HLA";
if($class =~ /^CYP/){$gene_class = "CYP"}
if($class =~ /^KIR/){$gene_class = "KIR"}
my $ref="$db/$gene_class/ref/split/$class.exon.fasta";
my (%hashs,%hash_max);
my %idks;
`ls $dir/$class.*.fasta|grep -v total >$dir/tfile.list`;
open TE, "$dir/tfile.list" or die "$!\n";
open OUT, ">$dir/$class.total.fasta";
open SOUT, ">$dir/result.$class.fasta";
while(my $file=<TE>){
	chomp $file;
	my $id = (split /\//,$file)[-1];
	$id =~ s/\.fasta//;
	print OUT ">$id\n";
	my @tes = (split /\./,$id); pop @tes;my $idk = join(".",@tes);$idks{$idk}=$idk;
	my $seq;
	open TT, "$file" or die "$!\n";
	while(my $line=<TT>){
		chomp $line;
		next if($line=~/^>/);
		$seq .= $line;
	}
	close TT;
	print OUT "$seq\n";
	$hashs{$id} = $seq;
}
close TE;
close OUT;
`$bin/blastn -query $dir/$class.total.fasta -out $dir/$class.blast.out -db $ref -outfmt 6 -num_threads 4 -max_target_seqs 100 -perc_identity 95 -qcov_hsp_perc 10`;
open BIN, "$dir/$class.blast.out" or die "$!\n";
my (%hash1, %hash2,%hash3, %hashas, %hashhk, %hash4);
while(<BIN>){
                chomp;
                my ($id, $hla, $t, $m, $i, $si) = (split)[0,1,3,4,5,11];
		$hash1{$id}{$hla} += $m + $i;
                $hash2{$id}{$hla} += $t - $i;
                $hash4{$id}{$hla} += $si;
}
close BIN;
my %hashh;
foreach my $id (sort keys %hash1){
                my (%hashs,$mscore);
		$mscore = 0 ;
                foreach my $hla(sort keys %{$hash1{$id}}){
                                my $mis = $hash1{$id}{$hla};
                                my $len = $hash2{$id}{$hla};
                                my $score = 100 * (1 - $mis/$len);
                              
                                my $kid = $hla;
    		        	my $nhla = $hla;
				my $fre=0;
                                my $tk = "$id\t$nhla";
                                $hashhk{$tk} = "$score\t$fre";
                                $hash3{$id}{$nhla} = $score;
                                if($score >= $mscore){$mscore = $score}
                                $hashs{$mscore} .= "\t$hla";
                                $hashas{$id} .= "$hla\t";
                }
                $hashh{$id} = "$mscore\t$hashs{$mscore}";
}
my %hashss;
my $mscore = 0;
foreach my $idk (sort keys %idks){
                 my $ts = 0;
                 for(my $i=1;$i<=$k;$i++){
                         my $idd = "$idk"."."."$i";
                         my $s = (split /\t/, $hashh{$idd})[0];
                         $ts += $s;
                 }
                 $ts = $ts / $k;
                  print "$idk\t$ts\n";
                 if($ts >= $mscore){$mscore = $ts}
                 $hashss{$ts} .= "$idk\t";
}
my @tids = (split /\t/, $hashss{$mscore});
my (%fhash,$ffre,%thash);
$ffre=0;
foreach my $tid(@tids){
                 my $tfre=0;
                 for(my $i=1;$i<=$k;$i++){
                         my $hkey = "$class"."_"."$i";
                         my $idd = "$tid"."."."$i";
                         my @arrs = (split /\t+/, $hashh{$idd});        
                         my $aa = shift @arrs;
                         my (%nhash,$mfre);
			 $mfre =0;
                         foreach $od(@arrs){
                                 my $id = $od;
				 #my $nd = $od;
				 my $fre=0;
				 #if(exists $hash{$id}){$fre=$hash{$id}}
                                 $nhash{$fre} .= "$id\t";
                                 #               print "$idd\t$od\t$nd\t$fre\n";
				 if($fre >= $mfre){$mfre = $fre}
                         }
                         $tfre += $mfre;
                         if($mfre == 0){$tfre = $tfre/2;}
                         $thash{$idd} = $nhash{$mfre};
                         #print "$idd\t$tfre\t$aa\n";
                 }
                 if($tfre >=$ffre){$ffre = $tfre}
                 $fhash{$tfre} = $tid;
}
my $select_id = $fhash{$ffre};
print "Select_id\t$select_id\n";
my ($hout,$scout,$nout);
for($i=1;$i<=$k;$i++){
                 $nout .= "\t$class"."_"."$i";
                 my $idd = "$select_id"."."."$i";
                 my $seq = $hashs{$idd};
                 print SOUT ">$idd\n$seq\n";
                 my ($hh,$mscorel) = ("",0);
                 $hh = (split /\t/,$thash{$idd})[0];
	 
}	 

close SOUT;
close OUT;

#my $tfa = "$dir/result.$class.fasta";
#`less $tfa |head -2 >$dir/allele.1.$class.fasta`;
#`less $tfa |tail -2 >$dir/allele.2.$class.fasta`;

