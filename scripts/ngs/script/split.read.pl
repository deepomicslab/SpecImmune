#!/usr/bin/perl -w
my ($infa,$samfile,$outfa) = @ARGV;
my (%hash, %hashl, %hashf, %hashq,%hashs);

open OUT, ">$outfa";
open FA, "$infa" or die "$infa\t$!\n";
while(<FA>){
	chomp;
	s/^\@//;
	my $id = $_;
	my $fa = <FA>;
	<FA>;
	my $qual = <FA>;
	chomp $fa; chomp $qual;
	$hashf{$id} = $fa;
	$hashq{$id} = $qual;
}
close FA;

open IN, "$samfile" or die "$!\n";
while(<IN>){
	chomp;
        next if(/^\@/);
        my ($id,$gene,$pos,$qual,$cigar,$nm,$sa) = (split)[0,2,3,4,5,11,20];
        my $seq = $hashf{$id};
	my $len = length($seq);
        if(!exists $hashl{$id}){$hashl{$id} = $len;}
        $nm =~ s/NM:i://;
        my $cigar1 = $cigar;
        my $match = 0;
        next if ($cigar eq "*");
        while($cigar1){
                $cigar1 =~ s/^(\d+)([M|S|H|I|D]+)//;
                if($2 eq "M"){$match += $1}
        }
        next if ( $match < 500);
        my $score = $nm / $match;
        next if($score >0.04);
	
	if($sa){
		next if($sa =~ m/$gene/);
		my $value = "$gene\t$pos\t$qual\t$cigar\t$nm\t$score\t$match";
		$hash{$id} .= "$value\n";
		#print "$id\t$gene\t$pos\t$qual\t$cigar\t$nm\t$sa\t$score\t$match\n";}
	}
}
close IN;

foreach my $key(sort keys %hash){
	my @lines = (split /\n/,$hash{$key});
	next if($#lines == 0);
	if($#lines > 1){$oo = join("\n",@lines); print "more\t$key\n$oo\n"}
	
	my %hasht;
	my @cigars;
	my ($p,$pos1,$pos2) = (0,0,0);
	my $len = $hashl{$key};
	if($#lines > 1){
		foreach my $line(@lines){
			my $cigar = (split /\t/,$line)[3];
			if($cigar =~ /^(\d+)[S|H]\S+[M|I|D](\d+)[S|H]$/){
			    $pos1 = $1; $pos2 = $len - $2;
			    if($1>500 && $2>500){$hashp{$key} = "$pos1\t$pos2";}
			}
		}
	}
	else{
	   foreach my $line(@lines){
		my ($s1, $s2) = (0,0);
		my ($cigar,$score) = (split /\t/, $line)[3,-2];
		if($cigar =~ /^(\d+)[S|H]\S+[M|I|D](\d+)[S|H]$/){
			if($1>$2){$s1 = $1; $pos1 = $1 - 100}
			if($1<$2){$s2 = $2; $pos2 = $len - $s2 + 100}
		}
		elsif($cigar =~ /^(\d+)[S|H]/){next if($1<500);$s1 = $1; $pos1 = $1 -100}
		elsif($cigar =~ /(\d+)[S|H]$/){next if($1<500);$s2 = $1; $pos2 = $len - $s2 + 100}
		my $pos = $pos2;
		if($pos1 > $pos2){$pos = $pos1}
		$p += $pos;
		#print "$key\t$pos\t$line\n";
	   }
	$p = $p /2;
        $hashp{$key} = int($p) ;
	#print "$key\t$len\t$p\n";
     }
}

foreach my $id(sort keys %hashf){
        my $fa = $hashf{$id};
	my $qual = $hashq{$id};
	if(!exists $hashp{$id}){print OUT "\@$id\n$fa\n+\n$qual\n"}
	else{
		#print "$id\t$hashl{$id}\t$hashp{$id}\n";
		my @ps = (split /\t/,$hashp{$id});
		my ($s,$i) = (0,0);
		foreach my $p(@ps){
		     my $read1 = "\@$id"."_$i";
		     my $fa1 = substr($fa,$s,$p);
		     my $qual1 = substr($qual,$s,$p);
		     $s = $p; $i += 1;
		     print OUT "$read1\n$fa1\n+\n$qual1\n";
		}
		my $read2 = "\@$id"."_$i";
                my $fa2 = substr($fa,$s);
		my $qual2 = substr($qual,$s);
		print OUT "$read2\n$fa2\n+\n$qual2\n";

	}
}

close OUT;
