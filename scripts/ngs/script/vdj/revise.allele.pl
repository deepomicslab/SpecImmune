#!/usr/bin/perl -w
my ($ifile,$ofile) = @ARGV;
my (%hashv,%hashj,%hashd,%hashr,%hashV,%hashD,%hashJ,%hashs);

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

open IN, "$ifile" or die "$!\n";
while(<IN>){
	chomp;
	my ($id,$v,$j,$d) = split;
        my $V = (split /\*/, $v)[0];
	my $J = (split /\*/, $j)[0];
	$hashj{$J} .= ";$j";
	my @vvs = split /;/,$v;
	my @nvv;
	foreach my $vv(@vvs){
		my $g = (split /\*/, $vv)[0];
		push @nvv, $g;
	}
        @nvv = uniq(@nvv);
	next if($#nvv > 0);
	$hashs{$V} += 1;
	$hashs{$J} += 1;
	$hashv{$V} .= ";$v";
	if($d){
		my $D = (split /\*/, $d)[0];
		$hashd{$D} .= ";$d";
		$hashs{$D} += 1;
	}
}

open OUT, ">$ofile";
foreach my $j(sort keys %hashj){
	my %hasht;
	my @arrs = (split /;/,$hashj{$j});
	shift @arrs;
	my $sum;
	foreach my $arr(@arrs){
		$hasht{$arr} += 1;
	}
	my @ass;
	foreach my $k(sort{$hasht{$b} <=> $hasht{$a}} keys %hasht){
		push @ass, $k;
		$sum += $hasht{$k};
	}
	if($#ass == 0){$hashJ{$j} = $ass[0];}
	else{
		my $a1 = $ass[0];
		my $a2 = $ass[1];
		my $c1 = $hasht{$a1};
		my $c2 = $hasht{$a2};
		if($c1/$sum >0.9){$hashJ{$j} = $a1}
		else{$hashJ{$j} = "$a1\t$a2";}
	}
}
foreach my $v(sort keys %hashv){
        my %hasht;
        my @arrs = (split /;/,$hashv{$j});
        shift @arrs;
        my $sum;
        foreach my $arr(@arrs){
                $hasht{$arr} += 1/($#arr + 1);
        }
        my @ass;
        foreach my $k(sort{$hasht{$b} <=> $hasht{$a}} keys %hasht){
                push @ass, $k;
                $sum += $hasht{$k};
        }
        if($#ass == 0){$hashV{$j} = $ass[0];}
        else{
                my $a1 = $ass[0];
                my $a2 = $ass[1];
                my $c1 = $hasht{$a1};
                my $c2 = $hasht{$a2};
                if($c1/$sum >0.9){$hashV{$j} = $a1}
		else{$hashV{$j} = "$a1\t$a2";}
        }
}

foreach my $d(sort keys %hashd){
        my %hasht;
        my @arrs = (split /;/,$hashd{$j});
        shift @arrs;
        my $sum;
        foreach my $arr(@arrs){
                $hasht{$arr} += 1/($#arr + 1);
        }
        my @ass;
        foreach my $k(sort{$hasht{$b} <=> $hasht{$a}} keys %hasht){
                push @ass, $k;
                $sum += $hasht{$k};
        }
        if($#ass == 0){$hashD{$j} = $ass[0];}
        else{
                my $a1 = $ass[0];
                my $a2 = $ass[1];
                my $c1 = $hasht{$a1};
                my $c2 = $hasht{$a2};
                if($c1/$sum >0.7){$hashD{$j} = $a1}
		else{$hashD{$j} = "$a1\t$a2";}
        }
}

open IN, "$ifile" or die "$!\n";
while(<IN>){
	chomp;
	my ($id,$v,$j,$d) = split;
        my $V = (split /\*/, $v)[0];
        my $J = (split /\*/, $j)[0];
	my @vvs = split /;/,$v;
        my @nvv;
        foreach my $vv(@vvs){
                my $g = (split /\*/, $vv)[0];
                push @nvv, $g;
        }
        @nvv = uniq(@nvv);
	my $ss = 0;
	if($#nvv > 0){
		foreach my $nv(@nvv){
			my $tt = (split /\*/,$nv)[0];
			if($hashs{$tt} > $ss){$V = $tt; $ss = $hashs{$tt}}
		}
		
	}
        my @vvs = split /\t/, $hashV{$V};
	my @jjs = split /\t/, $hashJ{$J};
	if($#vvs == 0){$v = $vvs[0];}
	else{
		my @vs = split /;/,$v;
		my $tag;
		foreach my $i(@vs){
			if($i eq $vvs[0]){$tag .= "a";}
			if($i eq $vvs[1]){$tag .= "b";}
		}
		if($tag =~ /a/){
			if($tag =~ /b/){$v = "$vvs[0];$vvs[1]";}
			else{$v = $vvs[0]}
		}elsif($tag =~ /b/){
			$v = $vvs[1];
		}
	}
	if($#jjs == 0){$j = $jjs[0];}
	else{
                my @js = split /;/,$j;
                my $tag;
                foreach my $i(@js){
                        if($i eq $jjs[0]){$tag .= "a";}
                        if($i eq $jjs[1]){$tag .= "b";}
                }
                if($tag =~ /a/){
                        if($tag =~ /b/){$j = "$jjs[0];$jjs[1]";}
                        else{$j = $jjs[0]}
                }elsif($tag =~ /b/){
                        $j = $jjs[1];
                }
        }
        if($d){
		my $D = (split /\*/, $d)[0];
                my @dds = split /\t/, $hashD{$D};
		if($#dds == 0){$d = $dds[0]}
		else{
			my @ds = split /;/, $d;
			my $tag;
                        foreach my $i(@ds){
                           if($i eq $vvs[0]){$tag .= "a";}
                           if($i eq $vvs[1]){$tag .= "b";}
                        }
                        if($tag =~ /a/){
                            if($tag =~ /b/){$v = "$vvs[0];$vvs[1]";}
                            else{$v = $vvs[0]}
                        }elsif($tag =~ /b/){
                            $v = $vvs[1];
                        }
		}
	}
}

		




}
close IN;

close OUT;


