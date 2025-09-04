use strict;use warnings;

my $cdna=shift;
my $align=shift;

my $TIS_codon="ATG,GTG,CTG,TTG,ACG,ATA,ATC,ATT";

my %cdna;my %start_codon;
print "sequence reading\n";
&READ_CDNA($cdna,\%cdna,\%start_codon);

my %read_site;
print "align file reading\n";
&READ_ALIGN($align,\%read_site);

my %uTIS;
print "uORF predicting\n";
&uORF(\%read_site,\%uTIS);

$align=~s/\S+\///;
open(F,">uTIS_predicted_$align");
foreach my $tr_name(keys %uTIS){
	foreach my $stop(sort{$a<=>$b} keys %{$uTIS{$tr_name}}){
		my $start=$uTIS{$tr_name}{$stop}{'start'};
		my $codon=$start_codon{$tr_name}{$stop}{$start};
		print F "$tr_name\t$uTIS{$tr_name}{$stop}{'start'}\t$uTIS{$tr_name}{$stop}{'stop'}\t$uTIS{$tr_name}{$stop}{'read'}\t$uTIS{$tr_name}{$stop}{'inframe'}\t$codon\t$uTIS{$tr_name}{$stop}{'pvalue'}\n";
	}
}
close(F);

### sub function ###

sub uORF{
	my $ref_read_site=shift;
	my $ref_uTIS=shift;
	
	my $count=0;
	open(F,"| gzip -c >potential_uTIS\.gz");
	foreach my $tr_name(keys %{$ref_read_site}){
		
		my($main_start,$main_stop,$tr_len)=$tr_name=~/(\d+):(\d+):(\d+)$/;
		
		foreach my $stop(sort {$a<=>$b} keys %{$$ref_read_site{$tr_name}}){
			foreach my $start(sort {$b<=>$a} keys %{$$ref_read_site{$tr_name}{$stop}}){
				my $uORF_len=$stop-$start;
				if(scalar(keys %{$$ref_read_site{$tr_name}{$stop}{$start}})/$uORF_len<0.1){
					next;
				}
				
				my $output="$tr_name\t$stop\t$start";
				my $total=0;my $inframe=0;my $uorf_site=0;
				for(my $i=$start;$i<$stop;$i++){
					$uorf_site++;
					last if($i>=$main_start);
					
					unless(exists $$ref_read_site{$tr_name}{$stop}{$start}{$i}){
						$output.="\t0";
					}
					else{
						$output.="\t$$ref_read_site{$tr_name}{$stop}{$start}{$i}";
						$total+=$$ref_read_site{$tr_name}{$stop}{$start}{$i};
						my $relative_pos=$i-$start;
						if($relative_pos%3==0){
							$inframe+=$$ref_read_site{$tr_name}{$stop}{$start}{$i};
						}
					}
				}
				if($total>6 and $uorf_site>3){
					$inframe/=$total;
					if($inframe>0.4){
						print F "$output\n";
						$count++;
						#print "$count\n";
					}
				}
			}
		}
	}
	close(F);
	
	print "calculate pvalues\n";
	if(-f "potential_uTIS_pvalue"){
		system("rm potential_uTIS_pvalue");
	}
	system("Rscript WILTEST1.R $count");
	
	my %potential_uTIS;
	open(F,"potential_uTIS_pvalue");
	while(my $line=<F>){
		$line=~s/\s+$//;
		$line=~s/^\s+//;
		$line=~s/\s+/\t/g;
	
		my @info=split(/\t/,$line);
	
		my $tr_name=$info[0];my $stop=$info[1];my $start=$info[2];
		my $pvalue=$info[4];
	
		next if($pvalue>0.1);
		
		$potential_uTIS{$tr_name}{$stop}{$start}=$pvalue;
	}
	close(F);
	
	foreach my $tr_name(keys %potential_uTIS){
		
		my($main_start,$main_stop,$tr_len)=$tr_name=~/(\d+):(\d+):(\d+)$/;
		
		foreach my $stop(keys %{$potential_uTIS{$tr_name}}){
			foreach my $start(sort{$b<=>$a} keys %{$potential_uTIS{$tr_name}{$stop}}){
				
				my $inframe=0;my $total=0;
				for(my $i=$start;$i<$stop;$i++){
					
					last if($i>=$main_start);
					
					if(exists $$ref_read_site{$tr_name}{$stop}{$start}{$i}){
						$total+=$$ref_read_site{$tr_name}{$stop}{$start}{$i};
						my $relative_pos=$i-$start;
						if($relative_pos%3==0){
							$inframe+=$$ref_read_site{$tr_name}{$stop}{$start}{$i};
						}
					}
				}

				my $len=$stop-$start;
				my $read_density=$total/$len;	
				unless(exists $$ref_uTIS{$tr_name}{$stop}){
					$$ref_uTIS{$tr_name}{$stop}{'stop'}=$stop;
					$$ref_uTIS{$tr_name}{$stop}{'start'}=$start;
					$$ref_uTIS{$tr_name}{$stop}{'read'}=$total;
					$$ref_uTIS{$tr_name}{$stop}{'inframe'}=$inframe;
					$$ref_uTIS{$tr_name}{$stop}{'read_density'}=$read_density;
					$$ref_uTIS{$tr_name}{$stop}{'pvalue'}=$potential_uTIS{$tr_name}{$stop}{$start};
				}
				else{
					my $extend_len=$start-$$ref_uTIS{$tr_name}{$stop}{'start'};
					my $extend_read=$total-$$ref_uTIS{$tr_name}{$stop}{'read'};	
					
					next if($extend_read<5);
					
					my $extend_density=$extend_read/$extend_len;
					my $extend_inframe=($inframe-$$ref_uTIS{$tr_name}{$stop}{'inframe'})/$extend_read;

					next if($extend_density/$$ref_uTIS{$tr_name}{$stop}{'read_density'}<0.5 or $extend_inframe<0.6);
					
					$$ref_uTIS{$tr_name}{$stop}{'stop'}=$stop;
					$$ref_uTIS{$tr_name}{$stop}{'start'}=$start;
					$$ref_uTIS{$tr_name}{$stop}{'read'}=$total;
					$$ref_uTIS{$tr_name}{$stop}{'inframe'}=$inframe;
					$$ref_uTIS{$tr_name}{$stop}{'read_density'}=$read_density;
					$$ref_uTIS{$tr_name}{$stop}{'pvalue'}=$potential_uTIS{$tr_name}{$stop}{$start};
				}
			}
		}
	}
}


sub READ_CDNA{
	my $cdna=shift;
	my $ref_cdna=shift;
	my $ref_start_codon=shift;
	
	my $proc=0;
	open(F,"gzip -dc $cdna|");
	while(my $line=<F>){
		$line=~s/\s+$//;$line=~s/>//;
		my $tr_name=$line;
		my($start,$stop,$tr_len)=$line=~/(\d+):(\d+):(\d+)$/;
		
		$line=<F>;
		$line=~s/\s+$//;
		
		$proc++;
		print STDERR "stop codons extracting, $proc genes\r";
		
		for(my $i=3;$i<$start+600;$i++){
			my $relative=$i-$start;
			next if($i>$start and $relative%3==0);   #exclude main cds in-frame";
			next if($i>$stop);   #not consider TIS in 3' UTR;
			
			my $codon=substr($line,$i,3);
			if($codon eq 'TAA' or $codon eq 'TAG' or $codon eq 'TGA'){
				for(my $j=$i-3;$j>0;$j-=3){
					my $codon1=substr($line,$j,3);
					if($codon1 eq "TAA" or $codon1 eq "TAG" or $codon1 eq "TGA"){
						last;
					}
					
					if($TIS_codon=~/$codon1/){
						push @{$$ref_cdna{$tr_name}{$i}},$j;
						$$ref_start_codon{$tr_name}{$i}{$j}=$codon1;
					
					}
				}
			}
		}
	}
	close(F);
	print "\n";
}

sub READ_ALIGN{
	my $align=shift;
	my $ref_reads_site=shift;
	
	open(F,"gzip -dc $align|");
	while(my $line=<F>){
		$line=~s/\s+$//;
		my @tps=split(/,/,$line);
		
		my $tr_name=shift @tps;
		my($start,$stop,$tr_len)=$tr_name=~/(\d+):(\d+):(\d+)$/;
		
		next unless(exists $cdna{$tr_name});
		
		for(my $psite=0;$psite<scalar(@tps);$psite++){
			
			last if($psite>=$start);
			next if($tps[$psite]==0);
			
			foreach my $stop(keys %{$cdna{$tr_name}}){
				foreach my $start(@{$cdna{$tr_name}{$stop}}){
					if($psite>=$start and $psite<=$stop){
						$$ref_reads_site{$tr_name}{$stop}{$start}{$psite}+=$tps[$psite];
					}
				}
			}
		}
	}
	close(F);
}
