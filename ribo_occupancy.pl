use strict;use warnings;

my $cdna=shift;
my $reads=shift;

open(F,"gzip -dc $cdna|");
my %cdna;
while(my $line=<F>){
	my($tr_name)=$line=~/(ENS\w{0,3}T\d+)/;
	$line=<F>;
	$line=~s/\s+$//;
	$cdna{$tr_name}=$line;
}
close(F);

my %codon_reads;my %codon_reads_count;
open(F,"gzip -dc $reads|");
while(my $line=<F>){
	$line=~s/\s+$//;
	my @tps=split(/,/,$line);
	my $name_info=shift @tps;
	my($tr_name,$gene_name,$start,$stop,$tr_len)=split(/:/,$name_info);

	my $total=0;my $total_count=0;
	for(my $pos=$start+99; $pos<=$stop-99;$pos+=3){
		$total+=$tps[$pos];
		if($tps[$pos]>0){
			$total_count++;
		}
	}
	
	
	next if($total<32 or $total_count<5);
    my $cds_len=$stop-$start;
	$total=$total/($cds_len/3);
	
	my %tmp;my %tmp1;
	for(my $pos=$start+99; $pos<=$stop-99;$pos+=3){
		my $codon=substr($cdna{$tr_name},$pos,3);
		$tmp{$codon}+=($tps[$pos]/$total);
		$tmp1{$codon}++;
	}
	
	foreach my $codon(keys %tmp){
		$tmp{$codon}/=$tmp1{$codon};
		$codon_reads{$codon}+=$tmp{$codon};
		$codon_reads_count{$codon}++;
	}
	
}
close(F);

$reads=~s/^\S+\///;
$reads=~s/\.gz//;
open(F,">ribo_occupancy_$reads");
foreach my $codon(sort keys %codon_reads){
    next if($codon eq "TAA" or $codon eq "TGA" or $codon eq "TAG");
	$codon_reads{$codon}/=$codon_reads_count{$codon};
	print F "$codon\t$codon_reads{$codon}\t$codon_reads_count{$codon}\n";
}
close(F);
