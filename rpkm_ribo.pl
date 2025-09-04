use strict;use warnings;

my $reads_to_tr=shift;

my $count=0;
open(F,"gzip -dc $reads_to_tr|");
my %read;my $total_reads=0;
while(my $line=<F>){
	$line=~s/\s+$//;
	my @tps=split(/,/,$line);

	my $name=shift @tps;
	my($tr_name,$gene_name,$start,$stop,$tr_len)=split(/:/,$name);
	for(my $pos=$start;$pos<$stop;$pos++){
		#		next if($tps[$pos]>10000);
		$total_reads+=$tps[$pos];
		$read{$name}+=$tps[$pos];
	}
}
close(F);

$total_reads/=1000000;

$reads_to_tr=~s/^\S+\///;
open(F,">rpkm_$reads_to_tr");
foreach my $name(keys %read){
	my($tr_name,$gene_name,$start,$stop,$tr_len)=split(/:/,$name);
	my $cds_len=$stop-$start;$cds_len/=1000;
	my $rpkm=$read{$name}/($total_reads*$cds_len);
	
	print F "$name\t$rpkm\n";
}
close(F);
