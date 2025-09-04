use strict;use warnings;

my $reads_aligned_transcriptom=shift;
my $soft_clip=shift;

unless(defined $soft_clip){
	$soft_clip=3;
}

open(F,"samtools view $reads_aligned_transcriptom|");
my %reads;
while(my $line=<F>){
	$line=~s/\s+$//;
	my @tps=split(/\t/,$line);
	my $CIGAR=$tps[5];

	my $soft_clip_left=0;my $soft_clip_right=0;
	$soft_clip_left=$1 if($CIGAR=~/^(\d+)S\d+M/);
	$soft_clip_right=$1 if($CIGAR=~/(\d+)S$/);

	next if($soft_clip_left>$soft_clip);
    #next if($soft_clip_right>$soft_clip);

	my $tr_name=$tps[2];

	my $pos_5end=$tps[3]-1;
	my $pos_psite=$pos_5end+12;
	my $pos_asite=$pos_5end+15;
	$reads{$tr_name}{'pos_5end'}{$pos_5end}++;
	$reads{$tr_name}{'pos_psite'}{$pos_psite}++;
	$reads{$tr_name}{'pos_asite'}{$pos_asite}++;
}
close(F);

$reads_aligned_transcriptom=~s/\S+\///;
open(F,"| gzip -c >reads_to_tr_5end_$reads_aligned_transcriptom\.gz");
open(F1,"| gzip -c >reads_to_tr_psite_$reads_aligned_transcriptom\.gz");
open(F2,"| gzip -c >reads_to_tr_asite_$reads_aligned_transcriptom\.gz");
foreach my $tr_name(keys %reads){
	
	my($name,$gene_name,$start,$stop,$tr_len)=split(/:/,$tr_name);
	
	foreach my $group(keys %{$reads{$tr_name}}){
		my $output=$tr_name;
		for(my $i=0;$i<$tr_len;$i++){
			unless(exists $reads{$tr_name}{$group}{$i}){
				$output.=",0";
			}
			else{
				$output.=",$reads{$tr_name}{$group}{$i}";
			}
		}
		if($group=~/5end/){
			print F $output,"\n";
		}
		elsif($group=~/psite/){
			print F1 $output,"\n";
		}
		elsif($group=~/asite/){
			print F2 $output,"\n";
		}
	}
}
close(F);close(F1);close(F2);
