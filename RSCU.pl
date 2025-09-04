use strict;use warnings;

my %code=();
$code{"TTT"}=$code{"TTC"}="F";
$code{"TCT"}=$code{"TCC"}=$code{"TCA"}=$code{"TCG"}="S";
$code{"TAT"}=$code{"TAC"}="Y";
$code{"TGT"}=$code{"TGC"}="C";
$code{"TGG"}="W";
$code{"CTT"}=$code{"CTC"}=$code{"CTA"}=$code{"CTG"}=$code{"TTA"}=$code{"TTG"}="L";
$code{"CCT"}=$code{"CCC"}=$code{"CCA"}=$code{"CCG"}="P";
$code{"CAT"}=$code{"CAC"}="H";
$code{"CAA"}=$code{"CAG"}="Q";
$code{"CGT"}=$code{"CGC"}=$code{"CGA"}=$code{"CGG"}=$code{"AGA"}=$code{"AGG"}="R";
$code{"ATT"}=$code{"ATC"}=$code{"ATA"}="I";
$code{"ATG"}="M";
$code{"ACT"}=$code{"ACC"}=$code{"ACA"}=$code{"ACG"}="T";
$code{"AAT"}=$code{"AAC"}="N";
$code{"AAA"}=$code{"AAG"}="K";
$code{"AGT"}=$code{"AGC"}="S";
$code{"GTT"}=$code{"GTC"}=$code{"GTA"}=$code{"GTG"}="V";
$code{"GCT"}=$code{"GCC"}=$code{"GCA"}=$code{"GCG"}="A";
$code{"GAT"}=$code{"GAC"}="D";
$code{"GAA"}=$code{"GAG"}="E";
$code{"GGT"}=$code{"GGC"}=$code{"GGA"}=$code{"GGG"}="G";
$code{"TAA"}=$code{"TAG"}=$code{"TGA"}="stop";

my %aa;
foreach my $codon(keys %code){
	push @{$aa{$code{$codon}}},$codon;
}

my $cdna=shift;
my $reads=shift;

my %cdna;my %start_pos;
open(F,"gzip -dc $cdna|");
while(my $header=<F>){
	my $seq=<F>;$seq=~s/\s+$//;
	
	$header=~s/\s+$//;
	$header=~s/^>//;
	
	my($start,$stop,$tr)=$header=~/(\d+):(\d+):(\d+)$/;
	$start_pos{$header}{'start'}=$1;
	$start_pos{$header}{'stop'}=$2;
	$start_pos{$header}{'cds_len'}=$2-$1;
	$cdna{$header}=$seq;
}
close(F);

my %codon_count;my %aa_count;
open(F,$reads);
<F>;
while(my $line=<F>){
	$line=~s/\s+$//;
	my @tps=split(/\t/,$line);
	my $name=shift @tps;
	my $mean=0;
	my $mean_count=0;
	foreach (@tps){
		$mean+=$_;
		$mean_count++;
	}
	next if($mean<1);

	my $cds=substr($cdna{$name},$start_pos{$name}{'start'},$start_pos{$name}{'cds_len'});
	my @codons=$cds=~/(\w{3})/g;
	foreach my $codon(@codons){
		unless(exists $code{$codon}){
			next;
		}

		$codon_count{$code{$codon}}{$codon}++;
		$aa_count{$code{$codon}}++;
	}
}
close(F);

open(F,">RSCU\_$reads");
foreach my $aa(keys %aa_count){
	my $aa_num=scalar(@{$aa{$aa}});
	my $exp_count=$aa_count{$aa}/$aa_num;
	foreach my $codon(@{$aa{$aa}}){
		unless(exists $codon_count{$aa}{$codon}){
			$codon_count{$aa}{$codon}=0;
		}
		
		my $rscu=$codon_count{$aa}{$codon}/$exp_count;
		print F "$aa\t$codon\t$rscu\t$codon_count{$aa}{$codon}\t$exp_count\t$aa_count{$aa}\n";
	}
}
close(F);

