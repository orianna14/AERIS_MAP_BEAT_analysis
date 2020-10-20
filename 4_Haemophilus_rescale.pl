#! /usr/bin/perl -w
# This script rescales the microbiome relative abundance data by downscaling Haemophilus to its average abundance.

# read in OTU table
open (IN, "L6.txt");
$header=<IN>;
chop $header;
@headers=split("\t",$header);
$sum=0;
my %hash=();
my %transhash=();
my %transhash2=();
while (<IN>) {
	chop;
	@a=split("\t",$_);
	for my $i (1..$#a) {
		$hash{$a[0]}{$headers[$i]}=$a[$i];
		$transhash{$headers[$i]}{$a[0]}=$a[$i];
	}
	if ($a[0] eq 'Haemophilus') {
		for my $j (1..$#a) {
			$sum += $transhash{$headers[$i]}{$a[0]};
		}
	}
	$average=$sum/$#a;
}
for my $key (sort keys %transhash) {
	if ($transhash{$key}{"Haemophilus"}>=0.4) {
		$transhash2{$key}{"Haemophilus"}=$average;
		$scale=(1-$average)/(1-$transhash{$key}{"Haemophilus"});
		for my $key2 (sort keys %hash) {
			next if ($key2 eq 'Haemophilus');
			$transhash2{$key}{$key2}=$transhash{$key}{$key2}*$scale;
		}
	}
	else {
		for my $key2 (sort keys %hash) {
			$transhash2{$key}{$key2}=$transhash{$key}{$key2};
		}
	}
}
print "SampleID";
for my $key (sort keys %transhash) {
	print "\t$key";
}
print "\n";
for my $key (sort keys %hash) {
	print $key;
	for my $key2 (sort keys %transhash) {
		print "\t$transhash2{$key2}{$key}";
	}
	print "\n";
}