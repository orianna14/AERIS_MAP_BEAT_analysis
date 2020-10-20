#! /usr/bin/perl -w
# This script calculates odds-ratio for the association between microbiome change-points and clinical events of interest (here being the endotypic switches).
# The script utilize an embedded R script (6_calcor.R) which should be placed in the same folder with this script.

unless (-e "tables") {
	system ("mkdir OR_table OR_results");
}
# get total number of samples/events for neutrophil and eosinophil subgroups
$totalneutro=71;
$totaleosino=87;
$total=549;
open (IN, $ARGV[0]); ## use the output from 5_changepoint_calc.pl here
$dump=<IN>;
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$total{$a[12]}++;
	if ($a[8] eq 'N' and $a[9] eq 'Y') {
		$neutro{$a[12]}++;
	}
	if ($a[10] eq 'N' and $a[11] eq 'Y') {
		$eosino{$a[12]}++;
	}
}
for my $key (keys %total) {
	open (OUT, ">OR_table/$key.txt");
	print OUT "\tCase\tNocase\n";
	if (exists $neutro{$key}) {
		$a=$neutro{$key};
	}
	else {
		$a=0;
	}
	$b=$total{$key}-$neutro{$key};
	$c=$totalneutro-$a;
	$d=($total-$totalneutro)-$b;
	print OUT "Change\t$a\t$b\n";
	print OUT "Nochange\t$c\t$d\n";
	system ("Rscript calcor.R OR_table/$key.txt OR_results/$key.txt");

}