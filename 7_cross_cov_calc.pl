#! /usr/bin/perl -w
# This script calculates cross-covariance between microbiome taxa and clinical measures of interest (here being eosinophilic/neutrophilic counts) for longitudinal samples of each patient.
# The script utilizes an embedded R script (7_ccf.R) which should be placed in the same folder with this script.

unless (-e "tables") {
	system ("mkdir tables ccf_results");
}

## read in OTU table ###
open (IN, $ARGV[0]); ##otu table
$header = <IN>;
chop $header;
@headers = split ("\t", $header);
open (REF, ">$ARGV[0]\_match");
$num = 1;
while (<IN>) {
	chop;
	@a = split ("\t", $_);
	my $taxon = "Taxon$num";
	for my $i (1..$#a) {
		$match{$taxon} = $a[0];
		$otu{$taxon}{$headers[$i]} = $a[$i];
		$transotu{$headers[$i]}{$taxon} = $a[$i];
	}
	print REF $taxon."\t".$match{$taxon}."\n";
	$num ++;
}

## read in metadata table which contains columns in order: sampleID\tSubjectID\tDate\tVisittype\tGroup\tNeutrophilic_perc\tEosinophilic_perc ###
open (IN, $ARGV[1]);
$dump=<IN>;
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$date{$a[1]}{$a[0]}=$a[2];
	$visit{$a[1]}{$a[0]}=$a[3];
	$class{$a[0]}=$a[4];
	$neutrop{$a[0]}=$a[5];
	$eosinop{$a[0]}=$a[6];
	$num{$a[1]} ++;
}

for my $key (sort keys %date) {
	next unless ($num{$key}>=5);
	for my $key2 (sort keys %otu) {
		open (OUT, ">tables/$key\_$key2.txt");
		print OUT "SampleID\tdate\tVisittype\tAbundance\tNeutro\tEosino\n";
		for my $key3 (sort {$date{$key}{$a}<=>$date{$key}{$b}} keys %{$date{$key}}) {
			print OUT $key3."\t".$date{$key}{$key3}."\t".$visit{$key}{$key3}."\t".$otu{$key2}{$key3}."\t".$neutrop{$key3}."\t".$eosinop{$key3}."\n";
		}
		system ("Rscript ccf.R tables/$key\_$key2.txt ccf_results/$key\_$key2.txt");
		open (IN2, "ccf_results/$key\_$key2.txt");
		my $count=1;
		while (<IN2>) {
			chop;
			if ($count==4) {
				s/^\s+//g;
				@a=split(/\s+/,$_);
			}
			if ($count==5) {
				@b=split(/\s+/,$_);
			}
			$count ++;
		}
		for my $val (0..$#a) {
			if ($a[$val]==0) {
				$cor{$key}{$key2}=$b[$val];
				$abscor{$key}{$key2}=abs($b[$val]);
				$lag{$key}{$key2}=$a[$val];
			}
		}
	}
}

for my $key (keys %cor) {
	$count = scalar keys %{$cor{$key}};
	my $j = 1;
	open (OUT3, ">ccf_results_zerolag/$key\-results");
	print OUT3 "Taxon\tCor\tLag\n";
	for my $key2 (sort {$abscor{$key}{$b} <=> $abscor{$key}{$a}} keys %{$cor{$key}}) {
		#$fdr = $pvalue{$key}{$key2}*$count/$j;
		print OUT3 $match{$key2}."\t".$cor{$key}{$key2}."\t".$lag{$key}{$key2}."\n";
		$j ++;
	}

}
