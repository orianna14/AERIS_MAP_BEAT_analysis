#! /usr/bin/perl -w
# This script calculates change-point for microbiome taxa for longitudinal microbiome datasets using PELT algorithm.
# The script utilize an embedded R script (5_PELT.R) which should be placed in the same folder with this script.

unless (-e "tables") {
	system ("mkdir tables results");
}

## read in OTU table ##
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
open (IN, "metadata.txt");
$header=<IN>;
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
		## run PELT algorithm in R ###
		system ("Rscript PELT.R tables/$key\_$key2.txt results/$key\_$key2.txt");
		open (IN1, "tables/$key\_$key2.txt");
		$dump=<IN1>;
		while (<IN1>) {
			chop;
			@a=split("\t",$_);
			$visittype{$count}=$a[2];
			$abundance{$count}=$a[3];
			$count ++;
		}
		### parse results and obtain changepoint for each taxa and its associated visittype ###
		open (IN2, "results/$key\_$key2.txt");
		while (<IN2>) {
			chop;
			if (/Changepoint Locations\s+\:\s+(.+)/) {
				last if ($1 eq "");
				@b=split(/\s+/,$1);
				for my $val (@b) {
					if ($visit{$val-1} eq $visit{$val}) {
						$flag=0;
					}
					else {
						$flag=1;
					}
					$diff=$abundance{$val}-$abundance{$val-1};
					print $sub."\t".$val."\t".$visit{$val-1}."\t".$visit{$val}."\t".$abundance{$val-1}."\t".$abundance{$val}."\t".$diff."\t".$flag."\t".$taxa{$taxa}."\n";
				}
			}
		}	
	}
}
