#! /usr/bin/perl -w
# This script performs taxonomic check for potential contamination based on list from Salter et al. BMC Biology paper.
# This script calls upon a text file 'negative_salter_paper.txt' which should be placed in the same folder with the script.

# read in the contamination taxa list from salter paper
open (IN, "3_negative_salter_paper.txt");
while (<IN>) {
	chop;
	$neg{$_} = 1;
}

# read in OTU table 
open (IN1, "L6.txt");
my $header = <IN1>;
@header = split ("\t", $header);

while (<IN1>) {
	chop;
	@a = split ("\t", $_);
	for my $i (1..$#a) {
		$hash{$a[0]}{$header[$i]} = $a[$i];
		#print $header[$i]."\t".$a[$i]."\n";
	}
}

for my $key (keys %hash) {
	if ($key =~ /g__(\S+)/) {
		if (exists $neg{$1}) {
			($tmp) = ($key =~ /g__(\S+)/);
			my $i = 0;
			my $j = 0;
			my $k = 0;
			my $total = 0;
			for my $key2 (keys %{$hash{$key}}) {
				#print $key2."\n";
				$i ++;
				if ($hash{$key}{$key2} > 0) {
					$j ++;
				}
				if ($hash{$key}{$key2} > 0.1) {
					$k ++;
				}
				$total += $hash{$key}{$key2};
			}
			my $average = $total/$i;
			my $occur = $j/$i;
			my $occur2 = $k/$i;
			$print{$tmp} = $occur."\t".$occur2."\t".$average;
			#print $tmp."\t".$occur."\t".$occur2."\t".$average."\n";
		}
	}
}

open (IN2, "3_negative_salter_paper.txt");
while (<IN2>) {
	chop;
	if (exists $print{$_}) {
		print $_."\t".$print{$_}."\n";
	}
	else {
		print $_."\t0\t0\t0\n";
	}
}
