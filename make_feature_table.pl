#!/usr/bin/env perl

#    (C) Copyright 2017 Sur Herrera Paredes
#
#    This file is part of Gene enrichment analysis.
#
#    Gene enrichment analysis is free software: you can redistribute itand/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Gene enrichment analysis is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Gene enrichment analysis. If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;

my $usage = "make_feature_table.pl <infile> <outfile>\n";
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $genome_col = 1;
my $gene_col = 3;

die $usage unless -f $infile;

open(IN,$infile) or die $!;
$genome_col--;
$gene_col--;
my (@line,$genome,$gene,%table,%genes,%genomes);
#my $gene_n = 0;
#my $genome_n = 0;
my $i = 0;
while(<IN>){
	chomp;
	@line = split(/\t/,$_);
	$genome = $line[$genome_col];
	$gene = $line[$gene_col];
	$gene =~ s/\s+$//;
	#print "==$gene==\n";

	$genes{$gene} = 1;
	$genomes{$genome} = 1;
	if(exists($table{$gene})){
		if(exists($table{$gene}->{$genome})){
			$table{$gene}->{$genome}++;
		}else{
			$table{$gene}->{$genome} = 1;
		}
	}else{
		$table{$gene}->{$genome} = 1
	}

	$i++;
	#last if $i > 10000;
}
close IN;

open(OUT,'>',$outfile) or die $!;
my @gene_list= keys %genes;
my @genome_list = keys %genomes;
print OUT "\t" . join ("\t",@genome_list) . "\n";
for $gene (@gene_list){
	print OUT "$gene";
	for $genome (@genome_list){
		if(exists($table{$gene}->{$genome})){
			print OUT "\t$table{$gene}->{$genome}";
		}else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close OUT;

print "$i genes processed\n";

