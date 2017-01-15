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

use warnings;
use strict;

my $usage = "\$ phylo_glm_all_datasets.pl <indir> <outdir> <fundir>\n";
my $indir = $ARGV[0];
my $outdir = $ARGV[1];
my $fundir = $ARGV[2];
my $logdir = 'logs/';

die $usage unless -d $indir && -d $outdir && -d $fundir;
die "$logdir does not exist\n" unless -d $logdir;

opendir(IN,$indir) or die $!;
my @files = grep{-f "$indir/$_" && /\w+_dataset\.rdat/} readdir IN;

for(@files){
	my ($group,@trash) = split(/_/,$_);
	my $command = "bsub -o $logdir/$group.phyloglm.log -e $logdir/$group.phyloglm.err -J $group.phyloglm Rscript $fundir/2.1.test_phyloglm.r $indir/$_ $outdir/$group.phyloglm.res.txt $fundir";
	print "Executing:\n\t>$command\n\n";
	system($command)
}

