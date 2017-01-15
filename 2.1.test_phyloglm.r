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
library(ape)
library(phylolm)

opts <- commandArgs(trailingOnly = TRUE)

# Read parameters
# First always load extra functions
fundir <- opts[3]
if(is.na(fundir)){
  fundir <- "~/rhizogenomics/src/trunk/compare_genomes/gene_enrichments/"
}
source(paste(fundir,"/functions.r",sep = ""))

# Then read input files
datfile <- opts[1]
if(is.na(datfile)){
  datfile <- "Burkholderiales_dataset.rdat"
}

outfile <- opts[2]
if(is.na(outfile)){
  outfile <- "Burkholderiales_phyloglm.res.txt"
}

# Read input
load(file = datfile)

#### Phylogenetic logistic regression ####
date()
phyloglm.res <- allgenes_phyloglm(tree = Dat$tree,genes = Dat$genes,Map = Dat$Map)
phyloglm.res$BH <- p.adjust(phyloglm.res$p.value,method = 'fdr')
phyloglm.res$Bonferroni <- p.adjust(phyloglm.res$p.value,method = 'bonferroni')
write.table(phyloglm.res,file = outfile, sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)
date()

sessionInfo()
