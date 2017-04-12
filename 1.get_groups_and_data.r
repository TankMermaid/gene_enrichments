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

# Read parameters
# First always load extra functions
library(ape)
fundir <- opts[5]
if(is.na(fundir)){
  fundir <- "~/rhizogenomics/src/trunk/compare_genomes/gene_enrichments/"
}
source(paste(fundir,"/functions.r",sep = ""))

# Then read input files
mapfile <- opts[1]
if(is.na(mapfile)){
  mapfile <- "~/rhizogenomics/data/compgen_asaf/metadata_3812genomes.tsv"
}

treefile <- opts[2]
if(is.na(treefile)){
  treefile <- "~/rhizogenomics/data/compgen_asaf/fastree_3812_genomes_31_markers.newick"
}

genesfile <- opts[3]
if(is.na(genesfile)){
  genesfile <- "~/rhizogenomics/data/compgen_asaf/cogs_matrix_isolates.tsv"
}

outdir <- opts[4]
if(is.na(outdir)){
  outdir <- "./"
}

class1 <- opts[6]
if(is.na(outdir)){
  class1 <- "PA"
}

class2 <- opts[7]
if(is.na(outdir)){
  class2 <- "NPA"
}

# Read data
Map <- read.table(mapfile,sep = "\t", header = TRUE, comment.char = "")
tree <- read.tree(treefile)
genes <- read.table(genesfile, sep = "\t",
                    header = TRUE, row.names = 1, comment.char = "")
#genes <- t(genes)
row.names(genes) <- sub(pattern = "^X",replacement = "",x = row.names(genes))

# Homogenize samples in all sets
shared <- intersect(row.names(genes),tree$tip.label)
length(shared)
head(shared)
shared <- intersect(shared, as.character(Map$taxon_oid))
length(shared)
head(shared)
Map <- subset(Map,taxon_oid %in% shared)
genes <- genes[shared,]
to_remove <- tree$tip.label[ !(tree$tip.label %in% shared) ] 
tree <- drop.tip(phy = tree,tip = to_remove)

#plot(tree)
ftable(Classification ~ Taxonomic_Assignment, Map)

# Choose a group
for(group in levels(Map$Taxonomic_Assignment)){
  #group <- "Burkholderiales"
  Dat <- get_group(group = group,tree = tree,genes = genes, Map = Map,
                   classes = c(class2,class1))
  filename <- paste(outdir,"/",group,"_dataset.rdat", sep = "")
  save(Dat, file = filename)
  rm(Dat,filename)
}
