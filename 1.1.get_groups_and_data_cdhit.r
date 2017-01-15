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

genesdir <- opts[3]
if(is.na(genesdir)){
  genesdir <- "./data/"
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

genesfiles <- list.files(genesdir)
genesfiles
#pattern <- paste("_",class1,"_",class2,"[.]faa[.]cdhit[.]clstr$",sep = "")
pattern <- paste("_all[.]cdhit[.]clstr$",sep = "")
genesfiles <- genesfiles[grep(pattern = pattern, x = genesfiles) ]
genesfiles

# genes <- read.table(genesfile, sep = "\t",
#                     header = TRUE, row.names = 1)
# genes <- t(genes)
# row.names(genes) <- sub(pattern = "^X",replacement = "",x = row.names(genes))

#plot(tree)
ftable(Classification ~ Taxonomic_Assignment, Map)

# Choose a group
for(file in genesfiles){
  #group <- "Burkholderiales"
  #file <- genesfiles[1]
  
  group <- strsplit(file,split = "_")[[1]][3]
  genes <- read.table(paste(genesdir,"/",file, sep = ""),
                      sep = "\t", header = TRUE, row.names = 1)
  genes <- as.matrix(genes)
  
  Dat <- get_group(group = group,tree = tree,genes = genes, Map = Map,
                   classes = c(class2,class1))
  filename <- paste(outdir,"/",group,"_dataset.rdat", sep = "")
  save(Dat, file = filename)
  rm(Dat,filename)
}
