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
library(ggtree)

ladder <- TRUE
genes <- c("pfam00356","pfam13377","pfam00248")


# plot_gene(gene = "pfam00356",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes, which = 1)

load("pfam_datasets/Actinobacteria-Corynebacterineae_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Actinobacteria-Corynebacterineae")
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
p1 <- gheatmap(p1,as.data.frame(Dat$genes[,genes]),
               colnames = TRUE, color = NA,colnames_position = "top") +
  scale_fill_gradientn(colours = topo.colors(10),trans = "log2",na.value = "white")
p1
ggsave("pfams_tree_actinocoryne.png",p1, width = 8,height = 10)


load("pfam_datasets/Alphaproteobacteria_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Alphaproteobacteria")
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
p1 <- gheatmap(p1,as.data.frame(Dat$genes[,genes]),
               colnames = TRUE, color = NA,colnames_position = "top") +
  scale_fill_gradientn(colours = topo.colors(10),trans = "log2",na.value = "white")
p1
ggsave("pfams_tree_alphaproteo.png",p1, width = 8,height = 10)


library(phylolm)

f1 <- Classification ~ gene
Map <- Dat$Map
Map$gene <- Dat$genes[,"pfam00248"]
ftable(Classification ~ gene, Map)
Map$Classification <- as.numeric(Map$Classification) - 1

m1 <- phyloglm(f1, data = Map, phy = Dat$tree, method = "logistic_IG10")
