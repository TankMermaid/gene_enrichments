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

library(AMOR)
library(ape)

# dataset_file <- "datasets/Burkholderiales_dataset.rdat"
# enrichment_file <- "enrichments/Burkholderiales.phyloglm.res.txt"
dataset_file <- "datasets/Pseudomonas_dataset.rdat"
enrichment_file <- "enrichments/Pseudomonas.phyloglm.res.txt"

load(dataset_file)
Enrich <- read.table(enrichment_file,sep = "\t", header = TRUE)

Dat.pca <- PCA(create_dataset(t(Dat$genes),Dat$Map), cor = TRUE)
summary(Dat.pca)
plotgg(Dat.pca, col = "Classification")

siggenes <- subset(Enrich, !is.na(BH) & BH < 0.05)$gene.id
siggenes <- subset(Enrich, !is.na(BH) & BH < 0.05 & Estimate > 0)$gene.id
length(siggenes)
Dat.pca <- PCA(remove_taxons(create_dataset(t(Dat$genes),Dat$Map),
                             taxons = colnames(Dat$genes)[ !(colnames(Dat$genes) %in% siggenes) ]),
               cor = TRUE)
summary(Dat.pca)
plotgg(Dat.pca, col = "IMG_Genus", shape = "Classification", point_size = 3) +
  scale_shape_manual(values = c(19,21))
