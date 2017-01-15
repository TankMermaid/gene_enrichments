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
#library(ggplot2)
library(glmnet)
library(phylolm)
#library(phytools)
source("~/rhizogenomics/src/trunk/compare_genomes/gene_enrichments/functions.r")

load(file = "Burkholderiales_dataset.rdat")

#### Phylogenetic Independent contrasts ####
date()
pic.res <- allgenes_pic(tree = Dat$tree, genes = Dat$genes, class = Dat$Map$Classification)
write.table(pic.res,file = "pic.res.txt", sep = "\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE)
date()

# head(pic.res[ order(pic.res$p.value),])
# gene <- "COG3046"
# plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 1)
# plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 2)
# plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 3)

#### Hypergeometric test ####
date()
hyp.res <- allgenes_hyp(tree = Dat$tree,genes = 1*(Dat$genes >= 1), Map = Dat$Map)
write.table(hyp.res,file = "hyp.res.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)
date()

# ggplot(hyp.res,aes(x = p.value)) + geom_histogram(bins = 20)
# ggplot(hyp.res,aes(x = full.p.value)) + geom_histogram(bins = 20)
# head(hyp.res[ order(hyp.res$p.value),])
# head(hyp.res[ order(hyp.res$full.p.value),])
# gene <- "COG0467"
# subset(hyp.res,gene.id == gene)
# plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 1)
# plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 2)
# plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 3)
# ftable(Dat$Map$Classification ~ Dat$genes[,gene])

#### Phylogenetic logistic regression ####
date()
phyloglm.res <- allgenes_phyloglm(tree = Dat$tree,genes = Dat$genes,Map = Dat$Map)
write.table(phyloglm.res,file = "phyloglm.res.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)
date()

#### Binary Phylogenetic logistic regression ####
date()
phylobinglm.res <- allgenes_phyloglm(tree = Dat$tree,
                                  genes = 1*(Dat$genes >= 1),
                                  Map = Dat$Map)
write.table(phylobinglm.res,file = "phylobinglm.res.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)
date()

#### Ridge logistic regression with PCA controlling for phylogeny ####
# Get distances and make PCA
phendists <- cophenetic(Dat$tree)
pcadist.sum <- summary(prcomp(phendists))

# Get design matrix
X <- as.matrix(Dat$genes)
pcols <- pcadist.sum$x[,1:3]
all(row.names(pcols) == row.names(X))
X <- cbind(X,pcols)
date()
set.seed(4299)
logis.ridge.pca <- cv.glmnet(x = X,y = Dat$Map$Classification, family = "binomial",
                             alpha = 0,nfolds = 10, type.measure = "deviance",
                             lambda.min.ratio = 0.000001,standardize = FALSE,
                             nlambda = 200)
date()
plot(logis.ridge.pca)
logis.ridge.pca.res <- as.matrix(coef(logis.ridge.pca,s = "lambda.min"))
logis.ridge.pca.res <- data.frame(logis.ridge.pca.res)
logis.ridge.pca.res$gene.id <- row.names(logis.ridge.pca.res)
logis.ridge.pca.res$Estimate <- logis.ridge.pca.res$X1
logis.ridge.pca.res$X1 <- NULL
logis.ridge.pca.res <- logis.ridge.pca.res[ -1, ]
logis.ridge.pca.res <- logis.ridge.pca.res[ !(logis.ridge.pca.res$gene.id %in% colnames(pcols)), ]
write.table(logis.ridge.pca.res,file = "logisridgepca.res.txt", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)

qqnorm(logis.ridge.pca.res$Estimate)
qqline(logis.ridge.pca.res$Estimate)

# #### Pagels test #### MOVED TO KURE
# date()
# pagel.res <- allgenes_pagel(tree = Dat$tree, Map = Dat$Map, genes = Dat$genes)
# write.table(pagel.res, file = "pagel.res.txt", sep = "\t", quote = FALSE, row.names = FALSE,
#             col.names = TRUE)
# date()

# #### Hidden rates models #### NOT BEING USED UNTIL I FIGURE OUT CORRECT LRT TEST
# library(corHMM)
# # Model for correlation
# dat <- data.frame(ID = as.character(Dat$Map$taxon_oid),
#                   Classification = as.numeric(Dat$Map$Classification) - 1,
#                   Gene = 1*(Dat$genes[,"COG0467"] >= 1))
# date()
# m1 <- corDISC(phy = Dat$tree,data = dat, ntraits = 2,model = "ARD")
# m2 <- corDISC(phy = Dat$tree,data = dat, ntraits = 2,model = "ER")
# m3 <- corDISC(phy = Dat$tree,data = dat, ntraits = 2,model = "SYM",
#               node.states = "marginal", diagn = FALSE)
# date()
# 1 - pchisq(2*(m1$loglik - m3$loglik),df = 4)

sessionInfo()
