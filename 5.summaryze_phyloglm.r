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
library(ggplot2)
#library(qvalue)
library(clusterProfiler)
source("~/rhizogenomics/src/trunk/compare_genomes/gene_enrichments/functions.r")
library(ape)
library(ggtree)
########## FUNCTIONS ###############
summarize_files <- function(dir, prefix = ""){
  files <- list.files(dir)
  
  Res <- NULL
  DA <- NULL
  for(file in files){
    #file <- files[1]
    #file <- files[5]
    taxonomic_group <- strsplit(file,split = "[.]")[[1]][1]
    
    enrich <- read.table(paste(dir,"/",file,sep = ""),header = TRUE)
    ngenes <- nrow(enrich)
    ntests <- ngenes - sum(is.na(enrich$p.value))
    enrich$q.value <- p.adjust(enrich$p.value, 'fdr')
    #qvals <- qvalue(enrich$p.value)
    
    enrichments <- as.character(enrich$gene.id[enrich$q.value < 0.05 & enrich$Estimate > 0])
    enrichments <- enrichments[!is.na(enrichments)]
    n.enrich <- length(enrichments)
    n.dep <- sum(enrich$q.value < 0.05 & enrich$Estimate < 0, na.rm = TRUE)
    
    if(n.enrich > 0){
      da <- data.frame(Taxon = taxonomic_group,
                       gene.id = enrichments,
                       type = "Enrichment")
      DA <- rbind(DA,da)
      rm(da)
    }
        
    p1 <- ggplot(enrich,aes(x = p.value)) +
      ggtitle(taxonomic_group) +
      geom_histogram(bins = 20)
    #p1
    #print(p1)
    ggsave(paste(prefix,".",taxonomic_group,".pvals.hist.svg",sep=""),p1, width = 4, height = 3)
    res <- data.frame(Taxon = taxonomic_group,
                      ngenes = ngenes,
                      ntests = ntests,
                      n.enrich = n.enrich,
                      n.dep = n.dep) 
    Res <- rbind(Res,res)
  }
  
  return(list(Tab = Res,DA = DA))
}

###################################

dir <- "~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/cdhit_enrichments/enrichments/"
#dir <- opts[1]
res.cdhit <- summarize_files(dir = dir,prefix = "cdhit")
res.cdhit$Tab
#head(res.cdhit$DA)
sort(table(res.cdhit$DA$gene.id),decreasing = TRUE)[1:50]

dir <- "~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/cog_enrichments/enrichments/"
res.cog <- summarize_files(dir = dir, prefix = "cog")
res.cog$Tab
#head(res.cog$DA)
sort(table(res.cog$DA$gene.id),decreasing = TRUE)[1:50]

dir <- "~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/enrichments/"
res.ko <- summarize_files(dir = dir, prefix = "ko")
res.ko$Tab
#head(res.ko$DA)
sort(table(res.ko$DA$gene.id),decreasing = TRUE)[1:50]

dir <- "~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/pfam_enrichments/enrichments/"
res.pfam <- summarize_files(dir = dir, prefix = "pfam")
res.pfam$Tab
#head(res.pfam$DA)
sort(table(res.pfam$DA$gene.id),decreasing = TRUE)[1:50]

dir <- "~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/tigrfam_enrichments/enrichments/"
res.tigrfam <- summarize_files(dir = dir, prefix = "tifgrfam")
res.tigrfam$Tab
#head(res.tigrfam$DA)
sort(table(res.tigrfam$DA$gene.id),decreasing = TRUE)[1:50]

# Compare signal accross gene definitions
res <- data.frame(cog = 100*res.cog$Tab$n.enrich / res.cog$Tab$ntests,
                  ko = 100*res.ko$Tab$n.enrich / res.ko$Tab$ntests,
                  cdhit = 100*res.cdhit$Tab$n.enrich / res.cdhit$Tab$ntests,
                  pfam = 100*res.pfam$Tab$n.enrich / res.pfam$Tab$ntests,
                  tigrfam = 100*res.tigrfam$Tab$n.enrich / res.tigrfam$Tab$ntests)
res <- res[ rowSums(res) != 0, ]
#gplots::heatmap.2(log2(as.matrix(res) + 1), trace = "none")
#gplots::heatmap.2(sqrt(as.matrix(res)), trace = "none")
#gplots::heatmap.2(as.matrix(res), trace = "none")


# Use kegg to look for module enrichments
# Using clusterProfiler package to get KEGG modules
# Ignore p-values at this point. I am not sure they are correct
res.ko$Tab
dir <- "~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/enrichments/"
outdir <- "summaries/"
dir.create(outdir)
files <- list.files(dir)

for(file in files){
  #file <- files[1]
  taxonomic_group <- strsplit(file,split = "[.]")[[1]][1]
  
  x <- read.table(paste(dir,"/",file,sep=""), header = TRUE,sep = "\t", quote = "")
  x$q.value <- p.adjust(x$p.value, method = 'fdr')
  x$gene.id <- sub(pattern = "^KO:", replacement = "", x = as.character(x$gene.id))
  universe <- x[ !is.na(x$Estimate), ]
  enriched <- universe$gene.id[ universe$Estimate > 0 & universe$q.value < 0.05 ]
  depleted <- universe$gene.id[ universe$Estimate < 0 & universe$q.value < 0.05 ]
  universe <- universe$gene.id
  xmod <- enrichMKEGG(gene = enriched,
                      organism = "ko",
                      keyType = 'kegg',
                      pvalueCutoff = 1,
                      pAdjustMethod = "BH",
                      universe = universe,
                      minGSSize = 1,
                      maxGSSize = 1000,
                      qvalueCutoff = 1)
  if(!is.null(xmod)){
    print(dim(xmod@result))
    print(head(xmod@result,10))
    #hist(xmod@result$pvalue)
    write.table(xmod@result,file = paste(outdir,"/",taxonomic_group,"_enrichmods.txt",sep = ""),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  xmod <- enrichKEGG(gene = enriched,
                      organism = "ko",
                      keyType = 'kegg',
                      pvalueCutoff = 1,
                      pAdjustMethod = "BH",
                      universe = universe,
                      minGSSize = 1,
                      maxGSSize = 1000,
                      qvalueCutoff = 1)
  if(!is.null(xmod)){
    print(dim(xmod@result))
    print(head(xmod@result,10))
    #hist(xmod@result$pvalue)
    write.table(xmod@result,file = paste(outdir,"/",taxonomic_group,"_enrichpath.txt",sep = ""),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  
  
  # Depletions
  xmod <- enrichMKEGG(gene = depleted,
                      organism = "ko",
                      keyType = 'kegg',
                      pvalueCutoff = 1,
                      pAdjustMethod = "BH",
                      universe = universe,
                      minGSSize = 1,
                      maxGSSize = 1000,
                      qvalueCutoff = 1)
  if(!is.null(xmod)){
    print(dim(xmod@result))
    print(head(xmod@result,10))
    #hist(xmod@result$pvalue)
    write.table(xmod@result,file = paste(outdir,"/",taxonomic_group,"_depletemods.txt",sep = ""),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  xmod <- enrichKEGG(gene = depleted,
                      organism = "ko",
                      keyType = 'kegg',
                      pvalueCutoff = 1,
                      pAdjustMethod = "BH",
                      universe = universe,
                      minGSSize = 1,
                      maxGSSize = 1000,
                      qvalueCutoff = 1)
  if(!is.null(xmod)){
    print(dim(xmod@result))
    print(head(xmod@result,10))
    #hist(xmod@result$pvalue)
    write.table(xmod@result,file = paste(outdir,"/",taxonomic_group,"_depletepath.txt",sep = ""),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
}


# Plot each phylogeny
library(ggtree)
dir.create("pa_npa_trees/")
ladder <- TRUE

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Acinetobacter_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Acinetobacter") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Acinetobacter.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.13, gene.pos = 0.135)
# title("Acinetobacter")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Bacillus_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Bacillus") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Bacillus.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.47, gene.pos = 0.5)
# title("Bacillus")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Bacteroidetes_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Bacteroidetes") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Bacteroidetes.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.8)
# title("Bacteroidetes")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Burkholderiales_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Burkholderiales") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Burkholderiales.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.53)
# title("Burkholderiales")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Caulobacteraceae_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Caulobacteraceae") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Caulobacteraceae.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.22)
# title("Caulobacteraceae")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Microbacteriaceae_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Microbacteriaceae") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Microbacteriaceae.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.33)
# title("Microbacteriaceae")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Micrococcaceae_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Micrococcaceae") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Micrococcaceae.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.24)
# title("Micrococcaceae")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Mycobacterium_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Mycobacterium") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Mycobacterium.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.15)
# title("Mycobacterium")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Nocardiaceae_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Nocardiaceae") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Nocardiaceae.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K10543",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.17)
# title("Nocardiaceae")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Nocardioidaceae_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Nocardioidaceae") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Nocardioidaceae.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.3)
# title("Nocardioidaceae")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Paenibacillus_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Paenibacillus") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Paenibacillus.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.2)
# title("Paenibacillus")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Pseudomonas_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Pseudomonas") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Pseudomonas.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.26, gene.pos = 0.27)
# title("Pseudomonas")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Rhizobiales_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Rhizobiales") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Rhizobiales.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.5)
# title("Rhizobiales")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Sphingomonas_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Sphingomonas") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Sphingomonas.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.25, gene.pos = 0.26)
# title("Sphingomonas")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Streptomyces_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Streptomyces") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Streptomyces.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.19, gene.pos = 0.20)
# title("Streptomyces")

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/ko_enrichments/datasets/Xanthomonadaceae_dataset.rdat")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Xanthomonadaceae") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
#p1
ggsave("pa_npa_trees/Xanthomonadaceae.tree.svg", p1, width = 5, height = 6)
# plot_gene(gene = "KO:K00001",tree = Dat$tree,
#           Map = Dat$Map,genes = Dat$genes,which = 1,
#           meta.pos = 0.49)
# title("Xanthomonadaceae")

# Plot Asaf's genes
load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/cdhit_enrichments/datasets2/Rhizobiales_dataset.rdat")
clusters <- c(80592,68111,131108, 139834, 41371, 132279,
              21014, 73034, 120818, 67240,61121, 100143,
              26740, 67653, 37724, 60512, 28411, 48892)
clusters <- paste("Cluster", clusters,sep = "")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Rhizobiales") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
p1 <- gheatmap(p1,as.data.frame(Dat$genes[,clusters]),
               colnames = FALSE, color = NA) +
  scale_fill_gradientn(colours = topo.colors(10))
#p1
ggsave("pa_npa_trees/asaf_rhizobiales.tree.svg", p1, width = 7, height = 6)

load("~/rhizogenomics/experiments/2016/2016-04-28.gene_enrichments/cdhit_enrichments/datasets/Burkholderiales_dataset.rdat")
clusters <- c(177294,177921)
clusters <- paste("Cluster", clusters,sep = "")
p1 <- ggtree(Dat$tree, ladderize = ladder) +
  ggtitle("Burkholderiales") +
  geom_treescale() +
  theme_tree()
p1 <- p1 %<+% Dat$Map + geom_tippoint(aes(color = Classification)) +
  scale_color_manual(values = c("brown","green"))
p1 <- gheatmap(p1,as.data.frame(Dat$genes[,clusters]),
               colnames = FALSE, color = NA) +
  scale_fill_gradientn(colours = topo.colors(10))
#p1
ggsave("pa_npa_trees/asaf_burkholderiales.tree.svg", p1, width = 5, height = 6)



