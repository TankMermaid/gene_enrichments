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
# Load libraries
library(ape)
library(ggplot2)
library(reshape2)
source("~/rhizogenomics/src/trunk/compare_genomes/gene_enrichments/functions.r")

# Load data
dir <- "~/rhizogenomics/experiments/2016/2016-04-20.gene_enrichments/"
load(file = paste(dir,"/Burkholderiales_dataset.rdat",sep = ""))
hyp.res <- read.table(paste(dir,"/hyp.res.txt",sep=""), sep = "\t", header = TRUE)
logis.res <- read.table(paste(dir,"/logisridgepca.res.txt",sep = ""), sep = "\t", header = TRUE)
pagel.res <- read.table(paste(dir,"/pagel.res.txt",sep =""), sep = "\t", header = TRUE)
phylobinglm.res <- read.table(paste(dir,"/phylobinglm.res.txt",sep=""), sep = "\t", header = TRUE)
phyloglm.res <- read.table(paste(dir,"/phylobinglm.res.txt",sep=""), sep = "\t", header = TRUE)
pic.res <- read.table(paste(dir,"/pic.res.txt",sep=""), sep = "\t", header = TRUE)

# Get log fold change per gene between classes
logfc <- NULL
for(gene in colnames(Dat$genes)){
  #gene <- colnames(Dat$genes)[1]
  Map <- Dat$Map
  Map$gene <- Dat$genes[,gene]
  mean1 <- (sum(Map$gene[ Map$Classification == "PA" ]) + 1) / (sum(Map$Classification == "PA") + 1)
  mean2 <- (sum(Map$gene[ Map$Classification == "NPA" ]) + 1) / (sum(Map$Classification == "NPA") + 1)
  logfc <- c(logfc,log2(mean1/mean2))
}
rm(Map,gene,mean1,mean2)

# Get p-values for logistic ridge
logis.res$p.value <- sapply(logis.res$Estimate,
                            FUN = function(x,mu = mean(logis.res[,2]), sdev = sd(logis.res[,2])){
                              pval <- 2*(1-pnorm(q = abs(x), mean = mu, sd = sdev))})

# Create master dataset
dat <- data.frame(gene.id = colnames(Dat$genes),
                  gene.mean = log2(1000*apply(Dat$genes,2,mean)),
                  gene.fc = logfc,
                  hypbin.score = hyp.res$score,
                  hypbin.pval = hyp.res$p.value,
                  hypbin.qval = p.adjust(hyp.res$p.value, method = "fdr"),
                  hyp.score = hyp.res$full.score,
                  hyp.pval = hyp.res$full.p.value,
                  hyp.qval = p.adjust(hyp.res$full.p.value, method = "fdr"),
                  ridgelogis.score = logis.res$Estimate,
                  ridgelogis.pval = logis.res$p.value,
                  ridgelogis.qval = p.adjust(logis.res$p.value, method = "fdr"),
                  pagel.score = pagel.res$LR,
                  pagel.pval = pagel.res$p.value,
                  pagel.qval = p.adjust(pagel.res$p.value, method = "fdr"),
                  #phylobinglm.score = phylobinglm.res$Estimate,
                  #phylobinglm.pval = phylobinglm.res$p.value,
                  #phylobinglm.qval = p.adjust(phylobinglm.res$p.value,method = "fdr"),
                  phyloglm.score = phyloglm.res$Estimate,
                  phyloglm.pval = phyloglm.res$p.value,
                  phyloglm.qval = p.adjust(phyloglm.res$p.value,method = "fdr"),
                  pic.score = pic.res$Estimate,
                  pic.pval = pic.res$p.value,
                  pic.qval = p.adjust(pic.res$p.value,method = "fdr"),
                  row.names = NULL)

head(dat)

#Make p-value histograms
p1 <- ggplot(dat,aes(x = hypbin.pval)) + 
  geom_histogram(bins = 20) +
  xlab("p-value") +
  ggtitle("hypbin") +
  theme_classic() 
#p1
ggsave("hypbin.pvalue.hist.svg",p1,width = 4,height = 3)
p1 <- ggplot(dat,aes(x = hyp.pval)) + 
  geom_histogram(bins = 20) +
  xlab("p-value") +
  ggtitle("hyp") +
  theme_classic() 
#p1
ggsave("hyp.pvalue.hist.svg",p1,width = 4,height = 3)
p1 <- ggplot(dat,aes(x = ridgelogis.pval)) + 
  geom_histogram(bins = 20) +
  xlab("p-value") +
  ggtitle("ridgelogis") +
  theme_classic() 
#p1
ggsave("ridgelogis.pvalue.hist.svg",p1,width = 4,height = 3)
p1 <- ggplot(dat,aes(x = pagel.pval)) + 
  geom_histogram(bins = 20) +
  xlab("p-value") +
  ggtitle("pagel") +
  theme_classic() 
#p1
ggsave("pagel.pvalue.hist.svg",p1,width = 4,height = 3)
p1 <- ggplot(dat,aes(x = phyloglm.pval)) + 
  geom_histogram(bins = 20) +
  xlab("p-value") +
  ggtitle("phyloglm") +
  theme_classic() 
#p1
ggsave("phyloglm.pvalue.hist.svg",p1,width = 4,height = 3)
p1 <- ggplot(dat,aes(x = pic.pval)) + 
  geom_histogram(bins = 20) +
  xlab("p-value") +
  ggtitle("pic") +
  theme_classic() 
#p1
ggsave("pic.pvalue.hist.svg",p1,width = 4,height = 3)


# Check number of significant hits
d <- apply(dat[ , grep(pattern = "qval$",x = colnames(dat)) ] < 0.05,2,table)
d
d <- melt(d)
colnames(d) <- c("Significant","Test","N")
p1 <- ggplot(d,aes(x = Test,y = N, fill = Significant)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion of significant") +
  theme(axis.text.x = element_text(angle = 90, color = "black"),
        axis.text.y = element_text(color = "black"))
p1
ggsave("propsignificant.svg",p1,width = 4, height = 5)
rm(d)

# Plot significant as a function of prevalence and fold change
p1 <- ggplot(dat, aes(x = gene.mean, y = gene.fc)) +
  geom_point(aes(color = hyp.qval < 0.05), alpha = 0.2) +
  geom_hline(yintercept = 0) +
  ylab("log2(PA_IDR/NPA_IDR)") +
  xlab("log2(# occurrences / 1000 genomes)") +
  AMOR::theme_blackbox
#p1
ggsave("hyp.volcano.svg",p1,width = 5, height = 4)
p1 <- ggplot(dat, aes(x = gene.mean, y = gene.fc)) +
  geom_point(aes(color = hypbin.qval < 0.05), alpha = 0.2) +
  geom_hline(yintercept = 0) +
  ylab("log2(PA_IDR/NPA_IDR)") +
  xlab("log2(# occurrences / 1000 genomes)") +
  AMOR::theme_blackbox
#p1
ggsave("hypbin.volcano.svg",p1,width = 5, height = 4)
p1 <- ggplot(dat, aes(x = gene.mean, y = gene.fc)) +
  geom_point(aes(color = ridgelogis.qval < 0.05), alpha = 0.2) +
  geom_hline(yintercept = 0) +
  ylab("log2(PA_IDR/NPA_IDR)") +
  xlab("log2(# occurrences / 1000 genomes)") +
  AMOR::theme_blackbox
#p1
ggsave("ridgelogis.volcano.svg",p1,width = 5, height = 4)
p1 <- ggplot(dat, aes(x = gene.mean, y = gene.fc)) +
  geom_point(aes(color = pagel.qval < 0.05), alpha = 0.2) +
  geom_hline(yintercept = 0) +
  ylab("log2(PA_IDR/NPA_IDR)") +
  xlab("log2(# occurrences / 1000 genomes)") +
  AMOR::theme_blackbox
#p1
ggsave("pagel.volcano.svg",p1,width = 5, height = 4)
p1 <- ggplot(dat, aes(x = gene.mean, y = gene.fc)) +
  geom_point(aes(color = phyloglm.qval < 0.05), alpha = 0.2) +
  geom_hline(yintercept = 0) +
  ylab("log2(PA_IDR/NPA_IDR)") +
  xlab("log2(# occurrences / 1000 genomes)") +
  AMOR::theme_blackbox
#p1
ggsave("phyloglm.volcano.svg",p1,width = 5, height = 4)
p1 <- ggplot(dat, aes(x = gene.mean, y = gene.fc)) +
  geom_point(aes(color = pic.qval < 0.05), alpha = 0.2) +
  geom_hline(yintercept = 0) +
  ylab("log2(PA_IDR/NPA_IDR)") +
  xlab("log2(# occurrences / 1000 genomes)") +
  AMOR::theme_blackbox
#p1
ggsave("pic.volcano.svg",p1,width = 5, height = 4)

# Plot PPV
d <- dat[ , grep(pattern = "qval$",x = colnames(dat)) ]
d <- d < 0.05
d2 <- dat[ , grep(pattern = "score$",x = colnames(dat)) ]
d2 <- d2 > 0
res <- matrix(NA,nrow = ncol(d), ncol = ncol(d))
row.names(res) <- colnames(d)
colnames(res) <- colnames(d)
for(i in 1:ncol(d)){
  for(j in 1:ncol(d)){
    #i <- 2
    #j <- 5
    x <- table(d[,i] & d2[,i],d[,j] & d2[,j])
    res[i,j] <- 100 * x[2,2] / sum(x[2,])
    #res[i,j] <- 100*table(d[,i],d[,j])[2,2] / sum(d[,i],na.rm = TRUE)
  }
}
rm(d,d2,i,j,x)
round(res,2)

#heatmap(res,margins = c(10,10))

res <- melt(res)
colnames(res) <- c("Target","Test","PPV")
p1 <- ggplot(res,aes(x = Target, y = Test, fill = PPV)) +
  geom_tile() +
  geom_text(aes(label = round(PPV,2)), col = "black") +
  scale_fill_gradient2(low = "#c51b7d", mid = "#ffffbf",
                       high = "#4d9221", midpoint = 50) +
  theme(axis.text.x = element_text(angle = 90, color = "black"),
        axis.text.y = element_text(color = "black"))
#p1
ggsave("ppv_heatmap.svg",p1, width = 6, height = 6)
rm(res)

# Plot an interesting COG Catched by hypergeometric, but not by phyloglm
d <- subset(dat,hyp.qval < 0.05 & phyloglm.qval > 0.05)
nrow(d)
d[ order(d$hyp.qval,decreasing = FALSE), c("gene.id","hyp.qval","phyloglm.score","phyloglm.qval")]

gene <- "COG2509"
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 1)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 2)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 3)

gene <- "COG4706"
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 1)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 2)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 3)

# Plot an interesting COG Catched by phyloglm but not hypergeometric
d <- subset(dat,hyp.qval > 0.05 & phyloglm.qval < 0.05 & phyloglm.score > 0)
nrow(d)
d[ order(d$hyp.qval,decreasing = FALSE), c("gene.id","hyp.qval","phyloglm.score","phyloglm.qval")]
#head(d[ order(d$hyp.qval,decreasing = FALSE), c("gene.id")])

gene <- "COG2947"
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 1)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 2)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 3)

gene <- "COG0364"
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 1)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 2)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 3)

# Plot an interesting COG phyloglm negative
d <- subset(dat,phyloglm.qval < 0.05 & phyloglm.score < 0)
nrow(d)
d[ order(d$hyp.qval,decreasing = FALSE), c("gene.id","hyp.qval","phyloglm.score","phyloglm.qval")]
head(d[ order(d$hyp.qval,decreasing = FALSE), c("gene.id")])

gene <- "COG0047"
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 1)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 2)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 3)

gene <- "COG0050"
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 1)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 2)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 3)

# Other
gene <- "COG3946"
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 1)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 2)
plot_gene(gene = gene,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 3)

rm(d)

