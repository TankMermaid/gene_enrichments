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

#' Apply Pagel's test
allgenes_pagel <- function(tree, Map, genes){
  #tree <- Dat$tree
  #Map <- Dat$Map
  #genes <- Dat$genes
  
  x <- Map$Classification
  #x <- as.numeric(Map$Classification) - 1
  names(x) <- row.names(Map)
  Res <- NULL
  for(gene in 1:ncol(genes)){
    #gene <- 3
    
    gene.id <- colnames(genes)[gene]
    Map$gene <- genes[,gene]
    #ftable(gene ~ Classification, Map)
    y <- 1*(Map$gene >= 1)
    
    if(length(unique(y)) == 1){
      m1 <- list(lik.ratio = NA, P = NA)
    }else{
      names(y) <- row.names(Map)
      m1 <- tryCatch(fitPagel(tree = tree,x = x, y = y),
                     error = function(e) list(lik.ratio = NA, P = NA))
    }
    
    res <- data.frame(gene.id = gene.id, LR = m1$lik.ratio,
                      p.value = m1$P)
    Res <- rbind(Res,res)
    rm(m1,res)
    if((gene %% 100) == 0 ){
      cat(gene," genes processed\n")
    }
    #cat(gene,"\n")
  }
  
  return(Res)
}

#' Phylogenetic logistic regression
allgenes_phyloglm <- function(tree, genes, Map){
  #tree <- Dat$tree
  #genes <- Dat$genes
  #Map <- Dat$Map
  
  
  f1 <- Classification ~ gene
  Map$Classification <- as.numeric(Map$Classification) - 1
  Res <- NULL
  for(gene in 1:ncol(genes)){
    #gene <- 11
    
    gene.id <- colnames(genes)[gene]
    # plot_gene(gene = gene.id,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 1)
    # plot_gene(gene = gene.id,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 2)
    # plot_gene(gene = gene.id,tree = Dat$tree,Map = Dat$Map,genes = Dat$genes,which = 3)
    Map$gene <- genes[,gene]
    #ftable(gene ~ Classification , Map)
    
    # dat <- data.frame(ID = as.character(Dat$Map$taxon_oid),
    #                   Classification = as.numeric(Dat$Map$Classification) - 1,
    #                   Gene = 1*(Dat$genes[,"COG0467"] >= 1))
    
    
    m1 <- tryCatch(phyloglm(f1, data = Map, phy = tree,
                            method = "logistic_IG10"),
                   error = function(e) list(coefficients = NA))
    
    if(is.na(coef(m1)[1])){
      res <- data.frame(gene.id = gene.id,
                        Estimate = NA,
                        SE = NA,
                        z.value = NA,
                        p.value = NA)
    }else{
      m1.sum <- summary(m1)
      res <- data.frame(gene.id = gene.id,
                        Estimate = m1.sum$coefficients["gene",1],
                        SE = m1.sum$coefficients["gene",2],
                        z.value = m1.sum$coefficients["gene",3],
                        p.value = m1.sum$coefficients["gene",4])
    }
    
    rm(m1,m1.sum)
    Map$gene <- NULL
    Res <- rbind(Res,res)
  }
  return(Res)
}


#' Hypergeometric test on gene counts
allgenes_hyp <- function(tree,genes,Map){
  #tree <- Dat$tree
  #Map <- Dat$Map
  #genes <- Dat$genes >=1
  
  N.pa <- sum(genes[ Map$Classification == "PA", ])
  N.npa <- sum(genes[ Map$Classification == "NPA", ])
  #N.genomes <- nrow(Map)
  Res <- NULL
  for(gene in 1:ncol(genes)){
    #gene <- 1
    
    gene.id <- colnames(genes)[gene]
    Map$gene <- genes[,gene]
    
    #genes.per.class <- aggregate(gene ~ Classification, data = Map,FUN = sum)
    #genomes.per.class <- aggregate(gene ~ Classification, data = Map,FUN = length)
    #N.genes <- sum(Map$gene)
    
    #ftable(gene ~ Classification , Map)
    
    #dhyper(x = 0, m = 2, n = 2, k = 2)
    #dhyper(x = 1, m = 2, n = 2, k = 2)
    #dhyper(x = 2, m = 2, n = 2, k = 2)
    #phyper(q = 0, m = 2, n = 2, k = 2)
    #phyper(q = 1, m = 2, n = 2, k = 2)
    #phyper(q = 2, m = 2, n = 2, k = 2)
    # Following is wrong, counts genomes for zero and genes for 1+
    # pval <- phyper(q = sum(subset(Map, Classification == "PA" & gene > 0)$gene) - 1,
    #                m = sum(Map$gene),
    #                n = nrow(subset(Map,gene == 0)),
    #                k = sum(subset(Map, Classification == "PA")$gene))
    #pval <- 1 - pval
    
    # Binary version
    pval <- phyper(q = nrow(subset(Map, Classification == "PA" & gene > 0)) - 1,
                   m = sum(Map$gene > 0),
                   n = nrow(subset(Map,gene == 0)),
                   k = nrow(subset(Map, Classification == "PA")))
    pval <- 1 - pval
    pval[ pval == 0 ] <- 1e-16
    score <- -log10(pval)
    
    pval2 <- phyper(q = sum(subset(Map, Classification == "PA")$gene) - 1, 
                    m = N.pa, n = N.npa, k = sum(Map$gene))
    pval2 <- 1 - pval2
    pval2[ pval2 == 0 ] <- 1e-16 
    
    res <- data.frame(gene.id = gene.id, score = score, p.value = pval,
                      full.score = -log10(pval2), full.p.value = pval2)
    Map$gene <- NULL
    Res <- rbind(Res,res)
  }
  
  Res <- data.frame(gene.id = Res$gene.id, score = Res$score,
                    z.score = (Res$score - mean(Res$score)) / sd(Res$score),
                    p.value = Res$p.value,
                    full.score = Res$full.score,
                    full.z.score = (Res$full.score - mean(Res$full.score)) / sd(Res$full.score),
                    full.p.value = Res$full.p.value)
  return(Res)
}

#' Plot a gene occurrence
#' 
#' Makes several types of plots
plot_gene <- function(gene,tree,Map,genes,which = 1,meta.pos = 0.52,
                      gene.pos = meta.pos + 0.03){
  
  #gene <- "COG3946"
  
  #ftable(genes.sub[,"COG3946"] ~ Dat$Classification)
  #tips <- row.names(Dat$genes)[ Dat$genes[ , gene] > 0 ]
  #subtree <- extract.clade(phy = Dat$tree, node = getMRCA(Dat$tree,tip = tips))
  
  if(which == 1){
    ntaxa <- nrow(genes)
    plot(tree,show.tip.label = FALSE,x.lim = gene.pos + 0.01)
    points(rep(meta.pos,ntaxa),1:ntaxa, pch = 19,
           col = c("brown","green")[as.numeric(Map$Classification)],
           cex = 0.5)
    points(rep(gene.pos,ntaxa),1:ntaxa, pch = 19,
           col = topo.colors(max(genes[,gene])+1)[ genes[,gene]+1 ],
           cex = 0.5)
  }else if(which == 2){
    dat <- Map
    dat$gene <- genes[,gene]
    p1 <- ggplot(dat,aes(x = Classification, y = ..count.., fill = factor(gene))) +
      geom_bar(position = "fill")
    p1
  }else if(which == 3){
    dat <- Map
    dat$gene <- genes[,gene]
    p1 <- ggplot(dat,aes(x = factor(gene), y = ..count.., fill = Classification)) +
      geom_bar(position = "fill") +
      scale_fill_manual(values = c('brown',"green"))
    p1
  }else{
    stop("ERROR: Not supported")
  }
}

#' Test all genes with phylogenetic independent contrasts
allgenes_pic <- function(tree,genes,class){
  #tree <- subtree
  #genes <- genes.sub
  #class <- Dat$Classification
  
  pic.class <- pic(x = as.numeric(class),phy = tree,scaled = TRUE)
  #hist(pic.class)
  
  Res <- NULL
  for(gene in 1:ncol(genes)){
    #gene <- 1
    #ftable(genes.sub[,gene] ~ Dat$Classification)
    #boxplot(genes.sub[,gene] ~ Dat$Classification)
    #plot(as.numeric(Dat$Classification),genes.sub[,gene])
    
    gene.id <- colnames(genes)[gene]
    pic.gene <- pic(x = genes[,gene], phy = tree, scaled = TRUE)
    
    #plot(pic.gene,pic.class, col = c("brown","green")[as.numeric(Dat$Classification)])
    m1 <- lm(pic.class ~ pic.gene)
    #m1.anova <- anova(m1)
    m1.sum <- summary(m1)
    
    
    if(is.na(coef(m1)["pic.gene"])){
      res <- data.frame(gene.id = gene.id,
                        Estimate = NA,
                        SE = NA,
                        t.value = NA,
                        p.value = NA)
    }else{
      res <- data.frame(gene.id = gene.id,
                        Estimate = m1.sum$coefficients["pic.gene",1],
                        SE = m1.sum$coefficients["pic.gene",2],
                        t.value = m1.sum$coefficients["pic.gene",3],
                        p.value = m1.sum$coefficients["pic.gene",4])
    }
    
    Res <- rbind(Res,res)
  }
  
  return(Res)
}

#' Get group data
#' 
#' Get subtree, informative genes and metadatada from group
get_group <- function(group, tree, genes, Map, classes = c("NPA","PA")){
  #classes <- c("NPA","PA")
  
  to_remove <- subset(Map,!(Classification %in% classes & Taxonomic_Assignment == group))
  to_remove <- as.character(to_remove$taxon_oid)
  # Get tree from group (only PA & NPA)
  #to_remove <- as.character(Map$taxon_oid[ Map$Taxonomic_Assignment != group | Map$Classification == "soil" ])
  subtree <- drop.tip(phy = tree,tip = to_remove)
  #plot(subtree)
  
  # Get Classification from group
  # Check if all tree labels are in the map
  if(all(tree$tip.label == intersect(tree$tip.label,Map$taxon_oid))){
    Dat <- subset(Map,taxon_oid %in% subtree$tip.label)
    row.names(Dat) <- Dat$taxon_oid
    Dat <- Dat[subtree$tip.label,]
    Dat <- droplevels(Dat)
    Dat$Classification <- factor(Dat$Classification, levels = classes)
    #plot(subtree,tip.color = c("brown","green")[as.numeric(Dat$Classification)])
  }else{
    stop("ERROR: Map and tree don't have the same taxons")
  }
  
  # Get genes (features) per genome
  if(all(subtree$tip.label == intersect(subtree$tip.label,row.names(genes)))){
    # Remove constant and non informative genes
    genes.sub <- genes[ subtree$tip.label, ]
    genes.sub <- genes.sub[ , !apply(genes.sub,2,function(x) all(x[1] == x)) ]
    #genes.sub <- genes.sub[ , apply(genes.sub,2,function(x) min(table(x)) > 1) ]
    genes.sub <- genes.sub[ ,colSums(genes.sub) >= 5 ]
    #dim(genes.sub)
    #gc()
  }else{
    stop("ERROR: Genes and subtree don't have the same taxons")
  }
  
  return(list(tree = subtree, Map = Dat, genes = genes.sub))
}
