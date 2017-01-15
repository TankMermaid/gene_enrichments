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
#library(glmnet)
#library(phylolm)
library(phytools)
source("functions.r")

load(file = "Burkholderiales_dataset.rdat")
#### Pagels test ####
date()
pagel.res <- allgenes_pagel(tree = Dat$tree, Map = Dat$Map, genes = Dat$genes)
write.table(pagel.res, file = "pagel.res.txt", sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = TRUE)
date()
