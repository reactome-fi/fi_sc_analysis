# The following R script is copied from the tutorial.ipynb hosted at https://github.com/theislab/dca/blob/master/tutorial.ipynb.
# This R script is used to generate some artificial file for software development.

setwd("/Users/wug/git/reactome-fi/fi_sc_analysis/RSrc") # replace with getwd() or path param

# make sure that splatter is installed: https://github.com/Oshlack/splatter
library(splatter)
# --------------------------------------------------------

simulate <- function(nGroups=2, nGenes=200, batchCells=2000, dropout=3)
{
    if (nGroups > 1) method <- 'groups'
    else             method <- 'single'
    
    group.prob <- rep(1, nGroups) / nGroups
    
    # new splatter requires dropout.type
    if ('dropout.type' %in% slotNames(newSplatParams())) {
        if (dropout)
            dropout.type <- 'experiment'
        else
            dropout.type <- 'none'
        cat(paste("With dropout.type: ", dropout.type, "\n", sep = ""))
        sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells,
                             dropout.type=dropout.type, method=method,
                             seed=42, dropout.shape=-1, dropout.mid=dropout)
        
    } else {
        sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells,
                             dropout.present=!dropout, method=method,
                             seed=42, dropout.shape=-1, dropout.mid=dropout)        
    }
    
    counts     <- as.data.frame(counts(sim))
    truecounts <- as.data.frame(assays(sim)$TrueCounts)
    
    dropout    <- assays(sim)$Dropout
    mode(dropout) <- 'integer'
    
    cellinfo   <- as.data.frame(colData(sim))
    geneinfo   <- as.data.frame(rowData(sim))
    
    list(counts=counts,
         cellinfo=cellinfo,
         geneinfo=geneinfo,
         truecounts=truecounts)
}

transpose <- function(in.file, out.file) {
    df <- read.delim(in.file, sep = "\t", header = TRUE)
    col.names <- names(df)
    row.names <- df[, 1]
    df.t <- t(df[, 2 : length(df)])
    colnames(df.t) <- row.names
    rownames(df.t) <- col.names[2:length(col.names)]
    write.csv(df.t, out.file, quote = FALSE)
}

# For two groups
group <- 2
dropout <- 3

# For 6 groups
group <- 6
dropout <- 1

data.dir <- "../data/simulated"
cellinfo.file <- paste(data.dir, paste(group, "cell_info.csv", sep = "_"), sep = "/")
geneinfo.file <- paste(data.dir, paste(group, "gene_info.csv", sep = "_"), sep = "/")
drop.file <- paste(data.dir, paste(group, "group_drop.csv", sep = "_"), sep = "/")
true.file <- paste(data.dir, paste(group, "group_true.csv", sep = "_"), sep = "/")
drop.file.t <- paste(data.dir, paste(group, "group_drop_t.csv", sep = "_"), sep = "/")
true.file.t <- paste(data.dir, paste(group, "group_true_t.csv", sep = "_"), sep = "/")

sim <- simulate(nGroups = group, dropout = dropout)
counts <- sim$counts
geneinfo <- sim$geneinfo
cellinfo <- sim$cellinfo
truecounts <- sim$truecounts
# Write out the data
cat("Export the data\n")
write.csv(counts, file = drop.file, quote = FALSE)
write.csv(truecounts, file = true.file, quote = FALSE)
write.csv(t(counts), file = drop.file.t, quote = FALSE)
write.csv(t(truecounts), file = true.file.t, quote = FALSE)
write.csv(cellinfo, file = cellinfo.file, quote = FALSE)
write.csv(geneinfo, file = geneinfo.file, quote = FALSE)