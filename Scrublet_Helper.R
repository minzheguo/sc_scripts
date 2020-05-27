
# export UMIs matrix from RNA assay of Seurat3 object for Scrublet
exportCounts <- function(object, group.by, rm.nonexpressed=T, min.exp=0) {
  groups <- sort(as.character(unique(object@meta.data[, group.by])))
  cat("Exporting count matrix of ", length(groups), " groups: ", paste0(groups, collapse = ","), "...\n",   sep="")
  for (gi in groups) {
    cat(gi, "...", sep="")
    gi.cells <- rownames(object@meta.data)[which(object@meta.data[, group.by]==gi)]
    gi.counts <- as.matrix(object@assays$RNA@counts[, gi.cells])
    if (rm.nonexpressed == T) {
      gi.counts <- gi.counts[which(rowSums(gi.counts>min.exp)>0), ]
    }
    cat(dim(gi.counts)[1], " features X ", dim(gi.counts)[2], " cells ... ", sep="")
    gi.counts <- data.frame(Gene=rownames(gi.counts), gi.counts, check.names = F)
    write.table(gi.counts, file=paste(gi,".counts.txt", sep=""), sep = "\t", col.names = T, row.names = F, quote=FALSE)
    cat("done\n")
  }
}
