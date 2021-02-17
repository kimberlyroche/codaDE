# Parse yeast and Ciona experiment data from Athanasiadou et al. paper
# This involves normalizing counts as per their code (mostly copied here)

file_dir <- "C:/Users/kim/Documents/codaDE/data/Athanasiadou_2021/S1CodeandData"

parse_dataset <- function(SI_file, RNA_file, ERCC_annot, groups) {
  # Spike-in counts and RNA counts
  YmatSI <- as.matrix(read.table(file = file.path(file_dir, SI_file)))
  YmatRNA <- as.matrix(read.table(file = file.path(file_dir, RNA_file)))
  # head(YmatRNA)
  
  # Spike-in amounts added (in attomoles
  ERCC <- read.table(file = file.path(file_dir, ERCC_annot))
  ERCC <- ERCC[rownames(YmatSI),]
  
  # Spike-in library sizes
  LibSizeSI <- colSums(YmatSI)
  
  # RNA library sizes
  LibSizeRNA <- colSums(YmatRNA)
  
  Ratio <- LibSizeSI/LibSizeRNA
  # data.frame(LibSizeSI, LibSizeRNA, Ratio) # print
  
  # Perform the normalization as per RMarkdown file
  # Total amol of each spike-in
  Nvec <- ERCC[,"attomoles"]
  # In sample to be sequenced
  names(Nvec) <- rownames(ERCC)
  AttomoleToMoleculesPerCell <- (1.e-18)*(6.22e23)*(1.e-7)
  MoleculesPerCell <- Nvec*AttomoleToMoleculesPerCell
  names(MoleculesPerCell) <- names(Nvec)
  
  # From RMarkdown comments: "Compute spike-in proportions. Choose as the reference spike-in,
  #   the one with the largest proportion; record its proportion and attomoles; and compute the
  # library-specific $\nu_j$ calibration factors"
  f.vec <- rowSums(YmatSI)/sum(YmatSI) # fraction of total spike-in
  # Counts accounted for by each spike-in
  INDmax <- which(f.vec == max(f.vec)) # index of spike-in
  f.ref <- f.vec[INDmax] # fraction (empirical proportion) for reference
  n.ref <- Nvec[names(f.vec[INDmax])] # corresponding amol
  nu <- (f.ref/n.ref)*LibSizeSI # for counts to nominal amol conversion
  Zmat <- YmatRNA%*%diag(1/nu) # normalization giving amol
  colnames(Zmat) <- colnames(YmatRNA)
  return(Zmat)
}

counts <- parse_dataset(SI_file = "YmatSIyeast.txt",
                        RNA_file = "YmatRNAyeast.txt",
                        ERCC_annot = "ERCCyeast.txt",
                        groups = as.factor(rep(c("C12", "C20", "C30"), rep(3,3))))
saveRDS(list(counts = counts, groups = groups), file = file.path(file_dir, "yeast_parsed.rds"))

counts <- parse_dataset(SI_file = "YmatSIciona.txt",
                        RNA_file = "YmatRNAciona.txt",
                        ERCC_annot = "ERCCciona.txt",
                        groups = as.factor(rep(c("lacz", "dnfgfr", "camras"), rep(3,3))))
saveRDS(list(counts = counts, groups = groups), file = file.path(file_dir, "ciona_parsed.rds"))
