library(ggplot2)

# Using genes identified as stable during iPS reprogramming in this paper
# https://www.nature.com/articles/s41598-018-26707-8

hk1 <- "Atp5f1"
hk2 <- "Pgk1"
hk3 <- "Gapdh"

counts <- readRDS("C:/Users/kim/Documents/codaDE/data/GTEx_data/parsed_GTEx.rds")
labels <- readRDS("C:/Users/kim/Documents/codaDE/data/GTEx_data/parsed_GTEx_org_by_tissue.rds")

hk1_idx <- which(counts$Description == toupper(hk1))
hk2_idx <- which(counts$Description == toupper(hk2))
hk3_idx <- which(counts$Description == toupper(hk3))

selected_tissue_types <- c(1, 3, 4, 18, 19, 38, 49, 55)
hk_idx <- hk3_idx
points <- c()
point_labels <- c()
for(tissue in selected_tissue_types) {
  s1 <- labels[[tissue]]
  subset_counts <- as.matrix(counts[,colnames(counts) %in% s1])
  if(ncol(subset_counts) > 0) {
    CPM_counts <- apply(subset_counts, 2, function(x) (x/sum(x))*1e6)
    points <- c(points, unname(unlist(CPM_counts[hk_idx,])))
    point_labels <- c(point_labels, rep(names(labels)[tissue], ncol(CPM_counts)))
  }
}

ggplot(data.frame(abundance = points, tissue = as.factor(point_labels)), aes(x = tissue, y = abundance)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
