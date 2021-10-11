# Plot supplemental variable importance figure for the RF predictive models.

source("path_fix.R")

library(caret)

use_baseline <- "self"
partial_flag <- "nopartial"

# ------------------------------------------------------------------------------
#   Define feature lookup
# ------------------------------------------------------------------------------

mapping <- list(METHOD = "method",
                P = "number of features",
                PARTIAL = "totals expected to be partially informative (e.g. spike-in normalized)",
                TOTALS_C_FC = "absolute fold change in mean totals (A vs. B)",
                TOTALS_C_D = "absolute change in mean totals",
                TOTALS_C_MAX_D = "max delta in totals",
                TOTALS_C_MED_D = "median delta in totals",
                TOTALS_C_SD_D = "SD in totals",
                CORR_RA_MED = "median correlation of relative abundances",
                CORR_RA_SD = "SD correlation of relative abundances",
                CORR_RA_SKEW = "skew correlation of relative abundances",
                CORR_LOG_MED = "median correlation of log + PC counts",
                CORR_LOG_SD = "SD correlation of log + PC counts",
                CORR_LOG_SKEW = "skew correlation of log + PC counts",
                CORR_CLR_MED = "median correlation of CLR features",
                CORR_CLR_SD = "SD correlation of CLR features",
                CORR_CLR_SKEW = "skew correlation of CLR features",
                COMP_C_P0_A = "percent features == 0 in A",
                COMP_C_P0_B = "percent features == 0 in B",
                COMP_C_P1_A = "percent features == 1 in A",
                COMP_C_P1_B = "percent features == 1 in B",
                COMP_C_P5_A = "percent features <= 5 in A",
                COMP_C_P5_B = "percent features <= 5 in B",
                COMP_RA_P01_A = "percent features < 0.1% relative abundance in A",
                COMP_RA_P01_B = "percent features < 0.1% relative abundance in B",
                COMP_RA_P1_A = "percent features < 1% relative abundance in A",
                COMP_RA_P1_B = "percent features < 1% relative abundance in B",
                COMP_RA_P5_A = "percent features < 5% relative abundance in A",
                COMP_RA_P5_B = "percent features < 5% relative abundance in B",
                COMP_RA_MAX_A = "max relative abundance in A",
                COMP_RA_MED_A = "median relative abundance in A",
                COMP_RA_SD_A = "SD relative abundance in A",
                COMP_RA_SKEW_A = "skew relative abundance in A",
                COMP_RA_MAX_B = "max relative abundance in B",
                COMP_RA_MED_B = "median relative abundance in B",
                COMP_RA_SD_B = "SD relative abundance in B",
                COMP_RA_SKEW_B = "skew relative abundance in B",
                COMP_C_ENT_A = "entropy in A",
                COMP_C_ENT_B = "entropy in B",
                FW_RA_MAX_D = "max change in relative abundance",
                FW_RA_MED_D = "median change in relative abundance",
                FW_RA_SD_D = "SD change in relative abundance",
                FW_RA_PPOS_D = "percent features with + change in relative abundances",
                FW_RA_PNEG_D = "percent features with - change in relative abundances",
                FW_RA_PFC05_D = "percent features with < 0.5 FC in relative abundance",
                FW_RA_PFC1_D = "percent features with < 1 FC in relative abundance",
                FW_RA_PFC2_D = "percent features with < 2 FC in relative abundance",
                FW_LOG_MAX_D = "max change in log + PC counts",
                FW_LOG_MED_D = "median change in log + PC counts",
                FW_LOG_SD_D = "SD change in log + PC counts",
                FW_LOG_PPOS_D = "percent features with + change in log + PC counts",
                FW_LOG_PNEG_D = "percent features with - change in log + PC counts",
                FW_LOG_PFC05_D = "percent features with < 0.5 FC in log + PC counts",
                FW_LOG_PFC1_D = "percent features with < 1 FC in log + PC counts",
                FW_LOG_PFC2_D = "percent features with < 2 FC in log + PC counts",
                FW_CLR_MAX_D = "max change in CLR",
                FW_CLR_MED_D = "median change in CLR",
                FW_CLR_SD_D = "SD change in CLR",
                FW_CLR_PPOS_D = "percent features with + change in CLR",
                FW_CLR_PNEG_D = "percent features with - change in CLR",
                FW_CLR_PFC05_D = "percent features with < 0.5 FC in CLR",
                FW_CLR_PFC1_D = "percent features with < 1 FC in CLR",
                FW_CLR_PFC2_D = "percent features with < 2 FC in CLR")

for(task in c("classification", "regression")) {
  for(type in c("TPR", "FPR")) {
    for(DE_method in c("ALDEx2", "DESeq2", "scran")) {
      save_fn <- file.path("output",
                           "predictive_fits",
                           paste0(use_baseline, "_", partial_flag),
                           task,
                           paste0(use_baseline, "_", type, "_", DE_method, "_variable-importance.rds"))
      if(file.exists(save_fn)) {
        plot_df <- readRDS(save_fn)
      } else {
        stop("Variable importance results not found!")
      }
      
      # Fix data.frame
      plot_df$feature <- rownames(plot_df)
      rownames(plot_df) <- NULL
      colnames(plot_df) <- c("weight", "feature")
      
      # Map features to readable names
      relabel <- suppressMessages(data.frame(feature = plot_df$feature) %>%
        left_join(data.frame(feature = names(mapping), string = unname(unlist(mapping)))))
      
      plot_df$label <- relabel$string
      
      # Truncate by relative feature weight
      plot_df <- plot_df %>%
        arrange(desc(weight)) %>%
        mutate(rel_weight = weight/sum(plot_df$weight)) #%>%
        # filter(rel_weight > 0.01)
      plot_df <- plot_df[1:10,]
      
      ggplot(plot_df, aes(x = weight, y = reorder(label, weight))) +
        geom_bar(stat = "identity") +
        theme_bw() +
        labs(x = "feature weight",
             y = paste0("feature in ", ifelse(type == "TPR", "sensitivity", "specificity"), " model")) +
        theme(axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
              axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))
      ggsave(file.path("output",
                       "images",
                       paste0(task, "_varImp_", type, "-", DE_method, ".png")),
             units = "in",
             dpi = 100,
             # height = ifelse(type == "TPR", 3, 6),
             height = 3.5,
             width = 7)
    }
  }
}
