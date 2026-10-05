library(tidyverse)
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(ggrepel)
library(ggpubr)

source("plot_parameters_TRiX.R")

# Read in ASM and ASE data frames. Make TRiXi data frame.
# start with ASM.
suppl_table_4 <- fread("Supplementary_tables/Supplementary_Table4.tsv")

# calculate skewing function
compute_folded_skew <- function(df) {
  
  k <- df$H1_Xa + df$H2_Xi
  n <- df$H1_Xa + df$H2_Xa + df$H1_Xi + df$H2_Xi
  
  logLik_fun <- function(P) {
    
    if(P <= 0 || P >= 0.5) return(-Inf)
    
    log_term1 <- dbinom(k, n, P, log=TRUE)
    log_term2 <- dbinom(k, n, 1-P, log=TRUE)
    
    sum(log(exp(log_term1) + exp(log_term2)))
  }
  
  # MLE
  opt <- optimize(logLik_fun, interval=c(1e-8, 0.5-1e-8), maximum=TRUE)
  P_hat <- opt$maximum
  logLik_max <- opt$objective
  
  # Likelihood cutoff
  cutoff <- logLik_max - 0.5 * qchisq(0.95, df=1)
  
  # CI bounds
  #lower <- uniroot(function(P) logLik_fun(P) - cutoff,
  #                 interval=c(1e-8, P_hat))$root
  
  #upper <- uniroot(function(P) logLik_fun(P) - cutoff,
  #                 interval=c(P_hat, 0.5-1e-8))$root
  
  # Scaling
  skew_scaled      <- 100 - 100 * P_hat
  #skew_scaled_low  <- 100 - 100 * upper
  #skew_scaled_high <- 100 - 100 * lower
  
  # Return as data frame
  data.frame(
    P_hat = P_hat,
    #CI_lower = lower,
    #CI_upper = upper,
    skew_scaled = skew_scaled#,
    #skew_scaled_CI_lower = skew_scaled_low,
    #skew_scaled_CI_upper = skew_scaled_high
  )
}

# calc skewing
ASM_df <- data.frame()

for (i in unique(suppl_table_4$sample)){
  
  skew <- suppl_table_4[suppl_table_4$sample == i,]
  
P_folded <- compute_folded_skew(skew)

temp_df <- data.frame(P_folded, i)

ASM_df <- rbind(ASM_df, temp_df)

}

# change colnames
colnames(ASM_df) <- c("P_hat", "value",      "sample")

# make TRiXi df
trixi_results <- fread("data/revision_new_experiment/TRiXi/method_comparison_samples.tsv")

# Read in Supplementary Table 5
ASE_df_filt <- fread("Supplementary_tables/supplementary_table5.tsv")

# summarize (essentially calculate skewing per individual (the median is skewing)).
ASE_df_filt_calc <- ASE_df_filt %>% dplyr::group_by(participant) %>% rstatix::get_summary_stats(scaled, type="common")

# merge data frames.
# prepare TRiXi, ASE and ASM data frames
trixi_results <- trixi_results[,c("sample", "mean_skewing")]
colnames(trixi_results) <- c("sample", "value")

ASE_df_filt_calc_short <- ASE_df_filt_calc[,c("participant", "median")]
colnames(ASE_df_filt_calc_short) <- c("sample", "value")

# add variable
trixi_results$variable <- "TRiXi"
ASE_df_filt_calc_short$variable <- "ASE"
ASM_df$variable <- "ASM"

# Make final df
combined_df <- rbind(ASM_df[ASM_df$sample %in% c("xxh006", "xxh014", "xxh020", "xxh027", "xxh034", "xxh176"), -1], ASE_df_filt_calc_short, trixi_results)

# export source data
write.table(combined_df, "source_data/figure_1H_1I.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

# read in source data
figure_1H_1I_source_data <- fread("source_data/figure_1H_1I.tsv")

######## Figure 1H ###########
ggsave2(filename = "01_plots/Figure_1H.pdf",
ggplot(figure_1H_1I_source_data,
       aes(x=factor(sample,
                    levels=c("xxh020", "xxh176", "xxh014", "xxh006", "xxh034", "xxh027")),
           y=value)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(col=variable))+
  geom_text_repel(aes(label=variable), min.segment.length = 0.0000001) +
  theme_AL_box(legend.position="none") +
  labs(x="", y="skewing (%)")
)


########### Figure 1i ############
# long to wide
figure_1I_source_data_wide <- dcast(figure_1H_1I_source_data, sample ~ variable, value.var = "value")

# ASM vs TRiXi
p1 <- ggscatter(figure_1I_source_data_wide, x = "TRiXi", y = "ASM",
                add = "reg.line", conf.int = TRUE) +
  stat_cor(method = "pearson", label.x = min(figure_1I_source_data_wide$ASM), label.y = max(figure_1I_source_data_wide$TRiXi)) +
  ggtitle("TRiXi vs ASM")+ coord_cartesian(xlim=c(50,90), ylim=c(50,100)) + geom_hline(yintercept = c(75, 90), lty=2) + geom_vline(xintercept = c(75, 90), lty=2)

# ASE vs TRiXi
p2 <- ggscatter(figure_1I_source_data_wide, x = "TRiXi", y = "ASE",
                add = "reg.line", conf.int = TRUE) +
  stat_cor(method = "pearson", label.x = min(figure_1I_source_data_wide$ASE), label.y = max(figure_1I_source_data_wide$TRiXi)) +
  ggtitle("TRiXi vs ASE") + coord_cartesian(xlim=c(50,90), ylim=c(50,100)) + geom_hline(yintercept = c(75, 90), lty=2) + geom_vline(xintercept = c(75, 90), lty=2)

# ASM vs ASE
p3 <- ggscatter(figure_1I_source_data_wide, x = "ASM", y = "ASE",
                add = "reg.line", conf.int = TRUE) +
  stat_cor(method = "pearson", label.x = min(figure_1I_source_data_wide$ASM), label.y = max(figure_1I_source_data_wide$ASE)) +
  ggtitle("ASM vs ASE")+ coord_cartesian(xlim=c(50,90), ylim=c(50,100)) + geom_hline(yintercept = c(75, 90), lty=2) + geom_vline(xintercept = c(75, 90), lty=2)


# combine and plot
ggsave2(filename = "01_plots/Figure_1I.pdf",
        plot_grid(p1, p2, p3, ncol = 3)
)



# get metrics for method comparisons
figure_1H_1I_source_data_diff_vs_ASM <- merge(figure_1H_1I_source_data[figure_1H_1I_source_data$variable == "TRiXi",], figure_1H_1I_source_data[figure_1H_1I_source_data$variable == "ASM",], by = "sample")
figure_1H_1I_source_data_diff_vs_ASM$diff <- abs(figure_1H_1I_source_data_diff_vs_ASM$value.x - figure_1H_1I_source_data_diff_vs_ASM$value.y)

figure_1H_1I_source_data_diff_vs_ASE <- merge(figure_1H_1I_source_data[figure_1H_1I_source_data$variable == "TRiXi",], figure_1H_1I_source_data[figure_1H_1I_source_data$variable == "ASE",], by = "sample")
figure_1H_1I_source_data_diff_vs_ASE$diff <- abs(figure_1H_1I_source_data_diff_vs_ASE$value.x - figure_1H_1I_source_data_diff_vs_ASE$value.y)

figure_1H_1I_source_data_diff_ASE_vs_ASM <- merge(figure_1H_1I_source_data[figure_1H_1I_source_data$variable == "ASE",], figure_1H_1I_source_data[figure_1H_1I_source_data$variable == "ASM",], by = "sample")
figure_1H_1I_source_data_diff_ASE_vs_ASM$diff <- abs(figure_1H_1I_source_data_diff_ASE_vs_ASM$value.x - figure_1H_1I_source_data_diff_ASE_vs_ASM$value.y)

mean(figure_1H_1I_source_data_diff_vs_ASM$diff)
mean(figure_1H_1I_source_data_diff_vs_ASE$diff)
mean(figure_1H_1I_source_data_diff_ASE_vs_ASM$diff)


