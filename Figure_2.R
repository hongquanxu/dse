## Figure 2
## Ying and Xu (2026+). Efficient Representation and Construction of Space-Filling Designs via a Double Stratification Enumerator. 

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

source("generate_soa.R")
source("dse_tables.R")

#############
## generate_table(): from dse_tables.R that generate the y and z values using Theorem 6
## get.soa.new(): from generate_soa.R that generates SOA
## main_p2(): from dse_tables.R that takes a design and return SOA, OA and Algorithm
#############
final_res <- generate_table(s = 2, p = 2) 
y <- as.numeric(final_res[1,4]) # y for CD2, updated on 9/17/26
z <- as.numeric(final_res[1,5]) # z for CD2
s <- 2

gsoa32x9 <- get.soa.new(N = 32, k = 9) %/% 2 # GSOA(32, 9, 2^2, 3)
design <- gsoa32x9
index <- c(34,49,11,84)
comb <- combinat::combn(9,6)

col_names <- paste0("D", rep(1:4, each = 2), c("_sp", "_all"))
CD = MP = setNames(as.data.frame(matrix(0, nrow = 10000, ncol = 8)), col_names)

set.seed(19080499)
for (i in 1:4) {
  design_new <- design[, comb[, index[i]]]
  cols <- (2 * i - 1):(2 * i)
  
  CD[, cols] <- main_p2(design_new, s, crit = "CD2")[,c(1,3)]
  MP[, cols] <- main_p2(design_new, s, crit = "MaxPro")[, c(1, 3)]
}

## ggplot theme
jrssb_theme <- function() {
  theme_bw(base_size = 11) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "white", color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

## boxplot for Figure 2
make_box_plot <- function(data, ylabel) {
  data %>%
    pivot_longer(cols = everything(), names_to = "Group", values_to = "Value") %>%
    separate(Group, into = c("Design", "Permute"), sep = "_", remove = FALSE) %>%
    mutate(
      Design = as.factor(Design),
      Permute = as.factor(ifelse(Permute == "sp", "(s,2)-level permutation", "all level permutation"))
    ) %>%
    ggplot(aes(x = Design, y = Value, fill = Permute)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.6) +
    xlab("Design") +
    ylab(ylabel) +
    scale_fill_manual(values = c(
      "(s,2)-level permutation" = "#F8766D",
      "all level permutation"   = "#00BA38"
    )) +
    labs(fill = NULL) +
    jrssb_theme()
}


plot_fig2 <- list(
  make_box_plot(CD, "CD"),
  make_box_plot(MP, "MaxPro")
)


fig_2 <- wrap_plots(plot_fig2, nrow = 1, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
  
## display Figure 2
print(fig_2)

## Examples.R : Table 2/Example 4; Tables S1, Example S1-S2
## Fig S1 (dse-simulation.R)
