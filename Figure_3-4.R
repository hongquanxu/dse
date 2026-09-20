## Figures 3 and 4
## Ying and Xu (2026+). Efficient Representation and Construction of Space-Filling Designs via a Double Stratification Enumerator. 

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(UniDOE)   
library(MaxPro)
source("dse_functions.R")
source("generate_soa.R")
source("dse_tables.R")

###################
## generate_table(): from dse_tables.R that generate the y and z values using Theorem 6
## get.soa.new(): from generate_soa.R that generates SOA
## main_p2(): from dse_tables.R that takes a design and return SOA, Alg, OA
## dse(): from dse_functions.R that calculates the double stratification enumerator value
###################

## GSOAs used in Figures 3 and 4
gsoa32x9 <- get.soa.new(N = 32, k = 9) %/% 2 # GSOA(32, 9, 2^2, 3)
gsoa81x9 <- read.csv("soa81x9.csv")[,-1] # GSOA(81,9, 3^2, 3)

##

make_panel_plot <- function(k, design, s, p, y, z, crit, n, alg_name = "GenUD") {
  if (k != 9) {
    comb <- combinat::combn(9, k)
    values <- unlist(parallel::mclapply(
      seq_len(ncol(comb)),
      function(i) dse(design[, comb[, i]], s, p, y, z),
      mc.cores = parallel::detectCores()
    ))
    design_new <- design[, comb[, which.min(values)]]
  } else {
    design_new <- design
  }
  
  res <- main_p2(design_new, s, crit = crit)
  df <- data.frame(S2LP = res[,1], ALG = res[,2], ALP = res[,3])
  names(df)[2] <- alg_name
  
  df %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
    mutate(Variable = ifelse(Variable == "S2LP", " S2LP", Variable)) %>%
    ggplot(aes(x = Variable, y = Value, fill = Variable)) +
    geom_boxplot() +
    ggtitle(paste0(n, "x", k)) +
    xlab("") +
    ylab("") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = "white")
    )
}

make_panel_grid <- function(design, crit, s, p, n, y, z, seed, alg_name = "GenUD") {
  set.seed(seed)
  plot <- lapply(3:9, make_panel_plot, design = design, s = s, p = p, y = y, z = z, 
                 crit = crit, n = n, alg_name = alg_name)
  wrap_plots(plot, nrow = 2, ncol = 4) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    labs(fill = NULL)
}

## CD config
final_res_cd <- generate_table(s = 2, p = 2)
plot_cd <- make_panel_grid(
  design = gsoa32x9,
  crit = "CD2", s = 2, p = 2, n = 32,
  y = as.numeric(final_res_cd[1, 4]), ## 9/17 changed 
  z = as.numeric(final_res_cd[1, 5]), ## 9/17 changed
  seed = 99804
)

## MaxPro config
final_res_mp <- generate_table(s = 3, p = 2)
plot_maxpro <- make_panel_grid(
  design = gsoa81x9,
  crit = "MaxPro", s = 3, p = 2, n = 81,
  y = as.numeric(final_res_mp[5, 4]), ## 9/17 changed this
  z = as.numeric(final_res_mp[5, 5]), ## 9/17 changed this
  seed = 804,
  alg_name = "MaxProQQ"
)

## make Figure 3
fig_3 <- plot_cd
print(fig_3)

## make Figure 4
fig_4 <- plot_maxpro
print(fig_4)
