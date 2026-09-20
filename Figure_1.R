## Figure 1 
## Ying and Xu (2026+). Efficient Representation and Construction of Space-Filling Designs via a Double Stratification Enumerator. 

library(ggplot2)
library(dplyr)


######################
## 2D Projection Plots
######################
df_soa <- as.data.frame(matrix(c(
  7,7, 3,6, 6,3, 2,2, 5,5, 1,4, 4,1, 0,0), ncol = 2, byrow = TRUE))
df_oa <- as.data.frame(matrix(c(
  7,7, 3,4, 5,3, 1,0, 6,6, 2,5, 4,2, 0,1), ncol = 2, byrow = TRUE))

## Setting up the grids
grid_specs <- list(
  "2x2" = list(x = seq(-0.5, 7.5, by = 4), y = seq(-0.5, 7.5, by = 4)),
  "2x4" = list(x = seq(-0.5, 7.5, by = 4), y = seq(-0.5, 7.5, by = 2)),
  "4x2" = list(x = seq(-0.5, 7.5, by = 2), y = seq(-0.5, 7.5, by = 4)))

panel_levels <- c("OA_2x2","OA_2x4","OA_4x2","SOA_2x2","SOA_2x4","SOA_4x2")

## Build point data: OA-based design in top row (3 panels), SOA in bottom row (3 panels)
plot_df <- bind_rows(df_oa  %>% mutate(panel = "OA_2x2"),
  df_oa  %>% mutate(panel = "OA_2x4"),
  df_oa  %>% mutate(panel = "OA_4x2"),
  df_soa %>% mutate(panel = "SOA_2x2"),
  df_soa %>% mutate(panel = "SOA_2x4"),
  df_soa %>% mutate(panel = "SOA_4x2")) %>%
  mutate(panel = factor(panel, levels = panel_levels))

## Grid lines for each panel
vline_df <- bind_rows(
  data.frame(panel = "OA_2x2",  xintercept = grid_specs[["2x2"]]$x),
  data.frame(panel = "OA_2x4",  xintercept = grid_specs[["2x4"]]$x),
  data.frame(panel = "OA_4x2",  xintercept = grid_specs[["4x2"]]$x),
  data.frame(panel = "SOA_2x2", xintercept = grid_specs[["2x2"]]$x),
  data.frame(panel = "SOA_2x4", xintercept = grid_specs[["2x4"]]$x),
  data.frame(panel = "SOA_4x2", xintercept = grid_specs[["4x2"]]$x)) %>%
  mutate(panel = factor(panel, levels = panel_levels))

hline_df <- bind_rows(
  data.frame(panel = "OA_2x2",  yintercept = grid_specs[["2x2"]]$y),
  data.frame(panel = "OA_2x4",  yintercept = grid_specs[["2x4"]]$y),
  data.frame(panel = "OA_4x2",  yintercept = grid_specs[["4x2"]]$y),
  data.frame(panel = "SOA_2x2", yintercept = grid_specs[["2x2"]]$y),
  data.frame(panel = "SOA_2x4", yintercept = grid_specs[["2x4"]]$y),
  data.frame(panel = "SOA_4x2", yintercept = grid_specs[["4x2"]]$y)) %>%
  mutate(panel = factor(panel, levels = panel_levels))

## checkerboard shadowing for OA-based design projection
make_checkerboard <- function(panel_id, xbreaks, ybreaks) {
  xmins <- xbreaks[-length(xbreaks)]
  xmaxs <- xbreaks[-1]
  ymins <- ybreaks[-length(ybreaks)]
  ymaxs <- ybreaks[-1]
  grid <- expand.grid(col = seq_along(xmins), row = seq_along(ymins))
  grid$xmin <- xmins[grid$col]
  grid$xmax <- xmaxs[grid$col]
  grid$ymin <- ymins[grid$row]
  grid$ymax <- ymaxs[grid$row]
  grid <- grid[(grid$row + grid$col) %% 2 == 1, ]
  grid$panel <- panel_id
  grid
}

## Checkerboard shading: only for OA's 2x4 and 4x2 panels (where OA fails stratification)
rect_df <- bind_rows(
  make_checkerboard("OA_2x4", grid_specs[["2x4"]]$x, grid_specs[["2x4"]]$y),
  make_checkerboard("OA_4x2", grid_specs[["4x2"]]$x, grid_specs[["4x2"]]$y)) %>%
  mutate(panel = factor(panel, levels = panel_levels))

fig_1 <- ggplot(plot_df, aes(V1, V2)) +
  geom_rect(
    data = rect_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey80") +
  geom_point(size = 1) +
  geom_vline(data = vline_df, aes(xintercept = xintercept), linetype = "dashed", color = "grey50") +
  geom_hline(data = hline_df, aes(yintercept = yintercept), linetype = "dashed", color = "grey50") +
  facet_wrap(~panel, nrow = 2, ncol = 3) +
  coord_fixed(xlim = c(-0.5, 7.5), ylim = c(-0.5, 7.5)) +
  scale_x_continuous(breaks = 0:7) +
  scale_y_continuous(breaks = 0:7) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = NULL, y = NULL)

## display Figure 1
print(fig_1)

