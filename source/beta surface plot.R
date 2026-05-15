

library(ggplot2)
library(patchwork)  # for side-by-side layout

# ---- Build a single long data frame ----
df_plot <- data.frame(
  s          = grid_df$s_grid,
  age        = grid_df$age,
  true       = beta_true_grid,
  est        = surface_mean,
  bias       = surface_mean - beta_true_grid,
  in_ci      = (beta_true_grid >= surface_q025) & (beta_true_grid <= surface_q975)
)

# Shared colour limits for true vs estimated panels
zlim <- range(c(df_plot$true, df_plot$est))

# ---- True surface ----
p_true <- ggplot(df_plot, aes(x = s, y = age, fill = true)) +
  geom_tile() +
  scale_fill_viridis_c(limits = zlim, name = "β(s, a)") +
  labs(title = "True surface", x = "s", y = "Age") +
  theme_bw()

# ---- Posterior mean surface ----
p_est <- ggplot(df_plot, aes(x = s, y = age, fill = est)) +
  geom_tile() +
  scale_fill_viridis_c(limits = zlim, name = "β(s, a)") +
  labs(title = "Estimated surface (posterior mean)", x = "s", y = "Age") +
  theme_bw()

# ---- Bias surface ----
bias_abs <- max(abs(df_plot$bias))
p_bias <- ggplot(df_plot, aes(x = s, y = age, fill = bias)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 0, limits = c(-bias_abs, bias_abs),
    name = "Bias"
  ) +
  labs(title = "Bias (est − true)", x = "s", y = "Age") +
  theme_bw()

# ---- Pointwise coverage ----
p_cover <- ggplot(df_plot, aes(x = s, y = age, fill = in_ci)) +
  geom_tile() +
  scale_fill_manual(values = c("TRUE" = "seagreen3", "FALSE" = "tomato2"),
                    name = "In 95% CrI") +
  labs(title = "Pointwise 95% CrI coverage", x = "s", y = "Age") +
  theme_bw()

# ---- Combine ----
(p_true | p_est) / (p_bias | p_cover)