# =============================================================
# Visualizer for:  Y = b0 + b1*X + b2*(X - c)*W + e
#
# Key behaviour:
#   The term b2*(X - c)*W is ZERO whenever X = c, regardless of W.
#   => Every line (one per value of W) PIVOTS through the same point
#      at X = c, where Y = b0 + b1*c. This is the "inflection" point.
#   For X > c, larger W tilts the slope up (if b2>0); for X < c it
#   tilts the other way -> a fan / pivot pattern centred on X = c.
# =============================================================

# install.packages(c("ggplot2", "dplyr", "tidyr"))  # if needed
library(ggplot2)
library(dplyr)
library(tidyr)

# ---- The model (deterministic part; e is noise) -------------
model_Y <- function(X, c, W, b0 = 0, b1 = 1, b2 = 0.5, e = 0) {
  b0 + b1 * X + b2 * (X - c) * W + e
}

# ---- Parameters you can play with ---------------------------
b0 <- 0      # intercept
b1 <- 1      # main effect of X
b2 <- 0.6    # interaction strength of (X - c) and W

X_seq  <- seq(-5, 5, length.out = 400)
W_vals <- c(-2, -1, 0, 1, 2)     # moderator levels (lines)
c_vals <- c(-2, 0, 2)            # pivot locations (panels)

# ---- Build a grid and evaluate ------------------------------
grid <- expand.grid(X = X_seq, W = W_vals, c = c_vals) |>
  mutate(
    Y      = model_Y(X, c, W, b0, b1, b2),
    W      = factor(W),
    c_lab  = factor(paste0("c = ", c), levels = paste0("c = ", c_vals))
  )

# Pivot points: where X = c -> Y = b0 + b1*c (same for all W)
pivots <- data.frame(c = c_vals) |>
  mutate(
    X     = c,
    Y     = b0 + b1 * c,
    c_lab = factor(paste0("c = ", c), levels = paste0("c = ", c_vals))
  )

# =============================================================
# PLOT 1: Y vs X, lines = W, panels = c. Vertical line at X = c.
# =============================================================
p1 <- ggplot(grid, aes(X, Y, colour = W, group = W)) +
  geom_line(linewidth = 1) +
  # the inflection: vertical guide at X = c
  geom_vline(data = pivots, aes(xintercept = X),
             linetype = "dashed", colour = "grey40") +
  # the single point through which all W-lines pass
  geom_point(data = pivots, aes(X, Y),
             inherit.aes = FALSE, size = 3, colour = "black") +
  geom_text(data = pivots, aes(X, Y, label = "X = c"),
            inherit.aes = FALSE, vjust = -1, hjust = -0.1, size = 3.5) +
  facet_wrap(~ c_lab) +
  scale_colour_viridis_d(name = "W") +
  labs(
    title    = expression(Y == b[0] + b[1]*X + b[2]*(X - c)*W + e),
    subtitle = sprintf("b0 = %.1f, b1 = %.1f, b2 = %.1f. Lines pivot at X = c.",
                       b0, b1, b2),
    x = "X", y = "Y"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

print(p1)
ggsave("interaction_pivot_by_c.png", p1, width = 11, height = 4, dpi = 150)

# =============================================================
# PLOT 2: Effect of varying b2 (interaction strength), fixed c = 0
# Shows how the "fan" opens up around the pivot.
# =============================================================
b2_vals <- c(0, 0.3, 0.6, 1.0)
grid2 <- expand.grid(X = X_seq, W = W_vals, b2 = b2_vals) |>
  mutate(
    Y     = model_Y(X, c = 0, W = W, b0 = b0, b1 = b1, b2 = b2),
    W     = factor(W),
    b2lab = factor(paste0("b2 = ", b2), levels = paste0("b2 = ", b2_vals))
  )

p2 <- ggplot(grid2, aes(X, Y, colour = W, group = W)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ b2lab) +
  scale_colour_viridis_d(name = "W") +
  labs(
    title    = "Interaction strength b2 controls the fan around X = c (c = 0)",
    subtitle = "b2 = 0 -> all W lines collapse onto one line (no moderation)",
    x = "X", y = "Y"
  ) +
  theme_minimal(base_size = 13)

print(p2)
ggsave("interaction_vary_b2.png", p2, width = 10, height = 7, dpi = 150)

# =============================================================
# PLOT 3 (optional): a single scattered realisation with noise e
# =============================================================
set.seed(1)
n   <- 600
dat <- data.frame(
  X = runif(n, -5, 5),
  W = sample(W_vals, n, replace = TRUE)
) |>
  mutate(
    c = 0,
    Y = model_Y(X, c, W, b0, b1, b2, e = rnorm(n, 0, 2)),
    W = factor(W)
  )

p3 <- ggplot(dat, aes(X, Y, colour = W)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_viridis_d(name = "W") +
  labs(title = "Simulated data with noise e ~ N(0, 2), c = 0",
       x = "X", y = "Y") +
  theme_minimal(base_size = 13)

print(p3)
ggsave("interaction_simulated.png", p3, width = 8, height = 5, dpi = 150)

cat("Saved: interaction_pivot_by_c.png, interaction_vary_b2.png, interaction_simulated.png\n")
