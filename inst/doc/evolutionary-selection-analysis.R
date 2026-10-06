## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  eval = TRUE, echo = TRUE, message = FALSE, warning = FALSE,
  collapse = TRUE, comment = "#>",
  fig.width = 6.5, fig.height = 4.5, dpi = 72, out.width = "100%"
)
set.seed(1)
library(lande)

## ----echo=FALSE, out.width="100%", fig.cap="Function dependencies. Drawn from workflow.dot with Graphviz."----
knitr::include_graphics("workflow.png")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(lande)

## ----message=FALSE, warning=FALSE---------------------------------------------
set.seed(42)

raw_data <- data.frame(
  survival = rbinom(100, 1, 0.6),
  size = rnorm(100, 10, 2),
  colour = rnorm(100, 5, 1)
)

prepared <- prepare_selection_data(
  data = raw_data,
  fitness_col = "survival",
  trait_cols = c("size", "colour"),
  standardize = TRUE,
  add_relative = TRUE,
  na_action = "drop"
)

cat("size:", nrow(prepared), "\n")
head(prepared)

## ----message=FALSE, warning=FALSE---------------------------------------------
# Binary fitness
fitness_binary <- c(1, 0, 1, 1, 0, 1)
detect_family(fitness_binary)

# Continuous fitness
fitness_cont <- c(2.3, 4.1, 3.2, 5.6, 1.8)
detect_family(fitness_cont)

## ----message=FALSE, warning=FALSE---------------------------------------------
# Calculate selection differential for size
S_size <- selection_differential(
  data = prepared,
  fitness_col = "relative_fitness",
  trait_col = "size",
  standardized = TRUE,
  use_relative = TRUE
)

# Calculate for colour
S_colour <- selection_differential(
  data = prepared,
  fitness_col = "relative_fitness",
  trait_col = "colour",
  standardized = TRUE,
  use_relative = TRUE
)

cat("Selection Differentials (S):\n")
cat("Size :", round(S_size, 4), "\n")
cat("Colour:", round(S_colour, 4), "\n")

## ----message=FALSE, warning=FALSE---------------------------------------------
linear_results <- analyze_linear_selection(
  data = prepared,
  fitness_col = "relative_fitness",
  trait_cols = c("size", "colour"),
  fitness_type = "continuous"
)

# View coefficients
summary(linear_results$model)

## ----message=FALSE, warning=FALSE---------------------------------------------
nonlinear_results <- analyze_nonlinear_selection(
  data = prepared,
  fitness_col = "relative_fitness",
  trait_cols = c("size", "colour"),
  fitness_type = "continuous"
)

# Extract quadratic term for size (multiply by 2 for gamma)
coef(summary(nonlinear_results$model))["I(size^2)", ]

## ----message=FALSE, warning=FALSE---------------------------------------------
disruptive_test <- analyze_disruptive_selection(
  data = prepared,
  fitness_col = "relative_fitness",
  trait_col = "size",
  fitness_type = "continuous"
)

print(disruptive_test)

## ----message=FALSE, warning=FALSE, cache=FALSE--------------------------------
# k is capped from the number of distinct trait values; the default is fine here
spline_fit <- univariate_spline(
  data = prepared,
  fitness_col = "survival",
  trait_col = "size",
  fitness_type = "binary",
  bootstrap = TRUE
)

head(spline_fit$grid)

## ----message=FALSE, warning=FALSE, fig.align="center", out.width="80%", fig.cap="Figure 1. Correlated fitness function for body size"----
p <- plot_univariate_fitness(
  spline_fit,
  "size",
  title = "Correlated Fitness Function: Body Size"
)
print(p)

## ----message=FALSE, warning=FALSE---------------------------------------------
surface <- correlated_fitness_surface(
  data = prepared,
  fitness_col = "survival",
  trait_cols = c("size", "colour"),
  method = "gam"
)

# Find optimum (maximum fitness)
optimum <- surface$grid[which.max(surface$grid$.fit), ]
print(optimum)

## ----message=FALSE, warning=FALSE, fig.width=8, fig.height=6, fig.cap="Correlated fitness surface for size and colour"----
p <- plot_correlated_fitness(
  tps = surface,
  trait_cols = c("size", "colour"),
  bins = 12
)

print(p)

## ----message=FALSE, warning=FALSE, fig.width=8, fig.height=6, fig.cap="Correlated fitness surface with original data points"----
p <- plot_correlated_fitness_enhanced(
  tps = surface,
  trait_cols = c("size", "colour"),
  original_data = prepared,
  fitness_col = "survival",
  bins = 12
)

print(p)

## ----dev="png", fig.show="hold"-----------------------------------------------
landscape <- adaptive_landscape(
  data = prepared,
  fitness_model = surface$model,
  trait_cols = c("size", "colour"),
  grid_n = 50,
  simulation_n = 1000
)

landscape$optimum

## ----message=FALSE, warning=FALSE, fig.width=8, fig.height=6, fig.cap="Adaptive landscape showing mean fitness as function of population mean size and colour"----
p <- plot_adaptive_landscape(
  landscape = landscape,
  trait_cols = c("size", "colour"),
  bins = 12
)

print(p)

## ----message=FALSE, warning=FALSE, fig.width=10, fig.height=7, fig.cap="3D adaptive landscape showing the shape of mean fitness"----
plot_adaptive_landscape_3d(
  landscape = landscape,
  trait_cols = c("size", "colour"),
  theta = -30,
  phi = 30
)

## ----message=FALSE, warning=FALSE---------------------------------------------
results_binary <- selection_coefficients(
  data = prepared,
  fitness_col = "survival", # Binary 0/1
  trait_cols = c("size", "colour"), # Only available traits
  fitness_type = "binary",
  standardize = TRUE
)

print(results_binary)

## ----message=FALSE, warning=FALSE---------------------------------------------
selection_report(prepared, "survival", c("size", "colour"), fitness_type = "binary")

## ----message=FALSE, warning=FALSE---------------------------------------------
comparison <- compare_fitness_surfaces_data(
  correlated_surface = surface,
  adaptive_landscape = landscape,
  trait_cols = c("size", "colour")
)

comparison$summary_stats

## ----message=FALSE, warning=FALSE, fig.width=12, fig.height=5, fig.cap="Overlay comparison"----
plots <- plot_fitness_surfaces_comparison(
  comparison_data = comparison,
  bins = 15
)

print(plots$overlay)

## ----message=FALSE, warning=FALSE, fig.width=12, fig.height=5, fig.cap="Side-by-side comparison"----
plots <- plot_fitness_surfaces_comparison(
  comparison_data = comparison,
  bins = 15
)

print(plots$side_by_side)

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
finch_prep <- prepare_selection_data(finch_yearly, "survived", "beak_pc1")
years <- temporal_landscape(finch_prep, "survived", "beak_pc1", "year", smoothing = "REML", simulation_n = 200)

## -----------------------------------------------------------------------------
years$summary[, c("time", "n", "mean_fitness", "edf", "optimum_beak_pc1", "peaks")]

## ----message=FALSE, warning=FALSE, fig.width=10, fig.height=5.5, fig.cap="Survival against beak size, one panel per year. Solid: fitness function. Dashed: adaptive landscape. Diamond: highest survival."----
plot_temporal_landscape(years, ncol = 4)

## ----message=FALSE, warning=FALSE, fig.width=7, fig.height=4.5, fig.cap="The same fitness functions as a heat map, blank where a year has no birds."----
plot_temporal_landscape(years, type = "heatmap")

## ----message=FALSE, warning=FALSE---------------------------------------------
set.seed(1)
ca <- canonical_analysis(bumpus, "survival", c("total_length", "weight", "humerus"))
ca

## ----message=FALSE, warning=FALSE, fig.width=7, fig.height=5, fig.cap="Fitness surface along the first and last canonical axes"----
plot_canonical_axes(ca, which = c(1, 3), grid_n = 30)

## ----message=FALSE, warning=FALSE---------------------------------------------
check_selection_assumptions(bumpus, "survival", c("total_length", "weight", "humerus"))

