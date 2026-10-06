#============================================================
# Run TMB REM Poisson model
#
# Model:
#   y_m ~ Poisson(mu_m)
#
#   log(mu_m) =
#       X_m beta
#     + g[giver_m]
#     + r[receiver_m]
#     + gpe[giver_m]
#     + rpe[receiver_m]
#     + sg[group_m]
#     + offset_m
#
# Random effects:
#
#   [g, r] ~ N(0, G ⊗ Sigma_g)
#
#   Sigma_g =
#     [ var_g    cov_gr ]
#     [ cov_gr   var_r  ]
#
#   [gpe_i, rpe_i] ~ N(0, Sigma_pe)
#
#   Sigma_pe =
#     [ var_gpe       cov_gpe_rpe ]
#     [ cov_gpe_rpe   var_rpe     ]
#
#   sg_k ~ N(0, var_sg)
#============================================================

library(TMB)
library(Matrix)

#============================================================
# 1. Read data
#============================================================

rem_table_all <- readRDS("rem_table_all.rds")
G <- readRDS("G.rds")

# make sure G is matrix-like
G <- as.matrix(G)

#============================================================
# 2. Response
#============================================================

y <- rem_table_all$y

#============================================================
# 3. Fixed-effect matrix X
#============================================================

# Your current fixed effects:
# columns 7:20 and 24:29
fixed_df <- rem_table_all[, c(7:23)]

X <- model.matrix(~ ., data = fixed_df)

#============================================================
# 4. Offset
#============================================================

offset <- rem_table_all$offset

#============================================================
# 5. Animal levels
#============================================================

animal_levels <- rownames(G)

#============================================================
# 6. giver / receiver index
#============================================================

giver_id <- match(as.character(rem_table_all$giver), as.character(animal_levels)) - 1
receiver_id <- match(as.character(rem_table_all$receiver), as.character(animal_levels)) - 1

n_animal <- length(animal_levels)

#============================================================
# 7. group index
#============================================================

group_factor <- factor(rem_table_all$group)
group_id <- as.integer(group_factor) - 1
n_group <- nlevels(group_factor)

#============================================================
# 8. G inverse and log determinant of G
#============================================================

Ginv <- solve(G)
Ginv <- as(Ginv, "dgCMatrix")

logdet_G <- as.numeric(determinant(G, logarithm = TRUE)$modulus)

#============================================================
# 9. TMB data list
#============================================================

tmb_data <- list(
  y           = as.numeric(y),
  X           = X,
  offset      = as.numeric(offset),
  
  giver_id    = as.integer(giver_id),
  receiver_id = as.integer(receiver_id),
  
  group_id    = as.integer(group_id),
  
  Ginv        = Ginv,
  logdet_G    = as.numeric(logdet_G),
  
  n_animal    = as.integer(n_animal),
  n_group     = as.integer(n_group)
)

#============================================================
# 10. Initial parameters
#============================================================

params <- list(
  beta = rep(0, ncol(X)),
  
  # genetic random effects
  g    = rep(0, n_animal),  # giver genetic effect
  r    = rep(0, n_animal),  # receiver genetic effect
  
  # permanent environmental random effects
  gpe  = rep(0, n_animal),  # giver PE
  rpe  = rep(0, n_animal),  # receiver PE
  
  # social group effect
  sg   = rep(0, n_group),
  
  # genetic variance components
  log_var_g = log(0.1),
  log_var_r = log(0.1),
  
  # PE variance components
  log_var_gpe = log(0.1),
  log_var_rpe = log(0.1),
  
  # social group variance
  log_var_sg = log(0.1),
  
  # correlations on unconstrained scale
  theta_gr = 0,
  theta_gpe_rpe = 0
)

#============================================================
# 11. Compile and load TMB model
#============================================================

compile("rem_poisson_gr_cor.cpp")
dyn.load(dynlib("rem_poisson_gr_cor"))

#============================================================
# 12. Build TMB objective function
#============================================================

obj <- MakeADFun(
  data       = tmb_data,
  parameters = params,
  random     = c("g", "r", "gpe", "rpe", "sg"),
  DLL        = "rem_poisson_gr_cor",
  silent     = TRUE
)

#============================================================
# 13. Optimize fixed effects + variance components
#============================================================
cat("Starting optimization:", Sys.time(), "\n")
opt <- nlminb(
  start     = obj$par,
  objective = obj$fn,
  gradient  = obj$gr,
  control   = list(
    eval.max = 1000,
    iter.max = 1000,
    trace    = 1
  )
)

cat("Convergence code:", opt$convergence, "\n")
cat("Message:", opt$message, "\n")
cat("Objective:", opt$objective, "\n")

#============================================================
# 14. Check gradient
#============================================================

gr <- obj$gr(opt$par)
cat("Max absolute gradient:", max(abs(gr)), "\n")

#============================================================
# 15. Report estimates and standard errors
#============================================================

rep <- sdreport(obj, par.fixed = opt$par)

cat("\nFixed parameters:\n")
print(summary(rep, "fixed"))

cat("\nReported parameters:\n")
print(summary(rep, "report"))

#============================================================
# 16. Extract variance components
#============================================================

report_summary <- summary(rep, "report")
print(report_summary)

#============================================================
# 17. AIC-like value
#
# Note:
# This uses the marginal negative log likelihood from TMB.
# Random effects are integrated out by Laplace approximation.
# k is the number of non-random parameters optimized by nlminb.
#============================================================

k <- length(opt$par)
AIC_tmb <- 2 * opt$objective + 2 * k

cat("\nNumber of non-random parameters:", k, "\n")
cat("AIC:", AIC_tmb, "\n")

#============================================================
# 18. Save output
#============================================================

saveRDS(
  list(
    opt = opt,
    rep = rep,
    report_summary = report_summary,
    AIC = AIC_tmb,
    parameter = opt$par,
    gradient = gr,
    max_abs_gradient = max(abs(gr))
  ),
  file = "rem_poisson_gr_cor_fit.rds"
)
