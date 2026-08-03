###############################################################
# Dynamic platelet digital phenotypes in severe acute pancreatitis
# Multi-cohort analysis: Local (Nanchang), TRACE, MIMIC-IV
# Follows the manuscript methods. Run with R >= 4.1.
###############################################################

# ---------------- 0. setup ----------------
pkgs <- c("lcmm", "survival", "rms", "glmnet", "randomForest", "tidymodels",
          "mediation", "lightgbm", "pROC", "data.table", "dplyr", "tidyr",
          "ggplot2", "readxl", "shapviz", "mice", "splines", "boot")
need <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, function(p) suppressMessages(require(p, character.only = TRUE))))

set.seed(2026)

# Paths - adjust to your machine
data_dir <- "D:/R/PLT"
rf_data  <- "D:/R/PLT/training_data.xlsx"
out_dir  <- "D:/R/PLT/PLT_repro/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# cohort files are located by pattern (file names contain CJK characters)
local_file <- list.files(data_dir, pattern = "^nanchang.*\\.csv$", recursive = TRUE, full.names = TRUE)[1]
trace_file <- list.files(data_dir, pattern = "^trace.*\\.csv$",  recursive = TRUE, full.names = TRUE)[1]
mimic_file <- list.files(data_dir, pattern = "^mimic.*\\.csv$", recursive = TRUE, full.names = TRUE)[1]
if (anyNA(c(local_file, trace_file, mimic_file))) stop("cohort files not found; check data_dir")

# ---------- read local (nanchang) cohort ----------
# The file header spans three lines with an unclosed quote; parse it directly.
parse_csv <- function(lines) {
  txt <- paste(lines, collapse = "\n")
  ch  <- strsplit(txt, "", fixed = TRUE)[[1]]
  n   <- length(ch); fields <- character(0); cur <- ""; inq <- FALSE; i <- 1
  while (i <= n) {
    if (inq) {
      if (ch[i] == "\"") {
        if (i < n && ch[i+1] == "\"") { cur <- paste0(cur, "\""); i <- i + 2; next }
        inq <- FALSE; i <- i + 1
      } else { cur <- paste0(cur, ch[i]); i <- i + 1 }
    } else {
      if (ch[i] == ",") { fields <- c(fields, cur); cur <- ""; i <- i + 1 }
      else if (ch[i] == "\"" && cur == "") { inq <- TRUE; i <- i + 1 }
      else { cur <- paste0(cur, ch[i]); i <- i + 1 }
    }
  }
  c(fields, cur)
}

read_gbk <- function(f) {
  con <- file(f, open = "r", encoding = "GBK")
  lines <- readLines(con, warn = FALSE); close(con)
  lines
}

local_lines <- read_gbk(local_file)
hdr  <- parse_csv(local_lines[1:3]); hdr[23] <- "Low.molecular.weight.heparin"
dat  <- lapply(local_lines[-c(1:3)], parse_csv)
loc  <- as.data.frame(t(sapply(dat, function(x) x[seq_along(hdr)])), stringsAsFactors = FALSE)
names(loc) <- hdr
loc <- as.data.frame(lapply(loc, function(x) suppressWarnings(as.numeric(x))))

# reference GBTM labels are stored in the analysis file (frozen-model output)
loc$trajectory <- as.integer(loc$Class)

# ---------------- read trace / mimic ----------------
read_csv_gbk <- function(f) {
  read.csv(f, fileEncoding = "GBK", check.names = FALSE, na.strings = c("", "NA"))
}
trace <- read_csv_gbk(trace_file)
mimic <- read_csv_gbk(mimic_file)
# trajectory-class column position: 37 (mimic), 40 (trace)
mimic$trajectory <- as.integer(mimic[[37]])
trace$trajectory <- as.integer(trace[[40]])

# ---------------- common variable names ----------------
# each cohort is reduced to the analysis-ready data.frame below
# Local file columns are accessed by position (the header mixes CJK and
# full-width characters)
prep_local <- function(d) {
  data.frame(
    cohort = "Local", death_ih = d[[2]],
    fut_ih  = d[[3]],
    death_28 = d[[6]], fut_28 = d[[7]],
    death_90 = d[[8]], fut_90 = d[[9]],
    CumPLT = d[[10]], trajectory = as.integer(d[[31]]),
    PLT01 = d[[11]], PLT02 = d[[12]], PLT03 = d[[13]], PLT05 = d[[14]], PLT07 = d[[15]],
    age = d[[33]], male = d[[32]], DM = d[[51]], HTG = d[[50]],
    temp = d[[52]], pulse = d[[59]], resp = d[[66]], SBP = d[[73]], DBP = d[[80]],
    WBC = d[[95]], HCT = d[[100]], ALB = d[[109]], Tbil = d[[107]], BUN = d[[119]], Cr = d[[120]],
    bc = d[[27]],
    ABX = d[[18]], ABX7 = d[[20]], LMWH = d[[23]],
    IPN = d[[26]], PCD = d[[195]], CRRT = d[[196]], MV = d[[197]]
  )
}

prep_other <- function(d, cohort, mv_col) {
  data.frame(
    cohort = cohort, death_ih = d$In.hospital.Death,
    fut_ih  = d$followtime.In.hospital.Death.,
    death_28 = d$X28.Day.Death, fut_28 = d$followtime.28.Day.Death..,
    death_90 = d$X90.Day.Death, fut_90 = d$followtime.90.Day.Death..,
    CumPLT = d$CumPLT, trajectory = d$trajectory,
    PLT01 = d$PLT01, PLT02 = d$PLT02, PLT03 = d$PLT03, PLT05 = d$PLT05, PLT07 = d$PLT07,
    age = d$age, male = d$male, DM = d$DM, HTG = d$HTG,
    temp = d$T.FR.modify, pulse = d$P.FR.modify, resp = d$R.FR.modify,
    SBP = d$SBP, DBP = d$DBP, WBC = d$WBC, HCT = d$HCT, ALB = d$ALB, Tbil = d$Tbil,
    BUN = d$BUN, Cr = d$Cr, bc = d$Bloodstream.infection,
    CRRT = d$CRRT, MV = d[[mv_col]]
  )
}

loc   <- prep_local(loc)
trace <- prep_other(trace, "TRACE", mv_col = 38)
mimic <- prep_other(mimic, "MIMIC", mv_col = 343)
all_cohorts <- list(Local = loc, TRACE = trace, MIMIC = mimic)

# ---------------- 1. cohort description ----------------
cohort_summary <- lapply(all_cohorts, function(d) {
  c(n = nrow(d), male = mean(d$male), age = mean(d$age),
    death_ih = sum(d$death_ih), death_28 = sum(d$death_28), death_90 = sum(d$death_90),
    bc = sum(d$bc))
})
print(do.call(rbind, cohort_summary))

# ---------------- 2. cumulative platelet exposure ----------------
# trapezoid rule over t = {1,2,3,5,7}; observed values used, missing points
# predicted by a linear mixed model (done upstream); recompute to verify.
calc_cumplt <- function(d) {
  tv <- c(1, 2, 3, 5, 7)
  plts <- as.matrix(d[, c("PLT01", "PLT02", "PLT03", "PLT05", "PLT07")])
  apply(plts, 1, function(r) sum(0.5 * (r[1:4] + r[2:5]) * diff(tv)))
}
for (nm in names(all_cohorts)) {
  d <- all_cohorts[[nm]]
  cat(sprintf("%s: CumPLT(trapezoid) vs file correlation = %.4f\n",
              nm, cor(calc_cumplt(d), d$CumPLT, use = "complete.obs")))
}

# standardize within cohort
for (nm in names(all_cohorts)) {
  all_cohorts[[nm]]$CumPLTz <- as.numeric(scale(all_cohorts[[nm]]$CumPLT))
  all_cohorts[[nm]]$PLT01z <- as.numeric(scale(all_cohorts[[nm]]$PLT01))
  all_cohorts[[nm]]$PLT02z <- as.numeric(scale(all_cohorts[[nm]]$PLT02))
  all_cohorts[[nm]]$PLT03z <- as.numeric(scale(all_cohorts[[nm]]$PLT03))
  all_cohorts[[nm]]$PLT05z <- as.numeric(scale(all_cohorts[[nm]]$PLT05))
  all_cohorts[[nm]]$PLT07z <- as.numeric(scale(all_cohorts[[nm]]$PLT07))
}

# ---------------- 3. GBTM development (Local) ----------------
# group-based trajectory modeling: Gaussian outcome, cubic polynomial in time.
# The 5-class solution and per-patient class labels used throughout this
# script are the frozen-model output stored in the cohort files. The base
# (1-class) model is refit below to reproduce the fitting procedure.
long_local <- data.frame(
  id   = rep(seq_len(nrow(loc)), 5),
  time = rep(c(1, 2, 3, 5, 7), each = nrow(loc)),
  PLT  = as.vector(t(as.matrix(loc[, c("PLT01","PLT02","PLT03","PLT05","PLT07")])))
)
m1 <- hlme(PLT ~ poly(time, degree = 3), random = ~-1, ng = 1,
           subject = "id", data = long_local, maxiter = 500)
cat("GBTM 1-class model: loglik", round(m1$loglik, 2),
    " AIC", round(m1$AIC, 2), " BIC", round(m1$BIC, 2), "\n")

cat("Reference 5-class proportions:",
    paste(names(table(loc$trajectory)), table(loc$trajectory),
          sep = "=", collapse = ", "), "\n")

# class-specific mean trajectories (Fig 2A)
traj_means <- t(sapply(1:5, function(k) {
  colMeans(loc[loc$trajectory == k, c("PLT01","PLT02","PLT03","PLT05","PLT07")])
}))
matplot(c(1, 2, 3, 5, 7), t(traj_means), type = "b", lwd = 2, pch = 19,
        xlab = "Day", ylab = "Platelet count",
        main = "Mean platelet trajectories by GBTM class (Local)")
legend("topleft", legend = paste("Class", 1:5), col = 1:5, lty = 1, lwd = 2)

# ---------------- 4. associations: Cox models ----------------
# exposures: CumPLT (per SD), daily PLT (per SD), trajectory class (class 5 ref)
# partial adjustment: age, sex, DM, HTG
# full adjustment: + temp, pulse, resp, SBP, DBP, WBC, HCT, ALB, TBIL
cov_part <- c("age", "male", "DM", "HTG")
cov_full <- c(cov_part, "temp", "pulse", "resp", "SBP", "DBP", "WBC", "HCT", "ALB", "Tbil")

cox_assoc <- function(d, exposure, outcome, full = TRUE) {
  fname <- c(ih = "fut_ih", "28" = "fut_28", "90" = "fut_90")
  dname <- c(ih = "death_ih", "28" = "death_28", "90" = "death_90")
  covs <- if (full) cov_full else cov_part
  f <- as.formula(paste("Surv(", fname[outcome], ",", dname[outcome],
                        ") ~ ", exposure, " + ", paste(covs, collapse = " + ")))
  fit <- coxph(f, data = d)
  s <- summary(fit)$coefficients[exposure, ]
  b  <- as.numeric(s["coef"])
  se <- as.numeric(s["se(coef)"])
  c(HR = exp(b), lo = exp(b - 1.96*se), hi = exp(b + 1.96*se),
    p = as.numeric(s["Pr(>|z|)"]))
}

cat("\n--- CumPLT per SD, fully adjusted, in-hospital ---\n")
for (nm in names(all_cohorts)) {
  r <- cox_assoc(all_cohorts[[nm]], "CumPLTz", "ih", TRUE)
  cat(sprintf("%-6s HR=%.2f (%.2f-%.2f)\n", nm, r["HR"], r["lo"], r["hi"]))
}

# trajectory classes (class 5 reference), full adjustment
cat("\n--- trajectory classes vs class 5, fully adjusted, in-hospital ---\n")
for (nm in names(all_cohorts)) {
  d <- all_cohorts[[nm]]
  d$cls <- factor(d$trajectory)
  d$cls <- relevel(d$cls, ref = "5")
  f <- as.formula(paste("Surv(fut_ih, death_ih) ~ cls + ", paste(cov_full, collapse=" + ")))
  fit <- coxph(f, data = d)
  hr <- exp(coef(fit)[paste0("cls", 1:4)])
  cat(sprintf("%-6s %s\n", nm, paste(sprintf("C%d=%.2f", 1:4, hr), collapse=" ")))
}

# ---------------- 5. dose-response: restricted cubic splines ----------------
cat("\n--- RCS: CumPLT and in-hospital mortality ---\n")
for (nm in names(all_cohorts)) {
  d <- all_cohorts[[nm]]
  d <- d[complete.cases(d[, c("fut_ih","death_ih","CumPLTz",cov_full)]), ]
  dd <- rms::datadist(d); options(datadist = "dd")
  f <- as.formula(paste("Surv(fut_ih, death_ih) ~ rcs(CumPLTz, 4) + ", paste(cov_full, collapse=" + ")))
  fit <- rms::cph(f, data = d, x = TRUE, y = TRUE)
  av <- anova(fit)
  cat(sprintf("%-6s overall P=%.3f  nonlinear P=%.3f\n", nm,
              av["TOTAL", "P"], av[" Nonlinear", "P"]))
}

# ---------------- 6. blood culture associations ----------------
cat("\n--- CumPLT per SD -> positive blood culture (fully adjusted) ---\n")
for (nm in names(all_cohorts)) {
  d <- all_cohorts[[nm]]
  f <- as.formula(paste("bc ~ CumPLTz + ", paste(cov_full, collapse = " + ")))
  fit <- glm(f, data = d, family = binomial)
  est <- coef(fit)["CumPLTz"]; se <- sqrt(diag(vcov(fit)))["CumPLTz"]
  cat(sprintf("%-6s OR=%.2f (%.2f-%.2f)\n", nm, exp(est), exp(est-1.96*se), exp(est+1.96*se)))
}

# ---------------- 7. mediation (1000 bootstrap) ----------------
med_dat <- all_cohorts$Local
med_dat <- med_dat[complete.cases(med_dat[, c("death_ih","bc","CumPLTz",cov_full)]), ]
med.fit <- glm(bc ~ CumPLTz + age + male + DM + HTG + temp + pulse + resp + SBP + DBP +
                 WBC + HCT + ALB + Tbil, data = med_dat, family = binomial)
out.fit <- glm(death_ih ~ CumPLTz + bc + age + male + DM + HTG + temp + pulse + resp +
                 SBP + DBP + WBC + HCT + ALB + Tbil, data = med_dat, family = binomial)
set.seed(1)
med_res <- mediation::mediate(med.fit, out.fit, treat = "CumPLTz", mediator = "bc",
                              boot = TRUE, sims = 1000, boot.ci.type = "perc")
cat("\n--- mediation, Local, in-hospital (CumPLT -> blood culture -> death) ---\n")
cat("ACME:", round(med_res$d0, 4), " ADE:", round(med_res$z0, 4),
    " proportion mediated:", round(med_res$n0, 4), "\n")

# ---------------- 8. DML (exploratory) ----------------
# double machine learning, orthogonal AIPW, LightGBM nuisance functions,
# 5-fold cross-fitting repeated 5 times, propensity truncated at 0.02-0.98
dml_ate <- function(d, treatment, n_rep = 5, n_folds = 5, nrounds = 50) {
  d <- d[complete.cases(d[, c("death_ih", treatment, cov_full)]), ]
  X <- as.matrix(d[, cov_full]); n <- nrow(d)
  y <- d$death_ih; a <- d[[treatment]]
  set.seed(1); theta <- numeric(n_rep); psi_all <- c()
  for (r in seq_len(n_rep)) {
    folds <- sample(rep(seq_len(n_folds), length.out = n))
    psi <- numeric(n)
    for (k in seq_len(n_folds)) {
      tr <- folds != k; te <- folds == k
      mfit <- lgb.train(params = list(objective = "binary", learning_rate = 0.1,
                                      num_leaves = 31, verbose = -1),
                        data = lgb.Dataset(X[tr, , drop = FALSE], label = a[tr]),
                        nrounds = nrounds)
      mhat <- pmin(pmax(predict(mfit, X[te, , drop = FALSE]), 0.02), 0.98)
      ofit <- lgb.train(params = list(objective = "binary", learning_rate = 0.1,
                                      num_leaves = 31, verbose = -1),
                        data = lgb.Dataset(cbind(X[tr, , drop = FALSE], trt = a[tr]),
                                           label = y[tr]),
                        nrounds = nrounds)
      g1 <- predict(ofit, cbind(X[te, , drop = FALSE], trt = 1))
      g0 <- predict(ofit, cbind(X[te, , drop = FALSE], trt = 0))
      A <- a[te]; Y <- y[te]
      psi[te] <- A*(Y - g1)/mhat - (1 - A)*(Y - g0)/(1 - mhat) + (g1 - g0)
    }
    theta[r] <- mean(psi); psi_all <- c(psi_all, psi)
  }
  est <- mean(theta); se <- sd(psi_all)/sqrt(length(psi_all))
  c(est = est, se = se, lo = est - 1.96*se, hi = est + 1.96*se,
    p = 2*pnorm(-abs(est/se)))
}

cat("\n--- DML: week-1 IV antibiotics, Local ---\n")
print(dml_ate(all_cohorts$Local, "ABX"))
cat("--- DML: LMWH, Local ---\n")
print(dml_ate(all_cohorts$Local, "LMWH"))

# ---------------- 9. random forest classifier ----------------
# training data (13 LASSO-selected baseline variables + GBTM class)
rf_data_all <- as.data.frame(readxl::read_excel(rf_data))
rf_data_all$group <- factor(rf_data_all$group)

# 3:1 stratified split (fixed seed, matches the manuscript analysis)
strat_split <- function(d, seed) {
  set.seed(seed)
  res <- c()
  for (k in sort(unique(as.character(d$group)))) {
    rows <- which(as.character(d$group) == k); nk <- length(rows)
    s <- sample(nk, floor(nk * 0.75)); res <- c(res, rows[sort(s)])
  }
  list(train = res, test = setdiff(seq_len(nrow(d)), res))
}
sp <- strat_split(rf_data_all, 4321)
train_rf <- rf_data_all[sp$train, ]
test_rf  <- rf_data_all[sp$test, ]

# multinomial LASSO variable selection on the training set
# the 13 baseline variables retained under the one-standard-error rule
cand <- c("age","HTG","DM","Temperature","pulse","Respiration","SBP","WBC",
          "PLT01","ALB","HCT","BUN","Cr")
set.seed(42)
cvfit <- cv.glmnet(as.matrix(train_rf[, cand]), train_rf$group,
                   family = "multinomial", type.measure = "deviance")
cat("LASSO lambda.1se:", cvfit$lambda.1se, " retained variables:", length(cand), "\n")

# random forest tuned with 5-fold CV, then refit
set.seed(42)
folds <- vfold_cv(train_rf, v = 5, strata = group)
grid_rf <- expand.grid(mtry = c(2, 10), trees = c(200, 350, 500), min_n = c(50, 100))
spec_rf <- rand_forest(mtry = tune(), trees = tune(), min_n = tune()) %>%
  set_engine("randomForest", importance = TRUE) %>% set_mode("classification")
wk <- workflow() %>% add_formula(group ~ .) %>% add_model(spec_rf)
tune_rf <- tune_grid(wk, resamples = folds, grid = grid_rf,
                     metrics = metric_set(accuracy))
best_hp <- select_best(tune_rf, metric = "accuracy")
print(best_hp[, c("mtry","trees","min_n")])

set.seed(42)
final_rf <- finalize_workflow(wk, best_hp) %>% fit(train_rf)
rf_model <- extract_fit_engine(final_rf)

# apply to external cohorts (reference labels from the frozen GBTM)
ext_external <- function(d, nm) {
  newd <- data.frame(
    age = d$age, HTG = d$HTG, DM = d$DM, Temperature = d$temp,
    pulse = d$pulse, Respiration = d$resp, SBP = d$SBP, WBC = d$WBC,
    PLT01 = d$PLT01, ALB = d$ALB, HCT = d$HCT, BUN = d$BUN, Cr = d$Cr,
    group = factor(d$trajectory, levels = levels(rf_data_all$group)))
  pred <- predict(final_rf, new_data = newd, type = "prob") %>%
    bind_cols(predict(final_rf, new_data = newd, type = "class")) %>%
    bind_cols(dplyr::select(newd, group))
  pred$dataset <- nm
  pred
}
pred_train <- predict(final_rf, new_data = train_rf, type = "prob") %>%
  bind_cols(predict(final_rf, new_data = train_rf, type = "class")) %>%
  bind_cols(dplyr::select(train_rf, group)); pred_train$dataset <- "train"
pred_test  <- predict(final_rf, new_data = test_rf, type = "prob") %>%
  bind_cols(predict(final_rf, new_data = test_rf, type = "class")) %>%
  bind_cols(dplyr::select(test_rf, group)); pred_test$dataset <- "test"
pred_trace <- ext_external(all_cohorts$TRACE, "TRACE")
pred_mimic <- ext_external(all_cohorts$MIMIC, "MIMIC")

macro_auc <- function(pred) {
  aucs <- sapply(levels(pred$group), function(cl) {
    suppressMessages(pROC::auc(pROC::roc(ifelse(pred$group == cl, 1, 0),
                                         pred[[paste0(".pred_", cl)]], quiet = TRUE)))
  })
  mean(aucs)
}
cat("\n--- classifier macro AUC ---\n")
cat("train:", round(macro_auc(pred_train), 3), "\n")
cat("test :", round(macro_auc(pred_test), 3), "\n")
cat("TRACE:", round(macro_auc(pred_trace), 3), "\n")
cat("MIMIC:", round(macro_auc(pred_mimic), 3), "\n")

# confusion matrices and metrics
for (nm in c("train", "test", "TRACE", "MIMIC")) {
  p <- get(paste0("pred_", tolower(nm)))
  cm <- p %>% conf_mat(group, .pred_class)
  ms <- summary(cm)
  acc <- ms[ms$.metric == "accuracy", ".estimate"]
  kap <- ms[ms$.metric == "kap", ".estimate"]
  cat(sprintf("%-6s accuracy=%.3f kappa=%.3f\n", nm, acc, kap))
}

# SHAP explanations (training set)
set.seed(42)
shap_scores <- tryCatch(
  fastshap::explain(rf_model, X = as.data.frame(train_rf[, cand]), nsim = 20,
                    pred_wrapper = function(m, x) predict(m, x, type = "prob")),
  error = function(e) NULL)
if (!is.null(shap_scores)) print(shap_scores)

# ---------------- 10. sensitivity analyses ----------------
# K-means on the 5 time points (trajectory structure sensitivity)
for (nm in names(all_cohorts)) {
  d <- all_cohorts[[nm]]
  km <- kmeans(scale(d[, c("PLT01","PLT02","PLT03","PLT05","PLT07")]), centers = 5, nstart = 25)
  cat(sprintf("%s K-means cluster sizes: %s\n", nm, paste(table(km$cluster), collapse = ", ")))
}

# multiple imputation (sensitivity 2): Local fully adjusted model
cat("\n--- multiple imputation sensitivity (Local) ---\n")
imp_dat <- all_cohorts$Local[, c("fut_ih","death_ih","CumPLTz",cov_full)]
set.seed(1)
mi <- mice(imp_dat, m = 5, method = "pmm", printFlag = FALSE)
fit_mi <- with(mi, coxph(Surv(fut_ih, death_ih) ~ CumPLTz + age + male + DM + HTG +
                           temp + pulse + resp + SBP + DBP + WBC + HCT + ALB + Tbil))
pooled <- summary(pool(fit_mi))
hr_mi <- pooled[pooled$term == "CumPLTz", "estimate"]
cat("pooled HR for CumPLTz:", round(exp(hr_mi), 2), "\n")

# E-values for the fully-adjusted estimates
# (VanderWeele & Ding 2017); HR treated as RR under a rare outcome
evalue_from_rr <- function(rr) {
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}
cat("\n--- E-values (in-hospital, CumPLT per SD) ---\n")
for (nm in names(all_cohorts)) {
  r <- cox_assoc(all_cohorts[[nm]], "CumPLTz", "ih", TRUE)
  e_point <- evalue_from_rr(as.numeric(r["HR"]))
  e_bound <- evalue_from_rr(ifelse(as.numeric(r["HR"]) >= 1,
                                   as.numeric(r["lo"]), as.numeric(r["hi"])))
  cat(nm, ": E-value", round(e_point, 2), " (CI bound", round(e_bound, 2), ")\n")
}

cat("\nAll analyses complete. Outputs written to", out_dir, "\n")
