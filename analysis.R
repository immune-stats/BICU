# FeLV and FIV: reproducible analyses for the accompanying Data article.
# Run from this folder: Rscript analysis.R
# Required packages: pROC and ResourceSelection. See README.md for cohort rules.
# Output CSVs are numerical source tables; the manuscript uses formatted tables.

# 1. Load data and record the R environment

packages <- c("pROC", "ResourceSelection")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages: ", paste(missing_packages, collapse = ", "))
options(contrasts = c("contr.treatment", "contr.poly"))
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args)) args[1] else "dataset.csv"
output_dir <- if (length(args) > 1) args[2] else "results"
if (!file.exists(input_file)) stop("Run from the repository root, or supply the CSV path as the first argument.")
input <- read.csv(input_file, na.strings = "", stringsAsFactors = FALSE, check.names = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
write.csv(data.frame(package = packages, version = vapply(packages, function(p) as.character(packageVersion(p)), "")),
          file.path(output_dir, "package_versions.csv"), row.names = FALSE)

# 2. Define factor reference categories and model specifications

levels_map <- list(breed = c("y", "n"), sex = c("f", "m"),
  sex_neuter = c("if", "im", "nf", "nm"), neuter_status = c("i", "n"),
  cohabitantcats = c("single cat", "multicat"), lifestyle = c("indoor", "outdoor"),
  concomitant_disorder = c("n", "y"), c_hct = c("normal", "low", "high"),
  c_leuco = c("normal", "low", "high"))
required <- c("period", "felv", "fiv", "age_group", "origin", names(levels_map))
stopifnot(setequal(required, names(input)))
stopifnot(all(is.na(input$age_group) | input$age_group %in% c("kitten", "adult", "senior")))
input$age_group <- factor(input$age_group, levels = c("kitten", "adult", "senior"))
stopifnot(all(is.na(input$origin) | input$origin %in% c("stray", "breeder", "shelter", "private_individual")))
for (v in names(levels_map)) {
  if (any(!is.na(input[[v]]) & !input[[v]] %in% levels_map[[v]])) stop("Unknown category in ", v)
  input[[v]] <- factor(input[[v]], levels = levels_map[[v]])
}
stopifnot(all(input$period %in% c("2013-2022", "2023-2025")))
for (v in c("felv", "fiv")) stopifnot(all(is.na(input[[v]]) | input[[v]] %in% 0:1))
input$period <- factor(input$period, levels = c("2013-2022", "2023-2025"))

original_vars <- list(FeLV = c("breed", "c_hct", "concomitant_disorder", "age_group"),
                     FIV = c("breed", "sex", "age_group", "lifestyle", "concomitant_disorder"))
new_vars <- list(FeLV = c(original_vars$FeLV, "c_leuco", "neuter_status"), FIV = original_vars$FIV)
candidate_vars <- c("breed", "sex", "neuter_status", "sex_neuter", "age_group",
                    "cohabitantcats", "lifestyle", "concomitant_disorder", "c_leuco", "c_hct")
# Full-precision coefficients and thresholds from the saved development exports.
# These constants are never fitted to the later admissions.
frozen_beta <- list(
  FeLV = c("(Intercept)" = -3.42952720397276, breedn = 1.22148308991807,
           c_hctlow = 0.947603557505741, c_hcthigh = 0.0503774598591411,
           concomitant_disordery = 0.521518249911755,
           age_groupadult = 0.408478355189583, age_groupsenior = -0.617429075880352),
  FIV = c("(Intercept)" = -3.96664949226748, breedn = 0.620349311566591,
          sexm = 0.552171497080684, age_groupadult = 1.04089728826907,
          age_groupsenior = 0.812493536296807, lifestyleoutdoor = 0.872084056537284,
          concomitant_disordery = 0.637734428949291))
frozen_threshold <- c(FeLV = 0.211358351099469, FIV = 0.211843372222884)

# 3. Statistical calculations shared by the analyses

model_formula <- function(vars) reformulate(vars, response = "outcome")
fit_model <- function(data, vars, link = "logit") {
  glm(model_formula(vars), data = data, family = binomial(link), na.action = na.fail)
}
roc_object <- function(y, p) pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
roc01 <- function(r) {
  t <- as.numeric(unlist(pROC::coords(r, "best", best.method = "closest.topleft", ret = "threshold")))
  t <- t[is.finite(t)]
  if (!length(t)) stop("No finite ROC01 threshold.")
  t[1]
}
hosmer_lemeshow <- function(y, p, g = 10L) {
  h <- ResourceSelection::hoslem.test(y, p, g = g)
  occupied <- rowSums(h$observed) > 0
  observed <- h$observed[occupied, , drop = FALSE]
  expected <- h$expected[occupied, , drop = FALSE]
  df <- sum(occupied) - 2L
  if (df <= 0 || any(!is.finite(expected) | expected <= 0)) {
    stop("Hosmer-Lemeshow is undefined for the occupied groups.")
  }
  # Remove empty factor levels only; all observations remain in the test.
  h$statistic <- c("X-squared" = sum((observed - expected)^2 / expected))
  h$parameter <- c(df = df)
  h$p.value <- pchisq(h$statistic, df, lower.tail = FALSE)
  h$observed <- observed; h$expected <- expected
  h
}
performance <- function(y, p, threshold = NULL) {
  r <- roc_object(y, p)
  if (is.null(threshold)) threshold <- roc01(r)
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  positive <- p >= threshold
  tp <- sum(positive & y == 1); tn <- sum(!positive & y == 0)
  fp <- sum(positive & y == 0); fn <- sum(!positive & y == 1)
  hl_error <- ""
  hl <- tryCatch(hosmer_lemeshow(y, p, g = 10), error = function(e) {
    hl_error <<- conditionMessage(e); NULL
  })
  se <- tp / (tp + fn); sp <- tn / (tn + fp)
  var_se <- se * (1 - se) / (tp + fn)
  var_sp <- sp * (1 - sp) / (tn + fp)
  z <- qnorm(0.975)
  ba <- (se + sp) / 2
  data.frame(sensitivity_lower_95 = se - z * sqrt(var_se),
    sensitivity_upper_95 = se + z * sqrt(var_se),
    specificity_lower_95 = sp - z * sqrt(var_sp),
    specificity_upper_95 = sp + z * sqrt(var_sp),
    balanced_accuracy_lower_95 = ba - z * sqrt((var_se + var_sp) / 4),
    balanced_accuracy_upper_95 = ba + z * sqrt((var_se + var_sp) / 4),
    n = length(y), cases = sum(y == 1), non_cases = sum(y == 0),
    auc = as.numeric(pROC::auc(r)), lower_95 = ci[1], upper_95 = ci[3], threshold = threshold,
    tp = tp, fn = fn, tn = tn, fp = fp, sensitivity = tp / (tp + fn), specificity = tn / (tn + fp),
    balanced_accuracy = (tp / (tp + fn) + tn / (tn + fp)) / 2,
    ppv = tp / (tp + fp), npv = tn / (tn + fn), accuracy = (tp + tn) / length(y),
    brier = mean((y - p)^2), hl_statistic = if (is.null(hl)) NA_real_ else unname(hl$statistic),
    hl_df = if (is.null(hl)) NA_real_ else unname(hl$parameter),
    hl_p = if (is.null(hl)) NA_real_ else hl$p.value, hl_error = hl_error)
}
coefficients_table <- function(fit) {
  s <- coef(summary(fit)); z <- qnorm(0.975)
  data.frame(term = rownames(s), beta = s[, 1], se = s[, 2], p = s[, 4],
             lower_beta = s[, 1] - z * s[, 2], upper_beta = s[, 1] + z * s[, 2],
             OR = exp(s[, 1]), lower_OR = exp(s[, 1] - z * s[, 2]),
             upper_OR = exp(s[, 1] + z * s[, 2]), row.names = NULL)
}
lrt <- function(base, augmented) {
  stopifnot(nobs(base) == nobs(augmented), identical(rownames(model.frame(base)), rownames(model.frame(augmented))))
  df <- attr(logLik(augmented), "df") - attr(logLik(base), "df")
  statistic <- max(0, 2 * as.numeric(logLik(augmented) - logLik(base)))
  data.frame(df = df, statistic = statistic, p = if (df > 0) pchisq(statistic, df, lower.tail = FALSE) else NA_real_)
}
complete_sample <- function(data, vars) data[complete.cases(data[, c("outcome", vars), drop = FALSE]), , drop = FALSE]
frozen_prediction <- function(d, infection) {
  mm <- model.matrix(reformulate(original_vars[[infection]]), d)
  b <- frozen_beta[[infection]]
  stopifnot(setequal(colnames(mm), names(b)))
  as.numeric(plogis(mm[, names(b), drop = FALSE] %*% b))
}

# 4. Construct cohorts and reproduce the article tables
# FeLV and FIV share the same workflow, using the specifications above.

run_analysis <- function() {
  folder <- output_dir
  tables <- list(); fits <- list(); curves <- list(); cohort_list <- list()
  add <- function(name, value) {
    tables[[name]] <<- if (is.null(tables[[name]])) value else rbind(tables[[name]], value)
  }
  d <- input
  for (infection in c("FeLV", "FIV")) {
    x <- d; x$outcome <- x[[tolower(infection)]]
    x <- x[!is.na(x$outcome), , drop = FALSE]
    dev_required <- candidate_vars
    # Preserve the historical FIV test-cohort eligibility requirements.
    test_required <- if (infection == "FeLV") original_vars$FeLV else
      unique(c(original_vars$FIV, "neuter_status", "sex_neuter"))
    dev <- complete_sample(x[x$period == "2013-2022", ], dev_required)
    test <- complete_sample(x[x$period == "2023-2025", ], test_required)
    combined <- rbind(dev, test)
    cohort_list[[infection]] <- combined
    # Ensure the historical development cohort recreates the original fixed coefficients.
    expected_n <- if (infection == "FeLV") 630 else 638
    expected_cases <- if (infection == "FeLV") 126 else 134
    check_fit <- fit_model(dev, original_vars[[infection]])
    delta <- max(abs(coef(check_fit)[names(frozen_beta[[infection]])] - frozen_beta[[infection]]))
    if (nrow(dev) != expected_n || sum(dev$outcome) != expected_cases || !is.finite(delta) || delta > 1e-7) {
      stop("Development reconstruction differs from the saved model. Review the input before proceeding.")
    }
    # Table 3: evaluate frozen equations at their original cut-offs.
    for (stage in c("Development", "Testing")) {
      sample <- if (stage == "Development") dev else test
      p <- frozen_prediction(sample, infection)
      perf <- performance(sample$outcome, p, frozen_threshold[infection])
      add("Table_3_Fixed_models", cbind(infection, stage, perf))
      curves[[paste(infection, stage)]] <- list(roc = roc_object(sample$outcome, p), performance = perf)
    }
    # Tables 1 and 2: unadjusted Pearson tests, separately by outcome.
    population_vars <- if (infection == "FeLV") new_vars$FeLV else original_vars$FIV
    for (outcome in 0:1) for (v in population_vars) {
      s <- combined[combined$outcome == outcome, ]
      tab <- table(s[[v]], s$period)
      active <- tab[rowSums(tab) > 0, , drop = FALSE]
      ct <- chisq.test(active, correct = FALSE)
      for (level in rownames(tab)) for (period in colnames(tab)) {
        add(paste0("Table_", if (infection == "FeLV") "1" else "2", "_Population"),
          data.frame(infection, outcome, variable = v, category = level, period,
            n = unname(tab[level, period]), denominator = sum(tab[, period]),
            percent = 100 * tab[level, period] / sum(tab[, period]), p = ct$p.value,
            minimum_expected_count = min(ct$expected)))
      }
    }
    # S1 compares the joint sex/neuter specification with the sex-only model.
    if (infection == "FIV") for (version in c("Published_FIV", "Original_FIV")) {
      vars <- if (version == "Published_FIV") c("breed", "sex_neuter", "age_group", "lifestyle", "concomitant_disorder") else original_vars$FIV
      f <- fit_model(dev, vars)
      add("Table_S1_Coefficients", cbind(model = version, coefficients_table(f)))
      dev_threshold <- roc01(roc_object(dev$outcome, fitted(f)))
      add("Table_S1_Performance", cbind(model = version, stage = "Development", AIC = AIC(f),
        performance(dev$outcome, fitted(f), dev_threshold)))
      add("Table_S1_Performance", cbind(model = version, stage = "Testing", AIC = NA_real_,
        performance(test$outcome, predict(f, newdata = test, type = "response"), dev_threshold)))
    }
    # Each univariable test uses its own complete-case sample.
    for (v in candidate_vars) {
      s <- complete_sample(combined, v)
      f <- fit_model(s, v); null <- glm(outcome ~ 1, data = s, family = binomial(), na.action = na.fail)
      r_uv <- roc_object(s$outcome, fitted(f))
      ci_uv <- as.numeric(pROC::ci.auc(r_uv, method = "delong"))
      add("Table_S3_Global_screening", cbind(infection, variable = v, n = nrow(s), cases = sum(s$outcome), lrt(null, f),
        auc = as.numeric(pROC::auc(r_uv)), lower_95 = ci_uv[1], upper_95 = ci_uv[3]))
      add("Table_S3_Univariable_coefficients", cbind(infection, variable = v, n = nrow(s), coefficients_table(f)))
    }
    # Fit the specified original and updated models to combined data.
    specifications <- if (infection == "FeLV") list(Original_FeLV = original_vars$FeLV, New_FeLV = new_vars$FeLV) else list(New_FIV = new_vars$FIV)
    for (model in names(specifications)) {
      vars <- specifications[[model]]
      s <- complete_sample(combined, vars)
      f <- fit_model(s, vars); fits[[model]] <- f
      add(if (model == "Original_FeLV") "Table_S2_Coefficients" else "Table_S4_Coefficients",
          cbind(model, n = nrow(s), coefficients_table(f)))
      add("Combined_model_performance", cbind(model, AIC = AIC(f), BIC = BIC(f),
        logLik = as.numeric(logLik(f)), parameters = attr(logLik(f), "df"),
        performance(s$outcome, fitted(f))))
      if (model != "Original_FeLV") curves[[paste(infection, "Combined")]] <-
        list(roc = roc_object(s$outcome, fitted(f)), performance = performance(s$outcome, fitted(f)))
    }
    add("Cohort_counts", data.frame(infection, stage = c("Development", "Testing", "Combined_before_additional_requirements"),
      n = c(nrow(dev), nrow(test), nrow(combined)), cases = c(sum(dev$outcome), sum(test$outcome), sum(combined$outcome))))
    for (v in c(candidate_vars, "origin")) {
      for (period in levels(combined$period)) {
        s <- combined[combined$period == period, ]
        add("Variable_availability", data.frame(infection, period, variable = v, n = nrow(s),
          available = sum(!is.na(s[[v]])), missing = sum(is.na(s[[v]]))))
      }
    }
  }
  # Main Tables 4 and 5 include reference rows and category-level Wald tests.
  for (infection in c("FeLV", "FIV")) {
    model <- paste0("New_", infection); s <- cohort_list[[infection]]
    fit <- fits[[model]]; adjusted <- coefficients_table(fit)
    for (v in new_vars[[infection]]) {
      uv <- complete_sample(s, v); fu <- fit_model(uv, v); simple <- coefficients_table(fu)
      term_names <- names(coef(fu))[-1]
      lv <- levels(s[[v]])
      for (j in seq_along(lv)) {
        if (j == 1) {
          values <- data.frame(OR_simple = 1, lower_simple = NA, upper_simple = NA, p_simple = NA,
                               OR_multiple = 1, lower_multiple = NA, upper_multiple = NA, p_multiple = NA)
        } else {
          term <- term_names[j - 1]; u <- simple[simple$term == term, ]; a <- adjusted[adjusted$term == term, ]
          values <- data.frame(OR_simple = u$OR, lower_simple = u$lower_OR, upper_simple = u$upper_OR, p_simple = u$p,
                               OR_multiple = a$OR, lower_multiple = a$lower_OR, upper_multiple = a$upper_OR, p_multiple = a$p)
        }
        add(paste0("Table_", if (infection == "FeLV") "4" else "5", "_Regression"),
          cbind(variable = v, category = lv[j], reference = j == 1, n_simple = nrow(uv), n_multiple = nobs(fit), values))
      }
    }
  }
  for (name in names(tables)) write.csv(tables[[name]], file.path(folder, paste0(name, ".csv")), row.names = FALSE, na = "")
  list(tables = tables, curves = curves)
}

# 5. Execute the analyses and retain statistical warnings

warning_log <- character()
results <- withCallingHandlers(run_analysis(), warning = function(w) {
  warning_log <<- c(warning_log, conditionMessage(w))
  invokeRestart("muffleWarning")
})
writeLines(if (length(warning_log)) warning_log else "No warnings recorded.",
           file.path(output_dir, "Model_and_test_warnings.txt"))

# 6. Recreate Figure 1 (ROC curves)

draw_figure <- function() {
  old_par <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.8, 0.8), mgp = c(2.4, 0.7, 0), xaxs = "i", yaxs = "i")
  on.exit(par(old_par))
  colors <- c("#0072B2", "#E69F00", "#009E73")
  for (infection in c("FeLV", "FIV")) {
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "1 - Specificity", ylab = "Sensitivity", main = paste(infection, "clinical decision support models"), cex.main = 0.95)
    abline(h = seq(0, 1, 0.1), v = seq(0, 1, 0.1), col = "grey88", lty = 2)
    abline(0, 1, col = "grey65")
    axis(1, at = seq(0, 1, .1), cex.axis = .8); axis(2, at = seq(0, 1, .1), las = 1, cex.axis = .8); box()
    labels <- character()
    for (j in 1:3) {
      stage <- c("Development", "Testing", "Combined")[j]
      obj <- results$curves[[paste(infection, stage)]]; r <- obj$roc; m <- obj$performance
      # Retain empirical coordinates, including ties; reverse pROC's order.
      lines(rev(1 - r$specificities), rev(r$sensitivities), type = "s", col = colors[j], lwd = 2)
      points(1 - m$specificity, m$sensitivity, pch = 21, bg = colors[j], cex = 1.1)
      label <- c("Original model: 2013-2022", "Original model: 2023-2025", "New model: 2013-2025")[j]
      labels <- c(labels, paste0(label, "\n", sprintf("AUC = %.2f; 95%% CI %.2f-%.2f", m$auc, m$lower_95, m$upper_95)))
    }
    legend("bottomright", legend = labels, col = colors, lwd = 2, cex = .66,
           y.intersp = 1.8, bg = "white", box.col = "grey60", inset = .015)
    mtext(if (infection == "FeLV") "A" else "B", side = 3, adj = -0.09, line = 1, font = 2, cex = 1.2)
  }
}
pdf(file.path(output_dir, "Figure_1_ROC.pdf"), width = 12, height = 6, useDingbats = FALSE)
draw_figure(); dev.off()
png(file.path(output_dir, "Figure_1_ROC.png"), width = 3600, height = 1800, res = 300, bg = "white")
draw_figure(); dev.off()
message("Completed. Results: ", normalizePath(output_dir))
