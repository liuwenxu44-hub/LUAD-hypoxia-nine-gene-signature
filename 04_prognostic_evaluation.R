# ============================================================================
# 04_prognostic_evaluation.R
# Prognostic model evaluation: nomogram, calibration, C-index,
# clinical baseline table, univariate/multivariate Cox (Figure 3)
# Requires: TCGA_LUAD_RiskScore_Clinical.csv from 03_LASSO_Cox_model.R
# ============================================================================

library(survival)
library(rms)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(tableone)
library(survminer)
library(timeROC)

dat <- read.csv("TCGA_LUAD_RiskScore_Clinical.csv", row.names = 1)

# ============================================================================
# 1. Data cleaning
# ============================================================================

# Standardize column names
names(dat)[names(dat) == "ajcc_pathologic_stage"] <- "Pathological_stage"
names(dat)[names(dat) == "gender"]                <- "Gender"
names(dat)[names(dat) == "age_at_index"]          <- "Age"

# Clean pathological stage
clean_stage <- function(x) {
  x <- str_to_upper(str_trim(as.character(x)))
  case_when(
    is.na(x) | x %in% c("", "NA", "NOT REPORTED", "NOT AVAILABLE", "--") ~ NA_character_,
    str_detect(x, "STAGE IV")  ~ "Stage IV",
    str_detect(x, "STAGE III") ~ "Stage III",
    str_detect(x, "STAGE II")  ~ "Stage II",
    str_detect(x, "STAGE I")   ~ "Stage I",
    TRUE ~ x
  )
}
dat$Pathological_stage <- clean_stage(dat$Pathological_stage)

# Clean T/N/M if present
clean_tnm <- function(x, prefix) {
  x <- str_to_upper(str_trim(as.character(x)))
  case_when(
    is.na(x) | x %in% c("", "NA", "NX", "MX", "TX",
                         "NOT REPORTED", "--") ~ NA_character_,
    str_detect(x, paste0(prefix, "4")) ~ paste0(prefix, "4"),
    str_detect(x, paste0(prefix, "3")) ~ paste0(prefix, "3"),
    str_detect(x, paste0(prefix, "2")) ~ paste0(prefix, "2"),
    str_detect(x, paste0(prefix, "1")) ~ paste0(prefix, "1"),
    TRUE ~ x
  )
}

tnm_cols <- c("ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m")
tnm_new  <- c("Pathologic_T", "Pathologic_N", "Pathologic_M")
for (i in seq_along(tnm_cols)) {
  if (tnm_cols[i] %in% names(dat)) {
    dat[[tnm_new[i]]] <- clean_tnm(dat[[tnm_cols[i]]],
                                    substring(tnm_new[i], 11, 11))
  }
}

# ============================================================================
# 2. Univariate Cox forest plot (Figure 3A)
# ============================================================================

plot_df <- dat %>%
  mutate(
    Age = as.numeric(Age) / 365.25,
    Gender = as.factor(Gender),
    Risk_group = factor(Risk_Group, levels = c("Low", "High")),
    Pathological_stage = factor(
      gsub("Stage [IV]+A|Stage [IV]+B", "Stage \\1", Pathological_stage),
      levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
    )
  ) %>%
  filter(Pathological_stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV"))

covariates <- c("Risk_group", "Age", "Gender", "Pathological_stage")

univ_results <- lapply(covariates, function(var) {
  f <- as.formula(paste0("Surv(OS.time/365, OS) ~ ", var))
  s <- summary(coxph(f, data = plot_df))
  data.frame(
    Variable = rownames(s$coefficients),
    HR    = round(s$coefficients[, "exp(coef)"], 3),
    Lower = round(s$conf.int[, "lower .95"], 3),
    Upper = round(s$conf.int[, "upper .95"], 3),
    P_val = round(s$coefficients[, "Pr(>|z|)"], 3),
    stringsAsFactors = FALSE
  )
})

univ_df <- do.call(rbind, univ_results)
univ_df$Variable <- gsub("Risk_groupHigh", "Risk group (High vs Low)",
                          univ_df$Variable)
univ_df$Variable <- gsub("Gendermale|GenderMale", "Gender (Male vs Female)",
                          univ_df$Variable)
univ_df$Variable <- gsub("Pathological_stage", "", univ_df$Variable)

univ_df$P_show <- ifelse(univ_df$P_val < 0.001, "< 0.001",
                          sprintf("%.3f", univ_df$P_val))
univ_df$HR_CI  <- paste0(univ_df$HR, " (", univ_df$Lower, "-",
                          univ_df$Upper, ")")
write.csv(univ_df, "Univariate_Cox_Results.csv", row.names = FALSE)

univ_df$Variable <- factor(univ_df$Variable, levels = rev(univ_df$Variable))

p_univ <- ggplot(univ_df, aes(x = HR, y = Variable)) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.25, color = "black", linewidth = 0.8) +
  geom_point(shape = 15, size = 3.5, color = "#00468B") +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "#ED0000", linewidth = 0.8) +
  theme_classic() +
  labs(x = "Hazard Ratio (95% CI)", y = "",
       title = "Univariate Cox Regression Analysis") +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        plot.title  = element_text(hjust = 0.5, face = "bold", size = 14))

ggsave("Figure_3A_Univariate_Cox.pdf", p_univ, width = 10, height = 6)

# ============================================================================
# 3. Multivariate Cox forest plot (Figure 3B)
# ============================================================================

multi_cox <- coxph(Surv(OS.time/365, OS) ~ Risk_group + Age + Gender +
                    Pathological_stage, data = plot_df)

pdf("Figure_3B_Multivariate_Cox.pdf", width = 14, height = 10, onefile = FALSE)
print(ggforest(multi_cox, data = plot_df,
  main = "Multivariate Cox Regression Analysis",
  cpositions = c(0.02, 0.22, 0.4), fontsize = 1.0, noDigits = 3))
dev.off()

# ============================================================================
# 4. Nomogram (Figure 3C)
# ============================================================================

nomo_df <- plot_df[!is.na(plot_df$Risk_Score) & !is.na(plot_df$Age), ]
keep_cols <- c("OS.time", "OS", "Risk_Score", "Age", "Gender",
               "Pathological_stage")
nomo_df <- na.omit(nomo_df[, keep_cols])

dd <- datadist(nomo_df)
options(datadist = "dd")

f <- cph(Surv(OS.time/365, OS) ~ Risk_Score + Age + Gender + Pathological_stage,
         data = nomo_df, x = TRUE, y = TRUE, surv = TRUE)

surv <- Survival(f)
surv1 <- function(x) surv(1, x)
surv3 <- function(x) surv(3, x)
surv5 <- function(x) surv(5, x)

nomo <- nomogram(f,
  fun = list(surv1, surv3, surv5), lp = FALSE,
  funlabel = c("1-year OS Prob.", "3-year OS Prob.", "5-year OS Prob."),
  maxscale = 100,
  fun.at = c(0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05))

pdf("Figure_3C_Nomogram_Clinical_Model.pdf", width = 15, height = 5.5)
plot(nomo, xfrac = 0.3)
dev.off()

# ============================================================================
# 5. Calibration curves (Figure 3D)
# ============================================================================

nomo_df$OS_years <- nomo_df$OS.time / 365
m_val <- round(nrow(nomo_df) / 4)

cal1 <- calibrate(f, cmethod = "KM", method = "boot",
                   u = 1, m = m_val, B = 1000)
cal3 <- calibrate(f, cmethod = "KM", method = "boot",
                   u = 3, m = m_val, B = 1000)
cal5 <- calibrate(f, cmethod = "KM", method = "boot",
                   u = 5, m = m_val, B = 1000)

pdf("Calibration_Curve_MultiPanel.pdf", width = 15, height = 5.5)
par(mfrow = c(1, 3))

plot(cal1, lwd = 2, lty = 1, errbar.col = "#00468B", col = "#00468B",
     xlab = "Nomogram-Predicted 1-Year OS",
     ylab = "Actual 1-Year OS",
     xlim = c(0, 1), ylim = c(0, 1),
     subtitles = FALSE, main = "1-Year OS Calibration")
abline(0, 1, lty = 3, lwd = 2, col = "black")

plot(cal3, lwd = 2, lty = 1, errbar.col = "#ED0000", col = "#ED0000",
     xlab = "Nomogram-Predicted 3-Year OS",
     ylab = "Actual 3-Year OS",
     xlim = c(0, 1), ylim = c(0, 1),
     subtitles = FALSE, main = "3-Year OS Calibration")
abline(0, 1, lty = 3, lwd = 2, col = "black")

plot(cal5, lwd = 2, lty = 1, errbar.col = "#42B540", col = "#42B540",
     xlab = "Nomogram-Predicted 5-Year OS",
     ylab = "Actual 5-Year OS",
     xlim = c(0, 1), ylim = c(0, 1),
     subtitles = FALSE, main = "5-Year OS Calibration")
abline(0, 1, lty = 3, lwd = 2, col = "black")
dev.off()

# ============================================================================
# 6. Time-dependent ROC
# ============================================================================

time_points <- c(365, 1095, 1825)
ROC_res <- timeROC(
  T = clin_sub$OS.time, delta = clin_sub$OS,
  marker = clin_sub$Risk_Score, cause = 1,
  weighting = "marginal", times = time_points, iid = TRUE
)

roc_data <- data.frame(
  FPR  = c(ROC_res$FP[, 1], ROC_res$FP[, 2], ROC_res$FP[, 3]),
  TPR  = c(ROC_res$TP[, 1], ROC_res$TP[, 2], ROC_res$TP[, 3]),
  Time = factor(rep(c("1-Year", "3-Year", "5-Year"),
                     each = nrow(ROC_res$TP)))
)

auc_labels <- paste0(c("1-Year", "3-Year", "5-Year"),
                      " (AUC = ", round(ROC_res$AUC, 3), ")")
levels(roc_data$Time) <- auc_labels

p_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Time)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#E64B35", "#4DBBD5", "#00A087")) +
  theme_classic() +
  labs(title = "Time-dependent ROC of Hypoxia Signature",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       color = "Time Point") +
  theme(legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

ggsave("Time_dependent_ROC_Hypoxia_Signature.pdf", p_roc,
       width = 14, height = 10)

# ============================================================================
# 7. Model setup for C-index comparison
# ============================================================================

cindex_data <- dat %>%
  mutate(
    OS_year  = as.numeric(OS.time / 365),
    OS_event = as.numeric(OS),
    Age      = as.numeric(age_at_index),
    Gender   = as.factor(gender)
  ) %>%
  filter(!is.na(OS_year), OS_year > 0, !is.na(OS_event),
         !is.na(Risk_Score), !is.na(Age), !is.na(Gender))

cox_risk     <- coxph(Surv(OS_year, OS_event) ~ Risk_Score,
                       data = cindex_data, x = TRUE)
cox_clinical <- coxph(Surv(OS_year, OS_event) ~ Age + Gender + Pathological_stage,
                       data = cindex_data, x = TRUE)
cox_combined <- coxph(Surv(OS_year, OS_event) ~ Risk_Score + Age + Gender +
                       Pathological_stage, data = cindex_data, x = TRUE)

# ============================================================================
# 8. C-index comparison
# ============================================================================

get_cindex <- function(fit, model_name) {
  s <- summary(fit)$concordance
  data.frame(Model = model_name,
    C_index = s[1], SE = s[2],
    Lower_95CI = s[1] - 1.96 * s[2],
    Upper_95CI = s[1] + 1.96 * s[2])
}

cindex_table <- bind_rows(
  get_cindex(cox_risk, "Risk score model"),
  get_cindex(cox_clinical, "Clinical model"),
  get_cindex(cox_combined, "Combined model")
) %>%
  mutate(C_index_95CI = paste0(round(C_index, 3), " (",
         round(Lower_95CI, 3), "-", round(Upper_95CI, 3), ")"))

write.table(cindex_table[, c("Model", "C_index_95CI")],
  "Supplementary_Table_S3_C_index.txt", sep = "\t",
  quote = FALSE, row.names = FALSE)

# ============================================================================
# 9. Clinical baseline table (Supplementary Table S2)
# ============================================================================

baseline_data <- dat %>%
  mutate(
    Age = as.numeric(Age),
    Gender = factor(str_to_title(as.character(Gender))),
    Risk_Group = factor(Risk_Group, levels = c("Low", "High")),
    Pathological_stage = factor(Pathological_stage,
      levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
  )

cont_vars <- c("Age")
cat_vars  <- intersect(c("Gender", "Pathological_stage", "Pathologic_T",
                          "Pathologic_N", "Pathologic_M"),
                        names(baseline_data))
all_vars  <- c(cont_vars, cat_vars)

table_s2 <- CreateTableOne(
  vars = all_vars, strata = "Risk_Group",
  data = baseline_data, factorVars = cat_vars,
  addOverall = TRUE, includeNA = TRUE
)

table_s2_mat <- print(table_s2, nonnormal = cont_vars,
  exact = cat_vars, showAllLevels = TRUE,
  quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

table_s2_df <- data.frame(Characteristic = rownames(table_s2_mat),
                           table_s2_mat, check.names = FALSE)

write.csv(table_s2_df,
  "Supplementary_Table_S2_Clinical_characteristics_by_risk_group.csv",
  row.names = FALSE)

message("04_prognostic_evaluation.R done.")
