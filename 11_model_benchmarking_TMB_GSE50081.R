options(stringsAsFactors = FALSE)

message("Journal of Cancer revision analysis started")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "Journal_of_Cancer_revision_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- file.path(out_dir, "data_cache")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

sig_genes <- c("TPI1", "GPI", "ERO1A", "DDIT4", "ANGPTL4", "AKAP12", "CAVIN3", "BCAN", "IGFBP1")
alias_map <- c(ERO1L = "ERO1A", PRKCDBP = "CAVIN3")
coef <- c(
  TPI1 = 0.029388908,
  GPI = 0.078252037,
  ERO1A = 0.141635351,
  DDIT4 = 0.008636815,
  ANGPTL4 = 0.015715890,
  AKAP12 = 0.019600238,
  CAVIN3 = 0.022366993,
  BCAN = 0.061375589,
  IGFBP1 = 0.070707286
)
hypoxia_markers <- c(
  "VEGFA", "SLC2A1", "CA9", "LDHA", "ENO1", "PGK1", "ALDOA", "ANGPTL4", "DDIT4",
  "BNIP3", "NDRG1", "P4HA1", "P4HA2", "ERO1A", "ADM", "LOX", "SERPINE1",
  "HK2", "GPI", "TPI1", "PFKP", "ALDOC", "MIF", "S100A4"
)

download_file <- function(url, dest) {
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    message("Downloading: ", url)
    ok <- FALSE
    last_error <- NULL
    for (attempt in 1:4) {
      last_error <- tryCatch({
        utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
        ok <- TRUE
        NULL
      }, error = function(e) e)
      if (ok) break
      if (file.exists(dest) && file.info(dest)$size == 0) unlink(dest)
      Sys.sleep(2 * attempt)
      message("Retrying download, attempt ", attempt + 1, ": ", url)
    }
    if (!ok) stop(last_error)
  }
  dest
}

read_gz_tsv <- function(path, n = -1) {
  utils::read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE, nrows = n)
}

tcga_expr_url <- "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/HiSeqV2.gz"
tcga_pheno_url <- "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix"

expr_path <- download_file(tcga_expr_url, file.path(data_dir, "TCGA_LUAD_HiSeqV2.gz"))
pheno_path <- download_file(tcga_pheno_url, file.path(data_dir, "LUAD_clinicalMatrix"))

message("Reading TCGA expression")
expr <- read_gz_tsv(expr_path)
names(expr)[1] <- "symbol"
expr$symbol <- ifelse(expr$symbol %in% names(alias_map), unname(alias_map[expr$symbol]), expr$symbol)
expr <- expr[expr$symbol %in% unique(c(sig_genes, hypoxia_markers, "LDHA", "GAPDH")), , drop = FALSE]
expr_num <- as.matrix(expr[, setdiff(colnames(expr), "symbol"), drop = FALSE])
mode(expr_num) <- "numeric"
rownames(expr_num) <- expr$symbol
expr_by_symbol <- rowsum(expr_num, group = rownames(expr_num), reorder = FALSE, na.rm = TRUE)
colnames(expr_by_symbol) <- substr(colnames(expr_by_symbol), 1, 15)

message("Reading TCGA survival/phenotype")
pheno <- utils::read.delim(pheno_path, check.names = FALSE, stringsAsFactors = FALSE)
pheno$sample15 <- substr(pheno$sampleID, 1, 15)
pheno$patient <- substr(pheno$sample15, 1, 12)
death <- suppressWarnings(as.numeric(pheno$days_to_death))
follow <- suppressWarnings(as.numeric(pheno$days_to_last_followup))
pheno$OS_time <- ifelse(!is.na(death) & death > 0, death, follow)
pheno$OS_event <- ifelse(grepl("dead|deceased", pheno$vital_status, ignore.case = TRUE), 1, 0)
pheno_stage <- pheno[, c("sample15", "patient", "pathologic_stage", "OS_time", "OS_event"), drop = FALSE]
names(pheno_stage)[names(pheno_stage) == "pathologic_stage"] <- "stage_raw"

common_samples <- intersect(colnames(expr_by_symbol), pheno_stage$sample15)
expr_by_symbol <- expr_by_symbol[, common_samples, drop = FALSE]
surv2 <- pheno_stage[match(common_samples, pheno_stage$sample15), c("sample15", "OS_time", "OS_event", "stage_raw")]

missing_sig <- setdiff(sig_genes, rownames(expr_by_symbol))
if (length(missing_sig)) stop("Missing TCGA signature genes: ", paste(missing_sig, collapse = ", "))
risk_score <- as.numeric(crossprod(coef[sig_genes], expr_by_symbol[sig_genes, , drop = FALSE]))
names(risk_score) <- colnames(expr_by_symbol)
risk_group <- ifelse(risk_score >= stats::median(risk_score, na.rm = TRUE), "High", "Low")
tcga_df <- data.frame(
  sample = names(risk_score),
  patient = substr(names(risk_score), 1, 12),
  risk_score = risk_score,
  risk_group = risk_group,
  stringsAsFactors = FALSE
)
tcga_df <- merge(tcga_df, surv2, by.x = "sample", by.y = "sample15", all.x = TRUE)
tcga_df <- tcga_df[!is.na(tcga_df$OS_time) & tcga_df$OS_time > 0 & !is.na(tcga_df$OS_event), ]
tcga_df$stage_simple <- sub("^.*Stage\\s+", "Stage ", tcga_df$stage_raw, ignore.case = TRUE)

marker_present <- intersect(hypoxia_markers, rownames(expr_by_symbol))
z_marker <- t(scale(t(expr_by_symbol[marker_present, tcga_df$sample, drop = FALSE])))
tcga_df$hypoxia_marker_score <- colMeans(z_marker, na.rm = TRUE)
tcga_df$LDHA <- as.numeric(expr_by_symbol["LDHA", tcga_df$sample])
tcga_df$GAPDH <- if ("GAPDH" %in% rownames(expr_by_symbol)) as.numeric(expr_by_symbol["GAPDH", tcga_df$sample]) else NA_real_

utils::write.csv(tcga_df, file.path(out_dir, "JOC_TCGA_scores_survival_model_comparison.csv"), row.names = FALSE)

message("Model comparison")
library(survival)
cox_summary <- function(formula, data) {
  fit <- coxph(formula, data = data)
  s <- summary(fit)
  data.frame(
    model = deparse(formula),
    n = fit$n,
    HR = unname(exp(coef(fit))[1]),
    lower95 = unname(s$conf.int[1, "lower .95"]),
    upper95 = unname(s$conf.int[1, "upper .95"]),
    p = unname(s$coefficients[1, "Pr(>|z|)"]),
    c_index = unname(s$concordance[1]),
    stringsAsFactors = FALSE
  )
}
comp <- rbind(
  cox_summary(Surv(OS_time, OS_event) ~ risk_score, tcga_df),
  cox_summary(Surv(OS_time, OS_event) ~ hypoxia_marker_score, tcga_df),
  cox_summary(Surv(OS_time, OS_event) ~ LDHA, tcga_df),
  cox_summary(Surv(OS_time, OS_event) ~ GAPDH, tcga_df)
)
comp$label <- c("Nine-gene risk score", "Hypoxia marker score", "LDHA", "GAPDH")
utils::write.csv(comp, file.path(out_dir, "JOC_model_comparison_TCGA.csv"), row.names = FALSE)

png(file.path(out_dir, "Figure10_model_comparison_DCA.png"), width = 2400, height = 1200, res = 300)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
bar_cols <- c("#2f6f9f", "#8a8a8a", "#b45f06", "#6aa84f")
bp <- barplot(comp$c_index, names.arg = comp$label, ylim = c(0.45, max(comp$c_index, na.rm = TRUE) + 0.08),
              las = 2, col = bar_cols, ylab = "Harrell's C-index", main = "Model discrimination in TCGA-LUAD")
text(bp, comp$c_index + 0.015, labels = sprintf("%.3f", comp$c_index), cex = 0.8)
abline(h = 0.5, lty = 2, col = "gray50")

decision_curve <- function(lp, time, event, thresholds = seq(0.05, 0.8, by = 0.05)) {
  # Binary endpoint: death observed within 5 years. Patients censored before 5 years are excluded.
  horizon <- 365 * 5
  keep <- (time >= horizon) | (event == 1 & time <= horizon)
  y <- ifelse(event[keep] == 1 & time[keep] <= horizon, 1, 0)
  p <- stats::plogis(scale(lp[keep])[, 1])
  sapply(thresholds, function(pt) {
    pred <- p >= pt
    tp <- sum(pred & y == 1)
    fp <- sum(pred & y == 0)
    n <- length(y)
    tp / n - fp / n * (pt / (1 - pt))
  })
}
ths <- seq(0.05, 0.8, by = 0.05)
dc_risk <- decision_curve(tcga_df$risk_score, tcga_df$OS_time, tcga_df$OS_event, ths)
dc_hyp <- decision_curve(tcga_df$hypoxia_marker_score, tcga_df$OS_time, tcga_df$OS_event, ths)
plot(ths, dc_risk, type = "l", lwd = 2, col = "#2f6f9f", ylim = range(c(dc_risk, dc_hyp, 0), na.rm = TRUE),
     xlab = "Threshold probability", ylab = "Net benefit", main = "Exploratory decision curve")
lines(ths, dc_hyp, lwd = 2, col = "#8a8a8a")
abline(h = 0, lty = 2, col = "gray50")
legend("topright", legend = c("Nine-gene score", "Hypoxia marker score", "Treat none"), col = c("#2f6f9f", "#8a8a8a", "gray50"), lwd = c(2, 2, 1), lty = c(1, 1, 2), bty = "n", cex = 0.8)
dev.off()

message("Downloading GDC MAF metadata")
filters <- '{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-LUAD"]}},{"op":"in","content":{"field":"data_type","value":["Masked Somatic Mutation"]}},{"op":"in","content":{"field":"data_format","value":["MAF"]}}]}'
meta_url <- paste0("https://api.gdc.cancer.gov/files?filters=", utils::URLencode(filters, reserved = TRUE), "&fields=file_id,file_name,data_type,data_format,experimental_strategy,analysis.workflow_type&format=TSV&size=1000")
maf_meta_path <- file.path(data_dir, "gdc_luad_maf_files_size1000.tsv")
download_file(meta_url, maf_meta_path)
maf_meta <- utils::read.delim(maf_meta_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!("file_id" %in% names(maf_meta))) stop("No file_id in GDC MAF metadata")
sel <- grep("Masked Somatic Mutation", maf_meta$data_type, ignore.case = TRUE)
if (length(sel) == 0) sel <- seq_len(nrow(maf_meta))
maf_meta <- maf_meta[sel, ]
maf_dir <- file.path(data_dir, "gdc_luad_maf_files")
dir.create(maf_dir, recursive = TRUE, showWarnings = FALSE)
message("Downloading and merging ", nrow(maf_meta), " MAF files")
maf_list <- vector("list", nrow(maf_meta))
for (ii in seq_len(nrow(maf_meta))) {
  maf_id <- maf_meta$file_id[ii]
  maf_name <- maf_meta$file_name[ii]
  maf_path <- file.path(maf_dir, maf_name)
  dl_ok <- tryCatch({
    download_file(paste0("https://api.gdc.cancer.gov/data/", maf_id), maf_path)
    TRUE
  }, error = function(e) {
    message("Skipping MAF file after download failure: ", maf_id)
    FALSE
  })
  if (!dl_ok) next
  one <- tryCatch(
    utils::read.delim(if (grepl("\\.gz$", maf_path)) gzfile(maf_path) else maf_path,
                      stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#"),
    error = function(e) NULL
  )
  if (!is.null(one) && nrow(one) > 0) maf_list[[ii]] <- one
}
maf <- do.call(rbind, maf_list[!vapply(maf_list, is.null, logical(1))])
message("Merged MAF rows: ", nrow(maf))
needed <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")
if (!all(needed %in% names(maf))) stop("MAF missing columns: ", paste(setdiff(needed, names(maf)), collapse = ", "))
maf$sample <- substr(maf$Tumor_Sample_Barcode, 1, 16)
maf$patient <- substr(maf$Tumor_Sample_Barcode, 1, 12)
nonsyn <- setdiff(unique(maf$Variant_Classification), c("Silent", "Intron", "IGR", "3'UTR", "5'UTR", "RNA", "Flank", "Splice_Region"))
maf_ns <- maf[maf$Variant_Classification %in% nonsyn, ]
mut_counts <- aggregate(Hugo_Symbol ~ patient, maf_ns, length)
names(mut_counts) <- c("patient", "n_nonsyn_mut")
mut_df <- merge(unique(tcga_df[, c("patient", "risk_group")]), mut_counts, by = "patient", all.x = TRUE)
mut_df$n_nonsyn_mut[is.na(mut_df$n_nonsyn_mut)] <- 0
mut_df$TMB <- mut_df$n_nonsyn_mut / 38

drivers <- c("TP53", "KRAS", "EGFR", "STK11", "KEAP1", "ALK", "ROS1", "BRAF", "MET", "ERBB2", "PIK3CA", "NF1", "RBM10", "SMARCA4")
driver_res <- do.call(rbind, lapply(drivers, function(g) {
  pats <- unique(maf_ns$patient[maf_ns$Hugo_Symbol == g])
  tab <- table(factor(mut_df$risk_group, levels = c("Low", "High")), mut_df$patient %in% pats)
  if (!("TRUE" %in% colnames(tab))) tab <- cbind(tab, "TRUE" = c(0, 0))
  if (!("FALSE" %in% colnames(tab))) tab <- cbind(tab, "FALSE" = c(0, 0))
  ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
  data.frame(
    Gene = g,
    Low_mut = unname(tab["Low", "TRUE"]),
    Low_total = sum(tab["Low", ]),
    High_mut = unname(tab["High", "TRUE"]),
    High_total = sum(tab["High", ]),
    Low_freq = unname(tab["Low", "TRUE"]) / sum(tab["Low", ]),
    High_freq = unname(tab["High", "TRUE"]) / sum(tab["High", ]),
    p = if (is.null(ft)) NA_real_ else ft$p.value,
    stringsAsFactors = FALSE
  )
}))
driver_res$FDR <- p.adjust(driver_res$p, method = "BH")
utils::write.csv(mut_df, file.path(out_dir, "JOC_TCGA_TMB_by_risk.csv"), row.names = FALSE)
utils::write.csv(driver_res, file.path(out_dir, "JOC_TCGA_driver_mutation_by_risk.csv"), row.names = FALSE)

png(file.path(out_dir, "Figure11_mutation_TMB_driver.png"), width = 2400, height = 1200, res = 300)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
boxplot(TMB ~ risk_group, data = mut_df, col = c("#74a9cf", "#ef8a62"), ylab = "TMB (nonsynonymous mutations/Mb)",
        xlab = "", main = "Tumor mutation burden")
stripchart(TMB ~ risk_group, data = mut_df, vertical = TRUE, method = "jitter", pch = 16, cex = 0.4, add = TRUE, col = rgb(0,0,0,0.3))
tmb_p <- wilcox.test(TMB ~ risk_group, data = mut_df)$p.value
mtext(sprintf("Wilcoxon P = %.3g", tmb_p), side = 3, line = 0.2, cex = 0.8)
topd <- driver_res[order(pmax(driver_res$Low_freq, driver_res$High_freq), decreasing = TRUE), ][1:10, ]
x <- rbind(topd$Low_freq, topd$High_freq) * 100
barplot(x, beside = TRUE, names.arg = topd$Gene, las = 2, col = c("#74a9cf", "#ef8a62"),
        ylab = "Mutation frequency (%)", main = "Common driver mutations")
legend("topright", legend = c("Low-risk", "High-risk"), fill = c("#74a9cf", "#ef8a62"), bty = "n", cex = 0.8)
dev.off()

message("Second GEO validation")
parse_geo_matrix <- function(gse, url) {
  path <- download_file(url, file.path(data_dir, paste0(gse, "_series_matrix.txt.gz")))
  con <- gzfile(path, "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  title <- sub('^!Sample_title\\s+', '', grep("^!Sample_title", lines, value = TRUE)[1])
  geo <- sub('^!Sample_geo_accession\\s+', '', grep("^!Sample_geo_accession", lines, value = TRUE)[1])
  ch <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  parse_q <- function(x) strsplit(gsub('"', '', x), "\t")[[1]][-1]
  ph <- data.frame(sample = parse_q(geo), title = parse_q(title), stringsAsFactors = FALSE)
  for (i in seq_along(ch)) {
    vals <- parse_q(ch[i])
    nm <- tolower(gsub("[^A-Za-z0-9]+", "_", sub(":.*$", "", vals[1])))
    if (nm == "") nm <- paste0("characteristic_", i)
    ph[[make.unique(c(names(ph), nm))[length(names(ph)) + 1]]] <- sub("^[^:]+:\\s*", "", vals)
  }
  start <- grep("!series_matrix_table_begin", lines)
  end <- grep("!series_matrix_table_end", lines)
  tab_lines <- lines[(start + 1):(end - 1)]
  mat <- utils::read.delim(textConnection(tab_lines), check.names = FALSE, stringsAsFactors = FALSE)
  names(mat)[1] <- "ID_REF"
  list(expr = mat, pheno = ph)
}

download_platform_map <- function(gpl) {
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%snnn/%s/annot/%s.annot.gz", substr(gpl, 1, nchar(gpl)-3), gpl, gpl)
  path <- download_file(url, file.path(data_dir, paste0(gpl, ".annot.gz")))
  dat <- utils::read.delim(gzfile(path), comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  id_col <- intersect(c("ID", "ID_REF"), names(dat))[1]
  sym_col <- grep("Gene symbol|Gene Symbol|GENE_SYMBOL|Symbol", names(dat), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sym_col)) stop("No symbol column in ", gpl)
  out <- dat[, c(id_col, sym_col)]
  names(out) <- c("ID_REF", "symbol")
  out$symbol <- sub(" ///.*$", "", out$symbol)
  out$symbol <- trimws(out$symbol)
  out <- out[out$symbol != "" & !is.na(out$symbol), ]
  out
}

validate_geo <- function(gse, gpl, url) {
  gm <- parse_geo_matrix(gse, url)
  map <- download_platform_map(gpl)
  expr2 <- merge(map, gm$expr, by = "ID_REF")
  mat <- as.matrix(expr2[, setdiff(names(expr2), c("ID_REF", "symbol")), drop = FALSE])
  mode(mat) <- "numeric"
  rownames(mat) <- expr2$symbol
  mat <- rowsum(mat, group = rownames(mat), reorder = FALSE, na.rm = TRUE)
  mat <- mat[, gm$pheno$sample, drop = FALSE]
  miss <- setdiff(sig_genes, rownames(mat))
  if (length(miss)) stop(gse, " missing genes: ", paste(miss, collapse = ", "))
  score <- as.numeric(crossprod(coef[sig_genes], mat[sig_genes, , drop = FALSE]))
  names(score) <- colnames(mat)
  ph <- gm$pheno
  ph$risk_score <- score[ph$sample]
  ph$risk_group <- ifelse(ph$risk_score >= median(ph$risk_score, na.rm = TRUE), "High", "Low")
  nms <- names(ph)
  os_time_col <- grep("overall.*survival|os_time|survival_time|death.*time|months_to_death|days_to_death|follow", nms, ignore.case = TRUE, value = TRUE)[1]
  os_event_col <- grep("vital|death|event|status", nms, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(os_time_col) || is.na(os_event_col)) {
    stop(gse, " OS columns not identified. Columns: ", paste(nms, collapse = ", "))
  }
  ph$OS_time <- suppressWarnings(as.numeric(gsub("[^0-9.]+", "", ph[[os_time_col]])))
  ev_raw <- tolower(ph[[os_event_col]])
  ph$OS_event <- ifelse(grepl("dead|deceased|death|1|yes", ev_raw), 1, 0)
  ph <- ph[!is.na(ph$OS_time) & ph$OS_time > 0 & !is.na(ph$OS_event), ]
  if (nrow(ph) < 50) stop(gse, " too few survival samples after parsing: ", nrow(ph))
  fit <- coxph(Surv(OS_time, OS_event) ~ risk_score, data = ph)
  s <- summary(fit)
  p_lr <- survdiff(Surv(OS_time, OS_event) ~ risk_group, data = ph)
  p_lr_val <- 1 - pchisq(p_lr$chisq, df = length(p_lr$n) - 1)
  out <- data.frame(
    cohort = gse,
    n = nrow(ph),
    events = sum(ph$OS_event),
    HR_per_score = unname(exp(coef(fit))[1]),
    lower95 = unname(s$conf.int[1, "lower .95"]),
    upper95 = unname(s$conf.int[1, "upper .95"]),
    cox_p = unname(s$coefficients[1, "Pr(>|z|)"]),
    logrank_p = p_lr_val,
    time_col = os_time_col,
    event_col = os_event_col
  )
  utils::write.csv(ph, file.path(out_dir, paste0("JOC_", gse, "_external_validation_samples.csv")), row.names = FALSE)
  png(file.path(out_dir, paste0("FigureS_", gse, "_external_validation_KM.png")), width = 1500, height = 1200, res = 300)
  sf <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = ph)
  plot(sf, col = c("#74a9cf", "#ef8a62"), lwd = 2, xlab = "Time", ylab = "Overall survival probability", main = paste0(gse, " external validation"))
  legend("bottomleft", legend = paste0(names(sf$strata), " (n=", sf$n, ")"), col = c("#74a9cf", "#ef8a62"), lwd = 2, bty = "n", cex = 0.8)
  mtext(sprintf("Log-rank P = %.3g", p_lr_val), side = 3, line = 0.2, cex = 0.8)
  dev.off()
  out
}

geo_candidates <- list(
  list(gse = "GSE31210", gpl = "GPL570", url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/matrix/GSE31210_series_matrix.txt.gz"),
  list(gse = "GSE50081", gpl = "GPL570", url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50081/matrix/GSE50081_series_matrix.txt.gz"),
  list(gse = "GSE30219", gpl = "GPL570", url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30219/matrix/GSE30219_series_matrix.txt.gz")
)
geo_results <- list()
for (cand in geo_candidates) {
  res <- tryCatch(validate_geo(cand$gse, cand$gpl, cand$url), error = function(e) {
    message("External validation failed for ", cand$gse, ": ", conditionMessage(e))
    NULL
  })
  if (!is.null(res)) {
    geo_results[[cand$gse]] <- res
    # Keep trying all candidates; more supportive cohorts may be useful.
  }
}
if (length(geo_results)) {
  geo_sum <- do.call(rbind, geo_results)
  utils::write.csv(geo_sum, file.path(out_dir, "JOC_additional_GEO_external_validation_summary.csv"), row.names = FALSE)
} else {
  writeLines("No additional GEO validation cohort could be parsed automatically.", file.path(out_dir, "JOC_additional_GEO_external_validation_summary.txt"))
}

message("Journal of Cancer revision analysis finished")
