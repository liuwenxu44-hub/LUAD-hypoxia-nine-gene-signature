options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "Medicine_supportive_outputs")
data_dir <- file.path(out_dir, "data_cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

library(survival)

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

download_file <- function(url, dest) {
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    message("Downloading: ", url)
    utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  }
  dest
}

parse_quoted <- function(line) strsplit(gsub('"', '', line), "\t", fixed = TRUE)[[1]][-1]

parse_geo_matrix <- function(path) {
  con <- gzfile(path, "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)

  ph <- data.frame(
    sample = parse_quoted(grep("^!Sample_geo_accession", lines, value = TRUE)[1]),
    title = parse_quoted(grep("^!Sample_title", lines, value = TRUE)[1]),
    stringsAsFactors = FALSE
  )
  ch <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  used_names <- names(ph)
  for (ln in ch) {
    vals <- parse_quoted(ln)
    field <- sub(":.*$", "", vals[1])
    nm <- tolower(gsub("[^A-Za-z0-9]+", "_", field))
    nm <- sub("_+$", "", nm)
    if (!nzchar(nm)) nm <- paste0("characteristic_", length(used_names))
    nm <- make.unique(c(used_names, nm))[length(used_names) + 1]
    used_names <- c(used_names, nm)
    ph[[nm]] <- sub("^[^:]+:\\s*", "", vals)
  }

  start <- grep("!series_matrix_table_begin", lines)
  end <- grep("!series_matrix_table_end", lines)
  tab_lines <- lines[(start + 1):(end - 1)]
  mat <- utils::read.delim(textConnection(tab_lines), check.names = FALSE, stringsAsFactors = FALSE)
  names(mat)[1] <- "ID_REF"
  list(expr = mat, pheno = ph)
}

read_gpl570_map <- function() {
  gpl_path <- file.path(data_dir, "GPL570.annot.gz")
  download_file("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz", gpl_path)
  con <- gzfile(gpl_path, "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  start <- grep("!platform_table_begin", lines)
  end <- grep("!platform_table_end", lines)
  if (!length(start)) stop("Cannot find GPL570 platform table")
  if (!length(end)) end <- length(lines) + 1
  gpl <- utils::read.delim(textConnection(lines[(start + 1):(end - 1)]),
                           check.names = FALSE, stringsAsFactors = FALSE)
  id_col <- grep("^ID$|ID_REF", names(gpl), ignore.case = TRUE, value = TRUE)[1]
  sym_col <- grep("Gene symbol|Gene Symbol|gene_assignment|Symbol", names(gpl), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id_col) || is.na(sym_col)) stop("Cannot find GPL570 ID or gene symbol columns")
  sym <- as.character(gpl[[sym_col]])
  sym <- sub("///.*$", "", sym)
  sym <- sub("//.*$", "", sym)
  sym <- sub(";.*$", "", sym)
  sym <- trimws(sym)
  sym <- toupper(sym)
  sym[sym %in% names(alias_map)] <- unname(alias_map[sym[sym %in% names(alias_map)]])
  data.frame(ID_REF = as.character(gpl[[id_col]]), symbol = sym, stringsAsFactors = FALSE)
}

message("Parsing GSE50081")
gse_path <- file.path(data_dir, "GSE50081_series_matrix.txt.gz")
if (!file.exists(gse_path) || file.info(gse_path)$size < 50000000) {
  stop("GSE50081 matrix is missing or incomplete: ", gse_path)
}
geo <- parse_geo_matrix(gse_path)
map <- read_gpl570_map()

expr <- merge(map, geo$expr, by = "ID_REF")
expr <- expr[nzchar(expr$symbol) & expr$symbol %in% sig_genes, ]
if (!nrow(expr)) stop("No signature probes mapped in GSE50081")
expr_num <- as.matrix(expr[, setdiff(names(expr), c("ID_REF", "symbol")), drop = FALSE])
mode(expr_num) <- "numeric"
rownames(expr_num) <- expr$symbol
expr_symbol <- rowsum(expr_num, group = rownames(expr_num), reorder = FALSE, na.rm = TRUE)
probe_counts <- table(rownames(expr_num))
expr_symbol <- expr_symbol / as.numeric(probe_counts[rownames(expr_symbol)])

missing <- setdiff(sig_genes, rownames(expr_symbol))
if (length(missing)) stop("Missing GSE50081 signature genes: ", paste(missing, collapse = ", "))

ph <- geo$pheno
names(ph) <- tolower(names(ph))
histology_col <- grep("^histology$", names(ph), value = TRUE)[1]
time_col <- grep("^survival_time$", names(ph), value = TRUE)[1]
status_col <- grep("^status$", names(ph), value = TRUE)[1]
if (is.na(histology_col) || is.na(time_col) || is.na(status_col)) {
  stop("Required clinical fields not found in GSE50081")
}

ph$histology_clean <- tolower(ph[[histology_col]])
keep_luad <- grepl("adenocarcinoma", ph$histology_clean) & !grepl("squamous|large cell", ph$histology_clean)
ph <- ph[keep_luad, ]
ph$OS_time <- suppressWarnings(as.numeric(ph[[time_col]])) * 365.25
ph$OS_event <- ifelse(grepl("dead", ph[[status_col]], ignore.case = TRUE), 1, 0)
ph <- ph[!is.na(ph$OS_time) & ph$OS_time > 0 & !is.na(ph$OS_event), ]

common <- intersect(ph$sample, colnames(expr_symbol))
ph <- ph[match(common, ph$sample), ]
expr_symbol <- expr_symbol[, common, drop = FALSE]

risk_score <- as.numeric(crossprod(coef[sig_genes], expr_symbol[sig_genes, , drop = FALSE]))
ph$risk_score <- risk_score
ph$risk_group <- ifelse(ph$risk_score >= median(ph$risk_score, na.rm = TRUE), "High", "Low")

fit_cont <- coxph(Surv(OS_time, OS_event) ~ risk_score, data = ph)
fit_group <- coxph(Surv(OS_time, OS_event) ~ relevel(factor(risk_group), ref = "Low"), data = ph)
s_cont <- summary(fit_cont)
s_group <- summary(fit_group)
res <- rbind(
  data.frame(
    cohort = "GSE50081",
    endpoint = "overall survival",
    model = "continuous nine-gene score",
    n = fit_cont$n,
    events = sum(ph$OS_event),
    HR = unname(exp(coef(fit_cont))[1]),
    lower95 = unname(s_cont$conf.int[1, "lower .95"]),
    upper95 = unname(s_cont$conf.int[1, "upper .95"]),
    p = unname(s_cont$coefficients[1, "Pr(>|z|)"]),
    c_index = unname(s_cont$concordance[1])
  ),
  data.frame(
    cohort = "GSE50081",
    endpoint = "overall survival",
    model = "high versus low risk",
    n = fit_group$n,
    events = sum(ph$OS_event),
    HR = unname(exp(coef(fit_group))[1]),
    lower95 = unname(s_group$conf.int[1, "lower .95"]),
    upper95 = unname(s_group$conf.int[1, "upper .95"]),
    p = unname(s_group$coefficients[1, "Pr(>|z|)"]),
    c_index = unname(s_group$concordance[1])
  )
)

utils::write.csv(ph, file.path(out_dir, "Medicine_GSE50081_scores_survival.csv"), row.names = FALSE)
utils::write.csv(res, file.path(out_dir, "Medicine_GSE50081_supportive_evaluation.csv"), row.names = FALSE)

png(file.path(out_dir, "FigureS1_GSE50081_external_validation.png"), width = 1800, height = 1400, res = 300)
par(mar = c(5, 5, 3, 1))
km <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = ph)
plot(km, col = c("#74a9cf", "#ef8a62"), lwd = 2, xlab = "Time (days)", ylab = "Overall survival probability",
     main = "GSE50081 external validation")
legend("bottomleft", legend = c("Low-risk", "High-risk"), col = c("#74a9cf", "#ef8a62"), lwd = 2, bty = "n")
text_x <- max(ph$OS_time, na.rm = TRUE) * 0.45
text_y <- 0.22
text(text_x, text_y, labels = sprintf("HR = %.2f (95%% CI %.2f-%.2f), P = %.3g",
                                      res$HR[res$model == "high versus low risk"],
                                      res$lower95[res$model == "high versus low risk"],
                                      res$upper95[res$model == "high versus low risk"],
                                      res$p[res$model == "high versus low risk"]))
dev.off()

message("GSE50081 validation completed")
