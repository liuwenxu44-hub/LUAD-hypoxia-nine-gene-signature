options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "Medicine_supportive_outputs")

tcga <- read.csv(file.path(out_dir, "Medicine_TCGA_scores_survival_model_comparison.csv"), check.names = FALSE)
comp <- read.csv(file.path(out_dir, "Medicine_model_comparison_TCGA.csv"), check.names = FALSE)
mut <- read.csv(file.path(out_dir, "Medicine_TCGA_TMB_by_risk.csv"), check.names = FALSE)
driver <- read.csv(file.path(out_dir, "Medicine_TCGA_driver_mutation_by_risk.csv"), check.names = FALSE)
gse <- read.csv(file.path(out_dir, "Medicine_GSE50081_scores_survival.csv"), check.names = FALSE)
gse_res <- read.csv(file.path(out_dir, "Medicine_GSE50081_supportive_evaluation.csv"), check.names = FALSE)

library(survival)

draw_fig10 <- function() {
  par(mar = c(6, 5, 3, 1), oma = c(0, 0, 0, 0), cex = 0.9)
  comp$short <- c("Nine-gene", "Hypoxia\nmarker", "LDHA", "GAPDH")
  cols <- c("#2f6f9f", "#8a8a8a", "#b45f06", "#6aa84f")
  yy <- comp$c_index
  xx <- seq_along(yy)
  plot(xx, yy, type = "n", xaxt = "n", xlab = "", ylab = "Harrell's C-index",
       xlim = c(0.55, 4.45), ylim = c(0.49, 0.65), main = "TCGA-LUAD discrimination")
  abline(h = 0.5, lty = 2, col = "gray55")
  segments(xx, 0.5, xx, yy, col = cols, lwd = 5)
  points(xx, yy, pch = 21, bg = cols, col = "black", cex = 1.8)
  axis(1, at = xx, labels = c("Nine-gene", "Hypoxia", "LDHA", "GAPDH"), cex.axis = 0.78)
  text(xx, yy + 0.006, labels = sprintf("%.3f", yy), cex = 0.85)
  mtext("A", side = 3, line = 1, at = par("usr")[1], adj = 0, font = 2, cex = 1.2)
}

draw_fig11 <- function() {
  par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), oma = c(0, 0, 0, 0), cex = 0.9)
  mut$risk_group <- factor(mut$risk_group, levels = c("Low", "High"))
  boxplot(log2(TMB + 1) ~ risk_group, data = mut, col = c("#74a9cf", "#ef8a62"),
          ylab = "log2(TMB + 1)", xlab = "", main = "Tumor mutation burden")
  stripchart(log2(TMB + 1) ~ risk_group, data = mut, vertical = TRUE, method = "jitter",
             pch = 16, cex = 0.35, add = TRUE, col = rgb(0, 0, 0, 0.25))
  tmb_p <- wilcox.test(TMB ~ risk_group, data = mut)$p.value
  mtext(sprintf("Wilcoxon P = %.3g", tmb_p), side = 3, line = 0.2, cex = 0.8)
  mtext("A", side = 3, line = 1, at = par("usr")[1], adj = 0, font = 2, cex = 1.2)

  driver <- driver[order(pmax(driver$Low_freq, driver$High_freq), decreasing = TRUE), ][1:10, ]
  x <- rbind(driver$Low_freq, driver$High_freq) * 100
  bp <- barplot(x, beside = TRUE, names.arg = driver$Gene, las = 2,
                ylim = c(0, max(x, na.rm = TRUE) + 8),
                col = c("#74a9cf", "#ef8a62"), ylab = "Mutation frequency (%)",
                main = "Common driver mutations")
  legend("topright", legend = c("Low-risk", "High-risk"),
         fill = c("#74a9cf", "#ef8a62"), bty = "n", cex = 0.82)
  sig <- driver$FDR < 0.05
  if (any(sig, na.rm = TRUE)) {
    mid <- colMeans(bp)
    text(mid[sig], apply(x[, sig, drop = FALSE], 2, max) + 4, labels = "*", cex = 1.2)
  }
  mtext("B", side = 3, line = 1, at = par("usr")[1], adj = 0, font = 2, cex = 1.2)
}

draw_figs1 <- function() {
  par(mar = c(5, 5, 3, 1), cex = 0.9)
  gse$risk_group <- factor(gse$risk_group, levels = c("Low", "High"))
  km <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = gse)
  plot(km, col = c("#74a9cf", "#ef8a62"), lwd = 2,
       xlab = "Time (days)", ylab = "Overall survival probability",
       main = "GSE50081 supportive evaluation")
  legend("bottomleft", legend = c("Low-risk", "High-risk"),
         col = c("#74a9cf", "#ef8a62"), lwd = 2, bty = "n")
  gr <- gse_res[gse_res$model == "high versus low risk", ]
  text(max(gse$OS_time, na.rm = TRUE) * 0.46, 0.22,
       labels = sprintf("HR = %.2f (95%% CI %.2f-%.2f), P = %.3g",
                        gr$HR, gr$lower95, gr$upper95, gr$p))
}

save_all <- function(stem, draw_fun, width, height) {
  png(file.path(out_dir, paste0(stem, ".png")), width = width * 300, height = height * 300, res = 300)
  draw_fun()
  dev.off()
  tiff(file.path(out_dir, paste0(stem, ".tif")), width = width, height = height, units = "in", res = 600, compression = "lzw")
  draw_fun()
  dev.off()
  pdf(file.path(out_dir, paste0(stem, ".pdf")), width = width, height = height)
  draw_fun()
  dev.off()
}

save_all("Figure10_model_comparison", draw_fig10, 4.8, 3.8)
save_all("Figure11_mutation_TMB_driver", draw_fig11, 7.6, 3.8)
save_all("FigureS1_GSE50081_external_validation", draw_figs1, 4.8, 3.8)

message("Revision figures regenerated")
