library(readr)
library(irlba)
library(Matrix)
library(dplyr)
library(wordspace)
library(scales)
library(aricode)
library(magrittr)
library(glmnet)
library(pROC)
library(ggplot2)
library(ggpubr)

to_laplacian <- function(A, regulariser=0) {
  row_degrees <- Matrix::rowSums(A)
  col_degrees <- Matrix::colSums(A)
  row_D <- 1 / sqrt(row_degrees + regulariser)
  col_D <- 1 / sqrt(col_degrees + regulariser)
  row_D[!is.finite(row_D)] <- 0
  col_D[!is.finite(col_D)] <- 0
  L <- Diagonal(x = row_D) %*% A %*% Diagonal(x = col_D)
  return(L)
}

cls_exmt <- function(X, y, train_prop=0.8) {
  sm <- sample(1:nrow(X), size = floor(nrow(X)*train_prop), replace=FALSE)
  model <- cv.glmnet(X[sm, ], y[sm], alpha = 0, family = "binomial")
  preds <- predict(model, X[-sm,], s='lambda.min')
  return(
    list(
      preds = c(preds),
      labels = y[-sm]
    )
  )
}

rocs_to_curve <- function(rocs,std_errs=2) {
  mc <- length(rocs)
  spec <- NULL
  sens <- NULL
  for (iter in 1:mc) {
    roc_li <- approx(x = rocs[[iter]]$specificities,
                    y = rocs[[iter]]$sensitivities,
                    n = 1000)
    spec <- rbind(spec, roc_li$x)
    sens <- rbind(sens, roc_li$y)
  }
  sds <- apply(sens, 2, sd)
  ses <- sds / sqrt(mc)
  sens_mean <- apply(sens, 2, mean)
  sens_upper <- sens_mean + std_errs * ses
  sens_lower <- sens_mean - std_errs * ses
  sens_df <- data.frame(sensitivity = sens_mean,
                        lower = sens_lower,
                        upper = sens_upper,
                        specificity = spec[1,])
  return(sens_df)
}

exmts_to_roc_df <- function(r_exmt) {
  full_100_roc_df <- rocs_to_curve(r_exmt$full_100)
  full_30_roc_df <- rocs_to_curve(r_exmt$full_30)
  core_100_roc_df <- rocs_to_curve(r_exmt$core_100)
  core_30_roc_df <- rocs_to_curve(r_exmt$core_30)
  cp_100_roc_df <- rocs_to_curve(r_exmt$cp_100)
  cp_30_roc_df <- rocs_to_curve(r_exmt$cp_30)
  full_100_roc_df$g <- "full (d=100)"
  full_30_roc_df$g <- "full (d=30)"
  core_30_roc_df$g <- "core (d=30)"
  core_100_roc_df$g <- "core (d=100)"
  cp_30_roc_df$g <- "core-periphery (d=30)"
  cp_100_roc_df$g <- "core-periphery (d=100)"
  roc_df <- do.call(rbind, list(full_100_roc_df,
                                full_30_roc_df,
                                core_100_roc_df,
                                core_30_roc_df,
                                cp_100_roc_df,
                                cp_30_roc_df))
  return(roc_df)
}

f_edge <- './edges-subj2-scan1.csv'
f_vertex <- './vertices-subj2-scan1.csv'
f_cp_vs_full_plot <- './cp_vs_full_plot.pdf'
f_cp_vs_core_plot <- './cp_vs_core_plot.pdf'

brain_edge_data <- readr::read_delim(f_edge, delim=' ',
                                     col_names = c('V1','V2'),
                                     col_types = list(
                                       V1 = col_character(),
                                       V2 = col_character()
                                     ))
brain_vertex_data <- readr::read_csv(f_vertex,
                                     col_types = list(
                                       name = col_character(),
                                       region = col_character(),
                                       Y = col_character()
                                     ))

region <- as.factor(pull(brain_vertex_data, region))
hemisphere <- as.factor(pull(brain_vertex_data, hemisphere))
tissue <- as.factor(pull(brain_vertex_data, tissue))
names(region) <- names(hemisphere) <- names(tissue) <- pull(brain_vertex_data, name)

v <- unique(pull(brain_vertex_data, "name"))
id <- 1:length(v); names(id) <- v
n <- length(v)

A <- sparseMatrix(i = id[pull(brain_edge_data, 1)],
                  j = id[pull(brain_edge_data, 2)],
                  symmetric = TRUE,
                  dims = c(n, n))
rownames(A) <- colnames(A) <- v
degrees <- rowSums(A)

L <- to_laplacian(A)

e_full <- irlba(L, 100)
X_full_100 <- e_full$u %*% diag(sqrt(e_full$d))
X_full_30 <- e_full$u[,1:30] %*% diag(sqrt(e_full$d[1:30]))
y <- tissue

mc <- 100
train_prop <- 0.8
exmt <- list()
regions = unique(region)
top_regions <- names(sort(table(region), decreasing = TRUE))
for (r in top_regions[1:6]) {
  print(paste0('Region ', r))
  idx <- region == r
  m <- sum(idx)
  A_core <- A[idx, idx]
  A_cp <- A[idx, ]
  L_core <- to_laplacian(A_core)
  L_cp <- to_laplacian(A_cp)
  y_r <- y[idx]
  e_core <- irlba(L_core, 100)
  e_cp <- irlba(L_cp, 100)
  X_core_100 <- e_core$u %*% diag(sqrt(e_core$d))
  X_core_30 <- e_core$u[,1:30] %*% diag(sqrt(e_core$d[1:30]))
  X_cp_100 <- e_cp$u %*% diag(sqrt(e_cp$d))
  X_cp_30 <- e_cp$u[,1:30] %*% diag(sqrt(e_cp$d[1:30]))
  # 
  r_exmt <- list(
    full_100 = list(),
    full_30 = list(),
    core_100 = list(),
    core_30 = list(),
    cp_100 = list(),
    cp_30 = list()
  )
  for (iter in 1:mc) {
    print(paste0(iter,' / ', mc))
    cls_full_100 <- cls_exmt(X_full_100[idx,], y[idx], train_prop)
    cls_full_30 <- cls_exmt(X_full_30[idx,], y[idx], train_prop)
    cls_core_100 <- cls_exmt(X_core_100, y[idx], train_prop)
    cls_core_30 <- cls_exmt(X_core_30, y[idx], train_prop)
    cls_cp_100 <- cls_exmt(X_cp_100, y[idx], train_prop)
    cls_cp_30 <- cls_exmt(X_cp_30, y[idx], train_prop)
    r_exmt$full_100[[iter]] <- roc(predictor = cls_full_100$preds,
                                   response = cls_full_100$labels,
                                   levels = c('gray', 'white'),
                                   direction = "<")
    r_exmt$full_30[[iter]] <- roc(predictor = cls_full_30$preds,
                                  response = cls_full_30$labels,
                                  levels = c('gray', 'white'),
                                  direction = "<")
    r_exmt$core_100[[iter]] <- roc(predictor = cls_core_100$preds,
                                   response = cls_core_100$labels,
                                   levels = c('gray', 'white'),
                                   direction = "<")
    r_exmt$core_30[[iter]] <- roc(predictor = cls_core_30$preds,
                                  response = cls_core_30$labels,
                                  levels = c('gray', 'white'),
                                  direction = "<")
    r_exmt$cp_100[[iter]] <- roc(predictor = cls_cp_100$preds,
                                 response = cls_cp_100$labels,
                                 levels = c('gray', 'white'),
                                 direction = "<")
    r_exmt$cp_30[[iter]] <- roc(predictor = cls_cp_30$preds,
                                 response = cls_cp_30$labels,
                                 levels = c('gray', 'white'),
                                 direction = "<")
  }
  exmt[[r]] <- r_exmt
}

# Plot core-periphery vs full
exmt_rocs <- lapply(exmt, exmts_to_roc_df)
for (r in names(exmt_rocs)) {
  exmt_rocs[[r]] <- exmt_rocs[[r]][!(exmt_rocs[[r]]$g %in% c("core (d=100)", "core (d=30)")), ]
}
theme_set(theme_minimal())
main_lwd <- 0.5
sec_lwd <- 0.1 * main_lwd
plts <- list()
rs <- names(exmt_rocs)
for (r in rs) {
  plts[[r]] <- exmt_rocs[[r]] %>%
    ggplot(aes(x = specificity, y = sensitivity, group = g)) +
    scale_x_reverse() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = g), alpha = 0.2) +
    geom_line(aes(specificity, lower, color = g), size = sec_lwd) +
    geom_line(aes(specificity, upper, color = g), size = sec_lwd) +
    geom_line(aes(color = g), size = main_lwd) +
    geom_abline(slope = 1, intercept = 1, color = "gray90") +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      plot.margin = unit(c(1.5, 0.5, 1, 0.5), "lines"),
      axis.title.x = element_text(size=8),
      axis.title.y = element_text(size=8),
      axis.text.x = element_text(size=7),
      axis.text.y = element_text(size=7)
    )
}
p <- ggarrange(plts[[1]],
  plts[[2]],
  plts[[3]],
  plts[[4]],
  plts[[5]],
  plts[[6]],
  labels = paste0("Region ", rs),
  common.legend = TRUE,
  legend = "bottom",
  font.label = list(size = 9.5,
                    face = "plain")
  # vjust=-0.8
)
ggsave(f_cp_vs_full_plot, p, "pdf", width = 7, height = 5.5)

# Plot core-periphery vs core
exmt_rocs <- lapply(exmt, exmts_to_roc_df)
for (r in names(exmt_rocs)) {
  exmt_rocs[[r]] <- exmt_rocs[[r]][!(exmt_rocs[[r]]$g %in% c("full (d=100)", "full (d=30)")), ]
  exmt_rocs[[r]]$g <- factor(exmt_rocs[[r]]$g, levels = c("core-periphery (d=30)", "core-periphery (d=100)", "core (d=30)", "core (d=100)"))
}
theme_set(theme_minimal())
main_lwd <- 0.5
sec_lwd <- 0.1 * main_lwd
plts <- list()
rs <- names(exmt_rocs)
for (r in rs) {
  plts[[r]] <- exmt_rocs[[r]] %>%
    ggplot(aes(x = specificity, y = sensitivity, group = g)) +
    scale_x_reverse() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = g), alpha = 0.2) +
    geom_line(aes(specificity, lower, color = g), size = sec_lwd) +
    geom_line(aes(specificity, upper, color = g), size = sec_lwd) +
    geom_line(aes(color = g), size = main_lwd) +
    geom_abline(slope = 1, intercept = 1, color = "gray90") +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      plot.margin = unit(c(1.5, 0.5, 1, 0.5), "lines"),
      axis.title.x = element_text(size=8),
      axis.title.y = element_text(size=8),
      axis.text.x = element_text(size=7),
      axis.text.y = element_text(size=7)
    )
}
p <- ggarrange(plts[[1]],
  plts[[2]],
  plts[[3]],
  plts[[4]],
  plts[[5]],
  plts[[6]],
  labels = paste0("Region ", rs),
  common.legend = TRUE,
  legend = "bottom",
  font.label = list(size = 9.5,
                    face = "plain")
  # vjust=-0.8
)
ggsave(f_cp_vs_core_plot, p, "pdf", width = 7, height = 5.5)