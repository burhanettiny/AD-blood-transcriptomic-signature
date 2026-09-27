## =================================================================
## Blood-Based Transcriptomic Signature for Alzheimer's Disease
## Full reproducible pipeline (discovery: GSE63060, validation: GSE63061)
## Son güncelleme: bu script, makalenin FINAL halinde raporlanan
## TÜM sayıları üretir (nested CV, bootstrap, kalibrasyon, demografik
## karşılaştırma, mito-ribozomal korelasyon dahil).
## =================================================================

## ---- BLOCK 1: Packages ----
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs_bioc <- c("GEOquery", "limma", "illuminaHumanv4.db", "clusterProfiler", "org.Hs.eg.db", "preprocessCore")
pkgs_cran <- c("randomForest", "pROC", "e1071", "ggplot2", "dplyr", "pheatmap")

for (p in pkgs_bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
for (p in pkgs_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(GEOquery); library(limma); library(illuminaHumanv4.db)
library(clusterProfiler); library(org.Hs.eg.db); library(preprocessCore)
library(randomForest); library(pROC); library(e1071)
library(dplyr); library(pheatmap); library(ggplot2)

set.seed(42)

## Sonuçların/görsellerin kaydedileceği klasör -- ihtiyaca göre değiştirin
OUTDIR <- "C:/temp_r"
dir.create(OUTDIR, showWarnings = FALSE)
setwd(OUTDIR)
getwd()

## ---- BLOCK 2: Discovery cohort (GSE63060) ----
gse60 <- getGEO("GSE63060", GSEMatrix = TRUE, getGPL = FALSE)
gse60 <- gse60[[1]]
expr60 <- exprs(gse60)
pheno60 <- pData(gse60)
cat("GSE63060 dims:", dim(gse60), "\n")
cat("Expression range:", range(expr60, na.rm = TRUE), "\n")   # 7.20-15.05 -> zaten log2, tekrar log2 ALMAYIN

## ---- BLOCK 3: Quantile normalization (discovery only, ön-inceleme için) ----
expr60_norm <- normalizeQuantiles(expr60)

## ---- BLOCK 4: Diagnosis / age / sex ----
## KRİTİK: doğru kolon "status:ch1" -- "status" kolonuyla KARIŞTIRMAYIN
## (o kolon GEO'nun genel yayın durumunu tutar, tanıyı değil)
diag_col <- "status:ch1"
print(table(pheno60[[diag_col]]))   # AD=145 CTL=104 MCI=80 bekleniyor

group <- rep(NA, nrow(pheno60))
group[grepl("^AD$|Alzheimer", pheno60[[diag_col]], ignore.case = TRUE)] <- "AD"
group[grepl("^CTL$|control|healthy", pheno60[[diag_col]], ignore.case = TRUE)] <- "CTL"
group[grepl("MCI|mild cognitive", pheno60[[diag_col]], ignore.case = TRUE)] <- "MCI"
pheno60$group <- factor(group, levels = c("CTL", "AD", "MCI"))
print(table(pheno60$group))

pheno60$age <- as.numeric(gsub("[^0-9.]", "", pheno60[["age:ch1"]]))
pheno60$sex <- factor(gsub(".*:\\s*", "", pheno60[["gender:ch1"]]))

## ---- BLOCK 5: Probe -> Gene Symbol annotation ----
probe_ids <- rownames(expr60_norm)
sym_map <- AnnotationDbi::select(illuminaHumanv4.db, keys = probe_ids,
                                  columns = "SYMBOL", keytype = "PROBEID")
sym_map <- sym_map[!duplicated(sym_map$PROBEID), ]
sym_map <- sym_map[!is.na(sym_map$SYMBOL), ]
expr60_annot <- expr60_norm[sym_map$PROBEID, ]

## ---- BLOCK 6: Differential Expression (AD vs CTL) ----
adctl_idx <- pheno60$group %in% c("AD", "CTL")
expr_adctl <- expr60_annot[, adctl_idx]
pheno_adctl <- pheno60[adctl_idx, ]
pheno_adctl$group <- droplevels(pheno_adctl$group)

## Unadjusted (ana model, 100-gen seçimi için kullanılır)
design0 <- model.matrix(~ group, data = pheno_adctl)
fit0 <- eBayes(lmFit(expr_adctl, design0))
top0 <- topTable(fit0, coef = "groupAD", number = Inf, adjust.method = "BH")
top0$Probe_ID <- rownames(top0)
top0 <- merge(top0, sym_map, by.x = "Probe_ID", by.y = "PROBEID")
cat("Unadjusted logFC range:", range(top0$logFC), "\n")   # sanity check: -1 ile +1 civarı olmalı

## Age-adjusted
design_age <- model.matrix(~ group + age, data = pheno_adctl)
fit_age <- eBayes(lmFit(expr_adctl, design_age))
top_age <- topTable(fit_age, coef = "groupAD", number = Inf, adjust.method = "BH")
top_age$Probe_ID <- rownames(top_age)

## Sex-adjusted
design_sex <- model.matrix(~ group + sex, data = pheno_adctl)
fit_sex <- eBayes(lmFit(expr_adctl, design_sex))
top_sex <- topTable(fit_sex, coef = "groupAD", number = Inf, adjust.method = "BH")
top_sex$Probe_ID <- rownames(top_sex)

## ---- BLOCK 7: Top 100 gene signature ----
sig_genes <- top0 %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.5)
top100 <- sig_genes %>% arrange(adj.P.Val) %>% slice_head(n = 100)
top100_probes <- top100$Probe_ID
cat("Top100 genes:", nrow(top100), "\n")

## Age/sex-adjusted overlap (makalede raporlanan 94/100 ve 100/100)
sig_age <- top_age %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.5)
sig_sex <- top_sex %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.5)
cat("Age-adjusted overlap:", sum(top100_probes %in% sig_age$Probe_ID), "/100\n")
cat("Sex-adjusted overlap:", sum(top100_probes %in% sig_sex$Probe_ID), "/100\n")

## ---- BLOCK 8: Figure 1 - Heatmap ----
mat_top100 <- expr_adctl[top100_probes, ]
mat_scaled <- t(scale(t(mat_top100)))
annot_col <- data.frame(Group = pheno_adctl$group)
rownames(annot_col) <- colnames(mat_scaled)

png("Fig1_heatmap.png", width = 1800, height = 2200, res = 220)
pheatmap(mat_scaled, annotation_col = annot_col, show_colnames = FALSE,
         show_rownames = TRUE, fontsize_row = 6,
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         main = "Figure 1. Differential gene expression profile (AD vs Control)")
dev.off()

## ---- BLOCK 9: Random Forest training ----
rf_data <- as.data.frame(t(mat_top100))
rf_data$group <- pheno_adctl$group
rf_model <- randomForest(group ~ ., data = rf_data, ntree = 500,
                          mtry = floor(sqrt(ncol(rf_data) - 1)), importance = TRUE)
print(rf_model)

imp <- importance(rf_model, type = 1)
imp_df <- data.frame(Probe_ID = rownames(imp), MDA = imp[, 1]) %>%
  arrange(desc(MDA)) %>% mutate(Rank = row_number())
imp_df <- merge(imp_df, sym_map, by.x = "Probe_ID", by.y = "PROBEID") %>% arrange(Rank)
print(head(imp_df, 20))

png("Fig3_varimp.png", width = 1600, height = 1400, res = 200)
varImpPlot(rf_model, n.var = 20, main = "Top 20 Predictive Genes (MDA)")
dev.off()

## ---- BLOCK 10: Internal ROC (RF / SVM / GLM) ----
rf_votes <- predict(rf_model, type = "prob")[, "AD"]
roc_rf_internal <- roc(rf_data$group, rf_votes, levels = c("CTL", "AD"), direction = "<")

svm_model <- svm(group ~ ., data = rf_data, kernel = "linear", probability = TRUE)
svm_pred <- predict(svm_model, rf_data, probability = TRUE)
svm_prob <- attr(svm_pred, "probabilities")[, "AD"]
roc_svm_internal <- roc(rf_data$group, svm_prob, levels = c("CTL", "AD"), direction = "<")

glm_data <- rf_data
glm_data$group <- ifelse(glm_data$group == "AD", 1, 0)
glm_model <- suppressWarnings(glm(group ~ ., data = glm_data, family = "binomial"))
glm_prob <- predict(glm_model, type = "response")
roc_glm_internal <- roc(glm_data$group, glm_prob)

cat("Internal AUC - RF:", auc(roc_rf_internal), " SVM:", auc(roc_svm_internal),
    " GLM:", auc(roc_glm_internal), "\n")

png("Fig2_internal_ROC.png", width = 1400, height = 1400, res = 200)
plot(roc_rf_internal, col = "orange", main = "Figure 2. Internal ROC comparison")
lines(roc_svm_internal, col = "green"); lines(roc_glm_internal, col = "blue")
legend("bottomright", legend = c(
  paste0("RF (AUC=", round(auc(roc_rf_internal), 3), ")"),
  paste0("SVM (AUC=", round(auc(roc_svm_internal), 3), ")"),
  paste0("GLM (AUC=", round(auc(roc_glm_internal), 3), ")")),
  col = c("orange", "green", "blue"), lty = 1)
dev.off()

## ---- BLOCK 11: Nested cross-validation (data leakage kontrolü) ----
k <- 10
n <- ncol(expr_adctl)
folds <- sample(rep(1:k, length.out = n))
nested_probs <- rep(NA, n)
names(nested_probs) <- colnames(expr_adctl)

for (i in 1:k) {
  test_idx <- which(folds == i)
  train_idx <- which(folds != i)
  expr_train <- expr_adctl[, train_idx]
  pheno_train <- pheno_adctl[train_idx, ]
  expr_test  <- expr_adctl[, test_idx]

  design_fold <- model.matrix(~ group, data = pheno_train)
  fit_fold <- eBayes(lmFit(expr_train, design_fold))
  top_fold <- topTable(fit_fold, coef = "groupAD", number = Inf, adjust.method = "BH")
  top_fold$Probe_ID <- rownames(top_fold)
  sig_fold <- top_fold[top_fold$adj.P.Val < 0.05 & abs(top_fold$logFC) > 0.5, ]
  sig_fold <- sig_fold[order(sig_fold$adj.P.Val), ]
  fold_genes <- head(sig_fold$Probe_ID, 100)

  rf_train_data <- as.data.frame(t(expr_train[fold_genes, ]))
  rf_train_data$group <- pheno_train$group
  rf_fold <- randomForest(group ~ ., data = rf_train_data, ntree = 500,
                           mtry = floor(sqrt(length(fold_genes))))

  rf_test_data <- as.data.frame(t(expr_test[fold_genes, , drop = FALSE]))
  nested_probs[test_idx] <- predict(rf_fold, rf_test_data, type = "prob")[, "AD"]
}

roc_nested <- roc(pheno_adctl$group, nested_probs, levels = c("CTL", "AD"), direction = "<")
cat("Nested CV AUC:", auc(roc_nested), " 95% CI:", ci.auc(roc_nested)[1], "-", ci.auc(roc_nested)[3], "\n")

## ---- BLOCK 12: Bootstrap stability of top biomarkers ----
n_boot <- 100
boot_top10 <- vector("list", n_boot)
for (b in 1:n_boot) {
  boot_idx <- sample(1:ncol(mat_top100), replace = TRUE)
  boot_data <- as.data.frame(t(mat_top100[, boot_idx]))
  boot_data$group <- pheno_adctl$group[boot_idx]
  if (length(unique(boot_data$group)) < 2) next
  rf_boot <- randomForest(group ~ ., data = boot_data, ntree = 300,
                           mtry = floor(sqrt(ncol(boot_data) - 1)), importance = TRUE)
  imp_boot <- importance(rf_boot, type = 1)
  boot_top10[[b]] <- rownames(imp_boot)[order(-imp_boot[, 1])][1:10]
}
freq_table <- sort(table(unlist(boot_top10)), decreasing = TRUE)
freq_df <- data.frame(Probe_ID = names(freq_table), Times_in_Top10 = as.integer(freq_table),
                       Pct = round(100 * as.integer(freq_table) / n_boot, 1))
freq_df <- merge(freq_df, sym_map, by.x = "Probe_ID", by.y = "PROBEID") %>% arrange(desc(Times_in_Top10))
print(head(freq_df[, c("SYMBOL", "Times_in_Top10", "Pct")], 15))
write.csv(freq_df, "Bootstrap_stability_top10.csv", row.names = FALSE)

## ---- BLOCK 13: GO enrichment ----
gene_list_real <- unique(top100$SYMBOL)
ego_real <- enrichGO(gene = gene_list_real, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                      ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego_df <- as.data.frame(ego_real)
cat("Enriched BP terms (q<0.05):", nrow(ego_df), "\n")
write.csv(ego_df, "GO_enrichment.csv", row.names = FALSE)

p5 <- dotplot(ego_real, showCategory = 15, title = "Figure 5. GO enrichment (AD signature)") +
  theme(axis.text.y = element_text(size = 10))
png("Fig5_GO_dotplot.png", width = 2600, height = 1300, res = 220)
print(p5)
dev.off()

## ---- BLOCK 14: GSE63061 - external validation cohort ----
gse61 <- getGEO("GSE63061", GSEMatrix = TRUE, getGPL = FALSE)
gse61 <- gse61[[1]]
expr61 <- exprs(gse61)
pheno61 <- pData(gse61)
cat("Expression range GSE63061:", range(expr61, na.rm = TRUE), "\n")  # zaten log2

diag_col61 <- "status:ch1"
group61 <- rep(NA, nrow(pheno61))
group61[grepl("^AD$|Alzheimer", pheno61[[diag_col61]], ignore.case = TRUE)] <- "AD"
group61[grepl("^CTL$|control|healthy", pheno61[[diag_col61]], ignore.case = TRUE)] <- "CTL"
pheno61$group <- factor(group61, levels = c("CTL", "AD"))
valid_idx <- !is.na(pheno61$group)
pheno61_valid <- pheno61[valid_idx, ]
pheno61_valid$age <- as.numeric(gsub("[^0-9.]", "", pheno61_valid[["age:ch1"]]))
pheno61_valid$sex <- factor(gsub(".*:\\s*", "", pheno61_valid[["gender:ch1"]]))

## ---- BLOCK 15 (ESKI YONTEM -- karsilastirma/referans icin saklandi):
## JOINT quantile normalization. UYARI: bu yontem discovery+external kohortlarini
## BIRLIKTE normalize eder, dolayisiyla tek bir yeni hastaya "tasinabilir" sekilde
## uygulanamaz (external kohortun tamami elde olmadan calismaz). Editor elestirisi
## tam olarak bu blokla ilgiliydi. Asagida BLOCK 15b-19b'de frozen-reference
## (sadece discovery'den turetilmis, dondurulmus hedef dagilim) alternatifi var --
## makalede asil raporlanacak olan odur. Bu blok sadece iki yontemi karsilastirmak
## icin (sensitivity analysis) calistirilmaya devam ediyor. ----
common_all_probes <- intersect(rownames(expr60), rownames(expr61))
combined_raw <- cbind(expr60[common_all_probes, ], expr61[common_all_probes, ])
combined_norm <- normalizeQuantiles(combined_raw)
n60 <- ncol(expr60)
expr60_joint <- combined_norm[, 1:n60]
expr61_joint <- combined_norm[, (n60 + 1):ncol(combined_norm)]

common_probes <- intersect(top100_probes, rownames(expr61_joint))
cat("Ortak probe sayısı:", length(common_probes), "/100\n")

mat_common_joint <- expr60_joint[common_probes, adctl_idx]
rf_data_joint <- as.data.frame(t(mat_common_joint))
rf_data_joint$group <- pheno_adctl$group
rf_model_joint <- randomForest(group ~ ., data = rf_data_joint, ntree = 500,
                                mtry = floor(sqrt(length(common_probes))), importance = TRUE)

svm_model_joint <- svm(group ~ ., data = rf_data_joint, kernel = "linear", probability = TRUE)

glm_data_joint <- rf_data_joint
glm_data_joint$group <- ifelse(glm_data_joint$group == "AD", 1, 0)
glm_model_joint <- suppressWarnings(glm(group ~ ., data = glm_data_joint, family = "binomial"))

## ---- BLOCK 16: External predictions + ROC ----
ext_data_joint <- as.data.frame(t(expr61_joint[common_probes, valid_idx]))

ext_probs_rf_j <- predict(rf_model_joint, ext_data_joint, type = "prob")[, "AD"]
roc_rf_ext_j <- roc(pheno61_valid$group, ext_probs_rf_j, levels = c("CTL", "AD"), direction = "<")

svm_pred_ext_j <- predict(svm_model_joint, ext_data_joint, probability = TRUE)
svm_prob_ext_j <- attr(svm_pred_ext_j, "probabilities")[, "AD"]
roc_svm_ext_j <- roc(pheno61_valid$group, svm_prob_ext_j, levels = c("CTL", "AD"), direction = "<")

glm_prob_ext_j <- predict(glm_model_joint, newdata = ext_data_joint, type = "response")
roc_glm_ext_j <- roc(pheno61_valid$group, glm_prob_ext_j, levels = c("CTL", "AD"), direction = "<")

cat("External AUC - RF:", auc(roc_rf_ext_j), " SVM:", auc(roc_svm_ext_j), " GLM:", auc(roc_glm_ext_j), "\n")
cat("External RF 95% CI:", ci.auc(roc_rf_ext_j)[1], "-", ci.auc(roc_rf_ext_j)[3], "\n")
cat("External SVM 95% CI:", ci.auc(roc_svm_ext_j)[1], "-", ci.auc(roc_svm_ext_j)[3], "\n")

png("Fig6_external_ROC.png", width = 1400, height = 1400, res = 200)
plot(roc_rf_ext_j, col = "orange", main = "Figure 6. External validation ROC (GSE63061)")
lines(roc_svm_ext_j, col = "green"); lines(roc_glm_ext_j, col = "blue")
legend("bottomright", legend = c(
  paste0("RF (AUC=", round(auc(roc_rf_ext_j), 3), ")"),
  paste0("SVM (AUC=", round(auc(roc_svm_ext_j), 3), ")"),
  paste0("GLM (AUC=", round(auc(roc_glm_ext_j), 3), ")")),
  col = c("orange", "green", "blue"), lty = 1)
dev.off()

## ---- BLOCK 17: Calibration + Brier score ----
brier_rf <- mean((ext_probs_rf_j - as.numeric(pheno61_valid$group == "AD"))^2)
cat("External RF Brier score:", brier_rf, "\n")

cal_df <- data.frame(pred = ext_probs_rf_j, obs = as.numeric(pheno61_valid$group == "AD"))
cal_df$bin <- cut(cal_df$pred, breaks = seq(0, 1, 0.2), include.lowest = TRUE)
cal_summary <- aggregate(cbind(pred, obs) ~ bin, data = cal_df, mean)
print(cal_summary)

png("Fig_calibration.png", width = 1200, height = 1200, res = 200)
plot(cal_summary$pred, cal_summary$obs, xlim = c(0, 1), ylim = c(0, 1), pch = 19,
     xlab = "Mean predicted probability", ylab = "Observed AD frequency",
     main = "External calibration (RF, GSE63061)")
abline(0, 1, col = "gray", lty = 2)
lines(cal_summary$pred, cal_summary$obs, col = "orange")
dev.off()

## ---- BLOCK 18: Demografik karşılaştırma (discovery vs external) ----
age_test <- t.test(pheno_adctl$age, pheno61_valid$age)
cat("Age discovery vs external, t-test p =", age_test$p.value, "\n")
sex_table <- table(c(rep("discovery", nrow(pheno_adctl)), rep("external", nrow(pheno61_valid))),
                    c(as.character(pheno_adctl$sex), as.character(pheno61_valid$sex)))
print(chisq.test(sex_table))

## ---- BLOCK 19: MCI risk stratification ----
mci_idx <- pheno60$group == "MCI"
expr_mci_joint <- expr60_joint[common_probes, mci_idx]
mci_data_joint <- as.data.frame(t(expr_mci_joint))
mci_probs_joint <- predict(rf_model_joint, mci_data_joint, type = "prob")[, "AD"]

rf_votes_joint <- predict(rf_model_joint, type = "prob")[, "AD"]
roc_rf_internal_joint <- roc(rf_data_joint$group, rf_votes_joint, levels = c("CTL", "AD"), direction = "<")
opt_cut_joint <- coords(roc_rf_internal_joint, "best", ret = "threshold")$threshold

mci_results <- data.frame(
  SampleID = colnames(expr_mci_joint), RiskScore = mci_probs_joint,
  Age = pheno60$age[mci_idx], Sex = pheno60$sex[mci_idx]
)
mci_results$Category <- ifelse(mci_results$RiskScore >= opt_cut_joint, "High Risk", "Low Risk")
cat("High risk:", sum(mci_results$Category == "High Risk"),
    " Low risk:", sum(mci_results$Category == "Low Risk"), "\n")

## MCI risk skoru ile yaş arasındaki ilişki (age-independence kontrolü)
print(cor.test(mci_results$RiskScore, mci_results$Age, method = "pearson"))

png("Fig4_MCI_histogram.png", width = 1400, height = 1000, res = 200)
hist(mci_results$RiskScore, breaks = 15, col = "orange",
     main = "Figure 4. AD-Risk stratification in MCI patients",
     xlab = "Probability of AD (Risk Score)")
dev.off()

write.csv(mci_results, "Table_3_MCI_Risk.csv", row.names = FALSE)

## =================================================================
## BLOCK 15b-19b (YENI -- TASINABILIR PIPELINE): Frozen-reference normalizasyon
## =================================================================
## Mantik: hedef dagilim SADECE discovery kohortundan (expr60) hesaplanir ve
## "dondurulur" (frozen). External kohort -- veya tek bir yeni hasta -- bu SABIT
## hedefe göre normalize edilir; discovery ve external ARTIK BIRLIKTE normalize
## EDILMIYOR. normalize.quantiles.use.target() her ornegi (kolonu) bagimsiz
## isleyebildigi icin, bu pipeline tek bir yeni hastaya bile calisir --
## "transportable" iddiasinin teknik karsiligi budur.

## ---- BLOCK 15b: Frozen reference hedef dagilimini SADECE discovery'den hesapla ----
common_all_probes <- intersect(rownames(expr60), rownames(expr61))
cat("Discovery-External ortak probe sayisi:", length(common_all_probes), "\n")

disc_common_raw <- as.matrix(expr60[common_all_probes, ])
frozen_target <- normalize.quantiles.determine.target(disc_common_raw)
saveRDS(frozen_target, "frozen_reference_target.rds")
cat("Frozen reference target kaydedildi -- ileride tek bir yeni hasta icin de kullanilabilir.\n")

## Discovery'yi bu dondurulmus hedefe göre normalize et, modeli bunun uzerinde egit
expr60_frozen <- normalize.quantiles.use.target(disc_common_raw, frozen_target)
dimnames(expr60_frozen) <- dimnames(disc_common_raw)

common_probes_frozen <- intersect(top100_probes, common_all_probes)
cat("Imza icindeki ortak probe sayisi (frozen):", length(common_probes_frozen), "/100\n")

mat_common_frozen <- expr60_frozen[common_probes_frozen, adctl_idx]
rf_data_frozen <- as.data.frame(t(mat_common_frozen))
rf_data_frozen$group <- pheno_adctl$group
rf_model_frozen <- randomForest(group ~ ., data = rf_data_frozen, ntree = 500,
                                 mtry = floor(sqrt(length(common_probes_frozen))), importance = TRUE)

svm_model_frozen <- svm(group ~ ., data = rf_data_frozen, kernel = "linear", probability = TRUE)

glm_data_frozen <- rf_data_frozen
glm_data_frozen$group <- ifelse(glm_data_frozen$group == "AD", 1, 0)
glm_model_frozen <- suppressWarnings(glm(group ~ ., data = glm_data_frozen, family = "binomial"))

## ---- BLOCK 16b: External kohortu SADECE frozen hedefe göre normalize et ----
## KRITIK: burada expr60 hic kullanilmiyor -- external ornekler discovery'nin
## kendisinden degil, discovery'den daha once turetilmis SABIT bir hedeften
## normalize ediliyor. Bu yuzden n_external=1 olsa bile ayni sonuc elde edilir.
ext_common_raw <- as.matrix(expr61[common_all_probes, valid_idx])
expr61_frozen <- normalize.quantiles.use.target(ext_common_raw, frozen_target)
dimnames(expr61_frozen) <- dimnames(ext_common_raw)

ext_data_frozen <- as.data.frame(t(expr61_frozen[common_probes_frozen, ]))

ext_probs_rf_f <- predict(rf_model_frozen, ext_data_frozen, type = "prob")[, "AD"]
roc_rf_ext_f <- roc(pheno61_valid$group, ext_probs_rf_f, levels = c("CTL", "AD"), direction = "<")

svm_pred_ext_f <- predict(svm_model_frozen, ext_data_frozen, probability = TRUE)
svm_prob_ext_f <- attr(svm_pred_ext_f, "probabilities")[, "AD"]
roc_svm_ext_f <- roc(pheno61_valid$group, svm_prob_ext_f, levels = c("CTL", "AD"), direction = "<")

glm_prob_ext_f <- predict(glm_model_frozen, newdata = ext_data_frozen, type = "response")
roc_glm_ext_f <- roc(pheno61_valid$group, glm_prob_ext_f, levels = c("CTL", "AD"), direction = "<")

cat("\n[FROZEN-REFERENCE] External AUC - RF:", auc(roc_rf_ext_f), " SVM:", auc(roc_svm_ext_f),
    " GLM:", auc(roc_glm_ext_f), "\n")
cat("[FROZEN-REFERENCE] External RF 95% CI:", ci.auc(roc_rf_ext_f)[1], "-", ci.auc(roc_rf_ext_f)[3], "\n")
cat("[FROZEN-REFERENCE] External SVM 95% CI:", ci.auc(roc_svm_ext_f)[1], "-", ci.auc(roc_svm_ext_f)[3], "\n")

png("Fig6_external_ROC_frozen.png", width = 1400, height = 1400, res = 200)
plot(roc_rf_ext_f, col = "orange", main = "Figure 6. External validation ROC (frozen-reference)")
lines(roc_svm_ext_f, col = "green"); lines(roc_glm_ext_f, col = "blue")
legend("bottomright", legend = c(
  paste0("RF (AUC=", round(auc(roc_rf_ext_f), 3), ")"),
  paste0("SVM (AUC=", round(auc(roc_svm_ext_f), 3), ")"),
  paste0("GLM (AUC=", round(auc(roc_glm_ext_f), 3), ")")),
  col = c("orange", "green", "blue"), lty = 1)
dev.off()

## ---- BLOCK 17b: Kalibrasyon + Brier skoru (frozen-reference) ----
brier_rf_f <- mean((ext_probs_rf_f - as.numeric(pheno61_valid$group == "AD"))^2)
cat("[FROZEN-REFERENCE] External RF Brier score:", brier_rf_f, "\n")

cal_df_f <- data.frame(pred = ext_probs_rf_f, obs = as.numeric(pheno61_valid$group == "AD"))
cal_df_f$bin <- cut(cal_df_f$pred, breaks = seq(0, 1, 0.2), include.lowest = TRUE)
cal_summary_f <- aggregate(cbind(pred, obs) ~ bin, data = cal_df_f, mean)
print(cal_summary_f)

png("Fig_calibration_frozen.png", width = 1200, height = 1200, res = 200)
plot(cal_summary_f$pred, cal_summary_f$obs, xlim = c(0, 1), ylim = c(0, 1), pch = 19,
     xlab = "Mean predicted probability", ylab = "Observed AD frequency",
     main = "External calibration - frozen reference (RF, GSE63061)")
abline(0, 1, col = "gray", lty = 2)
lines(cal_summary_f$pred, cal_summary_f$obs, col = "orange")
dev.off()

## ---- BLOCK 19b: MCI risk stratification -- frozen-reference model ile ----
## NOT: artik MCI hastalarini kendi (discovery-ici) kohortuna gore degil, ayni
## dondurulmus hedefe gore normalize ediyoruz -- boylece MCI skorlari da
## GSE63061'in o an elde olup olmamasindan bagimsiz hale geliyor.
mci_common_raw <- as.matrix(expr60[common_all_probes, mci_idx])
mci_frozen <- normalize.quantiles.use.target(mci_common_raw, frozen_target)
dimnames(mci_frozen) <- dimnames(mci_common_raw)

mci_data_frozen <- as.data.frame(t(mci_frozen[common_probes_frozen, ]))
mci_probs_frozen <- predict(rf_model_frozen, mci_data_frozen, type = "prob")[, "AD"]

rf_votes_frozen <- predict(rf_model_frozen, type = "prob")[, "AD"]
roc_rf_internal_frozen <- roc(rf_data_frozen$group, rf_votes_frozen, levels = c("CTL", "AD"), direction = "<")
opt_cut_frozen <- coords(roc_rf_internal_frozen, "best", ret = "threshold")$threshold

mci_results_frozen <- data.frame(
  SampleID = colnames(mci_frozen), RiskScore = mci_probs_frozen,
  Age = pheno60$age[mci_idx], Sex = pheno60$sex[mci_idx]
)
mci_results_frozen$Category <- ifelse(mci_results_frozen$RiskScore >= opt_cut_frozen, "High Risk", "Low Risk")
cat("[FROZEN-REFERENCE] High risk:", sum(mci_results_frozen$Category == "High Risk"),
    " Low risk:", sum(mci_results_frozen$Category == "Low Risk"), "\n")

print(cor.test(mci_results_frozen$RiskScore, mci_results_frozen$Age, method = "pearson"))

png("Fig4_MCI_histogram_frozen.png", width = 1400, height = 1000, res = 200)
hist(mci_results_frozen$RiskScore, breaks = 15, col = "orange",
     main = "Figure 4. AD-Risk stratification in MCI patients (frozen-reference)",
     xlab = "Probability of AD (Risk Score)")
dev.off()

write.csv(mci_results_frozen, "Table_3_MCI_Risk_frozen.csv", row.names = FALSE)

## =================================================================
## BLOCK 20b-20c (YENI): Ek confounder ve batch effect analizi
## =================================================================
## Editor elestirisine yanit: yas/cinsiyet disinda etnisite ve teknik
## batch (hibridizasyon chip'i) de kontrol edilir; PCA ile ornek
## yapisi gorsellestirilir.

## ---- BLOCK 20b: Etnisite ve chip batch confounder testleri ----
library(ggplot2)

pheno_adctl$ethnicity <- pheno60$`ethnicity:ch1`[adctl_idx]
cat("--- Etnisite dagilimi (AD vs CTL) ---\n")
print(table(pheno_adctl$group, pheno_adctl$ethnicity))
print(fisher.test(table(pheno_adctl$group, pheno_adctl$ethnicity),
                   simulate.p.value = TRUE, B = 10000))

chip_batch_all <- sub("_.*", "", pheno60$title)
pheno_adctl$chip_batch <- chip_batch_all[adctl_idx]
cat("\n--- Chip batch ---\n")
cat("Toplam benzersiz chip sayisi (AD+CTL):", length(unique(pheno_adctl$chip_batch)), "\n")
print(fisher.test(table(pheno_adctl$group, pheno_adctl$chip_batch),
                   simulate.p.value = TRUE, B = 10000))

## ---- BLOCK 20c: Etnisiteye gore duzeltilmis DGE + PCA gorsellestirme ----
ethnicity_factor <- factor(pheno_adctl$ethnicity)
group_factor <- factor(pheno_adctl$group, levels = c("CTL", "AD"))
design_eth <- model.matrix(~ group_factor + ethnicity_factor)
mat_adctl_eth <- expr60_annot[, adctl_idx]

fit_eth <- lmFit(mat_adctl_eth, design_eth)
fit_eth <- eBayes(fit_eth)
tt_eth <- topTable(fit_eth, coef = "group_factorAD", number = Inf, adjust.method = "BH")

sig_eth <- tt_eth[tt_eth$adj.P.Val < 0.05 & abs(tt_eth$logFC) > 0.5, ]
cat("\nEtnisiteye gore duzeltme sonrasi anlamli gen sayisi:", nrow(sig_eth), "\n")
overlap_eth <- intersect(rownames(sig_eth), top100_probes)
cat("100-gen imzasindan etnisite-duzeltmesi sonrasi anlamli kalan:", length(overlap_eth), "/100\n")

# PCA -- ornekler hastalik durumuna gore mi, batch'e gore mi ayrisiyor?
mat_pca <- expr60_annot[top100_probes, adctl_idx]
pca_res <- prcomp(t(mat_pca), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2],
                      Group = pheno_adctl$group, Batch = pheno_adctl$chip_batch)
cat("\nPCA aciklanan varyans -- PC1:", round(summary(pca_res)$importance[2,1]*100,1),
    "%  PC2:", round(summary(pca_res)$importance[2,2]*100,1), "%\n")

p_pca_group <- ggplot(pca_df, aes(PC1, PC2, color = Group)) + geom_point(size = 2) +
  theme_minimal() + ggtitle("PCA - Hastalik durumuna gore renklendirme")
p_pca_batch <- ggplot(pca_df, aes(PC1, PC2, color = Batch)) + geom_point(size = 2) +
  theme_minimal() + theme(legend.position = "none") +
  ggtitle("PCA - Chip batch'e gore renklendirme (legend gizli)")

ggsave("Fig_PCA_group.png", p_pca_group, width = 6, height = 5, dpi = 200)
ggsave("Fig_PCA_batch.png", p_pca_batch, width = 6, height = 5, dpi = 200)
cat("\nFig_PCA_group.png ve Fig_PCA_batch.png kaydedildi.\n")

## ---- BLOCK 20: Mitokondriyal-ribozomal ko-regülasyon ----
mito_genes <- c("NDUFA1", "NDUFS5", "ATP5ME", "ATP5F1E", "UQCRHL", "COX17")
ribo_genes <- c("RPL22", "RPL36AL", "RPS23", "RPS25", "RPS27A")
mito_probes <- top100$Probe_ID[top100$SYMBOL %in% mito_genes]
ribo_probes <- top100$Probe_ID[top100$SYMBOL %in% ribo_genes]
mito_score <- colMeans(mat_scaled[rownames(mat_scaled) %in% mito_probes, , drop = FALSE])
ribo_score <- colMeans(mat_scaled[rownames(mat_scaled) %in% ribo_probes, , drop = FALSE])
print(cor.test(mito_score, ribo_score, method = "spearman"))

## ---- BLOCK 21: Supplementary Table 1 (gen imzası) dışa aktar ----
supp1 <- top100 %>%
  left_join(imp_df[, c("Probe_ID", "MDA", "Rank")], by = "Probe_ID") %>%
  transmute(Gene_Symbol = SYMBOL, Probe_ID, RF_Importance_MDA = MDA, MDA_Rank = Rank,
            Log2_Fold_Change = logFC, Adjusted_p_value_FDR = adj.P.Val) %>%
  arrange(MDA_Rank)
write.csv(supp1, "Supplementary_Table_1.csv", row.names = FALSE)

## ---- BLOCK 22: SUMMARY ----
cat("\n\n==================== FINAL SUMMARY ====================\n")
cat("Discovery: AD=", sum(pheno_adctl$group == "AD"), " CTL=", sum(pheno_adctl$group == "CTL"), "\n")
cat("Internal OOB AUC (discovery-only norm., Block 9/10 modeli):", round(auc(roc_rf_internal), 3), "\n")
cat("Nested CV AUC:", round(auc(roc_nested), 3), "\n")
cat("Age mismatch p-value:", format(age_test$p.value, scientific = TRUE), "\n")

cat("\n---- [ESKI] Joint normalization sonuclari ----\n")
cat("External AUC - RF:", round(auc(roc_rf_ext_j), 3), " SVM:", round(auc(roc_svm_ext_j), 3),
    " GLM:", round(auc(roc_glm_ext_j), 3), "\n")
cat("MCI high-risk:", sum(mci_results$Category == "High Risk"), "/ 80\n")
cat("Brier score:", round(brier_rf, 3), "\n")

cat("\n---- [YENI] Frozen-reference normalization sonuclari  ----\n")
cat("External AUC - RF:", round(auc(roc_rf_ext_f), 3), " SVM:", round(auc(roc_svm_ext_f), 3),
    " GLM:", round(auc(roc_glm_ext_f), 3), "\n")
cat("External RF 95% CI:", round(ci.auc(roc_rf_ext_f)[1], 3), "-", round(ci.auc(roc_rf_ext_f)[3], 3), "\n")
cat("MCI high-risk:", sum(mci_results_frozen$Category == "High Risk"), "/ 80\n")
cat("Brier score:", round(brier_rf_f, 3), "\n")
cat("Imza icindeki ortak probe sayisi:", length(common_probes_frozen), "/100\n")

## ---- BLOCK 23: Çalışma alanını (workspace) kaydet ----
save.image("AD_analysis_workspace.RData")
cat("\nWorkspace kaydedildi:", file.path(getwd(), "AD_analysis_workspace.RData"), "\n")
cat("Tekrar açmak için: load('", file.path(getwd(), "AD_analysis_workspace.RData"), "')\n", sep="")
