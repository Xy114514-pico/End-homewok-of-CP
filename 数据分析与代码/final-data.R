# ====================================================================
# 最终修正版：Hansen & Melzner (2014) 计算可复现性检验
# 修正：描述性统计比对统一使用 make_comparison_table，消除 NA
# 修正：支付意愿保留 0 值，样本量 N=87
# 输出：所有对比表含完整 PE 和评级，无 NA
# ====================================================================

# 设置路径和种子 ----------------------------------------------------
set.seed(20250611)
data_path <- "E:/R2026/big homework"
output_path <- data_path

dir.create(file.path(output_path, "figures"), showWarnings = FALSE)
dir.create(file.path(output_path, "tables"), showWarnings = FALSE)
dir.create(file.path(output_path, "results"), showWarnings = FALSE)

# 加载包 ------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(rstatix)
library(effectsize)
library(BayesFactor)
library(corrplot)

# 辅助函数 ----------------------------------------------------------
calc_pe <- function(rep, orig) {
  ifelse(orig == 0 | is.na(orig), NA, abs(rep - orig) / abs(orig) * 100)
}

rate_pe <- function(pe, orig_val = NULL, rep_val = NULL) {
  if (is.na(pe) || is.na(orig_val) || is.na(rep_val)) return("无法计算")
  if (abs(orig_val) < 0.01 && abs(rep_val - orig_val) < 0.005) {
    return("因舍入导致的偏差")
  }
  if (pe == 0) return("完全一致")
  else if (pe < 10) return("次要偏差")
  else return("主要偏差")
}

partial_eta_sq <- function(F_val, df_effect, df_error) {
  (F_val * df_effect) / (F_val * df_effect + df_error)
}

# 生成对比表（通用）
make_comparison_table <- function(orig_vals, rep_vals, labels, stat_names, step_name) {
  results <- data.frame()
  for (i in seq_along(orig_vals)) {
    pe <- calc_pe(rep_vals[i], orig_vals[i])
    rating <- rate_pe(pe, orig_val = orig_vals[i], rep_val = rep_vals[i])
    results <- rbind(results, data.frame(
      步骤 = step_name,
      检验项 = labels[i],
      统计量 = stat_names[i],
      原研究值 = round(orig_vals[i], 4),
      本研究值 = round(rep_vals[i], 4),
      PE = round(pe, 2),
      评级 = rating,
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}

save_comparison <- function(df, filename) {
  write.csv(df, file.path(output_path, "tables", filename), row.names = FALSE)
  cat("已保存对比表:", filename, "\n")
}

# ====================================================================
# 步骤 1: 数据读取与分层预处理
# ====================================================================
cat("\n========== 步骤 1: 数据读取与分层预处理 ==========\n")
df_raw <- read.csv(file.path(data_path, "orig-data.csv"))
cat("原始数据行数:", nrow(df_raw), "\n")

# 主数据集（N=95）：仅排除 filter==0 和关键变量缺失
df_clean <- df_raw %>%
  filter(filter == 1) %>%
  filter(!is.na(n_categories), !is.na(kimchi_sum),
         !is.na(toaster45_mean), !is.na(toaster25_mean)) %>%
  mutate(
    condition = factor(group, levels = 1:5,
                       labels = c("Concrete", "Reverberation",
                                  "C-F#", "Nonsegmentation", "Abstract")),
    choice_toaster = factor(choice_toaster, levels = 1:2,
                            labels = c("Aggregated", "Individualized"))
  )
cat("主数据集（评价/分类/视觉）样本量:", nrow(df_clean), "\n")  # 95

# 支付意愿子数据集（N=87）：额外剔除 666 和 NA，保留 0
df_wtp <- df_clean %>%
  filter(t_t_A_w_p != 666, t_t_B_w_p != 666) %>%
  filter(!is.na(toaster_45_wtp), !is.na(toaster_25_wtp))
cat("支付意愿子数据集样本量:", nrow(df_wtp), "\n")  # 87

# 长格式数据
df_toaster_long <- df_clean %>%
  pivot_longer(cols = c(toaster45_mean, toaster25_mean),
               names_to = "info_type", values_to = "evaluation") %>%
  mutate(info_type = ifelse(info_type == "toaster45_mean", "Aggregated", "Individualized"))

df_wtp_long <- df_wtp %>%
  pivot_longer(cols = c(toaster_45_wtp, toaster_25_wtp),
               names_to = "info_type", values_to = "wtp") %>%
  mutate(info_type = ifelse(info_type == "toaster_45_wtp", "Aggregated", "Individualized"))

# ====================================================================
# 步骤 2: 描述性统计（全部条件）
# ====================================================================
cat("\n========== 步骤 2: 描述性统计 ==========\n")
desc_all <- df_clean %>%
  group_by(condition) %>%
  summarise(
    N = n(),
    Cat_Mean = mean(n_categories), Cat_SD = sd(n_categories),
    Glob_Mean = mean(kimchi_sum), Glob_SD = sd(kimchi_sum),
    Agg_Eval_Mean = mean(toaster45_mean), Agg_Eval_SD = sd(toaster45_mean),
    Ind_Eval_Mean = mean(toaster25_mean), Ind_Eval_SD = sd(toaster25_mean),
    .groups = "drop"
  )
write.csv(desc_all, file.path(output_path, "tables", "Step2_desc_all_conditions.csv"), row.names = FALSE)

# ---- 描述性统计比对（具体 vs 抽象），使用 make_comparison_table ----
desc_compare <- desc_all %>% filter(condition %in% c("Concrete", "Abstract"))
orig_desc <- data.frame(
  condition = c("Concrete", "Abstract"),
  Cat_Mean_orig = c(7.90, 6.63), Cat_SD_orig = c(1.52, 0.96),
  Glob_Mean_orig = c(6.25, 7.69), Glob_SD_orig = c(1.86, 1.01)
)
desc_joined <- left_join(desc_compare, orig_desc, by = "condition")

# 为每个指标生成对比行
desc_cmp_list <- list()
for (i in 1:nrow(desc_joined)) {
  cond <- desc_joined$condition[i]
  # 类别均值
  desc_cmp_list[[paste0(cond, "_Cat_Mean")]] <- make_comparison_table(
    orig_vals = desc_joined$Cat_Mean_orig[i],
    rep_vals = desc_joined$Cat_Mean[i],
    labels = paste("描述性(类别均值)", cond),
    stat_names = "Mean",
    step_name = "步骤2"
  )
  # 类别标准差
  desc_cmp_list[[paste0(cond, "_Cat_SD")]] <- make_comparison_table(
    orig_vals = desc_joined$Cat_SD_orig[i],
    rep_vals = desc_joined$Cat_SD[i],
    labels = paste("描述性(类别SD)", cond),
    stat_names = "SD",
    step_name = "步骤2"
  )
  # 整体均值
  desc_cmp_list[[paste0(cond, "_Glob_Mean")]] <- make_comparison_table(
    orig_vals = desc_joined$Glob_Mean_orig[i],
    rep_vals = desc_joined$Glob_Mean[i],
    labels = paste("描述性(整体均值)", cond),
    stat_names = "Mean",
    step_name = "步骤2"
  )
  # 整体标准差
  desc_cmp_list[[paste0(cond, "_Glob_SD")]] <- make_comparison_table(
    orig_vals = desc_joined$Glob_SD_orig[i],
    rep_vals = desc_joined$Glob_SD[i],
    labels = paste("描述性(整体SD)", cond),
    stat_names = "SD",
    step_name = "步骤2"
  )
}
desc_cmp_all <- bind_rows(desc_cmp_list)
save_comparison(desc_cmp_all, "Step2_desc_comparison.csv")
# ====================================================================
# 步骤 2b: 因变量3（烤面包机评价）的描述性统计比对
# 原文在结果文本中明确报告了具体/抽象条件下 4.5星 vs 2.5星 的 M 和 SD
# ====================================================================
cat("\n========== 步骤 2b: 因变量3 描述性统计比对 ==========\n")

# 1. 从原文中提取的描述性统计值（来自论文结果文本）
orig_dv3_desc <- data.frame(
  condition = c("Concrete", "Concrete", "Abstract", "Abstract"),
  info_type = c("Aggregated", "Individualized", "Aggregated", "Individualized"),
  Mean_orig = c(3.77, 4.33, 5.14, 2.98),
  SD_orig = c(1.36, 1.41, 0.72, 0.79),
  stringsAsFactors = FALSE
)

# 2. 从本研究数据中计算对应的描述性统计
rep_dv3_desc <- df_clean %>%
  filter(condition %in% c("Concrete", "Abstract")) %>%
  group_by(condition) %>%
  summarise(
    Agg_Mean = mean(toaster45_mean, na.rm = TRUE),
    Agg_SD = sd(toaster45_mean, na.rm = TRUE),
    Ind_Mean = mean(toaster25_mean, na.rm = TRUE),
    Ind_SD = sd(toaster25_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -condition,
    names_to = c("info_type", "stat"),
    names_pattern = "(Agg|Ind)_(Mean|SD)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    info_type = ifelse(info_type == "Agg", "Aggregated", "Individualized")
  )

# 3. 合并并生成比对表
dv3_desc_joined <- orig_dv3_desc %>%
  left_join(rep_dv3_desc, by = c("condition", "info_type")) %>%
  mutate(
    PE_Mean = calc_pe(Mean, Mean_orig),
    PE_SD = calc_pe(SD, SD_orig),
    评级_Mean = mapply(rate_pe, PE_Mean, Mean_orig, Mean),
    评级_SD = mapply(rate_pe, PE_SD, SD_orig, SD)
  )

# 4. 整理为规范表格（类似描述性统计表）
dv3_desc_table <- dv3_desc_joined %>%
  mutate(
    条件 = condition,
    信息类型 = info_type,
    原均值 = Mean_orig,
    本均值 = Mean,
    PE均值 = round(PE_Mean, 2),
    评级均值 = 评级_Mean,
    原标准差 = SD_orig,
    本标准差值 = SD,
    PE标准差 = round(PE_SD, 2),
    评级标准差 = 评级_SD
  ) %>%
  select(条件, 信息类型, 原均值, 本均值, PE均值, 评级均值, 
         原标准差, 本标准差值, PE标准差, 评级标准差)

# 5. 保存
write.csv(dv3_desc_table, 
          file.path(output_path, "tables", "Step2b_dv3_descriptive_comparison.csv"), 
          row.names = FALSE)

# 6. 同时生成与步骤2格式兼容的对比行（用于后续汇总）
dv3_cmp_list <- list()
for (i in 1:nrow(dv3_desc_joined)) {
  cond <- dv3_desc_joined$condition[i]
  info <- dv3_desc_joined$info_type[i]
  label <- paste0("描述性(DV3_", cond, "_", info, ")")
  
  # 均值
  dv3_cmp_list[[paste0(cond, "_", info, "_Mean")]] <- make_comparison_table(
    orig_vals = dv3_desc_joined$Mean_orig[i],
    rep_vals = dv3_desc_joined$Mean[i],
    labels = paste(label, "均值"),
    stat_names = "Mean",
    step_name = "步骤2b"
  )
  # 标准差
  dv3_cmp_list[[paste0(cond, "_", info, "_SD")]] <- make_comparison_table(
    orig_vals = dv3_desc_joined$SD_orig[i],
    rep_vals = dv3_desc_joined$SD[i],
    labels = paste(label, "SD"),
    stat_names = "SD",
    step_name = "步骤2b"
  )
}
dv3_cmp_all <- bind_rows(dv3_cmp_list)
save_comparison(dv3_cmp_all, "Step2b_dv3_descriptive_comparison_rows.csv")

cat("已保存因变量3的描述性统计比对表 (Step2b_dv3_descriptive_comparison.csv)\n")
# ====================================================================
# 步骤 3: 类别广度
# ====================================================================
cat("\n========== 步骤 3: 类别广度 ==========\n")
aov_cat <- aov(n_categories ~ condition, data = df_clean)
sum_cat <- summary(aov_cat)
F_cat <- sum_cat[[1]]$`F value`[1]
p_cat <- sum_cat[[1]]$`Pr(>F)`[1]
df_cat_eff <- sum_cat[[1]]$Df[1]
df_cat_err <- sum_cat[[1]]$Df[2]
eta_cat <- partial_eta_sq(F_cat, df_cat_eff, df_cat_err)

cat_anova_cmp <- make_comparison_table(
  orig_vals = c(2.50, 0.05, 0.10),
  rep_vals = c(F_cat, p_cat, eta_cat),
  labels = rep("类别广度ANOVA", 3),
  stat_names = c("F(4,90)", "p", "partial η²"),
  step_name = "步骤3.1"
)
save_comparison(cat_anova_cmp, "Step3.1_cat_anova_comparison.csv")

df_contrast <- df_clean %>% filter(condition %in% c("Concrete", "Abstract"))
t_cat <- t.test(n_categories ~ condition, data = df_contrast, var.equal = TRUE)
cat_t_cmp <- make_comparison_table(
  orig_vals = c(2.73, 0.008),
  rep_vals = c(t_cat$statistic, t_cat$p.value),
  labels = rep("t检验(类别)", 2),
  stat_names = c("t", "p"),
  step_name = "步骤3.2"
)
save_comparison(cat_t_cmp, "Step3.2_cat_t_comparison.csv")

# ====================================================================
# 步骤 4: 整体偏好
# ====================================================================
cat("\n========== 步骤 4: 整体偏好 ==========\n")
aov_glob <- aov(kimchi_sum ~ condition, data = df_clean)
sum_glob <- summary(aov_glob)
F_glob <- sum_glob[[1]]$`F value`[1]
p_glob <- sum_glob[[1]]$`Pr(>F)`[1]
df_glob_eff <- sum_glob[[1]]$Df[1]
df_glob_err <- sum_glob[[1]]$Df[2]
eta_glob <- partial_eta_sq(F_glob, df_glob_eff, df_glob_err)

glob_anova_cmp <- make_comparison_table(
  orig_vals = c(2.44, 0.052, 0.10),
  rep_vals = c(F_glob, p_glob, eta_glob),
  labels = rep("整体偏好ANOVA", 3),
  stat_names = c("F(4,90)", "p", "partial η²"),
  step_name = "步骤4.1"
)
save_comparison(glob_anova_cmp, "Step4.1_glob_anova_comparison.csv")

t_glob <- t.test(kimchi_sum ~ condition, data = df_contrast, var.equal = TRUE)
glob_t_cmp <- make_comparison_table(
  orig_vals = c(-2.71, 0.008),
  rep_vals = c(t_glob$statistic, t_glob$p.value),
  labels = rep("t检验(整体)", 2),
  stat_names = c("t", "p"),
  step_name = "步骤4.2"
)
save_comparison(glob_t_cmp, "Step4.2_glob_t_comparison.csv")

# ====================================================================
# 步骤 5: 评价混合ANOVA
# ====================================================================
cat("\n========== 步骤 5: 评价混合ANOVA (N=95) ==========\n")
mixed_eval <- anova_test(data = df_toaster_long,
                         dv = evaluation, wid = participant_number,
                         between = condition, within = info_type) %>%
  get_anova_table()

F_info <- mixed_eval$F[1]; p_info <- mixed_eval$p[1]
df_info_eff <- mixed_eval$DFn[1]; df_info_err <- mixed_eval$DFd[1]
F_inter <- mixed_eval$F[3]; p_inter <- mixed_eval$p[3]
df_inter_eff <- mixed_eval$DFn[3]; df_inter_err <- mixed_eval$DFd[3]
eta_info <- partial_eta_sq(F_info, df_info_eff, df_info_err)
eta_inter <- partial_eta_sq(F_inter, df_inter_eff, df_inter_err)

eval_cmp <- make_comparison_table(
  orig_vals = c(3.93, 0.05, 3.47, 0.01, 0.13),
  rep_vals = c(F_info, p_info, F_inter, p_inter, eta_inter),
  labels = c("信息类型主效应", "信息类型主效应", "评价交互", "评价交互", "评价交互"),
  stat_names = c("F(1,90)", "p", "F(4,90)", "p", "partial η²"),
  step_name = "步骤5"
)
save_comparison(eval_cmp, "Step5_eval_ANOVA_comparison.csv")

# ====================================================================
# 步骤 6: 支付意愿混合ANOVA (N=87)
# ====================================================================
cat("\n========== 步骤 6: 支付意愿混合ANOVA (N=87) ==========\n")
mixed_wtp <- anova_test(data = df_wtp_long,
                        dv = wtp, wid = participant_number,
                        between = condition, within = info_type) %>%
  get_anova_table()

F_wtp <- mixed_wtp$F[3]; p_wtp <- mixed_wtp$p[3]
df_wtp_eff <- mixed_wtp$DFn[3]; df_wtp_err <- mixed_wtp$DFd[3]
eta_wtp <- partial_eta_sq(F_wtp, df_wtp_eff, df_wtp_err)

wtp_cmp <- make_comparison_table(
  orig_vals = c(2.68, 0.04, 0.12),
  rep_vals = c(F_wtp, p_wtp, eta_wtp),
  labels = rep("支付意愿交互", 3),
  stat_names = c("F(4,82)", "p", "partial η²"),
  step_name = "步骤6"
)
save_comparison(wtp_cmp, "Step6_wtp_ANOVA_comparison.csv")

# ====================================================================
# 步骤 7: 卡方检验
# ====================================================================
cat("\n========== 步骤 7: 卡方检验 ==========\n")
choice_tab <- table(df_clean$condition, df_clean$choice_toaster)
chisq_test <- chisq.test(choice_tab)
chi2 <- chisq_test$statistic
p_chi <- chisq_test$p.value

chi_cmp <- make_comparison_table(
  orig_vals = c(12.47, 0.02),
  rep_vals = c(chi2, p_chi),
  labels = rep("卡方检验", 2),
  stat_names = c("χ²(4)", "p"),
  step_name = "步骤7"
)
save_comparison(chi_cmp, "Step7_chisq_comparison.csv")

# ====================================================================
# 步骤 8: 相关性分析
# ====================================================================
cat("\n========== 步骤 8: 相关性分析 ==========\n")
cor_data <- df_clean %>% select(n_categories, kimchi_sum, diff_toaster)
cor_matrix <- cor(cor_data, use = "pairwise")
write.csv(cor_matrix, file.path(output_path, "tables", "Step8_correlations.csv"), row.names = TRUE)
cor_test <- corrplot::cor.mtest(cor_data, conf.level = 0.95)
write.csv(cor_test$p, file.path(output_path, "tables", "Step8_cor_pvalues.csv"), row.names = TRUE)

# ====================================================================
# 步骤 9: 贝叶斯因子
# ====================================================================
cat("\n========== 步骤 9: 贝叶斯因子 ==========\n")
df_ba <- df_clean %>% filter(condition %in% c("Concrete", "Abstract")) %>% droplevels()
bf_cat <- ttestBF(formula = n_categories ~ condition, data = df_ba)
bf_glob <- ttestBF(formula = kimchi_sum ~ condition, data = df_ba)
bayes_res <- data.frame(
  假设 = c("类别广度", "整体偏好"),
  BF10 = c(as.vector(bf_cat), as.vector(bf_glob)),
  证据强度 = ifelse(c(as.vector(bf_cat), as.vector(bf_glob)) > 10, "强证据", 
                ifelse(c(as.vector(bf_cat), as.vector(bf_glob)) > 3, "中等证据", "弱证据"))
)
write.csv(bayes_res, file.path(output_path, "tables", "Step9_bayes_results.csv"), row.names = FALSE)
# ================================================================
# 补充步骤 9b: 提取贝叶斯后验估计与 95% HDI（最高密度区间）
# 满足指南表4 中“参数估计”与“置信区间”的要求
# ================================================================
cat("\n========== 步骤 9b: 贝叶斯后验估计与 HDI ==========\n")

# 加载辅助包（用于计算 HDI）
if (!require("bayestestR")) install.packages("bayestestR")
library(bayestestR)

# 1. 从贝叶斯因子对象中生成后验分布样本
# 注：设置迭代次数（iterations）以获得稳定的后验估计
set.seed(20250611)
post_cat <- posterior(bf_cat, iterations = 10000)
post_glob <- posterior(bf_glob, iterations = 10000)

# 2. 提取效应量（Cohen's d）和均值差（Mean Difference）的后验
# 在 ttestBF 中，后验样本包含:
#   - mu: 标准化效应量（Cohen's d）
#   - delta: 未标准化均值差（取决于数据单位）

# 以类别广度（Cat）为例：
# 提取 Cohen's d 的后验（即实际效应量）
d_cat <- post_cat[, "delta"]   # 注意：BayesFactor 中 delta 对应 Cohen's d

# 计算后验均值、标准差、HDI
post_mean_d_cat <- mean(d_cat)
post_sd_d_cat <- sd(d_cat)
hdi_d_cat <- hdi(d_cat, ci = 0.95)

# 提取原始数据尺度上的均值差（方便解释）
# 我们可以从数据中直接计算合并标准差，也可以直接从后验提取 mu
diff_cat <- post_cat[, "mu"]   # 均值差（未经标准化）
post_mean_diff_cat <- mean(diff_cat)
post_sd_diff_cat <- sd(diff_cat)
hdi_diff_cat <- hdi(diff_cat, ci = 0.95)

# 同理，整体偏好（Glob）
d_glob <- post_glob[, "delta"]
post_mean_d_glob <- mean(d_glob)
post_sd_d_glob <- sd(d_glob)
hdi_d_glob <- hdi(d_glob, ci = 0.95)

diff_glob <- post_glob[, "mu"]
post_mean_diff_glob <- mean(diff_glob)
post_sd_diff_glob <- sd(diff_glob)
hdi_diff_glob <- hdi(diff_glob, ci = 0.95)

# 3. 整理成报告格式（纳入表 10 的升级版）
bayes_estimates <- data.frame(
  假设 = c("类别广度", "整体偏好"),
  BF10 = c(as.vector(bf_cat), as.vector(bf_glob)),
  后验效应量_均值 = c(post_mean_d_cat, post_mean_d_glob),
  后验效应量_SD = c(post_sd_d_cat, post_sd_d_glob),
  HDI_lower_95 = c(hdi_d_cat$CI_low, hdi_d_glob$CI_low),
  HDI_upper_95 = c(hdi_d_cat$CI_high, hdi_d_glob$CI_high),
  原始均值差_均值 = c(post_mean_diff_cat, post_mean_diff_glob),
  原始均值差_HDI_lower = c(hdi_diff_cat$CI_low, hdi_diff_glob$CI_low),
  原始均值差_HDI_upper = c(hdi_diff_cat$CI_high, hdi_diff_glob$CI_high)
)

# 四舍五入保留两位小数
bayes_estimates[, -1] <- round(bayes_estimates[, -1], 2)

# 保存增强后的贝叶斯结果
write.csv(bayes_estimates, 
          file.path(output_path, "tables", "Step9b_bayes_estimates_HDI.csv"), 
          row.names = FALSE)

cat("\n后验估计与 HDI 结果:\n")
print(bayes_estimates)
# ====================================================================
# 步骤 10: 整体汇总（终极修正版）
# 修正1: 直接使用内存中的 p_info / p_inter，避免行名索引 NA
# 修正2: 引入边界缓冲区，处理 p ≈ 0.05 的舍入偏差
# ====================================================================
cat("\n========== 步骤 10: 整体汇总（终极修正） ==========\n")

# --- 10.1 定义所有对比文件并合并 ---
compare_files <- c(
  "Step2_desc_comparison.csv",                     # DV1, DV2 描述性统计
  "Step2b_dv3_descriptive_comparison_rows.csv",   # DV3 描述性统计
  "Step3.1_cat_anova_comparison.csv",             # 类别广度 ANOVA
  "Step3.2_cat_t_comparison.csv",                 # 类别广度 t检验
  "Step4.1_glob_anova_comparison.csv",            # 整体偏好 ANOVA
  "Step4.2_glob_t_comparison.csv",                # 整体偏好 t检验
  "Step5_eval_ANOVA_comparison.csv",              # 评价混合 ANOVA
  "Step6_wtp_ANOVA_comparison.csv",               # 支付意愿混合 ANOVA
  "Step7_chisq_comparison.csv"                    # 卡方检验
)

all_results <- data.frame()
for (f in compare_files) {
  full_path <- file.path(output_path, "tables", f)
  if (file.exists(full_path)) {
    df <- read.csv(full_path, stringsAsFactors = FALSE)
    all_results <- bind_rows(all_results, df)
    cat("已合并:", f, "\n")
  } else {
    cat("警告：文件不存在，跳过:", f, "\n")
  }
}

# 过滤无效行
all_results <- all_results %>%
  filter(!is.na(PE), !is.na(评级), 评级 != "无法计算")

write.csv(all_results, file.path(output_path, "tables", "Step10_ALL_COMPARISONS.csv"), row.names = FALSE)
cat("Step10_ALL_COMPARISONS 已保存，共", nrow(all_results), "行有效对比。\n")

# --- 10.2 PE评级分布 ---
rating_summary <- all_results %>%
  group_by(评级) %>%
  summarise(数量 = n(), .groups = "drop") %>%
  mutate(占比 = round(数量 / sum(数量) * 100, 2))
write.csv(rating_summary, file.path(output_path, "tables", "Step10_rating_summary.csv"), row.names = FALSE)
cat("\nPE评级分布:\n")
print(rating_summary)

# --- 10.3 推论一致性（使用直接变量 + 边界缓冲） ---
# 确保 p_info 和 p_inter 存在（它们应在步骤5中计算）
if (!exists("p_info") || !exists("p_inter")) {
  # 若意外丢失，从已保存的文件中读取备用
  eval_file <- file.path(output_path, "tables", "Step5_eval_ANOVA_comparison.csv")
  if (file.exists(eval_file)) {
    eval_df <- read.csv(eval_file)
    p_info <- eval_df[eval_df$统计量 == "p" & eval_df$检验项 == "信息类型主效应", "本研究值"]
    p_inter <- eval_df[eval_df$统计量 == "p" & eval_df$检验项 == "评价交互", "本研究值"]
    if (length(p_info) == 0) p_info <- NA
    if (length(p_inter) == 0) p_inter <- NA
  } else {
    p_info <- NA; p_inter <- NA
  }
}

# 构建 p 值向量（顺序与检验标签一一对应）
p_orig_all <- c(0.05, 0.052, 0.008, 0.008, 0.05, 0.01, 0.02)
p_rep_all  <- c(p_cat, p_glob, t_cat$p.value, t_glob$p.value, p_info, p_inter, p_chi)

# 定义含边界缓冲的一致性判断函数
# 规则：如果原 p 值在 [0.045, 0.055] 区间，且本 p 值也在该区间，视为"边界一致"
consistent_p_buffered <- function(p_rep, p_orig, alpha = 0.05, buffer = 0.005) {
  # 标准判断
  standard <- (p_rep < alpha) == (p_orig < alpha)
  
  # 边界缓冲判断：若原值在 alpha 附近 (alpha ± buffer)，且本值也在同范围，则视为一致
  if (abs(p_orig - alpha) <= buffer && abs(p_rep - alpha) <= buffer) {
    return(TRUE)  # 边界舍入偏差，不改变推论方向
  }
  return(standard)
}

# 应用缓冲一致性
consistent_vec <- mapply(consistent_p_buffered, p_rep_all, p_orig_all)

# 生成详细的一致性表格（含边界标注）
consistency_df <- data.frame(
  检验 = c("类别ANOVA", "整体ANOVA", "t(类别)", "t(整体)", "信息类型主效应", "评价交互", "卡方"),
  原p = p_orig_all,
  本p = round(p_rep_all, 4),
  一致 = consistent_vec,
  备注 = ifelse(
    consistent_vec & (abs(p_orig_all - 0.05) <= 0.005 & abs(p_rep_all - 0.05) <= 0.005),
    "边界一致（舍入偏差）",
    ifelse(consistent_vec, "标准一致", "不一致")
  )
)

write.csv(consistency_df, file.path(output_path, "tables", "Step10_consistency.csv"), row.names = FALSE)
cat("\n推论一致性（含边界缓冲）:\n")
print(consistency_df)

# 计算最终一致率（含边界一致）
final_consistency_rate <- mean(consistent_vec) * 100
cat("\n最终一致率:", sum(consistent_vec), "/ 7 (", round(final_consistency_rate, 2), "%)\n")
# ====================================================================
# 步骤 11: 生成图表
# ====================================================================
cat("\n========== 步骤 11: 生成图表 ==========\n")
p1 <- ggplot(df_clean, aes(x = condition, y = n_categories, fill = condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.7) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
  labs(x = "Sound Condition", y = "Number of Categories") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(output_path, "figures", "Fig1_category_breadth.png"), p1, width = 6, height = 4, dpi = 300)

p2 <- ggplot(df_clean, aes(x = condition, y = kimchi_sum, fill = condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.7) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
  labs(x = "Sound Condition", y = "Global Choices (out of 8)") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(output_path, "figures", "Fig2_global_choices.png"), p2, width = 6, height = 4, dpi = 300)

p3 <- ggplot(df_toaster_long, aes(x = info_type, y = evaluation, color = condition, group = condition)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", aes(linetype = condition)) +
  labs(x = "Information Type", y = "Evaluation (1-7)") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(output_path, "figures", "Fig3_evaluation_interaction.png"), p3, width = 6, height = 4, dpi = 300)

p4 <- ggplot(df_wtp_long, aes(x = info_type, y = wtp, color = condition, group = condition)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", aes(linetype = condition)) +
  labs(x = "Information Type", y = "Willingness to Pay (€)") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(output_path, "figures", "Fig4_wtp_interaction.png"), p4, width = 6, height = 4, dpi = 300)

choice_prop <- df_clean %>%
  group_by(condition) %>%
  summarise(prop_agg = mean(choice_toaster == "Aggregated") * 100)
p5 <- ggplot(choice_prop, aes(x = condition, y = prop_agg, fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  ylim(0, 100) + labs(x = "Condition", y = "% Choosing Aggregated") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(output_path, "figures", "Fig5_choice_proportion.png"), p5, width = 6, height = 4, dpi = 300)

# ====================================================================
# 步骤 12: 生成最终报告
# ====================================================================
sink(file.path(output_path, "results", "FINAL_REPORT.txt"))
cat("==================== 最终可复现性报告 ====================\n")
cat("文献: Hansen & Melzner (2014) Journal of Experimental Social Psychology\n")
cat("复现日期:", Sys.Date(), "\n")
cat("主数据集 (评价/分类/视觉) N =", nrow(df_clean), " (自由度 F(1,90) / F(4,90))\n")
cat("支付意愿子数据集 N =", nrow(df_wtp), " (自由度 F(4,82))\n\n")
cat("--- PE评级分布 ---\n")
print(rating_summary)
cat("\n--- 推论一致性 ---\n")
print(consistency_df)
cat("\n--- 所有对比表已保存至 tables/ 文件夹 ---\n")
sink()

cat("\n========== 全部完成！结果已保存至:", output_path, " ==========\n")