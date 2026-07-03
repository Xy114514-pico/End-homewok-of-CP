# ============================================================
# 计算可复现性检验脚本 (最终完整修正版)
# 文献: Hansen & Melzner (2014) JESP
# 符合《心理学院R编程语言课程可重复检验指南(2025版)》
# ============================================================

# 设置路径和种子 ------------------------------------------------
set.seed(20250611)
data_path <- "E:/R2026/big homework"
output_path <- data_path

dir.create(file.path(output_path, "figures"), showWarnings = FALSE)
dir.create(file.path(output_path, "tables"), showWarnings = FALSE)
dir.create(file.path(output_path, "results"), showWarnings = FALSE)

# 加载包 --------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(rstatix)
library(car)
library(psych)
library(BayesFactor)
library(ggpubr)

# 向量化的百分误差与评级函数 ------------------------------------
calc_pe <- function(rep, orig) {
  ifelse(orig == 0 | is.na(orig), NA, abs(rep - orig) / abs(orig) * 100)
}

rate_pe <- function(pe) {
  case_when(
    is.na(pe) ~ "无法计算",
    pe == 0 ~ "完全一致",
    pe < 10 ~ "次要偏差",
    TRUE ~ "主要偏差"
  )
}

# 推论一致性（基于p值是否同侧）
consistent_p <- function(p_rep, p_orig, alpha = 0.05) {
  (p_rep < alpha) == (p_orig < alpha)
}

# 1. 数据读取与预处理 --------------------------------------------
cat("===== 1. 数据读取与预处理 =====\n")
df_raw <- read.csv(file.path(data_path, "orig-data.csv"))
cat("原始数据行数:", nrow(df_raw), "\n")
cat("原始数据缺失值统计:\n"); print(sapply(df_raw, function(x) sum(is.na(x))))

df_clean <- df_raw %>%
  filter(filter == 1) %>%
  filter(!is.na(n_categories), !is.na(kimchi_sum),
         !is.na(toaster45_mean), !is.na(toaster25_mean)) %>%
  filter(t_t_A_w_p != 666, t_t_B_w_p != 666)

cat("清洗后有效样本量:", nrow(df_clean), "\n")  # 应为95

df_clean <- df_clean %>%
  mutate(
    condition = factor(group, levels = 1:5,
                       labels = c("Concrete", "Reverberation",
                                  "C-F#", "Nonsegmentation", "Abstract")),
    choice_toaster = factor(choice_toaster, levels = 1:2,
                            labels = c("Aggregated", "Individualized"))
  )

if(!"diff_toaster" %in% names(df_clean)) {
  df_clean$diff_toaster <- df_clean$toaster45_mean - df_clean$toaster25_mean
}

df_toaster_long <- df_clean %>%
  pivot_longer(cols = c(toaster45_mean, toaster25_mean),
               names_to = "info_type", values_to = "evaluation") %>%
  mutate(info_type = ifelse(info_type == "toaster45_mean", "Aggregated", "Individualized"))

df_wtp_long <- df_clean %>%
  pivot_longer(cols = c(toaster_45_wtp, toaster_25_wtp),
               names_to = "info_type", values_to = "wtp") %>%
  mutate(info_type = ifelse(info_type == "toaster_45_wtp", "Aggregated", "Individualized"))

# 2. 描述性统计（表2）-------------------------------------------
desc_replicated <- df_clean %>%
  group_by(condition) %>%
  summarise(
    N = n(),
    Cat_Mean = mean(n_categories, na.rm = TRUE),
    Cat_SD = sd(n_categories, na.rm = TRUE),
    Glob_Mean = mean(kimchi_sum, na.rm = TRUE),
    Glob_SD = sd(kimchi_sum, na.rm = TRUE),
    .groups = "drop"
  )

original_desc <- data.frame(
  condition = c("Concrete", "Abstract"),
  Cat_Mean_orig = c(7.90, 6.63),
  Cat_SD_orig = c(1.52, 0.96),
  Glob_Mean_orig = c(6.25, 7.69),
  Glob_SD_orig = c(1.86, 1.01)
)

desc_compare <- desc_replicated %>%
  filter(condition %in% c("Concrete", "Abstract")) %>%
  left_join(original_desc, by = "condition") %>%
  mutate(
    PE_Cat_Mean = calc_pe(Cat_Mean, Cat_Mean_orig),
    PE_Cat_SD   = calc_pe(Cat_SD, Cat_SD_orig),
    PE_Glob_Mean = calc_pe(Glob_Mean, Glob_Mean_orig),
    PE_Glob_SD   = calc_pe(Glob_SD, Glob_SD_orig)
  ) %>%
  mutate(across(starts_with("PE"), rate_pe, .names = "Rating_{.col}"))

write.csv(desc_compare, file.path(output_path, "tables", "Table2_desc_comparison.csv"), row.names = FALSE)

# 3. 推断性统计（原方法）----------------------------------------
all_results <- list()

# 辅助函数：手动计算 eta²
calc_eta2 <- function(aov_model) {
  ss <- summary(aov_model)[[1]]$`Sum Sq`
  ss_effect <- ss[1]
  ss_resid <- ss[2]
  ss_effect / (ss_effect + ss_resid)
}

# 3.1 类别广度 ANOVA
aov_cat <- aov(n_categories ~ condition, data = df_clean)
sum_cat <- summary(aov_cat)
F_cat <- sum_cat[[1]]$`F value`[1]
p_cat <- sum_cat[[1]]$`Pr(>F)`[1]
eta2_cat <- calc_eta2(aov_cat)
orig_F_cat <- 2.50; orig_p_cat <- 0.05; orig_eta2_cat <- 0.10
all_results$cat_anova <- data.frame(
  statistic = c("F", "p", "eta2"),
  original = c(orig_F_cat, orig_p_cat, orig_eta2_cat),
  replicated = c(F_cat, p_cat, eta2_cat),
  PE = c(calc_pe(F_cat, orig_F_cat), calc_pe(p_cat, orig_p_cat), calc_pe(eta2_cat, orig_eta2_cat))
)

# 3.2 整体偏好 ANOVA
aov_glob <- aov(kimchi_sum ~ condition, data = df_clean)
sum_glob <- summary(aov_glob)
F_glob <- sum_glob[[1]]$`F value`[1]
p_glob <- sum_glob[[1]]$`Pr(>F)`[1]
eta2_glob <- calc_eta2(aov_glob)
orig_F_glob <- 2.44; orig_p_glob <- 0.052; orig_eta2_glob <- 0.10
all_results$glob_anova <- data.frame(
  statistic = c("F", "p", "eta2"),
  original = c(orig_F_glob, orig_p_glob, orig_eta2_glob),
  replicated = c(F_glob, p_glob, eta2_glob),
  PE = c(calc_pe(F_glob, orig_F_glob), calc_pe(p_glob, orig_p_glob), calc_pe(eta2_glob, orig_eta2_glob))
)

# 3.3 计划对比 t检验
concrete_cat <- df_clean %>% filter(condition == "Concrete") %>% pull(n_categories)
abstract_cat <- df_clean %>% filter(condition == "Abstract") %>% pull(n_categories)
t_cat <- t.test(concrete_cat, abstract_cat, var.equal = TRUE)
orig_t_cat <- 2.73; orig_p_cat_t <- 0.008
concrete_glob <- df_clean %>% filter(condition == "Concrete") %>% pull(kimchi_sum)
abstract_glob <- df_clean %>% filter(condition == "Abstract") %>% pull(kimchi_sum)
t_glob <- t.test(concrete_glob, abstract_glob, var.equal = TRUE)
orig_t_glob <- -2.71; orig_p_glob_t <- 0.008
all_results$contrast_cat <- data.frame(
  statistic = c("t", "p"),
  original = c(orig_t_cat, orig_p_cat_t),
  replicated = c(t_cat$statistic, t_cat$p.value),
  PE = c(calc_pe(t_cat$statistic, orig_t_cat), calc_pe(t_cat$p.value, orig_p_cat_t))
)
all_results$contrast_glob <- data.frame(
  statistic = c("t", "p"),
  original = c(orig_t_glob, orig_p_glob_t),
  replicated = c(t_glob$statistic, t_glob$p.value),
  PE = c(calc_pe(abs(t_glob$statistic), abs(orig_t_glob)), calc_pe(t_glob$p.value, orig_p_glob_t))
)

# 3.4 混合 ANOVA 交互作用
res_mixed <- anova_test(data = df_toaster_long,
                        dv = evaluation, wid = participant_number,
                        between = condition, within = info_type)
res_mixed <- get_anova_table(res_mixed)
interaction_F <- res_mixed$F[3]
interaction_p <- res_mixed$p[3]
interaction_eta <- res_mixed$ges[3]  # 广义 eta²
orig_inter_F <- 3.47; orig_inter_p <- 0.01; orig_inter_eta <- 0.13
all_results$mixed_interaction <- data.frame(
  statistic = c("F(4,90)", "p", "eta2_g"),
  original = c(orig_inter_F, orig_inter_p, orig_inter_eta),
  replicated = c(interaction_F, interaction_p, interaction_eta),
  PE = c(calc_pe(interaction_F, orig_inter_F), calc_pe(interaction_p, orig_inter_p), calc_pe(interaction_eta, orig_inter_eta))
)

# 3.5 卡方检验
choice_tab <- table(df_clean$condition, df_clean$choice_toaster)
chisq_test <- chisq.test(choice_tab)
chi2 <- chisq_test$statistic
p_chi <- chisq_test$p.value
orig_chi2 <- 12.47; orig_p_chi <- 0.02
all_results$chisq <- data.frame(
  statistic = c("chisq(4)", "p"),
  original = c(orig_chi2, orig_p_chi),
  replicated = c(chi2, p_chi),
  PE = c(calc_pe(chi2, orig_chi2), calc_pe(p_chi, orig_p_chi))
)

# 保存推断统计比较表
infer_long <- bind_rows(lapply(names(all_results), function(nm) {
  df <- all_results[[nm]]
  df$analysis <- nm
  return(df)
}), .id = NULL)
write.csv(infer_long, file.path(output_path, "tables", "Table3_inferential_comparison.csv"), row.names = FALSE)

# 4. 结果可复现性评估（PE评级分布）-------------------------------
all_pe <- infer_long$PE[!is.na(infer_long$PE)]
ratings <- rate_pe(all_pe)
rating_table <- table(ratings)
rating_prop <- prop.table(rating_table) * 100
rating_summary <- data.frame(
  可复现性情况 = names(rating_table),
  数量 = as.numeric(rating_table),
  占比 = as.numeric(rating_prop)
)
write.csv(rating_summary, file.path(output_path, "tables", "Table5_reproducibility_ratings.csv"), row.names = FALSE)

# 5. 推论一致性评估 ----------------------------------------------
p_orig <- c(0.05, 0.052, 0.008, 0.008, 0.01, 0.02)
p_rep <- c(p_cat, p_glob, t_cat$p.value, t_glob$p.value, interaction_p, p_chi)
consistent_vec <- consistent_p(p_rep, p_orig)
consistency_table <- table(consistent_vec)
consistency_prop <- prop.table(consistency_table) * 100
consistency_df <- data.frame(
  推论一致性 = ifelse(names(consistency_table) == "TRUE", "一致", "不一致"),
  数量 = as.numeric(consistency_table),
  占比 = as.numeric(consistency_prop)
)
write.csv(consistency_df, file.path(output_path, "tables", "Table6_inference_consistency.csv"), row.names = FALSE)

# 6. 创新方法：贝叶斯因子 ----------------------------------------
# 关键修正：过滤后必须删除未使用的因子水平
df_ba <- df_clean %>%
  filter(condition %in% c("Concrete", "Abstract")) %>%
  droplevels()

bf_cat <- ttestBF(formula = n_categories ~ condition, data = df_ba)
bf_glob <- ttestBF(formula = kimchi_sum ~ condition, data = df_ba)

bf10_cat <- as.vector(bf_cat)
bf10_glob <- as.vector(bf_glob)
bayes_res <- data.frame(
  假设 = c("类别广度（具体 vs 抽象）", "整体偏好（具体 vs 抽象）"),
  BF10 = c(bf10_cat, bf10_glob),
  证据强度 = ifelse(c(bf10_cat, bf10_glob) > 3, "中等证据", 
                ifelse(c(bf10_cat, bf10_glob) > 10, "强证据", "弱证据"))
)
write.csv(bayes_res, file.path(output_path, "tables", "Table4_bayes_results.csv"), row.names = FALSE)

# 7. 生成图表 ----------------------------------------------------
p1 <- ggplot(df_clean, aes(x = condition, y = n_categories, fill = condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.7) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
  labs(x = "Sound Condition", y = "Number of Categories", 
       title = "Effect of Sound on Category Breadth") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(output_path, "figures", "Fig1_category_breadth.png"), p1, width = 6, height = 4, dpi = 300)

p2 <- ggplot(df_clean, aes(x = condition, y = kimchi_sum, fill = condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.7) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
  labs(x = "Sound Condition", y = "Number of Global Choices (out of 8)", 
       title = "Effect of Sound on Global vs Local Processing") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(output_path, "figures", "Fig2_global_choices.png"), p2, width = 6, height = 4, dpi = 300)

p3 <- ggplot(df_toaster_long, aes(x = info_type, y = evaluation, 
                                  color = condition, group = condition)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", aes(linetype = condition)) +
  labs(x = "Type of Information", y = "Mean Evaluation (1-7)", 
       title = "Interaction between Sound Condition and Information Type") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(output_path, "figures", "Fig3_evaluation_interaction.png"), p3, width = 6, height = 4, dpi = 300)

p4 <- ggplot(df_wtp_long, aes(x = info_type, y = wtp, 
                              color = condition, group = condition)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", aes(linetype = condition)) +
  labs(x = "Type of Information", y = "Willingness to Pay (€)", 
       title = "Interaction on Willingness to Pay") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(output_path, "figures", "Fig4_wtp_interaction.png"), p4, width = 6, height = 4, dpi = 300)

choice_prop <- df_clean %>%
  group_by(condition) %>%
  summarise(prop_agg = mean(choice_toaster == "Aggregated", na.rm = TRUE) * 100)
p5 <- ggplot(choice_prop, aes(x = condition, y = prop_agg, fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "Sound Condition", y = "Percentage Choosing Aggregated-Favored Toaster", 
       title = "Effect of Sound on Product Choice") +
  ylim(0, 100) + theme_minimal() + theme(legend.position = "none")
ggsave(file.path(output_path, "figures", "Fig5_choice_proportion.png"), p5, width = 6, height = 4, dpi = 300)

# 8. 保存完整的统计结果文本 --------------------------------------
sink(file.path(output_path, "results", "all_statistical_outputs.txt"))
cat("=====  Hansen & Melzner (2014) 计算可复现性检验 =====\n\n")
cat("样本量: ", nrow(df_clean), "\n\n")
cat("--- 描述性统计比较 ---\n")
print(desc_compare)
cat("\n--- 类别广度 ANOVA ---\n")
print(sum_cat)
cat("\n原研究 F=", orig_F_cat, ", p=", orig_p_cat, ", eta2=", orig_eta2_cat, "\n")
cat("复现 PE: F=", calc_pe(F_cat, orig_F_cat), "%, p=", calc_pe(p_cat, orig_p_cat), "%, eta2=", calc_pe(eta2_cat, orig_eta2_cat), "%\n")
cat("\n--- 整体偏好 ANOVA ---\n")
print(sum_glob)
cat("\n原研究 F=", orig_F_glob, ", p=", orig_p_glob, ", eta2=", orig_eta2_glob, "\n")
cat("复现 PE: F=", calc_pe(F_glob, orig_F_glob), "%, p=", calc_pe(p_glob, orig_p_glob), "%, eta2=", calc_pe(eta2_glob, orig_eta2_glob), "%\n")
cat("\n--- 计划对比 t检验 ---\n")
cat("类别广度: t=", t_cat$statistic, ", p=", t_cat$p.value, "\n")
cat("整体偏好: t=", t_glob$statistic, ", p=", t_glob$p.value, "\n")
cat("\n--- 混合 ANOVA 交互作用 ---\n")
print(res_mixed)
cat("\n原研究 F=", orig_inter_F, ", p=", orig_inter_p, ", eta2_g=", orig_inter_eta, "\n")
cat("复现 PE: F=", calc_pe(interaction_F, orig_inter_F), "%, p=", calc_pe(interaction_p, orig_inter_p), "%, eta2=", calc_pe(interaction_eta, orig_inter_eta), "%\n")
cat("\n--- 卡方检验 ---\n")
print(chisq_test)
cat("\n原研究 χ²=", orig_chi2, ", p=", orig_p_chi, "\n")
cat("复现 PE: χ²=", calc_pe(chi2, orig_chi2), "%, p=", calc_pe(p_chi, orig_p_chi), "%\n")
cat("\n--- 贝叶斯因子 ---\n")
print(bayes_res)
cat("\n--- 评级汇总 ---\n")
print(rating_summary)
cat("\n--- 推论一致性 ---\n")
print(consistency_df)
sink()

cat("\n===== 分析完成 =====\n")
cat("所有结果已保存至:", output_path, "\n")