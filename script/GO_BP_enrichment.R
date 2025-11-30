
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(purrr)
})

## ========== 0) 参数 & 路径 ==========

base_dir <- "path_to_de_results"

out_dir <- file.path(
  base_dir,
  "GO_BP_all_tissues_bubble_per_TissueGender"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## DE 基因筛选阈值
padj_thr  <- 0.05
lfc_thr   <- log2(1.5)

## enrichGO 阈值（在 run_enrichGO_bp 中先用它筛一遍）
ego_padj_thr <- 0.1

## 需要遍历的 Tissue / Gender
tissues <- c("Liver", "Muscle")
genders <- c("Male", "Female")

## ========== 1) 小工具：读 DE 结果并返回 up/down 基因 ==========

load_de_and_split_ud <- function(file, padj_thr, lfc_thr) {
  stopifnot(file.exists(file))
  df <- read.table(
    file,
    header      = TRUE,
    sep         = "\t",
    quote       = "",
    check.names = FALSE
  )
  need_cols <- c("gene", "log2FoldChange", "padj")
  miss <- setdiff(need_cols, colnames(df))
  if (length(miss) > 0) {
    stop("DE 结果缺少列: ", paste(miss, collapse = ", "), " in ", file)
  }

  df <- df %>%
    filter(!is.na(gene),
           !is.na(log2FoldChange),
           !is.na(padj))

  up_genes <- df %>%
    filter(padj < padj_thr,
           log2FoldChange >= lfc_thr) %>%
    pull(gene) %>%
    unique()

  down_genes <- df %>%
    filter(padj < padj_thr,
           log2FoldChange <= -lfc_thr) %>%
    pull(gene) %>%
    unique()

  universe <- df$gene %>% unique()

  list(
    up_genes   = up_genes,
    down_genes = down_genes,
    universe   = universe
  )
}

## ========== 2) 小工具：对某个基因集做 enrichGO，并转为 long df ==========

run_enrichGO_bp <- function(genes, universe, tissue, gender, group, direction,
                            padj_cut = ego_padj_thr) {
  if (length(genes) < 10) {
    return(NULL)
  }

  ego <- tryCatch(
    enrichGO(
      gene          = genes,
      universe      = universe,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,      # 不在这里截断
      qvalueCutoff  = 1,
      readable      = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(NULL)
  }

  df <- as.data.frame(ego)

  ## 只保留 p.adjust < padj_cut 的 term
  df <- df %>%
    filter(!is.na(p.adjust),
           p.adjust < padj_cut)

  if (nrow(df) == 0) return(NULL)

  df <- df %>%
    mutate(
      Tissue    = tissue,
      Gender    = gender,
      group     = group,      # "WT" 或 "XL"
      direction = direction   # "Up" / "Down"
    )

  df
}

## ========== 3) 遍历所有 Tissue×Gender，收集 GO BP 结果 ==========

ego_all_list <- list()

for (t in tissues) {
  for (g in genders) {
    folder_name <- paste0(tolower(t), "_", tolower(g))  # e.g. liver_male
    de_dir <- file.path(base_dir, folder_name)

    ## WT Old vs Young
    f_wt <- file.path(de_dir, "B1_WT_Old_vs_Young.tsv")
    ## XL Old vs Young (Heter)
    f_xl <- file.path(de_dir, "B2_Heter_Old_vs_Young.tsv")

    if (!file.exists(f_wt) || !file.exists(f_xl)) {
      message("跳过: ", t, " ", g, "，找不到 B1/B2 文件")
      next
    }

    message("Processing: ", t, " - ", g)

    ## WT
    de_wt <- load_de_and_split_ud(f_wt, padj_thr = padj_thr, lfc_thr = lfc_thr)
    ego_wt_up   <- run_enrichGO_bp(de_wt$up_genes,   de_wt$universe,
                                   tissue = t, gender = g,
                                   group = "WT", direction = "Up")
    ego_wt_down <- run_enrichGO_bp(de_wt$down_genes, de_wt$universe,
                                   tissue = t, gender = g,
                                   group = "WT", direction = "Down")

    ## XL
    de_xl <- load_de_and_split_ud(f_xl, padj_thr = padj_thr, lfc_thr = lfc_thr)
    ego_xl_up   <- run_enrichGO_bp(de_xl$up_genes,   de_xl$universe,
                                   tissue = t, gender = g,
                                   group = "XL", direction = "Up")
    ego_xl_down <- run_enrichGO_bp(de_xl$down_genes, de_xl$universe,
                                   tissue = t, gender = g,
                                   group = "XL", direction = "Down")

    ego_all_list <- c(
      ego_all_list,
      list(ego_wt_up, ego_wt_down, ego_xl_up, ego_xl_down)
    )
  }
}

ego_all <- ego_all_list %>%
  purrr::compact() %>%
  dplyr::bind_rows()

if (nrow(ego_all) == 0) {
  stop("没有任何 GO BP term 通过当前阈值，请调高 q/p 阈值后重试。")
}

## 保存一份 long 表，方便以后复用
go_tsv <- file.path(out_dir, "GO_BP_enrich_all_tissues_long.tsv")
write.table(
  ego_all,
  file      = go_tsv,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
message("已写出 GO BP 富集结果: ", go_tsv)

## ego_all 典型列包括：
## ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count,
## Tissue, Gender, group, direction

## 添加一些通用列
ego_all <- ego_all %>%
  mutate(
    contrast = paste(Tissue, Gender, group, sep = "_"),   # e.g. Liver_Male_WT
    neg_log10_padj = -log10(p.adjust),
    Description_wrap = stringr::str_wrap(Description, width = 70)
  )

## ========== 4) “每个 Tissue×Gender×Group 各选3个通路”并画 Up / Down 图 ==========

## -------- Up：每个(Tissue,Gender,group) 选3条，上调 --------

ego_up <- ego_all %>%
  filter(direction == "Up")

if (nrow(ego_up) > 0) {

  ## 对每个对比单独选 top3（按 p.adjust）
  ego_up_top_by_contrast <- ego_up %>%
    group_by(Tissue, Gender, group, contrast, Description_wrap) %>%
    summarise(
      min_padj = min(p.adjust, na.rm = TRUE),
      .groups  = "drop_last"
    ) %>%
    arrange(min_padj) %>%
    slice_head(n = 3) %>%        # 每个 Tissue×Gender×group 各最多 3 个
    ungroup()

  ## 所有选中的通路（并集）
  up_selected_terms <- ego_up_top_by_contrast %>%
    pull(Description_wrap) %>%
    unique()

  ego_up_plot <- ego_up %>%
    filter(Description_wrap %in% up_selected_terms)

  ## 通路顺序：按“在全体中最显著”的顺序，方便阅读
  term_order_up <- ego_up_plot %>%
    group_by(Description_wrap) %>%
    summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(min_padj) %>%
    pull(Description_wrap)

  ego_up_plot$Description_wrap <- factor(
    ego_up_plot$Description_wrap,
    levels = rev(term_order_up)   # 最显著放上面
  )

  ## 对比顺序：按 Tissue, Gender, group 的组合排序
  contrast_order_up <- ego_up_plot %>%
    distinct(Tissue, Gender, group, contrast) %>%
    arrange(Tissue, Gender, group) %>%
    pull(contrast)

  ego_up_plot$contrast <- factor(
    ego_up_plot$contrast,
    levels = contrast_order_up
  )

  p_up <- ggplot(
    ego_up_plot,
    aes(x = contrast, y = Description_wrap)
  ) +
    geom_point(
      aes(size = Count, fill = neg_log10_padj),
      shape  = 21,
      colour = "black",
      alpha  = 0.8
    ) +
    scale_size_continuous(
      name = "Gene count",
      range = c(1.5, 6)
    ) +
    scale_fill_gradient(
      name = expression(-log[10]("padj")),
      low  = "white",
      high = "#E31A1C"
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y  = element_text(colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid   = element_blank(),
      legend.position = "right"
    ) +
    ggtitle("GO BP enrichment (Upregulated genes)\nEach Tissue×Gender×Group: top 3 terms")

  pdf_up <- file.path(out_dir, "GO_BP_Up_top3_per_TissueGenderGroup_all_tissues.pdf")
  pdf(pdf_up, width = 4, height = 3.6)
  print(p_up)
  dev.off()
  message("已保存上调基因 GO BP bubble plot: ", pdf_up)
} else {
  message("没有上调方向的显著通路，未绘制 Up 图。")
}

## -------- Down：每个(Tissue,Gender,group) 选3条，下调 --------

ego_down <- ego_all %>%
  filter(direction == "Down")

if (nrow(ego_down) > 0) {

  ego_down_top_by_contrast <- ego_down %>%
    group_by(Tissue, Gender, group, contrast, Description_wrap) %>%
    summarise(
      min_padj = min(p.adjust, na.rm = TRUE),
      .groups  = "drop_last"
    ) %>%
    arrange(min_padj) %>%
    slice_head(n = 3) %>%
    ungroup()

  down_selected_terms <- ego_down_top_by_contrast %>%
    pull(Description_wrap) %>%
    unique()

  ego_down_plot <- ego_down %>%
    filter(Description_wrap %in% down_selected_terms)

  term_order_down <- ego_down_plot %>%
    group_by(Description_wrap) %>%
    summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(min_padj) %>%
    pull(Description_wrap)

  ego_down_plot$Description_wrap <- factor(
    ego_down_plot$Description_wrap,
    levels = rev(term_order_down)
  )

  contrast_order_down <- ego_down_plot %>%
    distinct(Tissue, Gender, group, contrast) %>%
    arrange(Tissue, Gender, group) %>%
    pull(contrast)

  ego_down_plot$contrast <- factor(
    ego_down_plot$contrast,
    levels = contrast_order_down
  )

  p_down <- ggplot(
    ego_down_plot,
    aes(x = contrast, y = Description_wrap)
  ) +
    geom_point(
      aes(size = Count, fill = neg_log10_padj),
      shape  = 21,
      colour = "black",
      alpha  = 0.8
    ) +
    scale_size_continuous(
      name = "Gene count",
      range = c(1.5, 6)
    ) +
    scale_fill_gradient(
      name = expression(-log[10]("padj")),
      low  = "white",
      high = "#E31A1C"
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y  = element_text(colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid   = element_blank(),
      legend.position = "right"
    ) +
    ggtitle("GO BP enrichment (Downregulated genes)\nEach Tissue×Gender×Group: top 3 terms")

  pdf_down <- file.path(out_dir, "GO_BP_Down_top3_per_TissueGenderGroup_all_tissues.pdf")
  pdf(pdf_down, width = 4, height = 3.6)
  print(p_down)
  dev.off()
  message("已保存下调基因 GO BP bubble plot: ", pdf_down)
} else {
  message("没有下调方向的显著通路，未绘制 Down 图。")
}
