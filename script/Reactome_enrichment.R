suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(purrr)
})

## ========== 0) 参数 & 路径 ==========

base_dir <- "path_to_de_results"

out_dir <- file.path(
  base_dir,
  "Reactome_all_tissues_bubble_per_TissueGender"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## DE 基因筛选阈值
padj_thr  <- 0.05
lfc_thr   <- log2(1.5)

## Reactome 富集阈值（在 run_enrichReactome 中先用它筛一遍）
reactome_padj_thr <- 0.1

## 需要遍历的 Tissue / Gender
tissues <- c("Liver", "Muscle", "Skin")
genders <- c("Male", "Female")

## ========== 1) 小工具：读 DE 结果并返回 up/down 基因（SYMBOL） ==========

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
    dplyr::filter(
      !is.na(gene),
      !is.na(log2FoldChange),
      !is.na(padj)
    )

  up_genes <- df %>%
    dplyr::filter(padj < padj_thr,
                  log2FoldChange >= lfc_thr) %>%
    dplyr::pull(gene) %>%
    unique()

  down_genes <- df %>%
    dplyr::filter(padj < padj_thr,
                  log2FoldChange <= -lfc_thr) %>%
    dplyr::pull(gene) %>%
    unique()

  universe <- df$gene %>% unique()

  list(
    up_genes   = up_genes,
    down_genes = down_genes,
    universe   = universe
  )
}

## ========== 2) 小工具：对某个基因集做 Reactome 富集，并转为 long df ==========

run_enrichReactome <- function(genes, universe, tissue, gender, group, direction,
                               padj_cut = reactome_padj_thr) {
  ## genes / universe 为 SYMBOL，需要先转 Entrez
  if (length(genes) < 10) {
    return(NULL)
  }

  ## 目标基因 SYMBOL -> ENTREZ
  gene_map <- tryCatch(
    bitr(genes,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = org.Mm.eg.db),
    error = function(e) NULL
  )
  if (is.null(gene_map) || nrow(gene_map) < 10) {
    return(NULL)
  }
  gene_entrez <- unique(gene_map$ENTREZID)

  ## 背景集 universe SYMBOL -> ENTREZ
  uni_map <- tryCatch(
    bitr(universe,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = org.Mm.eg.db),
    error = function(e) NULL
  )
  if (is.null(uni_map) || nrow(uni_map) < 10) {
    return(NULL)
  }
  universe_entrez <- unique(uni_map$ENTREZID)

  ## Reactome 富集
  er <- tryCatch(
    enrichPathway(
      gene          = gene_entrez,
      universe      = universe_entrez,
      organism      = "mouse",
      pvalueCutoff  = 1,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1,
      readable      = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(er)) return(NULL)

  df <- as.data.frame(er)
  if (nrow(df) == 0) return(NULL)

  df <- df %>%
    dplyr::filter(!is.na(p.adjust),
                  p.adjust < padj_cut)
  if (nrow(df) == 0) return(NULL)

  df <- df %>%
    dplyr::mutate(
      Tissue    = tissue,
      Gender    = gender,
      group     = group,      # "WT" 或 "XL"
      direction = direction   # "Up" / "Down"
    )

  df
}

## ========== 3) 遍历所有 Tissue×Gender，收集 Reactome 结果 ==========

react_all_list <- list()

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
    er_wt_up   <- run_enrichReactome(de_wt$up_genes,   de_wt$universe,
                                     tissue = t, gender = g,
                                     group = "WT", direction = "Up")
    er_wt_down <- run_enrichReactome(de_wt$down_genes, de_wt$universe,
                                     tissue = t, gender = g,
                                     group = "WT", direction = "Down")

    ## XL
    de_xl <- load_de_and_split_ud(f_xl, padj_thr = padj_thr, lfc_thr = lfc_thr)
    er_xl_up   <- run_enrichReactome(de_xl$up_genes,   de_xl$universe,
                                     tissue = t, gender = g,
                                     group = "XL", direction = "Up")
    er_xl_down <- run_enrichReactome(de_xl$down_genes, de_xl$universe,
                                     tissue = t, gender = g,
                                     group = "XL", direction = "Down")

    react_all_list <- c(
      react_all_list,
      list(er_wt_up, er_wt_down, er_xl_up, er_xl_down)
    )
  }
}

react_all <- react_all_list %>%
  purrr::compact() %>%
  dplyr::bind_rows()

if (nrow(react_all) == 0) {
  stop("没有任何 Reactome term 通过当前阈值，请调高 q/p 阈值后重试。")
}

## 保存一份 long 表，方便以后复用
react_tsv <- file.path(out_dir, "Reactome_enrich_all_tissues_long.tsv")
write.table(
  react_all,
  file      = react_tsv,
  sep       = "\t",
  row.names = FALSE
)
message("已写出 Reactome 富集结果: ", react_tsv)

## 添加一些通用列
react_all <- react_all %>%
  dplyr::mutate(
    contrast = paste(Tissue, Gender, group, sep = "_"),   # e.g. Liver_Male_WT
    neg_log10_padj = -log10(p.adjust),
    Description_wrap = stringr::str_wrap(Description, width = 60)
  )

## ========== 4) 每个 Tissue×Gender×Group 各选3个通路，画 Up / Down 图 ==========

## -------- Up：每个(Tissue,Gender,group) 选3条，上调 --------

react_up <- react_all %>%
  dplyr::filter(direction == "Up")

if (nrow(react_up) > 0) {

  react_up_top_by_contrast <- react_up %>%
    dplyr::group_by(Tissue, Gender, group, contrast, Description_wrap) %>%
    dplyr::summarise(
      min_padj = min(p.adjust, na.rm = TRUE),
      .groups  = "drop_last"
    ) %>%
    dplyr::arrange(min_padj) %>%
    dplyr::slice_head(n = 3) %>%   # 每个 Tissue×Gender×group 各最多 3 个
    dplyr::ungroup()

  up_selected_terms <- react_up_top_by_contrast %>%
    dplyr::pull(Description_wrap) %>%
    unique()

  react_up_plot <- react_up %>%
    dplyr::filter(Description_wrap %in% up_selected_terms)

  term_order_up <- react_up_plot %>%
    dplyr::group_by(Description_wrap) %>%
    dplyr::summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(min_padj) %>%
    dplyr::pull(Description_wrap)

  react_up_plot$Description_wrap <- factor(
    react_up_plot$Description_wrap,
    levels = rev(term_order_up)   # 最显著放上面
  )

  contrast_order_up <- react_up_plot %>%
    dplyr::distinct(Tissue, Gender, group, contrast) %>%
    dplyr::arrange(Tissue, Gender, group) %>%
    dplyr::pull(contrast)

  react_up_plot$contrast <- factor(
    react_up_plot$contrast,
    levels = contrast_order_up
  )

  p_up <- ggplot(
    react_up_plot,
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
    ggtitle("Reactome enrichment (Upregulated genes)\nEach Tissue×Gender×Group: top 3 terms")

  pdf_up <- file.path(out_dir, "Reactome_Up_top3_per_TissueGenderGroup_all_tissues.pdf")
  pdf(pdf_up, width = 7, height = 5)
  print(p_up)
  dev.off()
  message("已保存上调基因 Reactome bubble plot: ", pdf_up)
} else {
  message("没有上调方向的显著 Reactome 通路，未绘制 Up 图。")
}

## -------- Down：每个(Tissue,Gender,group) 选3条，下调 --------

react_down <- react_all %>%
  dplyr::filter(direction == "Down")

if (nrow(react_down) > 0) {

  react_down_top_by_contrast <- react_down %>%
    dplyr::group_by(Tissue, Gender, group, contrast, Description_wrap) %>%
    dplyr::summarise(
      min_padj = min(p.adjust, na.rm = TRUE),
      .groups  = "drop_last"
    ) %>%
    dplyr::arrange(min_padj) %>%
    dplyr::slice_head(n = 3) %>%
    dplyr::ungroup()

  down_selected_terms <- react_down_top_by_contrast %>%
    dplyr::pull(Description_wrap) %>%
    unique()

  react_down_plot <- react_down %>%
    dplyr::filter(Description_wrap %in% down_selected_terms)

  term_order_down <- react_down_plot %>%
    dplyr::group_by(Description_wrap) %>%
    dplyr::summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(min_padj) %>%
    dplyr::pull(Description_wrap)

  react_down_plot$Description_wrap <- factor(
    react_down_plot$Description_wrap,
    levels = rev(term_order_down)
  )

  contrast_order_down <- react_down_plot %>%
    dplyr::distinct(Tissue, Gender, group, contrast) %>%
    dplyr::arrange(Tissue, Gender, group) %>%
    dplyr::pull(contrast)

  react_down_plot$contrast <- factor(
    react_down_plot$contrast,
    levels = contrast_order_down
  )

  p_down <- ggplot(
    react_down_plot,
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
    ggtitle("Reactome enrichment (Downregulated genes)\nEach Tissue×Gender×Group: top 3 terms")

  pdf_down <- file.path(out_dir, "Reactome_Down_top3_per_TissueGenderGroup_all_tissues.pdf")
  pdf(pdf_down, width = 7, height = 5)
  print(p_down)
  dev.off()
  message("已保存下调基因 Reactome bubble plot: ", pdf_down)
} else {
  message("没有下调方向的显著 Reactome 通路，未绘制 Down 图。")
}

