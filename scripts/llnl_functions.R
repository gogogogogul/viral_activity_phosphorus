# General functions used from the custom llnl library I created for this project

library(vegan)
library(ggpubr)
library(ggplot2) 
library(dplyr)
library(rlang) 

abundance.normal <-
function (mean, reads.df, all_reads_col, mean_substring_rmv = ".Trimmed.Mean", 
    cov_int = 0) 
{
    cov <- as.data.frame(lapply(mean, function(x) ifelse(x >= 
        cov_int, 1, 0)))
    mean.cov <- mean * cov
    colnames(mean.cov) <- gsub(mean_substring_rmv, "", colnames(mean.cov))
    norm.counts <- data.frame(rows = 1:nrow(mean.cov))
    for (i in rownames(reads.df)) {
        new_col <- mean.cov[i]/(as.numeric(reads.df[i, ][all_reads_col])) * 
            10^9
        norm.counts <- cbind.data.frame(norm.counts, new_col)
    }
    return(norm.counts[, -which(names(norm.counts) == "rows")])
}
add_sheet_to_excel <-
function (file, sheet_name, data, overwrite_sheet = TRUE) 
{
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop("Package 'openxlsx' is required. Install it with install.packages('openxlsx')")
    }
    wb <- openxlsx::loadWorkbook(file)
    if (overwrite_sheet && sheet_name %in% names(wb)) {
        openxlsx::removeWorksheet(wb, sheet_name)
    }
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, data)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    invisible(wb)
}
alpha.div.plot <-
function (df, stats.df, x, y, fill, color, xlab, ylab, legend.title.color, 
    legend.title.fill, box.alpha = 0, box.line.size = 0.6, 
    box.outliers = FALSE, jitter.width = 0.15, 
    dodge.width = 0.9, jitter.height = 0, jitter.alpha = 0.7, 
    jitter.size = 2, rotate = 0, hide.ns = TRUE, stat.size = 6, 
    stat.step = 0.1, stat.tip = 0, stat.backet.nudge = 0, 
    legend = "none") 
{
    ggpubr::ggboxplot(df, x = x, y = y, fill = fill, color = fill, 
        outliers = box.outliers, alpha = box.alpha, size = box.line.size) + 
        ggplot2::geom_jitter(position = position_jitterdodge(jitter.width = jitter.width, 
            dodge.width = dodge.width, jitter.height = jitter.height), 
            alpha = jitter.alpha, size = jitter.size, aes(color = !!sym(fill))) + 
        xlab(xlab) + ylab(ylab) + guides(x = guide_axis(angle = rotate)) + 
        stat_pvalue_manual(stats.df, label = "p.adj.signif", 
            hide.ns = hide.ns, size = stat.size, step.increase = stat.step, 
            tip.length = stat.tip, bracket.nudge.y = stat.backet.nudge) + 
        theme(legend.position = legend) + font("axis.text") + 
        font("xylab") + scale_color_manual(values = col.lnll) + 
        guides(color = guide_legend(title = legend.title.color), 
            fill = guide_legend(title = legend.title.fill))
}
col.lnll <-
c(T0 = "burlywood4", T1 = "darkblue", T2 = "steelblue2", T3 = "lightblue", 
Phosphorus = "#A71930", `No Phosphorus` = "gray70", Dry = "burlywood4", 
Virome = "#1D1160", `Transcriptionally Active` = "darkgreen", 
Metatranscriptome = "darkgreen")
dist.table <-
function (abundance.table, singleton, dist.method, data.transform = "hellinger") 
{
    if ({
        {
            singleton
        }
    } == "Y") {
        df_votu <- abundance.table[rowSums((abundance.table > 
            0) * 1) > 1, ]
    }
    else if ({
        {
            singleton
        }
    } == "N") {
        df_votu <- abundance.table
    }
    else {
        print("Singleton is not Y or N")
    }
    df_votu <- t(df_votu)
    print(paste("Samples before droping samples without votus", 
        nrow(df_votu)))
    dropped <- df_votu[apply(df_votu, 1, function(x) all(x == 
        0)), ]
    df_votu.keep <- df_votu[apply(df_votu, 1, function(x) !all(x == 
        0)), ]
    print(paste("Samples after droping samples without votus", 
        nrow(df_votu.keep)))
    print(paste("Samples without votus", nrow(dropped)))
    df_votu_xform <- decostand(df_votu.keep, method = data.transform)
    df_votu_dist <- as.matrix(vegdist(df_votu_xform, method = dist.method))
    return(df_votu_dist)
}
group.xdis <-
c(T0 = 1, T1 = 2, T2 = 3, T3 = 4, T0 = 1, T1 = 2, T2 = 3, T3 = 4, 
`T0:Dry` = 1, `T1:No Phosphorus` = 1.8, `T1:Phosphorus` = 2.2, 
`T2:No Phosphorus` = 2.8, `T2:Phosphorus` = 3.2, 
`T3:No Phosphorus` = 3.8, `T3:Phosphorus` = 4.2
)

meta.mod <-
function (metadata, votu.dist.table) 
{
    meta.ordered <- metadata[order(row.names(metadata)), ]
    meta_ <- meta.ordered[(row.names(meta.ordered) %in% row.names(votu.dist.table)), 
        ]
    meta.ordered <- meta_[match(row.names(votu.dist.table), row.names(meta_)), 
        ]
    print(paste("Are samples in the metadata and the distance table in the same order?", 
        all(rownames(meta.ordered) == rownames(votu.dist.table))))
    if (all(rownames(meta.ordered) == rownames(votu.dist.table)) == 
        TRUE) {
        return(meta.ordered)
    }
    else {
        print("Samples in the metadata and distance table are not in the same order")
    }
}
normal.tukey.alpha <-
function (df, x, y, group, xlab, ylab, rotate, legend, x.dist = group.xdis, 
    fill = NULL, y.position.multi = 0.95, step.inc = 0.18, 
    p.value.txt.sz = 4, ns.hidden = TRUE, axis.text = 9, bracket.nudge = 0, 
    xy.lab = 11, col.pallete = col.lnll) 
{
    df$interactions <- interaction(df[[x]], df[[group]], sep = ":")
    A <- as.formula(paste(y, "~", x, "*", group))
    print("check the data is normal")
    print("Shapiro test:")
    lm.mod_ <- lm(A, data = df)
    print(shapiro_test(lm.mod_$residuals))
    stat.test <- df %>% lm(A, data = .) %>% aov(.) %>% tukey_hsd(.) %>% 
        drop_na()
    stat.test$xmin <- sapply(stat.test$group1, function(x) {
        if (x %in% names(x.dist)) {
            return(x.dist[x])
        }
        else {
            return(NA)
        }
    })
    stat.test$xmax <- sapply(stat.test$group2, function(x) {
        if (x %in% names(x.dist)) {
            return(x.dist[x])
        }
        else {
            return(NA)
        }
    })
    y.positions <- df %>% group_by(interactions) %>% summarize(richness.mean = mean(!!sym(y))) %>% 
        column_to_rownames("interactions")
    stat.test$y.group1 <- sapply(stat.test$group1, function(x) {
        if (x %in% rownames(y.positions)) {
            return(y.positions[x, ])
        }
        else {
            return(NA)
        }
    })
    stat.test$y.group2 <- sapply(stat.test$group2, function(x) {
        if (x %in% rownames(y.positions)) {
            return(y.positions[x, ])
        }
        else {
            return(NA)
        }
    })
    stat.test$y.position <- pmin(stat.test$y.group1, stat.test$y.group2) * 
        y.position.multi
    r <- ggboxplot(df, x = x, y = y, fill = group) + xlab(xlab) + 
        ylab(ylab) + guides(x = guide_axis(angle = rotate)) + 
        stat_pvalue_manual(stat.test, label = "p.adj.signif", 
            xmin = "xmin", xmax = "xmax", bracket.nudge.y = bracket.nudge, 
            hide.ns = ns.hidden, size = p.value.txt.sz, step.increase = step.inc, 
            tip.length = 0) + theme(legend.position = legend) + 
        font("axis.text", size = axis.text) + font("xylab", size = xy.lab) + 
        scale_fill_manual(name = group, values = col.pallete)
    return(list(r, stat.test, df))
}
permanova.table <-
function (distance.matrix, metadata.frame, variable, permutations = 999) 
{
    distance.matrix <- as.dist(distance.matrix)
    adonis2.formula <- as.formula(paste("distance.matrix ~", 
        variable))
    p_value.table <- as.data.frame(adonis2(adonis2.formula, data = metadata.frame, 
        permutations = permutations))
    p_value <- p_value.table$`Pr(>F)`
    p_value.table <- p_value.table %>% mutate(Signif = case_when(is.na(p_value) ~ 
        "", p_value <= 0.001 ~ "***", p_value <= 0.01 ~ "**", 
        p_value <= 0.05 ~ "*", p_value <= 0.1 ~ 
            ".", TRUE ~ "")) %>% mutate(F = as.numeric(F)) %>% 
        mutate(`Pr(>F)` = as.numeric(`Pr(>F)`)) %>% mutate(across(where(is.numeric), 
        ~round(., digits = 3))) %>% replace(is.na(.), "")
    table <- ggtexttable(p_value.table, theme = ttheme("light"))
    table <- tab_add_title(table, paste("PERMANOVA -", variable))
    return(table)
}
rel_abundance <-
function (votu.table, num_samples) 
{
    if (ncol(votu.table) == num_samples) {
        votu.table %>% mutate(across(everything(), ~./sum(.)))
    }
    else {
        print("Number of columns in the votu table do not match the number of samples that were input.\n          Do you need to convert your samples from rows to columns?")
    }
}
shape.lnll <-
c(Phosphorus = 15, `No Phosphorus` = 17, Dry = 16)
