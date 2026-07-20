library(magrittr)
library(survival)

# get sgi valid clusters memberships 
fget_m <- function(hc, th = -1){
  # imitate sgi to use helper funtions
  obj = sgi:::get_uc_vc(hc, th = th)
  obj$hc = hc
  sgi::get_vcps(obj)
}

# get all possible clusters from the feature hc tree
fget_all_clusters <- function(hc, returnm = F){
  m = fget_m(hc)
  m = cbind(l1 = m[,1]*0+1, m)
  # just return membership matrix
  if(returnm) return(m)
  
  # create list of clusters from m 
  ff <- function(x){
    cs = sort(unique(na.omit(x)))
    lapply( setNames(cs, cs), function(i) which(x ==i)) 
  }
  
  do.call("c", apply(m, 2, ff)) 
}

# to get sgi vc membersip Ms efficiently given different fs spaces
fget_multiple_m <- function(X, flist, minsize = 10, linkg = 'ward.D2'){
  # calculate squared distaces according to each f 
  fss = unique(unlist(flist))
  d2 = vector('list', 3)
  print(dim(X))
  print(length(fss))
  d2[fss] = lapply(fss, function(i) dist(X[,i])^2)
  
  # calculate distance matrix for each feature space 
  ds = lapply(flist, function(fs) sqrt(Reduce(`+`, d2[fs]))); rm(d2)
  
  # get ms for each sgi
  lapply(ds, function(d) fget_m(hclust(d, method = linkg), th = minsize))
}

# function to collect and return all tests given fs hc tree 
f_M <- function(hc_f, minsize = 10, linkg = 'ward.D2'){
  fs = fget_all_clusters(hc_f) %>% {.[order(-sapply(.,length))]}
  ms = fget_multiple_m(x, fs, minsize, linkg)
  # one giant M to collect all tests
  M = do.call(cbind, ms)
  apply(M, 2, function(a) a-min(a,na.rm = T))
}

# Simes procedure extended with effective number of test calculation
# multiple testing correction in multiple SGI runs 
f_correct <- function(p_values, p_cor_matrix, k_min = 50, k_length = 5, use_approx = FALSE) {
  o = order(p_values)
  p_cor_matrix = p_cor_matrix[o,o]
  p_values = p_values[o]
  
  # this part is replacible, this is based on 
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5676236/pdf/883.pdf
  get_effective_n <- function(p_cor) {
    p_cor <- as.matrix(p_cor)
    N <- nrow(p_cor)
    if (N <= 1) return(1)
    # eigen values
    Ls <- eigen(p_cor, symmetric = TRUE, only.values = TRUE)$values
    if (sum(Ls < 0) > 0) {
      Ls <- Ls[Ls > 0]
      Ls <- Ls / sum(Ls) * N
    }
    Ls <- Ls[Ls > 1]
    N - sum(Ls - 1)
  }
  
  k <- length(p_values)
  
  if(use_approx && k_min < k){
    # approximate eff n by k 
    ks <- c(seq(k_min), unique(round(seq(k_min + 1, k, length.out = k_length))))
    # use this exact values to approximate for ks
    exact_eff_ns <- sapply(ks, function(i) get_effective_n(p_cor_matrix[1:i, 1:i]))
    print("EFFECTIVE NUMBER OF TESTS")
    print(rev(exact_eff_ns)[1])
    eff_ns <- approx(x = ks, y = exact_eff_ns, xout = seq(k))$y
  }else {
    # exact computation 
    eff_ns <- sapply(seq(k), function(i) get_effective_n(p_cor_matrix[1:i, 1:i]))
  }
  # enforce monotocity
  eff_ns <- cummax(eff_ns)
  
  
  global_eff_n <- eff_ns[k]
  corrected_p_values <- pmin((global_eff_n * p_values) / eff_ns, 1)
  cummax(corrected_p_values)[order(o)]
}

#outcome tests
ff_factor <- function(x, y){
  x = x[!is.na(y)]
  y = y[!is.na(y)]
  
  # no level
  if( length(levels(x) )<2 ) {
    return(1)
  }
  
  # too many level
  if(length(levels(x)) == length(x)){
    warning("# of factor levels equals number of samples: length(levels(x)) = length(x)")
    return(1)
  }
  
  # 2 level: fisher test
  if(length(levels(x)) == 2){
    fu=try(fisher.test(x, y=y))
    
    if(inherits(fu, "try-error")){
      fu=try(fisher.test(x, y=y, simulate.p.value = T))
      if(inherits(fu, "try-error")) return(1)
    }
    pval=fu$p.value
    stat=fu$estimate
  }else{ # more than 2 levels
    cht = try( tryCatch( chisq.test(x, y=y), warning=function(w) w) )
    if(inherits(cht, "try-error")) return(1)
    
    if(inherits(cht, "warning")){ # if warning for simulate p val
      cht = try( tryCatch( chisq.test(x, y=y, simulate.p.value = T), warning=function(w) w) )
      if(inherits(cht, "try-error") | inherits(cht, "warning")) return(1)
    }
    pval=cht$p.value
    stat=cht$statistic
  }
  
  if(is.null(pval)) {pval=1}
  if(is.null(stat)) {stat=NA}
  
  return(unname(pval))
}
# pvalues of all SGI runs 
ff_numeric <- function(y, m){
  y = y[!is.na(m)]
  m = m[!is.na(m)]
  e = try(t.test(y~m)$p.value)
  if(inherits(e, 'try-error')) return(1)
  e
}
library(rms)
library(stringr)
assignInNamespace('formatNP', function(x, ...) x, 'rms')
get_model_stats <- function(x, precision = 60) {
  # remember old number formatting function
  # (which would round and transforms p-values to formats like "<0.01")
  old_format_np <- rms::formatNP
  # substitute it with a function which will print out as many digits as we want
  assignInNamespace("formatNP", function(x, ...)
    formatC(x, format = "f", digits = precision), "rms")
  
  # remember old width setting
  old_width <- options('width')$width
  # substitute it with a setting making sure the table will not wrap
  options(width = old_width + 4 * precision)
  
  # actually print the data and capture it
  cap <- capture.output(print(x))
  
  # restore original settings
  options(width = old_width)
  assignInNamespace("formatNP", old_format_np, "rms")
  
  # model stats
  stats <- c()
  stats$R2.adj <-
    str_match(cap, "R2 adj\\s+ (\\d\\.\\d+)") %>% na.omit() %>% .[, 2] %>% as.numeric()
  
  # coef stats lines
  coef_lines <-
    cap[which(str_detect(cap, "Coef\\s+S\\.E\\.")):(length(cap) - 1)]
  
  # parse
  coef_lines_table <-
    
    suppressWarnings(readr::read_table(coef_lines %>% stringr::str_c(collapse = "\n")))
  colnames(coef_lines_table)[1] <- "Predictor"
  
  list(stats = stats,
       coefs = coef_lines_table)
}

wilcox_test <- function(v, y) {
  # try to run wilcox.test  
  res <- tryCatch({
    this_fit <- rms::orm(v~y, family = "probit")
    this_res <-
      get_model_stats(this_fit) %>% .$coefs %>% data.frame()
    # format results
    names(this_res) <-
      c('analyte',
        'std_error',
        'estimate',
        'statistic',
        'p_value',
        'other')
    # return results if test
    this_res$p_value[nrow(this_res)]
    
  }, error = function(e){
    print(e$message)
    # something didn't work with the wilcox.test (probably no samples in one of the groups); return NA vector
    1
  })
  # return results
  return(res)
}
simes_correct <- function(p_values) {
  o = order(p_values)
  pvals = p_values[o]
  for (i in 1:length(pvals)) {
    pvals[i] = length(pvals)/i * pvals[i]
  }
  pvals = cummax(pvals)
  return(pvals[order(o)])
}
bonf_correct <- function(p_values) {
  return(p_values * length(p_values))
}
get_rat <- function(v) {
  round(1000 * sum(v <= 0.05)/length(v))
}

get_vals <- function(v) {
  pv = get_rat(v)
  vals = c(seq(5, pv, length.out = 7), seq(pv + 1, 1000 - pv, length.out = 8))
  round(vals)
}

get_corrections <- function(ps, vals, Mc) {
  lapply(vals, function(i){
    set.seed(1234)
    brkpoint = which(ps[order(ps)] > 0.05)[1]
    idxs = sample(brkpoint:length(ps), 1000-i)
    before = sample(seq(brkpoint - 1), i)
    ps1 = ps[order(ps)][c(before, idxs)]
    
    M1 = Mc[order(ps)[c(before, idxs)],order(ps)[c(before, idxs)]]
    
    # calculate correlations between tests 
    xc = M1
    
    # multiple testing correction 
    ph = f_correct(ps1, xc, use_approx = T)
    o = order(ps1)
    #plot(seq(length(ph)), ph[o], col = "red", lwd = 1, type = "l", ylab = "Corrected p-values", xlab = "")
    #title(paste0("Approximation for Simes Correction in QMDiab Plasma, ", i, "<= 0.05"))
    #lines(seq(length(ph)), ph2[o], col = "blue", lwd = 1, type = "l")
    #lines(seq(length(ph)), ph[o], col = "red", lwd = 1, type = "l")
    #legend(x = "bottomright",legend=c("Approximated", "Exact"),  
    #fill = c("red","blue") 
    #)
    #abline(h = 0.05)
    #print("ORIGINAL")
    #print(length(ps1[ps1 <= 0.05]))
    #print("AFTER CORRECTION")
    #print(length(ph[ph <= 0.05]))
    df = data.frame(approx = ph[o], simes = simes_correct(ps1), bonf = bonf_correct(ps1))
    df
  })
}
stretch_pval <- function(p) {
  ifelse(p <= 0.05, p / 0.05 * 0.5, 0.5 + (p - 0.05) / (1 - 0.05) * 0.5)
}

inv_stretch <- function(y) {
  ifelse(y <= 0.5, y * 0.05 / 0.5, 0.05 + (y - 0.5) * (1 - 0.05) / 0.5)
}

get_comb_plots2 <- function(cors, vals, label = "A") {
  plots = lapply(seq_along(cors), function(i){
    df = cors[[i]]
    colnames(df) = c("approx", "simes", "bonf")
    df = df[,c(1,2)]
    df$x = 1:1000
    df$title = paste0(vals[i], " initial significant p-values")
    
    # Sort each method’s p-values
    df <- df %>%
      mutate(across(c(approx, exact), sort))
    
    # Long format for plotting
    df_long <- df %>%
      pivot_longer(cols = c(approx, exact),
                   names_to = "method", values_to = "pval")
    
    df_long$method <- factor(df_long$method,
                             levels = c("exact", "approx"),
                             labels = c("Extended Simes (exact)",
                                        "Extended Simes (approximated)"))
    
    # Apply nonlinear transformation for p-values
    df_long$stretched_pval <- stretch_pval(df_long$pval)
    
    ggplot(df_long, aes(x = x, y = stretched_pval, color = method)) +
      geom_line(size = 1.1, alpha = 0.95) +
      scale_color_manual(
        values = c("Extended Simes (exact)" = "darkgreen",
                   "Extended Simes (approximated)" = "pink")
      ) +
      facet_grid(~ title) +
      xlim(0, 1000) +
      scale_y_continuous(
        name = "p-value\n",
        limits = c(0, 1),
        breaks = stretch_pval(c(0, 0.01, 0.025, 0.05, 0.25, 0.5, 1)),
        labels = signif(inv_stretch(stretch_pval(c(0, 0.01, 0.025, 0.05, 0.25, 0.5, 1))), 3)
      ) +
      geom_hline(yintercept = stretch_pval(0.05), linetype = 'dotted', col = 'black') +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "wheat", colour = "black"),
        panel.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(size = 22),
        legend.text = element_text(size = 17),
        legend.key.size = unit(3, "line")
      )
  })
  
  # Combine all panels
  pg = plot_grid(plotlist = lapply(seq_along(plots), function(i){
    p = plots[[i]] + theme(legend.position = "none")
    
    # Remove x or y labels depending on position
    if (!(i %in% c(1, 4, 7))) {
      p = p + labs(y = "\n") +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    if (i <= 0) {
      p = p + theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
    }
    p
  }), labels = c(""), label_size = 20, nrow = 3)
  
  # Add legend at bottom
  pg = plot_grid(
    pg,
    get_plot_component(plots[[1]] + theme(legend.position = 'bottom'), 'guide-box-bottom', return_all = TRUE),
    ncol = 1,
    rel_heights = c(1, 0.1)
  )
  
  pg
}

get_comb_plots3 <- function(cors, vals, label = "A") {
  plots <- lapply(seq_along(cors), function(i){
    df <- cors[[i]]
    colnames(df) <- c("approx", "simes", "bonf")
    
    # Keep only approx, simes, bonf
    df <- df[, c("approx", "simes", "bonf")]
    df$x <- 1:1000
    df$title <- paste0(vals[i], " initial significant p-values")
    
    # Sort each method's p-values
    df <- df %>%
      mutate(across(c(approx, simes, bonf), sort))
    
    # Pivot longer
    df_long <- df %>%
      pivot_longer(cols = c("bonf", "simes", "approx"),
                   names_to = "method", values_to = "pval")
    
    # Factor and rename methods with desired legend order
    df_long$method <- factor(df_long$method,
                             levels = c("bonf", "simes", "approx"),
                             labels = c("Bonferroni",
                                        "Simes",
                                        "Extended Simes (approximated)"))
    
    # Apply nonlinear transformation
    df_long$stretched_pval <- stretch_pval(df_long$pval)
    
    ggplot(df_long, aes(x = x, y = stretched_pval, color = method)) +
      geom_line(size = 1.1, alpha = 0.95) +
      scale_color_manual(values = c(
        "Bonferroni" = "orange",
        "Simes" = "blue",
        "Extended Simes (approximated)" = "pink"
      )) +
      facet_grid(~ title) +
      xlim(0, 1000) +
      scale_y_continuous(
        name = "p-value\n",
        limits = c(0, 1),
        breaks = stretch_pval(c(0, 0.01, 0.025, 0.05, 0.25, 0.5, 1)),
        labels = signif(inv_stretch(stretch_pval(c(0, 0.01, 0.025, 0.05, 0.25, 0.5, 1))), 3)
      ) +
      geom_hline(yintercept = stretch_pval(0.05), linetype = 'dotted', col = 'black') +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "wheat", colour = "black"),
        panel.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(size = 22),
        legend.text = element_text(size = 17),
        legend.key.size = unit(3, "line")
      )
  })
  
  # Combine all panels: 5 rows for 15 plots
  pg <- plot_grid(plotlist = lapply(seq_along(plots), function(i){
    p <- plots[[i]] + theme(legend.position = "none")
    if (!(i %% 3 == 1)) {
      p <- p + labs(y = "\n") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    p
  }), labels = c(""), label_size = 20, nrow = 5)
  
  # Add legend at bottom
  pg <- plot_grid(
    pg,
    get_plot_component(plots[[1]] + theme(legend.position = 'bottom'), 'guide-box-bottom', return_all = TRUE),
    ncol = 1,
    rel_heights = c(1, 0.1)
  )
  
  pg
}

set.seed(42)

x <- rosmap_data %>% as.matrix()

idxs <- sample(
  seq_len(nrow(x)),
  min(250, nrow(x))
)

x <- x[idxs, , drop = FALSE]
x <- x %>% scale()

hc_f <- hclust(
  dist(t(x)),
  method = "ward.D2"
)

M <- f_M(hc_f)

Mc <- cor(M, use = "pairwise.complete.obs")
Mc[is.na(Mc)] <- 0
Mc <- abs(Mc)

ordinal_y <- rosmap_clins$braaksc[idxs]

ordinal_ps <- apply(
  M,
  2,
  wilcox_test,
  v = ordinal_y
)

ordinal_ps2 <- ordinal_ps[!is.na(ordinal_ps)]

ordinal_corrections <- get_corrections(
  ordinal_ps2,
  get_vals(ordinal_ps2),
  Mc
)

ordinal_corrections2 <- lapply(ordinal_corrections, function(a) {
  b <- a
  b$simes <- pmin(unlist(b$simes), 1)
  b$bonf  <- pmin(unlist(b$bonf), 1)
  b
})

panel_indices <- c(1, 3, 5, 7, 9, 11, 13, 14, 15)

pdf(
  "supplementary_figure_7.pdf",
  width = 14,
  height = 23.5 * 3 / 5
)

get_comb_plots(
  ordinal_corrections2[panel_indices],
  get_vals(ordinal_ps2)[panel_indices],
  "C"
)

dev.off()

pdf(
  "supplementary_figure_8.pdf",
  width = 14,
  height = 23.5 * 3 / 5
)

get_comb_plots3(
  ordinal_corrections2[panel_indices],
  get_vals(ordinal_ps2)[panel_indices],
  "C"
)


set.seed(42)

x <- adni_bl_lipids %>% as.matrix()

idxs <- sample(
  seq_len(nrow(x)),
  min(250, nrow(x))
)

x <- x[idxs, , drop = FALSE]
x <- x %>% scale()

hc_f <- hclust(
  dist(t(x)),
  method = "ward.D2"
)

M <- f_M(hc_f)

Mc <- cor(M, use = "pairwise.complete.obs")
Mc[is.na(Mc)] <- 0
Mc <- abs(Mc)

ordinal_y <- adni_clinsl0$bmi[idxs]

ordinal_ps <- apply(
  M,
  2,
  wilcox_test,
  v = ordinal_y
)

ordinal_ps2 <- ordinal_ps[!is.na(ordinal_ps)]

ordinal_corrections <- get_corrections(
  ordinal_ps2,
  get_vals(ordinal_ps2),
  Mc
)

ordinal_corrections2 <- lapply(ordinal_corrections, function(a) {
  b <- a
  b$simes <- pmin(unlist(b$simes), 1)
  b$bonf  <- pmin(unlist(b$bonf), 1)
  b
})

panel_indices <- c(1, 3, 5, 7, 9, 11, 13, 14, 15)

pdf(
  "supplementary_figure_9.pdf",
  width = 14,
  height = 23.5 * 3 / 5
)

get_comb_plots(
  ordinal_corrections2[panel_indices],
  get_vals(ordinal_ps2)[panel_indices],
  "C"
)

dev.off()


