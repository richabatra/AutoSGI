library(magrittr)
library(survival)

# get sgi valid clusters memberships 
#wrapper that mimics internal sgi package logic to extract valid cluster (vc) memberships from a hclust tree
fget_m <- function(hc, th = -1){
  # imitate sgi to use helper funtions
  #get upper/valid cluster split object from the sgi package's internal function
  obj = sgi:::get_uc_vc(hc, th = th)
  #attach the hclust tree to the object so downstream sgi functions can use it
  obj$hc = hc
  #return the valid cluster partition memberships
  sgi::get_vcps(obj)
}

# get all possible clusters from the feature hc tree
#build the full list of clusters (as sample index sets) implied by every level of the hclust tree
fget_all_clusters <- function(hc, returnm = F){
  #get the membership matrix from fget_m
  m = fget_m(hc)
  #add a first "level 1" column where every sample belongs to a single root cluster
  m = cbind(l1 = m[,1]*0+1, m)
  # just return membership matrix
  #optionally return the raw membership matrix instead of a cluster list
  if(returnm) return(m)
  
  # create list of clusters from m 
  #helper to convert one membership column into a named list of sample-index sets, one per cluster label
  ff <- function(x){
    #unique non-NA cluster labels in this column, sorted
    cs = sort(unique(na.omit(x)))
    #for each cluster label, get the indices of samples belonging to it
    lapply( setNames(cs, cs), function(i) which(x ==i)) 
  }
  
  #apply ff to every column of m and flatten all resulting cluster lists into one combined list
  do.call("c", apply(m, 2, ff)) 
}

# to get sgi vc membersip Ms efficiently given different fs spaces
#compute valid-cluster membership matrices for multiple feature subsets (flist), reusing precomputed squared distances
fget_multiple_m <- function(X, flist, minsize = 10, linkg = 'ward.D2'){
  # calculate squared distaces according to each f 
  #unique feature indices used across all feature subsets, so distances are computed only once per feature
  fss = unique(unlist(flist))
  #placeholder list to store per-feature squared distance matrices
  d2 = vector('list', 3)
  #debug print of data dimensions
  print(dim(X))
  #debug print of number of unique features
  print(length(fss))
  #for each single feature, compute the squared pairwise distance contribution across samples
  d2[fss] = lapply(fss, function(i) dist(X[,i])^2)
  
  # calculate distance matrix for each feature space 
  #for each feature subset, sum the squared per-feature distances and take sqrt to get the combined distance matrix
  ds = lapply(flist, function(fs) sqrt(Reduce(`+`, d2[fs]))); rm(d2)
  
  # get ms for each sgi
  #for each feature-subset distance matrix, run hierarchical clustering and extract valid cluster memberships
  lapply(ds, function(d) fget_m(hclust(d, method = linkg), th = minsize))
}

# function to collect and return all tests given fs hc tree 
#build the full test membership matrix M across all clusters derived from the feature hclust tree
f_M <- function(hc_f, minsize = 10, linkg = 'ward.D2'){
  #get all clusters implied by the feature tree, sorted from largest to smallest
  fs = fget_all_clusters(hc_f) %>% {.[order(-sapply(.,length))]}
  #compute sgi valid-cluster membership matrices for each of these feature-defined clusters
  ms = fget_multiple_m(x, fs, minsize, linkg)
  # one giant M to collect all tests
  #combine all membership matrices column-wise into one big matrix representing every test
  M = do.call(cbind, ms)
  #shift each column so its minimum value is 0 (normalizing cluster labels)
  apply(M, 2, function(a) a-min(a,na.rm = T))
}

# Simes procedure extended with effective number of test calculation
# multiple testing correction in multiple SGI runs 
#apply an extended Simes procedure that accounts for correlation between tests via an "effective number of tests"
f_correct <- function(p_values, p_cor_matrix, k_min = 50, k_length = 5, use_approx = FALSE) {
  #order p-values from smallest to largest
  o = order(p_values)
  #reorder the correlation matrix rows/columns to match sorted p-value order
  p_cor_matrix = p_cor_matrix[o,o]
  #sort the p-values themselves
  p_values = p_values[o]
  
  # this part is replacible, this is based on 
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5676236/pdf/883.pdf
  #compute the effective number of independent tests from a correlation matrix via eigenvalue decomposition
  get_effective_n <- function(p_cor) {
    #ensure matrix format
    p_cor <- as.matrix(p_cor)
    #number of tests in this block
    N <- nrow(p_cor)
    #trivial case: a single test has an effective number of 1
    if (N <= 1) return(1)
    # eigen values
    #compute eigenvalues of the correlation matrix (symmetric, values only)
    Ls <- eigen(p_cor, symmetric = TRUE, only.values = TRUE)$values
    #if any eigenvalues are negative (numerical issue), keep only positive ones and rescale to sum to N
    if (sum(Ls < 0) > 0) {
      Ls <- Ls[Ls > 0]
      Ls <- Ls / sum(Ls) * N
    }
    #keep only eigenvalues greater than 1 (contribute to redundancy)
    Ls <- Ls[Ls > 1]
    #effective number of tests: total minus the "excess" from correlated eigenvalues
    N - sum(Ls - 1)
  }
  
  #total number of tests
  k <- length(p_values)
  
  #if approximation requested and there are more tests than the minimum threshold
  if(use_approx && k_min < k){
    # approximate eff n by k 
    #choose a sparse set of block sizes to compute exactly, then interpolate the rest
    ks <- c(seq(k_min), unique(round(seq(k_min + 1, k, length.out = k_length))))
    # use this exact values to approximate for ks
    #compute the exact effective number of tests at each chosen block size
    exact_eff_ns <- sapply(ks, function(i) get_effective_n(p_cor_matrix[1:i, 1:i]))
    #debug print label
    print("EFFECTIVE NUMBER OF TESTS")
    #print the effective n at the largest block size
    print(rev(exact_eff_ns)[1])
    #linearly interpolate effective n for every block size from k_min to k
    eff_ns <- approx(x = ks, y = exact_eff_ns, xout = seq(k))$y
  }else {
    # exact computation 
    #compute the exact effective number of tests at every possible block size (expensive but exact)
    eff_ns <- sapply(seq(k), function(i) get_effective_n(p_cor_matrix[1:i, 1:i]))
  }
  # enforce monotocity
  #ensure the effective number of tests never decreases as more tests are added (cumulative max)
  eff_ns <- cummax(eff_ns)
  
  
  #the effective number of tests using the full test set
  global_eff_n <- eff_ns[k]
  #Simes-style correction scaled by the ratio of global to local effective test counts, capped at 1
  corrected_p_values <- pmin((global_eff_n * p_values) / eff_ns, 1)
  #enforce monotonicity on corrected p-values and restore original (unsorted) order
  cummax(corrected_p_values)[order(o)]
}

#outcome tests
#association test between a categorical predictor x and outcome y (Fisher's exact or chi-squared depending on levels)
ff_factor <- function(x, y){
  #remove samples where the outcome is missing
  x = x[!is.na(y)]
  y = y[!is.na(y)]
  
  # no level
  #if x has fewer than 2 levels, there's nothing to test, so return p = 1
  if( length(levels(x) )<2 ) {
    return(1)
  }
  
  # too many level
  #if every sample is its own level, the test isn't meaningful; warn and return p = 1
  if(length(levels(x)) == length(x)){
    warning("# of factor levels equals number of samples: length(levels(x)) = length(x)")
    return(1)
  }
  
  # 2 level: fisher test
  #for a binary factor, use Fisher's exact test
  if(length(levels(x)) == 2){
    #attempt standard Fisher's exact test
    fu=try(fisher.test(x, y=y))
    
    #if it fails (e.g., table too large), fall back to simulated p-value version
    if(inherits(fu, "try-error")){
      fu=try(fisher.test(x, y=y, simulate.p.value = T))
      #if that also fails, give up and return p = 1
      if(inherits(fu, "try-error")) return(1)
    }
    #extract p-value and effect size estimate
    pval=fu$p.value
    stat=fu$estimate
  }else{ # more than 2 levels
    #for more than 2 levels, use chi-squared test, catching warnings (e.g., low expected counts) as well as errors
    cht = try( tryCatch( chisq.test(x, y=y), warning=function(w) w) )
    #if the test throws an actual error, return p = 1
    if(inherits(cht, "try-error")) return(1)
    
    #if a warning was raised (e.g., approximation may be incorrect), retry with simulated p-value
    if(inherits(cht, "warning")){ # if warning for simulate p val
      cht = try( tryCatch( chisq.test(x, y=y, simulate.p.value = T), warning=function(w) w) )
      #if that also errors or warns, give up and return p = 1
      if(inherits(cht, "try-error") | inherits(cht, "warning")) return(1)
    }
    #extract p-value and test statistic
    pval=cht$p.value
    stat=cht$statistic
  }
  
  #safety fallback: if p-value came back NULL, default to 1
  if(is.null(pval)) {pval=1}
  #safety fallback: if statistic came back NULL, default to NA
  if(is.null(stat)) {stat=NA}
  
  #return the unnamed p-value
  return(unname(pval))
}
# pvalues of all SGI runs 
#simple two-group t-test between numeric outcome y and binary cluster membership m
ff_numeric <- function(y, m){
  #drop samples with missing cluster membership
  y = y[!is.na(m)]
  m = m[!is.na(m)]
  #run t-test of y by group m, catching errors
  e = try(t.test(y~m)$p.value)
  #if the test fails (e.g., one group empty), return p = 1
  if(inherits(e, 'try-error')) return(1)
  e
}
library(rms)
library(stringr)
#override rms's internal p-value formatting function so it stops rounding/truncating p-values in printed output
assignInNamespace('formatNP', function(x, ...) x, 'rms')
#extract regression model statistics (R2 and coefficient table) by parsing the printed model summary text
get_model_stats <- function(x, precision = 60) {
  # remember old number formatting function
  # (which would round and transforms p-values to formats like "<0.01")
  #save the current formatNP function so it can be restored later
  old_format_np <- rms::formatNP
  # substitute it with a function which will print out as many digits as we want
  #temporarily replace formatNP with one that prints full-precision decimal numbers
  assignInNamespace("formatNP", function(x, ...)
    formatC(x, format = "f", digits = precision), "rms")
  
  # remember old width setting
  #save current console width option
  old_width <- options('width')$width
  # substitute it with a setting making sure the table will not wrap
  #widen the console so the printed table doesn't line-wrap and break parsing
  options(width = old_width + 4 * precision)
  
  # actually print the data and capture it
  #print the model object and capture the text output as a character vector
  cap <- capture.output(print(x))
  
  # restore original settings
  #restore the original console width
  options(width = old_width)
  #restore the original formatNP function
  assignInNamespace("formatNP", old_format_np, "rms")
  
  # model stats
  #initialize container for extracted stats
  stats <- c()
  #extract adjusted R-squared value via regex from the captured text
  stats$R2.adj <-
    str_match(cap, "R2 adj\\s+ (\\d\\.\\d+)") %>% na.omit() %>% .[, 2] %>% as.numeric()
  
  # coef stats lines
  #find the block of text lines corresponding to the coefficient table (from header to end)
  coef_lines <-
    cap[which(str_detect(cap, "Coef\\s+S\\.E\\.")):(length(cap) - 1)]
  
  # parse
  #parse the coefficient table text block into a data frame
  coef_lines_table <-
    
    suppressWarnings(readr::read_table(coef_lines %>% stringr::str_c(collapse = "\n")))
  #rename first column (predictor names) for clarity
  colnames(coef_lines_table)[1] <- "Predictor"
  
  #return both the summary stats and the parsed coefficient table
  list(stats = stats,
       coefs = coef_lines_table)
}

#ordinal regression-based test (probit) comparing a variable v against grouping y, returning the p-value of the last coefficient
wilcox_test <- function(v, y) {
  # try to run wilcox.test  
  #wrap the model fitting in error handling since some clusters may fail (e.g., degenerate groups)
  res <- tryCatch({
    #fit an ordinal/semiparametric probit regression model of v on y
    this_fit <- rms::orm(v~y, family = "probit")
    #extract the parsed coefficient table from the fitted model
    this_res <-
      get_model_stats(this_fit) %>% .$coefs %>% data.frame()
    # format results
    #assign standardized column names to the coefficient table
    names(this_res) <-
      c('analyte',
        'std_error',
        'estimate',
        'statistic',
        'p_value',
        'other')
    # return results if test
    #return the p-value of the final row (the term of interest)
    this_res$p_value[nrow(this_res)]
    
  }, error = function(e){
    #print the error message for debugging
    print(e$message)
    # something didn't work with the wilcox.test (probably no samples in one of the groups); return NA vector
    #fallback: treat as non-significant (p = 1) if the model fit fails
    1
  })
  # return results
  return(res)
}
#classic Simes procedure for multiple testing correction
simes_correct <- function(p_values) {
  #order p-values ascending
  o = order(p_values)
  pvals = p_values[o]
  #Simes adjustment: scale each sorted p-value by (n / rank)
  for (i in 1:length(pvals)) {
    pvals[i] = length(pvals)/i * pvals[i]
  }
  #enforce monotonicity via cumulative max
  pvals = cummax(pvals)
  #restore original order
  return(pvals[order(o)])
}
#simple Bonferroni correction: multiply each p-value by the total number of tests
bonf_correct <- function(p_values) {
  return(p_values * length(p_values))
}
#compute the count (per-mille, out of 1000) of values in v that are significant at the 0.05 threshold
get_rat <- function(v) {
  round(1000 * sum(v <= 0.05)/length(v))
}

#build a sequence of "number of initially significant p-values" scenarios to evaluate correction methods at
get_vals <- function(v) {
  #baseline rate of significant p-values (per 1000)
  pv = get_rat(v)
  #combine a lower range (below baseline) and upper range (above baseline) of candidate significant counts
  vals = c(seq(5, pv, length.out = 7), seq(pv + 1, 1000 - pv, length.out = 8))
  round(vals)
}

#for each target count of "significant" p-values, build a resampled p-value set and apply all correction methods
get_corrections <- function(ps, vals, Mc) {
  lapply(vals, function(i){
    #fixed seed for reproducible resampling
    set.seed(1234)
    #find the rank position where sorted p-values first exceed 0.05 (breakpoint between "significant" and not)
    brkpoint = which(ps[order(ps)] > 0.05)[1]
    #randomly sample values from the non-significant region to fill out the rest of the 1000-test set
    idxs = sample(brkpoint:length(ps), 1000-i)
    #randomly sample i values from the significant region
    before = sample(seq(brkpoint - 1), i)
    #combine sampled significant and non-significant p-values into one resampled vector
    ps1 = ps[order(ps)][c(before, idxs)]
    
    #subset the correlation matrix to match the same resampled test indices
    M1 = Mc[order(ps)[c(before, idxs)],order(ps)[c(before, idxs)]]
    
    # calculate correlations between tests 
    xc = M1
    
    # multiple testing correction 
    #apply the extended/approximate Simes correction to the resampled p-values
    ph = f_correct(ps1, xc, use_approx = T)
    #order for plotting/comparison
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
    #build a comparison data frame of the approximate extended-Simes, plain Simes, and Bonferroni corrected p-values
    df = data.frame(approx = ph[o], simes = simes_correct(ps1), bonf = bonf_correct(ps1))
    df
  })
}
#nonlinear transform that "stretches" the visual scale so p-values below 0.05 occupy the bottom half and above occupy the top half (for plotting emphasis near the significance threshold)
stretch_pval <- function(p) {
  ifelse(p <= 0.05, p / 0.05 * 0.5, 0.5 + (p - 0.05) / (1 - 0.05) * 0.5)
}

#inverse of stretch_pval, used to relabel the stretched axis back into real p-value units
inv_stretch <- function(y) {
  ifelse(y <= 0.5, y * 0.05 / 0.5, 0.05 + (y - 0.5) * (1 - 0.05) / 0.5)
}

#build a grid of panels comparing "approx" vs "exact" extended-Simes correction across different significant-count scenarios (note: references an "exact" column that isn't actually created above, so this function may error if run as-is)
get_comb_plots2 <- function(cors, vals, label = "A") {
  plots = lapply(seq_along(cors), function(i){
    #grab the correction results for this panel
    df = cors[[i]]
    colnames(df) = c("approx", "simes", "bonf")
    #keep only the approx and simes columns for this plot
    df = df[,c(1,2)]
    #x-axis index (1 to 1000 tests)
    df$x = 1:1000
    #panel title showing how many initially significant p-values this scenario used
    df$title = paste0(vals[i], " initial significant p-values")
    
    # Sort each method’s p-values
    #sort each method's p-values ascending for a clean step/line plot
    df <- df %>%
      mutate(across(c(approx, exact), sort))
    
    # Long format for plotting
    #reshape to long format so ggplot can color lines by method
    df_long <- df %>%
      pivot_longer(cols = c(approx, exact),
                   names_to = "method", values_to = "pval")
    
    #relabel method factor levels with display-friendly names, controlling legend order
    df_long$method <- factor(df_long$method,
                             levels = c("exact", "approx"),
                             labels = c("Extended Simes (exact)",
                                        "Extended Simes (approximated)"))
    
    # Apply nonlinear transformation for p-values
    #apply the stretch transform so values near 0.05 are visually emphasized
    df_long$stretched_pval <- stretch_pval(df_long$pval)
    
    #build the line plot of stretched p-values vs rank, colored by correction method
    ggplot(df_long, aes(x = x, y = stretched_pval, color = method)) +
      geom_line(size = 1.1, alpha = 0.95) +
      scale_color_manual(
        values = c("Extended Simes (exact)" = "darkgreen",
                   "Extended Simes (approximated)" = "pink")
      ) +
      facet_grid(~ title) +
      xlim(0, 1000) +
      #custom y-axis breaks/labels that show real p-value units despite the stretched scale
      scale_y_continuous(
        name = "p-value\n",
        limits = c(0, 1),
        breaks = stretch_pval(c(0, 0.01, 0.025, 0.05, 0.25, 0.5, 1)),
        labels = signif(inv_stretch(stretch_pval(c(0, 0.01, 0.025, 0.05, 0.25, 0.5, 1))), 3)
      ) +
      #dotted reference line marking the 0.05 significance threshold
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
  #arrange all individual panels into a grid, hiding legends on each subplot and trimming y-axis labels on non-left-edge panels
  pg = plot_grid(plotlist = lapply(seq_along(plots), function(i){
    p = plots[[i]] + theme(legend.position = "none")
    
    # Remove x or y labels depending on position
    #for panels not in the first column (positions 1, 4, 7), blank out the y-axis text/ticks to avoid repetition
    if (!(i %in% c(1, 4, 7))) {
      p = p + labs(y = "\n") +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    #dead branch (i is never <= 0): intended to strip x-axis labels for non-bottom-row panels
    if (i <= 0) {
      p = p + theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
    }
    p
  }), labels = c(""), label_size = 20, nrow = 3)
  
  # Add legend at bottom
  #extract a single shared legend from the first panel and place it below the combined grid
  pg = plot_grid(
    pg,
    get_plot_component(plots[[1]] + theme(legend.position = 'bottom'), 'guide-box-bottom', return_all = TRUE),
    ncol = 1,
    rel_heights = c(1, 0.1)
  )
  
  pg
}

#build a grid of panels comparing Bonferroni, Simes, and approximate extended-Simes corrections across scenarios
get_comb_plots3 <- function(cors, vals, label = "A") {
  plots <- lapply(seq_along(cors), function(i){
    #grab correction results for this panel
    df <- cors[[i]]
    colnames(df) <- c("approx", "simes", "bonf")
    
    # Keep only approx, simes, bonf
    #keep all three correction method columns
    df <- df[, c("approx", "simes", "bonf")]
    #x-axis index (1 to 1000 tests)
    df$x <- 1:1000
    #panel title showing scenario's initial significant-count
    df$title <- paste0(vals[i], " initial significant p-values")
    
    # Sort each method's p-values
    #sort each correction method's p-values ascending
    df <- df %>%
      mutate(across(c(approx, simes, bonf), sort))
    
    # Pivot longer
    #reshape to long format for multi-line plotting
    df_long <- df %>%
      pivot_longer(cols = c("bonf", "simes", "approx"),
                   names_to = "method", values_to = "pval")
    
    # Factor and rename methods with desired legend order
    #relabel and order method factor levels for consistent legend display
    df_long$method <- factor(df_long$method,
                             levels = c("bonf", "simes", "approx"),
                             labels = c("Bonferroni",
                                        "Simes",
                                        "Extended Simes (approximated)"))
    
    # Apply nonlinear transformation
    #apply the stretch transform for visual emphasis near the 0.05 threshold
    df_long$stretched_pval <- stretch_pval(df_long$pval)
    
    #build the line plot comparing all three correction methods
    ggplot(df_long, aes(x = x, y = stretched_pval, color = method)) +
      geom_line(size = 1.1, alpha = 0.95) +
      scale_color_manual(values = c(
        "Bonferroni" = "orange",
        "Simes" = "blue",
        "Extended Simes (approximated)" = "pink"
      )) +
      facet_grid(~ title) +
      xlim(0, 1000) +
      #custom y-axis breaks/labels to display real p-value units on the stretched scale
      scale_y_continuous(
        name = "p-value\n",
        limits = c(0, 1),
        breaks = stretch_pval(c(0, 0.01, 0.025, 0.05, 0.25, 0.5, 1)),
        labels = signif(inv_stretch(stretch_pval(c(0, 0.01, 0.025, 0.05, 0.25, 0.5, 1))), 3)
      ) +
      #dotted reference line at the 0.05 significance threshold
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
  #arrange all panels into a 5-row grid, hiding legends and blanking y-axis labels except on the leftmost column (every third panel)
  pg <- plot_grid(plotlist = lapply(seq_along(plots), function(i){
    p <- plots[[i]] + theme(legend.position = "none")
    if (!(i %% 3 == 1)) {
      p <- p + labs(y = "\n") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    p
  }), labels = c(""), label_size = 20, nrow = 5)
  
  # Add legend at bottom
  #extract shared legend from the first panel and place it below the grid
  pg <- plot_grid(
    pg,
    get_plot_component(plots[[1]] + theme(legend.position = 'bottom'), 'guide-box-bottom', return_all = TRUE),
    ncol = 1,
    rel_heights = c(1, 0.1)
  )
  
  pg
}

#fix random seed for reproducible sampling/clustering downstream
set.seed(42)

#convert ROSMAP metabolomics data to a plain matrix for this analysis
x <- rosmap_data %>% as.matrix()

#randomly sample up to 250 subjects (or all, if fewer) to keep computation tractable
idxs <- sample(
  seq_len(nrow(x)),
  min(250, nrow(x))
)

#subset to the sampled subjects
x <- x[idxs, , drop = FALSE]
#standardize (z-score) each feature column
x <- x %>% scale()

#hierarchically cluster the features (transpose so features are rows) using Ward's method
hc_f <- hclust(
  dist(t(x)),
  method = "ward.D2"
)

#build the full test membership matrix M from all clusters implied by the feature tree
M <- f_M(hc_f)

#compute pairwise correlations between all test columns of M (used to model test dependence)
Mc <- cor(M, use = "pairwise.complete.obs")
#replace any missing correlations with 0 (treat as uncorrelated)
Mc[is.na(Mc)] <- 0
#take absolute value of correlations (only care about magnitude of dependence, not direction)
Mc <- abs(Mc)

#outcome variable: Braak staging score for the sampled subjects
ordinal_y <- rosmap_clins$braaksc[idxs]

#run the ordinal regression-based test for every cluster/test column in M against the Braak score outcome
ordinal_ps <- apply(
  M,
  2,
  wilcox_test,
  v = ordinal_y
)

#drop any tests that failed to produce a valid p-value
ordinal_ps2 <- ordinal_ps[!is.na(ordinal_ps)]

#compute correction comparisons (approx extended-Simes, Simes, Bonferroni) across a range of "significant count" scenarios
ordinal_corrections <- get_corrections(
  ordinal_ps2,
  get_vals(ordinal_ps2),
  Mc
)

#cap Simes and Bonferroni corrected values at 1 for every scenario's results
ordinal_corrections2 <- lapply(ordinal_corrections, function(a) {
  b <- a
  b$simes <- pmin(unlist(b$simes), 1)
  b$bonf  <- pmin(unlist(b$bonf), 1)
  b
})

#select a subset of scenario panels to actually display in the figure
panel_indices <- c(1, 3, 5, 7, 9, 11, 13, 14, 15)

#open PDF device for the two-method comparison figure (ROSMAP, Braak score outcome)
pdf(
  "supplementary_figure_7.pdf",
  width = 14,
  height = 23.5 * 3 / 5
)

#render the combined panel figure comparing approx vs exact extended-Simes (note: get_comb_plots is used here though only get_comb_plots2/3 are defined above)
get_comb_plots(
  ordinal_corrections2[panel_indices],
  get_vals(ordinal_ps2)[panel_indices],
  "C"
)

#close the PDF device, finalizing supplementary_figure_7.pdf
dev.off()

#open PDF device for the three-method comparison figure (ROSMAP, Braak score outcome)
pdf(
  "supplementary_figure_8.pdf",
  width = 14,
  height = 23.5 * 3 / 5
)

#render the combined panel figure comparing Bonferroni, Simes, and approximate extended-Simes
get_comb_plots3(
  ordinal_corrections2[panel_indices],
  get_vals(ordinal_ps2)[panel_indices],
  "C"
)

#note: no dev.off() call here before the next pdf() call below, so this figure may not be properly finalized/closed

#reset seed for reproducibility of the second (ADNI) analysis
set.seed(42)

#convert ADNI baseline lipidomics data to a plain matrix
x <- adni_bl_lipids %>% as.matrix()

#randomly sample up to 250 subjects (or all, if fewer)
idxs <- sample(
  seq_len(nrow(x)),
  min(250, nrow(x))
)

#subset to sampled subjects
x <- x[idxs, , drop = FALSE]
#standardize each lipid feature column
x <- x %>% scale()

#hierarchically cluster the lipid features using Ward's method
hc_f <- hclust(
  dist(t(x)),
  method = "ward.D2"
)

#build the full test membership matrix M from all clusters implied by the lipid feature tree
M <- f_M(hc_f)

#compute pairwise correlations between all test columns of M
Mc <- cor(M, use = "pairwise.complete.obs")
#replace missing correlations with 0
Mc[is.na(Mc)] <- 0
#take absolute value of correlations
Mc <- abs(Mc)

#outcome variable: BMI for the sampled ADNI subjects
ordinal_y <- adni_clinsl0$bmi[idxs]

#run the ordinal regression-based test for every cluster/test column in M against BMI
ordinal_ps <- apply(
  M,
  2,
  wilcox_test,
  v = ordinal_y
)

#drop any tests that failed to produce a valid p-value
ordinal_ps2 <- ordinal_ps[!is.na(ordinal_ps)]

#compute correction comparisons across a range of "significant count" scenarios for the ADNI/BMI analysis
ordinal_corrections <- get_corrections(
  ordinal_ps2,
  get_vals(ordinal_ps2),
  Mc
)

#cap Simes and Bonferroni corrected values at 1
ordinal_corrections2 <- lapply(ordinal_corrections, function(a) {
  b <- a
  b$simes <- pmin(unlist(b$simes), 1)
  b$bonf  <- pmin(unlist(b$bonf), 1)
  b
})

#select the same subset of scenario panels to display
panel_indices <- c(1, 3, 5, 7, 9, 11, 13, 14, 15)

#open PDF device for the ADNI/BMI comparison figure
pdf(
  "supplementary_figure_9.pdf",
  width = 14,
  height = 23.5 * 3 / 5
)

#render the combined panel figure comparing correction methods for the ADNI/BMI analysis
get_comb_plots(
  ordinal_corrections2[panel_indices],
  get_vals(ordinal_ps2)[panel_indices],
  "C"
)

dev.off()
