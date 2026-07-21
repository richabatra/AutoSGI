library(rms)
library(stringr)
library(readr)
#get rosmap data
rosmap_data = read_excel("tmp_rosmap_brain_metabolomics_processed_medcor_data_pmicor.xlsx", sheet = "assay")
#get clinical data
rosmap_clins = read_excel("tmp_rosmap_brain_metabolomics_processed_medcor_data_pmicor.xlsx", sheet = "colData")
#get pathways
rosmap_paths = read_excel("tmp_rosmap_brain_metabolomics_processed_medcor_data_pmicor.xlsx", sheet = "rowData")
#transpose rosmap data, correct columns
rosmap_data <- t(rosmap_data)
#use first row (originally header row of assay sheet) as new column names after transpose
colnames(rosmap_data) <- rosmap_data[1, ]
#drop that first row now that it's been used as the header
rosmap_data <- rosmap_data[-1, ]
#define the clinical variables of interest to keep from colData
clin_vars = c("age_death", "msex", "educ", "pmi", "bmi", "anye4", "cogdx", "braaksc", "ceradsc", "sqrt_amyloid", "tangles", "gpath", "niareagansc", "diagnosis", "cogn_global", "cogng_random_slope")
#subset clinical data to only the selected variables
rosmap_clins <- rosmap_clins[, which(colnames(rosmap_clins) %in% clin_vars)]
#coerce to plain data frame
rosmap_clins <- as.data.frame(rosmap_clins)
#convert all assay values to numeric (they were character after transpose)
rosmap_data <- data.frame(apply(rosmap_data, 2, function(x) as.numeric(as.character(x))))
#pull out super pathway annotations for each metabolite
rosmap_super_pathway = rosmap_paths[, 4:4]$SUPER_PATHWAY 
#pull out sub pathway annotations for each metabolite
rosmap_sub_pathway = rosmap_paths[, 5:5]$SUB_PATHWAY 
#align row names of assay data with clinical data row names (subject IDs)
rownames(rosmap_data) <- rownames(rosmap_clins)
#set classes
#APOE e4 carrier status as factor
rosmap_clins$anye4 = as.factor(rosmap_clins$anye4)
#sex as factor
rosmap_clins$msex = as.factor(rosmap_clins$msex)
#NIA-Reagan score as factor
rosmap_clins$niareagansc = as.factor(rosmap_clins$niareagansc)
#diagnosis as factor
rosmap_clins$diagnosis = as.factor(rosmap_clins$diagnosis)
#reverse-code CERAD score so higher values mean more severe pathology (consistent direction with other measures)
rosmap_clins$ceradsc = 5 - rosmap_clins$ceradsc
#mark cognitive diagnosis as ordinal for modeling purposes
class(rosmap_clins$cogdx) = "ordinal"
#mark Braak staging as ordinal
class(rosmap_clins$braaksc) = "ordinal"
#mark (reverse-coded) CERAD score as ordinal
class(rosmap_clins$ceradsc) = "ordinal"
#keep an unmodified copy of the original clinical data for reference
orig_clins = rosmap_clins
library(rms)
#parser from stackoverflow by CoderGuy123
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
#ordinal regression-based test (probit) comparing a variable v against grouping y, returning both p-value and statistic of the last coefficient
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
    #return both the p-value and test statistic of the final row (the term of interest)
    list(pval = this_res$p_value[nrow(this_res)], stat = this_res$statistic[nrow(this_res)])

  }, error = function(e){
    #print the error message for debugging
    print(e$message)
    # something didn't work with the wilcox.test (probably no samples in one of the groups); return NA vector
    #fallback: return NA for both p-value and statistic if the model fit fails
    list(pval = NA, stat = NA)
  })

  # return results
  return(res)
}
#subset clinical data down to just the three ordinal outcome variables used for association testing
rosmap_clins2 = rosmap_clins[,which(colnames(rosmap_clins)%in%c("cogdx", "braaksc", "ceradsc"))]
#initialize sgi parameters: scaled metabolomics data, the ordinal outcomes, a minimum cluster size (5% of samples), and the custom ordinal test function
rosmap_params = sgi_params_init(scale(rosmap_data), rosmap_clins2, minsize = nrow(rosmap_data)/20, user_defined_tests = c(ordinal = wilcox_test))
#run the sgi hierarchical selection procedure, testing clusters against outcomes with Simes correction, generating plots and exporting result tables
rosmap_metab_tree_results = hierarchical_selection(rule = rule_init("ALL", 2), rosmap_params, cluster_min = 2, correction_opt = "simes", plot = T, supp_plot = T, summary_plot = T, output_names = list(sgi_plots = "rosmap-hierarchical-selection-results.pdf",
                                                                                                                                                                                                        summary = "rosmap-hierarchical-selection-summary.pdf",                                                                                                                                                                       cluster_results = "rosmap-hierarchical-selection-labels.xlsx",
                                                                                                                                                                                                        cluster_indices = "rosmap-hierarchical-selection-cluster-indices.xlsx",
                                                                                                                                                                                                       sgi_as_results = "rosmap-hierarchical-selection-association-results.xlsx")) 
library(fpc)
#bootstrap cluster stability (Jaccard similarity) for the subtree-619 metabolite cluster, split into k=2 groups
c619stability = clusterboot(dist(scale(rosmap_data[, which(colnames(rosmap_data) %in% rosmap_metab_tree_results$subtrees[[619]])])), distance = TRUE, B=500, clustermethod = disthclustCBI, method = "ward.D2", bootmethod ="boot", cut = "number", k = 2)
#bootstrap cluster stability for the subtree-224 metabolite cluster, split into k=2 groups
c358stability_k2 = clusterboot(dist(scale(rosmap_data[, which(colnames(rosmap_data) %in% rosmap_metab_tree_results$subtrees[[224]])])), distance = TRUE, B=500, clustermethod = disthclustCBI, method = "ward.D2", bootmethod ="boot", cut = "number", k = 2)
#bootstrap cluster stability for the same subtree-224 cluster, but split into k=3 groups instead
c358stability_k3 = clusterboot(dist(scale(rosmap_data[, which(colnames(rosmap_data) %in% rosmap_metab_tree_results$subtrees[[224]])])), distance = TRUE, B=500, clustermethod = disthclustCBI, method = "ward.D2", bootmethod ="boot", cut = "number", k = 3)
#fixed seed for reproducibility of the plotting/labeling steps below
set.seed(1234)
#set transparent background and outer margins for the histogram panels
par(xpd = NA, 
    bg = "transparent", 
    oma = c(2, 2, 0, 0)) 
#select the first bootstrap group (subgroup) from the k=2 clustering
j = 1
#plot histogram of Jaccard similarity values across bootstrap replicates for this subgroup
hist(c358stability_k2$bootresult[j,],xlim=c(0, 1),breaks=seq(0, 1, by = 0.05), 
    xlab="Jaccard similarity", 
    main=paste("Subgroup", 2), col = 'grey', xaxt='n')
#add a dashed red vertical line at the mean Jaccard similarity
abline(v = round(mean(c358stability_k2$bootresult[j,]), 2), col = 'red', lwd = 2, lty = 'dashed')
#default candidate x-axis tick positions
xt = c(0, 0.2, 0.4, 0.6, 0.8, 1)
#remove any default ticks that would sit too close to the mean line (to avoid label overlap)
xt = xt[-c(which(abs(xt - mean(c358stability_k2$bootresult[j,])) <= 0.07))]
#safety fallback: if all ticks got removed, restore the default set
if (is.null(xt) | length(xt) == 0) {
  xt <<- c(0, 0.2, 0.4, 0.6, 0.8, 1)
}
#add the mean value itself as an additional tick mark
xt = c(xt, round(mean(c358stability_k2$bootresult[j,]), 2))
#debug print of the final tick positions
print(xt)
#draw the custom x-axis with sorted tick positions/labels
axis(side=1,at=sort(xt),labels=sort(xt))
#save this histogram plot for later combination into the figure grid
s2 = recordPlot()
#reset plotting parameters for the next histogram panel
par(xpd = NA, 
    bg = "transparent", 
    oma = c(2, 2, 0, 0)) 
#select the second bootstrap group from the k=2 clustering
j = 2
#plot histogram of Jaccard similarity values for this subgroup
hist(c358stability_k2$bootresult[j,],xlim=c(0, 1),breaks=seq(0, 1, by = 0.05), 
    xlab="Jaccard similarity", 
    main=paste("Subgroup", 3), col = 'grey', xaxt='n')
#add dashed red line at the mean
abline(v = round(mean(c358stability_k2$bootresult[j,]), 2), col = 'red', lwd = 2, lty = 'dashed')
#default candidate tick positions
xt = c(0, 0.2, 0.4, 0.6, 0.8, 1)
#drop ticks too close to the mean line
xt = xt[-c(which(abs(xt - mean(c358stability_k2$bootresult[j,])) <= 0.07))]
#fallback if all ticks were removed
if (is.null(xt) | length(xt) == 0) {
  xt <<- c(0, 0.2, 0.4, 0.6, 0.8, 1)
}
#add mean value as an extra tick
xt = c(xt, round(mean(c358stability_k2$bootresult[j,]), 2))
#draw custom x-axis
axis(side=1,at=sort(xt),labels=sort(xt))
#save this histogram plot
s3 = recordPlot()
#reset plotting parameters for the next histogram panel
par(xpd = NA, 
    bg = "transparent", 
    oma = c(2, 2, 0, 0)) 
#select the second bootstrap group from the k=3 clustering
j = 2
#plot histogram of Jaccard similarity values for this subgroup
hist(c358stability_k3$bootresult[j,],xlim=c(0, 1),breaks=seq(0, 1, by = 0.05), 
    xlab="Jaccard similarity", 
    main=paste("Subgroup", 4), col = 'grey', xaxt='n')
#add dashed red line at the mean
abline(v = round(mean(c358stability_k3$bootresult[j,]), 2), col = 'red', lwd = 2, lty = 'dashed')
#default candidate tick positions
xt = c(0, 0.2, 0.4, 0.6, 0.8, 1)
#drop ticks too close to the mean line
xt = xt[-c(which(abs(xt - mean(c358stability_k3$bootresult[j,])) <= 0.07))]
#fallback if all ticks were removed
if (is.null(xt) | length(xt) == 0) {
  xt <<- c(0, 0.2, 0.4, 0.6, 0.8, 1)
}
#add mean value as an extra tick
xt = c(xt, round(mean(c358stability_k3$bootresult[j,]), 2))
#draw custom x-axis
axis(side=1,at=sort(xt),labels=sort(xt))
#save this histogram plot
s4 = recordPlot()
#reset plotting parameters for the final histogram panel
par(xpd = NA, 
    bg = "transparent", 
    oma = c(2, 2, 0, 0)) 
#select the third bootstrap group from the k=3 clustering
j = 3
#plot histogram of Jaccard similarity values for this subgroup
hist(c358stability_k3$bootresult[j,],xlim=c(0, 1),breaks=seq(0, 1, by = 0.05), 
    xlab="Jaccard similarity", 
    main=paste("Subgroup", 5), col = 'grey', xaxt='n')
#add dashed red line at the mean
abline(v = round(mean(c358stability_k3$bootresult[j,]), 2), col = 'red', lwd = 2, lty = 'dashed')
#default candidate tick positions
xt = c(0, 0.2, 0.4, 0.6, 0.8, 1)
#drop ticks too close to the mean line
xt = xt[-c(which(abs(xt - mean(c358stability_k3$bootresult[j,])) <= 0.07))]
#fallback if all ticks were removed
if (is.null(xt) | length(xt) == 0) {
  xt <<- c(0, 0.2, 0.4, 0.6, 0.8, 1)
}
#add mean value as an extra tick
xt = c(xt, round(mean(c358stability_k3$bootresult[j,]), 2))
#draw custom x-axis
axis(side=1,at=sort(xt),labels=sort(xt))
#save this final histogram plot
s5 = recordPlot()
#open PDF device for the combined 4-panel stability figure
pdf("supp_figure_9.pdf", height = 11.2, width = 13)
#arrange all four saved histogram plots into a 2x2 grid with panel labels a-d
plot_grid(plotlist = lapply(list(s2, s3, s4, s5), ggdraw), ncol = 2, labels = c("a", "b", "c", "d"))

dev.off()
