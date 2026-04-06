#reformatted from sgi-associations
#extracts relevant subgroup pairs
extract <- function(r, as, p_value) {
  r = cbind(cluster_pair = rownames(as$cluster_pairs), r)
  r = r[padj <= p_value, ] #add p_value constraint
  if (nrow(r) > 0) {
    vcp = as$cluster_pairs[r$cluster_pair, ]
    attr(r, "hvcp") = as.matrix(vcp)
  }
  class(r) = c("assoc_table", class(r))
  r
}
#split a table across several pages
#given some annotations and results, as well as max columns on a page
#returns ceiling(#annotations/max_columns) separate dataframes
split_data <- function(annotations, results, max_per_page, name = "pathways") {
  
  cnt = 0
  curr_annotations <- vector()
  curr_results <- vector()
  
  for (i in 1:length(annotations)) {
    cnt = cnt + 1
    curr_annotations <- append(curr_annotations, annotations[i])
    curr_results <- append(curr_results, results[i])
    if (i == length(annotations) | cnt == max_per_page) {
      data_plot <- data.frame(name = curr_annotations,
                                 valid = curr_results)
      names(data_plot)[1] = name
      table <- tableGrob(data_plot, rows = NULL)
      grid.newpage()
      grid.draw(table)
      curr_annotations <- c()
      curr_results <- c()
      cnt = 0
    }
  }
}

#returns an empty page with the message in the middle
plot_message <- function(message) {
  obj = ggplot() + xlim(0, 10) + ylim(0, 10) + theme_void() + 
    annotate("text", x = 5, y = 5, label = message, size = 4, colour = "black")
  return(obj)
}

#returns an invalid plot with a specific error message
invalid_plot <- function(message) {
  obj = ggplot() + xlim(0, 10) + ylim(0, 10) + theme_void() + 
    annotate("text", x = 5, y = 7, label = "No Results\n", size = 6, colour = "red") + 
    annotate("text", x = 5, y = 6.5, label = paste0("Error message: ", message, "\n"), size = 3, colour = "black")
  return(obj)
}

#supplemental plot displaying all the features of an SGI cluster
display_features <- function(features, tree = NULL, maximum = 10) {
  #tree based selection
  if (!is.null(tree)) {
    if (length(features) <= maximum) {
      #add labels of all metabolites to right of tree_point
      cluster_labs = data.frame(biomolecules = features)
      names_table = tableGrob(cluster_labs, rows = NULL)
      return(grid.arrange(arrangeGrob(tree, names_table, nrow = 2)))
    }
    else {
      empty = ggplot() + theme_void()
      #basically just keeping plot centered in middle without scaling height up
      grid.arrange(arrangeGrob(empty, tree, empty, nrow = 3, heights = c(7/4, 7/2, 7/4))) #default width
    }
  }
  else { 
    if (length(features) > maximum) {
      return(NULL)
    }
    #set selection, add features to page with lower margins in between
    empty = ggplot() + theme_void()
    #basically just keeping plot centered in middle without scaling height up
    cluster_labs = data.frame(biomolecules = features)
    names_table = tableGrob(cluster_labs, rows = NULL)
    grid.arrange(arrangeGrob(empty, names_table, empty, nrow = 3, heights = c(7/8, 5.25, 7/8))) #default width
  }
}

#creates excel sheet given the as objects
combined_as_sheet <- function(asx) {
  lapply(asx, function(as){
    as_info = do.call(cbind, lapply(seq(length(as$results)), function(i) {
      d = as.data.frame(as$results[[i]])
      rownames(d) = paste0(d$cid1, "vs", d$cid2)
      d = d[c(3)]
      colnames(d) = names(as$results)[i]
      d 
    }))
    as_info = cbind(as.data.frame(rownames(as_info)), as_info)
    colnames(as_info)[1] = "CID"
    as_info
  })
}

#runs sgi plotting
#returns a plot containing an error message OR a regular sgi plot
#and true/false depending on if the cluster is invalid
safe_run_plot <- function(as, p_value, clins, metabolite_cluster, annotation) {
  run <- function() {
    tryCatch({
      gg_tree = plot(as, padj_th = p_value)
      if (!is.null(dim(clins[[1]]))) {
        #custom outcome -- if followed like sgi tutorial this will work
        cnms = colnames(clins)
        cls = sapply(seq(ncol(clins)), function(c){return(class(clins[[c]][1])[1])})
        clins <- as.data.frame(lapply(seq(ncol(clins)), function(c){unname(unlist(clins[[c]][1]))}))
        colnames(clins) <- cnms
        lapply(seq(ncol(clins)), function(c){
          if (cls[c] == "binary" | cls[c] == "factor") {
            clins[[c]] <<- as.factor(clins[[c]])
          }
          else {
            class(clins[[c]]) <<- cls[c]
          }
        })
      }
      obj = plot_overview(gg_tree = gg_tree, as = as, outcomes = clins, xdata = metabolite_cluster, data_title = annotation) #include custom annotation for data title
      return(list(result = obj, invalid = FALSE)) #cluster is valid
    },
    error = function(e) {
      print(e$message)
      if (e$message == "no significant clusters for given p-value threshold!\n         
          still want to plot the tree, try plotting sgi.object created by sgi::sgi()!") {
        e$message = paste0(e$message, "\n(Likely too many x or y observations, or data is essentially constant).")
      }
      obj = invalid_plot(e$message)
      return(list(result = obj, invalid = TRUE)) #cluster is invalid
    })
  }
  ret <- run()
  return(ret)
}

#runs sgi clustering + initialization
#returns whether the cluster is invalid (true) or false, and a corresponding plot/sg object
safe_run_sgi <- function(metabolite_cluster, minsize, user_defined_tests = c(), sgi_distance, sgi_linkage, clins) {
  run <- function() {
    tryCatch({
      hc = hclust(dist(metabolite_cluster, method = sgi_distance), method = sgi_linkage)
      sg = sgi_init(hc, minsize, outcomes = clins, user_defined_tests = user_defined_tests)
      return(list(result = sg, invalid = FALSE)) #this is probably bad practice (fix todo)
    },
    error = function(e) {
      print(e$message)
      obj = invalid_plot(e$message)
      return(list(result = obj, invalid = TRUE))
    })
  }
  ret <- run()
  return(ret)
}
av <- function(x, rm = F) {
  if (rm) {return(unlist(x)[!is.na(unlist(x))])}
  else {return(unlist(x))}
}
#given sgi objects, association objects, and outcomes correct the padj in asx
correct_as <- function(sgs, asx, clin, correction_opt) {
  
  # function to collect and return all tests given sgs
  f_M <- function(sgs) {
    ms = lapply(sgs, function(sg){sgi::get_vcps(sg)})
    # one giant M to collect all tests
    M = do.call(cbind, ms)
    apply(M, 2, function(a) a-min(a,na.rm = T))
  }
  
  # Simes procedure extended with effective number of test calculation
  # multiple testing correction in multiple SGI runs 
  f_correct <- function(p_values, p_cor_matrix, k_min = 50, k_length = 40, use_approx = FALSE) {
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
    
    k_min = min(k_min, length(p_values))
    k_length = min(k_length, length(p_values) - k_min)
    
    k <- length(p_values)
    
    if(use_approx && k_min < k){
      # approximate eff n by k 
      ks <- c(seq(k_min), unique(round(seq(k_min + 1, k, length.out = k_length))))
      # use this exact values to approximate for ks
      exact_eff_ns <- sapply(ks, function(i) get_effective_n(p_cor_matrix[1:i, 1:i]))
      eff_ns <- approx(x = ks, y = exact_eff_ns, xout = seq(k))$y
    }else {
      # exact computation 
      eff_ns <- sapply(seq(k), function(i) {
        get_effective_n(p_cor_matrix[1:i, 1:i])
      })
    }
    # enforce monotocity
    eff_ns <- cummax(eff_ns)
    
    global_eff_n <- eff_ns[k]
    corrected_p_values <- pmin((global_eff_n * p_values) / eff_ns, 1)
    cummax(corrected_p_values)[order(o)]
  }
  
  #get all the tests
  M = f_M(sgs)
  
  # calculate correlations between tests 
  xc = cor(apply(M,2,function(x){x[is.na(x)] = mean(x,na.rm = T);x}))
  xc = abs(xc)
  
  #for each outcome we need to correct so we get pvalues for each as object
  lvcps = sapply(sgs, function(sg){ncol(sgi::get_vcps(sg))})
  pv <- do.call(cbind, lapply(colnames(clin), function(cx){
    tsts <- do.call("c", lapply(asx, function(as){
      as.data.frame(as$results[which(colnames(clin) == cx)])[, 4]
    }))
    tsts <- setNames(tsts, do.call("c", lapply(seq(length(sgs)), function(i){
      rep(i, times = lvcps[i])
    })))
  }))
  
  if (correction_opt == "simes") {
    pvd = as.data.frame(pv)
    simes_corrected <- as.data.frame(do.call(cbind, lapply(pvd, function(col){
      v = as.vector(col)
      na_idxs = which(is.na(v))
      v[na_idxs] = 1
      rord = order(v)
      corrections = sum(lvcps)/rank(v, ties = "max")
      v = v * corrections
      v = v[rord] %>% cummax() %>% pmin(.,1)
      v = v[order(rord)]
      v[na_idxs] = NA
      v
    })))
  }
  colnames(pv) = colnames(clin)
  
  #now we correct for each outcome
  lapply(colnames(clin), function(cx){
    ps = pv[, which(colnames(pv) == cx)]
    #make all NA equal 1
    ps[is.na(ps)] = 1
    # multiple testing correction (bonferroni or simes)
    if (correction_opt == "bonferroni") {
      sapply(seq(length(sgs)), function(i){
        asx[[i]]$results[[cx]][, 3]$padj <<- sapply(asx[[i]]$results[[cx]][, 4]$pval, function(x){min(1, x * sum(lvcps))}) #multiply by total sgs
      })
    }
    if (correction_opt == "simes") {
      sc_clin = simes_corrected[, which(colnames(pv) == cx)]
      sapply(seq(length(sgs)), function(i){
        asx[[i]]$results[[cx]][, 3]$padj <<- sc_clin[which(rownames(pv) == i)]
      })
    }
    if (correction_opt == "extended_simes") {
      #only need exact for good candidates (we know there are at least lvcps tests since on same tree)
      ph = f_correct(ps, xc, k_min = min(50, sum(ps < 0.05/min(lvcps))), k_length = 5, use_approx = T)
      sapply(seq(length(sgs)), function(i){
        #we dont divide by lvcps since we use pvalue
        asx[[i]]$results[[cx]][, 3]$padj <<- ph[which(rownames(pv) == i)]
      })
    }
  })
  
  #return updated as objects
  return(asx)
}
#given a list of clusters, returns all appropriate data need for tree based selection and network based selection
#processes clusters for unsupervised subsetting
process_clusters <- function(clusters, rule, dataset, clins, minsize, user_defined_tests, sgi_distance, sgi_linkage, plot, correction_opt) {
  
  p_value = rule$p_value
  cluster_indices = vector(mode = "list", length = length(subtrees))
  sgi_plots = list()
  sgs <- lapply(seq(length(clusters)), function(cluster_id) {
    
    if (cluster_id %% 50 == 0) {
      cat(cluster_id)
      cat(" SGI objects out of ")
      cat(length(clusters))
      cat(" processed")
      cat("\n")
    }
    
    metabolite_indices = which(colnames(dataset) %in% clusters[[cluster_id]])
    metabolite_cluster <- as.data.frame(dataset[, metabolite_indices])
    
    #run sgi on cluster
    ret_sgi = safe_run_sgi(metabolite_cluster, minsize, user_defined_tests = user_defined_tests, sgi_distance, sgi_linkage, clins)
    ret_sgi$result
  })
  asx <- lapply(sgs, function(sg) {
    sgi_run(sg)
  })
  asx <- correct_as(sgs, asx, clins, correction_opt)
  res <- lapply(seq(length(clusters)), function(cluster_id){
    
    if (cluster_id %% 50 == 0) {
      cat(cluster_id)
      cat(" SGI objects out of ")
      cat(length(clusters))
      cat(" plotted")
      cat("\n")
    }
    
    metabolite_indices = which(colnames(dataset) %in% clusters[[cluster_id]])
    metabolite_cluster <- as.data.frame(dataset[, metabolite_indices])
    
    
    cid <- NA
    val_tree <- inv_clust <- F
    
    
    if (!inv_clust && tree_validation(asx[[cluster_id]], rule)) {
      #tree is valid
      val_tree <- T
    }
    
    if (!inv_clust & plot & val_tree) { #only add interesting trees to pdf
      
      #title is annotation
      #if annotation is too big just use cluster index (might need to fix on sgi specifically)
      annotation = paste("clst", as.character(cluster_id))
      
      #run all methods needed for plotting
      #use run_plot
      ret_plot = safe_run_plot(asx[[cluster_id]], p_value, clins, metabolite_cluster, annotation)
      if (ret_plot$invalid) {
        if (annotation == "")
        inv_clust <- T
      }
      else {
        sgi_plots[[length(sgi_plots) + 1]] <<- ret_plot$result
        cid <- cluster_id
      }

    }
    
    #if cluster is invalid, tree is also not of interest
    if (inv_clust) {
      label <- "Invalid"
      val_tree <- F
    }
    else {
      if (val_tree) {
        label <- "True"
      }
      else {
        label <- "False"
      }
    }
    
    #pad metabolite_indices with NA
    metabolite_indices = c(metabolite_indices, rep(0, ncol(dataset) - length(metabolite_indices)))
    cluster_indices[[cluster_id]] <<- metabolite_indices
    return(list(results = val_tree, cluster_labels = label, cluster_id = cid))
  })
  res <- as.data.frame(do.call(rbind, res))
  
  return(list(results = av(res$results), cluster_labels = av(res$cluster_labels), cluster_indices = cluster_indices, 
             cluster_ids = av(res$cluster_id, T), sgi_plots = sgi_plots, asx = asx))
}
