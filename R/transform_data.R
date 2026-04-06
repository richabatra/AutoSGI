#' Transform data
#'
#' This is a function which performs interpretable dimensionality reduction on the data given by ssgsea or pseudotime
#'
#' @param dataset A dataset with metabolite columns
#' @param annotations A dataframe with the first column consisting on the metabolite names and the second column consisting of the corresponding annotations
#' @param method The method to perform reduction, ssgsea or pseudotime
#' @return Always returns an updated dataframe with one column for each annotation
#' @export
#' @examples
#' \dontrun{
#' }
#'
#'
transform_data <- function(dataset, annotations, method = "ssgsea") {
  
  if (ncol(dataset) != nrow(annotations)) {
    stop("Number of rows in annotations should equal number of columns in dataset")
  }
  if (method != "ssgsea" & method != "pseudotime") {
    stop("Available methods for dimensionality reduction are ssgsea and pseudotime only")
  }
  if (ncol(annotations) != 2) {
    stop("There should be 2 columns in annotations: the first column containing metabolite names and the second containing metabolite annotations")
  }
  if (!identical(sort(annotations[[1]]),sort(colnames(dataset)))) {
    stop("Dataset column names should be identical to first column in annotations")
  }
  
  classes = unique(annotations[[2]])
  filtered_classes = unique(classes[sapply(classes, function(g){return(length(which(annotations[[2]] == g)) > 1)})])
  if (length(classes) != length(filtered_classes)) {
    warning("Removing annotations which appear only once")
  }
  colnames(annotations) = c("features", "annotation")
  
  if (method == "ssgsea") {
    pathway_enrichment <- function(dt, fa, a) {
      gs = lapply(fa, function(g){
        return(rownames(dt)[which(a==g)])
      })
      param = ssgseaParam(as.matrix(dt), gs)
      gsva(param)
    }
    
    transformed_dataset <- dataset %>% t() %>% pathway_enrichment(.,filtered_classes, annotations[[2]]) %>% as.data.frame() %>% t()
    colnames(transformed_dataset) = filtered_classes
  }
  else {
    dpt_as_ssgsea <- function(eset, grp_list, dist_method='euclidean'){
      tmp<- lapply(grp_list, FUN=function(grp){
        this_eset <- eset[which(eset@featureData@data$annotation == grp), ]
        # diffusion map
        dm_obj <- DiffusionMap(this_eset, distance = dist_method, knn_params = list(method = 'covertree'))
        # pseudotimeå
        dpt_obj <- dm_obj %>% DPT() %>% plot.DPT()
        # extract time
        dpt <- dpt_obj$data$Colour
      }) 
      tmp <- do.call(rbind, tmp)
      return(tmp)
    }
    
    dpt_from_dt <- function(data, classes) {
      assay_data = data %>% t()
      #colData is a placeholder here
      se = SummarizedExperiment(assay=as.matrix(assay_data), colData=as.matrix(data), rowData=as.data.frame(annotations))
      eset = ExpressionSet(assayData=(se %>% assay() %>% as.matrix()),
                           phenoData=(se %>% colData() %>% data.frame() %>% AnnotatedDataFrame()),
                           featureData=(se %>% rowData() %>% data.frame() %>% AnnotatedDataFrame()))
      dpt_as_ssgsea(eset, classes) %>% as.data.frame()
    }
    
    transformed_dataset <- dataset %>% dpt_from_dt(.,filtered_classes) %>% t()
    return(transformed_dataset)
  }

  transformed_dataset
}