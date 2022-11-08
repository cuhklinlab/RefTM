#' @export

RefTM.LDA <- function(sc_data, k1 = 5, k2 = NULL, workflow = "LDA", covariate = NULL, refBeta){
  docvoc2dtm = function(doc_voc){
    if(length(which(colSums(doc_voc) == 0))>0){
      doc_voc = doc_voc[,-which(colSums(doc_voc) == 0)]
    }
    library(slam)
    dtm = as.simple_triplet_matrix(doc_voc)
    dimnames(dtm) = list(Docs = 1:dtm$nrow, Terms = 1:dtm$ncol)
    return(dtm)
  }
  if(!is.null(k2)){
      sprintf("Number of topics is fixed to be %s.", k1+k2)
      model_sc = lda_svi(docvoc2dtm(t(sc_data)),batchsize = 100,K = k1+k2,K0 = k1,maxiter = 1000,refBeta = refBeta,passes = 2)
  }else{
    model_sc = lda_svi(docvoc2dtm(t(sc_data)),batchsize = 100,K = k1,K0 = k1,maxiter = 1000,refBeta = refBeta,passes = 2)
  }
  return(model_sc = model_sc)
}



