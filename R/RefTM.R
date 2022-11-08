RefTM <- function(sc_data, ref_data, k1 = 5, k2 = NULL, workflow = "LDA", covariate = NULL){
  docvoc2dtm = function(doc_voc){
    if(length(which(colSums(doc_voc) == 0))>0){
      doc_voc = doc_voc[,-which(colSums(doc_voc) == 0)]
    }
    library(slam)
    dtm = as.simple_triplet_matrix(doc_voc)
    dimnames(dtm) = list(Docs = 1:dtm$nrow, Terms = 1:dtm$ncol)
    return(dtm)
  }
  perplexity = function(X, topic_word_distribution, doc_topic_distribution) {
    EPS = 1e-16
    p = X@p
    j = X@i
    x = X@x
    ll = 0
    for(i in 1:nrow(X)) {
      p1 = p[[i]]
      p2 = p[[i + 1L]]
      pointer = p1 + seq_len(p2 - p1)
      word_indices = j[pointer] + 1L
      word_counds = x[pointer]
      dot_prod = doc_topic_distribution[i, , drop = FALSE] %*%
        topic_word_distribution[ , word_indices, drop = FALSE]
      ll = ll +  log(dot_prod + EPS) %*% word_counds
    }
    # drop dimensions
    ll = as.numeric(ll)
    exp(-ll / sum(X@x))
  }
##RefTM-LDA workflow
  if(workflow == "LDA"){
    library(topicmodels)
    model_ref <- LDA(ref_data, k = k1, control = list(em = list(iter.max = 1000,tol = 10^-5)))
    print("Finished modeling the reference data.")
    if(!is.null(k2)){
      print(k2)
      sprintf("Number of topics is fixed to be %s.", k1+k2)
      if(length(which(rowSums(sc_data) == 0)) != 0){
        model_sc = lda_svi(docvoc2dtm(t(sc_data)),batchsize = 100,K = k1+k2,K0 = k1,maxiter = 1000,refBeta = exp(model_ref@beta[,-which(rowSums(sc_data) == 0)]),passes = 2)
      }else{
        model_sc = lda_svi(docvoc2dtm(t(sc_data)),batchsize = 100,K = k1+k2,K0 = k1,maxiter = 1000,refBeta = exp(model_ref@beta),passes = 2)
      }
    }else{
      print("Select the best model.")
      k = k1 + c(0, 2, 5, 10, 15, 20)
      models <- list()
      perplex <- c()
      for (i in 1:6) {
        if(length(which(rowSums(sc_data) == 0)) != 0){
          models[[i]] = lda_svi(docvoc2dtm(t(sc_data)),batchsize = 100,K = k[i],K0 = k1,maxiter = 1000,refBeta = exp(model_ref@beta[,-which(rowSums(sc_data) == 0)]),passes = 2)
        }else{
          models[[i]] = lda_svi(docvoc2dtm(t(sc_data)),batchsize = 100,K = k[i],K0 = k1,maxiter = 1000,refBeta = exp(model_ref@beta),passes = 2)
        }
        perplex[i] <- perplexity(as(t(sc_data),"sparseMatrix"),models[[i]]$beta,models[[i]]$theta)

      }
      log.lik <- data.frame(topics = k, LL = perplex)
      log.lik$first_derivative <- c(-Inf, (diff(log.lik$LL) / diff(log.lik$topics)))
      log.lik$second_derivative <- c(-Inf, -Inf, diff(log.lik$first_derivative)[-1]/diff(log.lik$topics[-1]))
      par(bty = 'n')
      plot(log.lik$topics, log.lik$LL, xlab="Number of topics", ylab="perplexity", type='o', pch=16, col='black', main='Model selection')
      points(log.lik$topics[which(log.lik$second_derivative == max(log.lik$second_derivative))], log.lik$LL[which(log.lik$second_derivative == max(log.lik$second_derivative))], pch=4, col='red', lwd = 7)
      model_sc <- models[[which(log.lik$second_derivative == max(log.lik$second_derivative))]]
    }
  }else

##RefTM-STM workflow
  if(workflow == "STM"){
    if(is.null(covariate)){
      stop("Please use LDA workflow.")
    }
    ref_data  = ref_data[,which(colSums(ref_data) != 0)]
    sc_data = sc_data[which(colSums(ref_data) != 0),]
    getCorpus = function(mat){
      corpus = list()
      for (i in 1:dim(mat)[1]) {
        term = which(mat[i,] != 0)
        corpus[[i]] = rbind(term,mat[i,term])
      }
      return(corpus)
    }
    model_ref <- RefTM.STM(getCorpus(ref_data), as.character(1:dim(ref_data)[2]), K = k1)
    print("Finished modeling the reference data.")
    if(!is.null(k2)){
      sprintf("Number of topics is fixed to be %s.", k1+k2)
      model_sc <- RefTM.STM(getCorpus(t(sc_data)), as.character(1:dim(sc_data)[1]), K = k1+k2, prevalence = ~covariate, refBeta = exp(model_ref$beta$logbeta[[1]]))
    }else{
      print("Select the best model.")
      perp = function(model){
        exp(-model$convergence$bound[model$convergence$its]/sum(model$settings$dim$wcounts$x))
      }
      k <- k1 + c(0,2,5,10,15,20)
      models <- list()
      perplex <- c()
      for (i in 1:6) {
        models[[i]] <- RefTM.STM(getCorpus(t(sc_data)), as.character(1:dim(sc_data)[1]), K = k[i], prevalence = ~covariate, refBeta = exp(model_ref$beta$logbeta[[1]]))
        perplex[i] <- perp(models[[i]])
      }
      log.lik <- data.frame(topics = k, LL = perplex)
      log.lik$first_derivative <- c(-Inf, (diff(log.lik$LL) / diff(log.lik$topics)))
      log.lik$second_derivative <- c(-Inf, -Inf, diff(log.lik$first_derivative)[-1]/diff(log.lik$topics[-1]))
      par(bty = 'n')
      plot(log.lik$topics, log.lik$LL, xlab="Number of topics", ylab="perplexity", type='o', pch=16, col='black', main='Model selection')
      points(log.lik$topics[which(log.lik$second_derivative == max(log.lik$second_derivative))], log.lik$LL[which(log.lik$second_derivative == max(log.lik$second_derivative))], pch=4, col='red', lwd = 7)
      model_sc <- models[[which(log.lik$second_derivative == max(log.lik$second_derivative))]]
    }
  }
  return(model_sc = model_sc)
}



