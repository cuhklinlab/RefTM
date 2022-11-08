#' Fit a Latent Dirichlet Allocation model to a text corpus
#'
#' @param dtm This may be a DocumentTermMatrix from the tm package, a simple_triplet_matrix from the slam package, or anything coercible to a simple_triplet_matrix. This includes a data frame with the document names in the first column, the terms in the second column, and the number or appearances in the third column.
#' @param passes How many times we look at each document.
#' @param batchsize The size of the minibatches.
#' @param maxiter The maximum iterations for the "E step" for each document (the updating of the per-document parameters within each minibatch). The default of 100 follows the reference implementation in python by the authors.
#' @param alpha hyperparameter.
#' @param kappa learning rate parameter. Lower values give greater weight to later iterations. For guaranteed convergence to a local optimum, kappa must lie in the interval (0.5,1].
#' @param tau_0 learning rate parameter. Higher values reduce the influence of early iterations.
#' @param K The number of topics
#' @importFrom methods is
#' @export

lda_svi <- function(dtm,passes=1,batchsize=256,maxiter=100,K,K0 = 0,alpha=1/K,kappa=0.7,tau_0=1024,refBeta=NULL){

	if (is(dtm,"DocumentTermMatrix")){
		if (!any(attr(dtm,'weighting') %in% c('term frequency','tf'))){
			stop('The DocumentTermMatrix object must use term frequency weighting')
		}
	}

	doc_ids <- dtm$i - 1#the c++ code expects 0-indexed ids
	docs <- dtm$dimnames$Docs

	term_ids <- dtm$j - 1#the c++ code expect 0-indexed ids
	terms <- dtm$dimnames$Terms

	counts <- dtm$v

	tmp = matrix(1/dtm$ncol,K,dtm$ncol)

	if(K0 != 0){
	  tmp[1:K0,] = refBeta
	}
	refBeta = tmp

	res_list <- lda_online_cpp(doc_ids,term_ids,counts,K,K0,passes,batchsize,maxiter=maxiter,alpha=alpha,tau_0=tau_0,kappa=kappa,refBeta = refBeta)

	gamma <- res_list$Gamma
	beta <- res_list$Beta

	colnames(gamma) <- seq(1:ncol(gamma))#topic labels
	rownames(gamma) <- unique(docs)

	colnames(beta) <- unique(terms)
	rownames(beta) <- seq(1:nrow(beta))

	# convert variational parameters to model parameters
	# (this follows from equation 2 in the paper)
	# Noting that the expectation of a Dirichlet(a) rv is a/sum(a)

	theta <- gamma

	for (i in seq(1,nrow(gamma))){
	  theta[i,] <- theta[i,]/sum(theta[i,])
	}

	for (i in seq(1,nrow(beta))){
	  beta[i,] <- beta[i,]/sum(beta[i,])
	}



	list('theta'=theta,'beta'=beta,'gamma'=gamma)#TODO: tidy output
}
