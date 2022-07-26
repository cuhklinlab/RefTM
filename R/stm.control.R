#Workhorse Function for the STM model
#compared to the original we have more initializations,
# more explicit options, trimmed fat, memoization

stm.control <- function(documents, vocab, settings, model=NULL,refBeta = NULL) {

  globaltime <- proc.time()
  verbose <- settings$verbose
  ##########
  #Step 1: Initialize Parameters
  ##########
  ngroups <- settings$ngroups
  if(is.null(model)) {
    if(verbose) cat(switch(EXPR=settings$init$mode,
                           Spectral = "Beginning Spectral Initialization \n",
                           LDA = "Beginning LDA Initialization \n",
                           Random = "Beginning Random Initialization \n",
                           Custom = "Beginning Custom Initialization \n"))
    #initialize
    model <- stm.init(documents, settings,refBeta)
    #if we were using the Lee and Mimno method of setting K, update the settings
    if(settings$dim$K==0) settings$dim$K <- nrow(model$beta[[1]])
    #unpack
    mu <- list(mu=model$mu)
    sigma <- model$sigma
    beta <- list(beta=model$beta)
    if(!is.null(model$kappa)) beta$kappa <- model$kappa
    lambda <- model$lambda
    convergence <- NULL
    #discard the old object
    rm(model)
  } else {
    if(verbose) cat("Restarting Model...\n")
    #extract from a standard STM object so we can simply continue.
    mu <- model$mu
    beta <- list(beta=lapply(model$beta$logbeta, exp))
    if(!is.null(model$beta$kappa)) beta$kappa <- model$beta$kappa
    sigma <- model$sigma
    lambda <- model$eta
    convergence <- model$convergence
    #manually declare the model not converged or it will stop after the first iteration
    convergence$stopits <- FALSE
    convergence$converged <- FALSE
    #iterate by 1 as that would have happened otherwise
    convergence$its <- convergence$its + 1
  }

  #Pull out some book keeping elements
  ntokens <- sum(settings$dim$wcounts$x)
  betaindex <- settings$covariates$betaindex
  stopits <- FALSE
  if(ngroups!=1) {
    # randomly assign groups so that subsample are representative
    groups <- base::split(seq_len(length(documents)),
                    sample(rep(seq_len(ngroups), length=length(documents))))
  }
  suffstats <- vector(mode="list", length=ngroups)

  if(settings$convergence$max.em.its==0) {
    stopits <- TRUE
    if(verbose) cat("Returning Initialization.")
  }
  ############
  #Step 2: Run EM
  ############
  while(!stopits) {

    #one set of updates with groups, another without.
    if(ngroups!=1) {
      #####
      #Blocked Updates
      #####
      # ordering of groups should be randomized
      for(i in sample(seq_len(ngroups))) {
        t1 <- proc.time()
        #update the group id
        gindex <- groups[[i]]
        #construct the group specific sets
        gdocs <- documents[gindex]
        if(is.null(mu$gamma)) {
          gmu <- mu$mu
        } else {
          gmu <- mu$mu[,gindex]
        }
        gbetaindex <- betaindex[gindex]
        glambda <- lambda[gindex,]

        #run the model
        suffstats[[i]] <- estep(documents=gdocs, beta.index=gbetaindex,
                                update.mu=(!is.null(mu$gamma)),
                                beta$beta, refBeta, glambda, gmu, sigma,
                                verbose)
        if(verbose) {
          msg <- sprintf("Completed Group %i E-Step (%d seconds). \n", i, floor((proc.time()-t1)[3]))
          cat(msg)
        }
        t1 <- proc.time()

        #if all slots are full.  Combine and run M-step
        if(!any(unlist(lapply(suffstats, is.null)))) {
          #Combine the sufficient statistics
          #(note this is somewhat kludgier than I would prefer
          # but it isn't very costly in terms of time so its fine)
          sigma.ss <- suffstats[[1]]$sigma
          lambda <- suffstats[[1]]$lambda
          beta.ss <- suffstats[[1]]$beta
          bound.ss <- suffstats[[1]]$bound
          for(j in 2:ngroups) {
            sigma.ss <- sigma.ss + suffstats[[j]]$sigma
            lambda <- rbind(lambda, suffstats[[j]]$lambda)
            for(a in 1:length(beta.ss)) {
              beta.ss[[a]] <- beta.ss[[a]] + suffstats[[j]]$beta[[a]]
            }
            bound.ss <- c(bound.ss, suffstats[[j]]$bound)
          }
          # Now do the updates themselves
          mu <- stm:::opt.mu(lambda=lambda, mode=settings$gamma$mode,
                       covar=settings$covariates$X, enet=settings$gamma$enet, ic.k=settings$gamma$ic.k,
                       maxits=settings$gamma$maxits)
          sigma <- stm:::opt.sigma(nu=sigma.ss, lambda=lambda,
                             mu=mu$mu, sigprior=settings$sigma$prior)
          beta <- opt.beta(beta.ss, beta$kappa, settings)

          if(verbose) {
           #M-step message
            timer <- floor((proc.time()-t1)[3])
            msg <- ifelse(timer>1,
                          sprintf("Completed M-Step (%d seconds). \n", floor((proc.time()-t1)[3])),
                          "Completed M-Step. \n")
            cat(msg)
          }
        }
      }
    } else {
      #####
      # Non-Blocked Updates
      #####
      t1 <- proc.time()
      #run the model
      suffstats <- estep(documents=documents, beta.index=betaindex,
                              update.mu=(!is.null(mu$gamma)),
                              beta$beta, refBeta,lambda, mu$mu, sigma,
                              verbose)
      msg <- sprintf("Completed E-Step (%d seconds). \n", floor((proc.time()-t1)[3]))
      if(verbose) cat(msg)
      t1 <- proc.time()
      sigma.ss <- suffstats$sigma
      lambda <- suffstats$lambda
      beta.ss <- suffstats$beta
      bound.ss <- suffstats$bound
      #do the m-step
      mu <- stm:::opt.mu(lambda=lambda,mode=settings$gamma$mode,
                   covar=settings$covariates$X, enet=settings$gamma$enet, ic.k=settings$gamma$ic.k,
                   maxits=settings$gamma$maxits)
      sigma <- stm:::opt.sigma(nu=sigma.ss, lambda=lambda,
                         mu=mu$mu, sigprior=settings$sigma$prior)
      beta <- opt.beta(beta.ss, beta$kappa, settings,refBeta)
      if(verbose) {
        timer <- floor((proc.time()-t1)[3])
        msg <- ifelse(timer>1,
                      sprintf("Completed M-Step (%d seconds). \n", floor((proc.time()-t1)[3])),
                      "Completed M-Step. \n")
        cat(msg)
      }
    }
    #Convergence
    convergence <- stm:::convergence.check(bound.ss, convergence, settings)
    stopits <- convergence$stopits

    #Print Updates if we haven't yet converged
    if(!stopits & verbose) stm:::report(convergence, ntokens=ntokens, beta, vocab,
                                       settings$topicreportevery, verbose)
  }
  #######
  #Step 3: Construct Output
  #######
  time <- (proc.time() - globaltime)[3]
  #convert the beta back to log-space
  beta$logbeta <- beta$beta
  for(i in 1:length(beta$logbeta)) {
    beta$logbeta[[i]] <- stm:::safelog(beta$logbeta[[i]])
  }
  beta$beta <- NULL
  lambda <- cbind(lambda,0)
  model <- list(mu=mu, sigma=sigma, beta=beta, settings=settings,
                vocab=vocab, convergence=convergence,
                theta=exp(lambda - log(rowSums(exp(lambda)))),
                #note altered from row.lse above because of a
                #Windows specific bug that was happening with
                #matrixStats package and large matrices 8/27
                eta=lambda[,-ncol(lambda), drop=FALSE],
                invsigma=solve(sigma), time=time, version=utils::packageDescription("stm")$Version)

  class(model) <- "STM"
  return(model)
}






