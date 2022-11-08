#include<unordered_map>


#include<RcppArmadillo.h>

#include "lda_svi.h"


using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List lda_online_cpp(IntegerVector doc_ids,IntegerVector terms,IntegerVector counts,int K,int K0,int passes,int batchsize,int maxiter,double tau_0,double kappa,double alpha,NumericMatrix refBeta){
	
	int n = doc_ids.size();
	unordered_map<int,unordered_map<int,int>> dtm;

	for (int i = 0;i<n;i++){
		dtm[doc_ids[i]][terms[i]] = counts[i];
	}

	int D = sort_unique(doc_ids).size();

	int V = sort_unique(terms).size();

	LDA_State lda(D,V,K,K0, dtm,alpha,refBeta);

	lda.fit_model(passes,batchsize,maxiter,tau_0,kappa,refBeta);

	return List::create(Rcpp::Named("Beta")=lda.beta,Rcpp::Named("Gamma")=lda.gamma);

}
