#ifndef _lda_svi_guard
#define _lda_svi_guard

#include<unordered_map>
#include <RcppArmadillo.h>


arma::rowvec dirichlet_expectation(const arma::rowvec &alpha);

//Represent documents and topics as integers; handle conversion in R code.

//TODO: make this a class
struct LDA_State{
	
	arma::mat gamma;//document/topic
	arma::mat beta;//word/topic

	int D;//total number of documents
	int V;//vocabulary size
	int K;//The number of topics we're fitting
                int K0;//The number of reference topics 
	int t;//how many minibatches we've processed so far
                Rcpp::NumericMatrix  refBeta;

	double alpha;

	// Data
	std::unordered_map<int,std::unordered_map<int,int>> dtm;

	public:
		LDA_State(int D,int V,int K,int K0,std::unordered_map<int,std::unordered_map<int,int>> &dtm,double alpha, Rcpp::NumericMatrix refBeta);
		void update_minibatch(std::vector<int> documents);
		void fit_model(int passes,int batchsize,int maxiter,double tau_0,double kappa, Rcpp::NumericMatrix refBeta);
		void update_minibatch(std::vector<int> docs,int maxiter,double tau_0,double kappa, Rcpp::NumericMatrix refBeta);//this could probably be private
};

#endif
