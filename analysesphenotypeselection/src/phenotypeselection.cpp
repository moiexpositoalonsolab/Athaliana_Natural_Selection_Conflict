#include <stdlib.h>
#include <cstdio>
#include <stdio.h>
#include <random>
#include <math.h>
#include <vector>
#include <list>
#include <string>
#include <iostream>

#include <fstream>
#include <sstream>

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// https://barumpark.com/blog/2019/Sampling-Integers/
#define ARMA_64BIT_WORD 1 /// https://stackoverflow.com/questions/40592054/large-matrices-in-rcpparmadillo-via-the-arma-64bit-word-define
//// when armadillo is loaded, //// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>


// [[Rcpp::export]]
Rcpp::List linearselectionC(arma::vec w, arma::mat z){
  arma::mat P=arma::cov(z); //
  arma::mat Pinv=arma::pinv(P);
  arma::mat d = std::move(arma::join_rows( w, z));
  arma::mat sall = arma::cov(d);
  sall.shed_row(0);
  arma::vec s=sall.col(0);
  arma::vec B = Pinv * s;

  Rcpp::List mylist(2);
  mylist.names() = Rcpp::CharacterVector::create("beta","s");
  mylist(0) = B;
  mylist(1) = s;

  return mylist;
}

// [[Rcpp::export]]
Rcpp::List linearresponseC(arma::vec B, arma::mat g ){
  arma::vec  deltaz = g * B;
  // return deltaz;
  return Rcpp::List::create(
          Rcpp::Named("deltaZ") = deltaz);
}

// [[Rcpp::export]]
Rcpp::List quadraticselectionC(arma::vec w, arma::mat z){
  arma::mat P=arma::cov(z);
  arma::mat Pinv=arma::pinv(P);
  arma::mat c(z.n_cols,z.n_cols,arma::fill::zeros);
  arma::vec tmp(z.n_rows,arma::fill::zeros);
  arma::mat covval(1,1,arma::fill::zeros);
  for(int i=0; i< z.n_cols ; i++){
    for(int j=0; j< z.n_cols ; j++){
      tmp= z.col(i) % z.col(j); // % = element-wise prod
      covval = arma::cov(w, tmp);
      c(i,j) = covval.at(0,0);
    }
  }
  arma::mat gamma = Pinv * c * Pinv;
  return Rcpp::List::create(
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("c") = c,
    Rcpp::Named("P") = c);
}

// [[Rcpp::export]]
Rcpp::List quadraticresponseC(arma::mat P, arma::mat g , arma::mat gamma, arma::vec s, arma::vec B){
  arma::mat deltaP = (P*gamma*P) - s * s.t();
  arma::mat deltaG = g * (gamma-(B *B.t())) * g;
  return Rcpp::List::create(
              Rcpp::Named("deltaP") = deltaP,
              Rcpp::Named("deltaG") = deltaG);
}

// Version 1 implementation, inefficient
// [[Rcpp::export]]
Rcpp::List PHENOSELECTIONc(arma::vec w, arma::mat z, arma::mat g, int replicates){

  // create a list where results will be stored. One per bootstrap replicate
  Rcpp::List toreport(replicates);
  // create an index for the bootstrap
  arma::uvec index = arma::linspace<arma::uvec>(0L, w.n_elem - 1L, w.n_elem);

  // iterate through replicates
  for(int i=0; i<replicates; i++){
    // bootstrap
    arma::uvec myindex = Rcpp::RcppArmadillo::sample(index, index.n_elem, true);

    arma::mat ztmp = z.rows(myindex);
    arma::vec wtmp = w.rows(myindex);

    Rcpp::List lins, quas, linr, quar;  // initialize lists
    lins = linearselectionC(w=wtmp,z=ztmp);
    quas = quadraticselectionC(w=wtmp,z=ztmp);
      arma::vec B = lins["beta"];
      arma::vec s = lins["s"];
      arma::mat gamma = quas["gamma"];
      arma::mat c = quas["c"];
      arma::mat P = quas["P"];
    linr = linearresponseC(B,g);
      arma::vec deltaZ = linr["deltaZ"];
    quar = quadraticresponseC(P,g,gamma,s,B);
      arma::mat deltaP = quar["deltaP"];
      arma::mat deltaG = quar["deltaG"];

    Rcpp::List mylist(7);
    mylist.names() = Rcpp::CharacterVector::create("beta","s","gamma","c","deltaZ","deltaP","deltaG");
    mylist(0) = B;
    mylist(1) = s;
    mylist(2) = gamma;
    mylist(3) = c;
    mylist(4) = deltaZ;
    mylist(5) = deltaP;
    mylist(6) = deltaG;

    toreport(i) = mylist;
  }

  return toreport;
}

// Version 2 implementation, only summary statistics
// [[Rcpp::export]]
Rcpp::List PHENOSELECTIONcboot(arma::vec w, arma::mat z, arma::mat g, int r, bool debug){
  // initialize list objects and create an index for the bootstrap
  arma::uvec index = arma::linspace<arma::uvec>(0L, w.n_elem - 1L, w.n_elem);
  Rcpp::List lins, quas, linr, quar;  // intermediate
  //  First run
    lins = linearselectionC(w=w,z=z);
    quas = quadraticselectionC(w=w,z=z);
    arma::vec B = lins["beta"];
    arma::vec s = lins["s"];
    arma::mat gamma = quas["gamma"];
    arma::mat c = quas["c"];
    arma::mat P = quas["P"];
    linr = linearresponseC(B,g);
    arma::vec deltaZ = linr["deltaZ"];
    quar = quadraticresponseC(P,g,gamma,s,B);
    arma::mat deltaP = quar["deltaP"];
    arma::mat deltaG = quar["deltaG"];
    //
    arma::vec B2 = pow(B,2); // sum squared values just using element-wise multip.
    arma::vec s2 = pow(s,2);
    arma::mat gamma2 = pow(gamma,2);
    arma::mat c2 = pow(c,2);
    arma::vec deltaZ2 = pow(deltaZ,2);
    arma::mat deltaP2 = pow(deltaP,2);
    arma::mat deltaG2 = pow(deltaG,2);

  for(int i=1; i<r; i++){ // starting after 1st iteration
    arma::uvec myindex = Rcpp::RcppArmadillo::sample(index, index.n_elem, true);
    arma::mat ztmp = z.rows(myindex);
    arma::vec wtmp = w.rows(myindex);
    lins = linearselectionC(w=wtmp,z=ztmp);
    quas = quadraticselectionC(w=wtmp,z=ztmp);
    arma::vec B_ = lins["beta"];
      // if(debug) std::cout << B_ << std::endl;
    arma::vec s_ = lins["s"];
      // if(debug) std::cout << s_ << std::endl;
    arma::mat gamma_ = quas["gamma"];
      // if(debug) std::cout << gamma_ << std::endl;
    arma::mat c_ = quas["c"];
      // if(debug) std::cout << c_ << std::endl;
    arma::mat P_ = quas["P"];
      // if(debug) std::cout << P_ << std::endl;
    linr = linearresponseC(B_,g);
    arma::vec deltaZ_ = linr["deltaZ"];
      // if(debug) std::cout << deltaZ_ << std::endl;
    quar = quadraticresponseC(P_,g,gamma_,s_,B_);
    arma::mat deltaP_ = quar["deltaP"];
      // if(debug) std::cout << deltaP_ << std::endl;
    arma::mat deltaG_ = quar["deltaG"];
      // if(debug) std::cout << deltaG_ << std::endl;
    // sum to previous set
    B += B_;
    s += s_;
    gamma += gamma_;
    c += c_;
    deltaZ += deltaZ_;
    deltaP += deltaP_;
    deltaG += deltaG_;
      // if(debug) std::cout << "error here" << std::endl;
    B2 += pow(B_,2); // sum squared values just using element-wise multip.
    s2 += pow(s_,2);
    gamma2 += pow(gamma_,2);
    c2 += pow(c_,2);
    deltaZ2 += pow(deltaZ_,2);
    deltaP2 += pow(deltaP_,2);
    deltaG2 += pow(deltaG_,2);

  }
  // Compute estimators
  // mean E[x] = sum x /n
  B = B / r;
  B2 = B2 / r - pow(B,2);
  s = s / r;
  s2 = s2 / r - pow(s,2);
  gamma = gamma / r;
  gamma2 = gamma2 / r - pow(gamma,2);
  c = c / r;
  c2 = c2 / r - pow(c,2);
  deltaZ = deltaZ / r;
  deltaZ2 = deltaZ2 / r - pow(deltaZ,2);
  deltaP = deltaP / r;
  deltaP2 = deltaP2 / r - pow(deltaP,2);
  deltaG = deltaG / r;
  deltaG2 = deltaG2 / r - pow(deltaG,2);
  //
  // lists to report
  //
  Rcpp::List meanboot(7);
  Rcpp::List seboot(7);
  meanboot.names() = Rcpp::CharacterVector::create("beta","s","gamma","c","deltaZ","deltaP","deltaG");
  seboot.names() = Rcpp::CharacterVector::create("beta","s","gamma","c","deltaZ","deltaP","deltaG");
  meanboot(0) = B;
  meanboot(1) = s;
  meanboot(2) = gamma;
  meanboot(3) = c;
  meanboot(4) = deltaZ;
  meanboot(5) = deltaP;
  meanboot(6) = deltaG;
  seboot(0) = B2;
  seboot(1) = s2;
  seboot(2) = gamma2;
  seboot(3) = c2;
  seboot(4) = deltaZ2;
  seboot(5) = deltaP2;
  seboot(6) = deltaG2;
  //// final lists
  return Rcpp::List::create(
                Rcpp::Named("mean") = meanboot,
                Rcpp::Named("se") = seboot);
  //// final lists
  // return Rcpp::List::create(
  //               Rcpp::Named("b") = B,
  //               Rcpp::Named("s") = s);
}

