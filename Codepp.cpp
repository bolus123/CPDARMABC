#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]



int InvertQ(const arma::colvec& coef){
  int out;
  arma::colvec coef1 = coef;
  int k = coef1.n_elem;
  arma::mat blockmat(k, k, arma::fill::zeros);
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::vec mod;
  
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++){
      if (i == 0) {
        blockmat(i, j) = coef1(j);
      } else if (i > 0) {
        if ((i - 1) == j) {
          blockmat(i, j) = 1;
        }
      }
    }
  }
  
  
  arma::eig_gen(eigval, eigvec, blockmat);
  
  mod = arma::pow(arma::pow(arma::real(eigval), 2) + arma::pow(arma::imag(eigval), 2), 0.5);
  
  if (mod.max() >= 1) {
    out = 0;
  } else {
    out = 1;
  }
  
  return out;
  
}

arma::mat parsMat(const int& n, const arma::colvec& coef){
  arma::mat out(n, n, arma::fill::zeros);
  arma::vec coef1 = coef;
  int k = coef1.n_elem;
  int i;
  int j;
  int r;
  
  int checkInv = InvertQ(coef1);
  
  if (checkInv == 1) {
  
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        for (r = 0; r < k; r++) {
          if (i == j) {
            out(i, j) = 1;
          } else if ((i - 1 - r) == j) {
            out(i, j) = coef1(r);
          }
        }
      }
    }
    
  } else {
    //Rcpp::stop("The parameter vector is not invertible/stationary");
    Rcpp::warning("The parameter vector is not invertible/stationary");
  }
  
  return out;
  
}


// [[Rcpp::export]]
arma::mat DifMat(const int& n, const arma::colvec& order) {
  arma::mat Dif(n, n, arma::fill::eye);
  arma::mat Dif1(n, n, arma::fill::eye);
  arma::mat Dif2;
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if ((i - 1) == j) {
        Dif1(i, j) = -1;
      }
    }
  }
  
  for (int r = 0; r < order(1); r++) {
    Dif2 = Dif1.submat(r + 1, r, n - 1, n - 1);
    Dif = Dif2 * Dif;
  }
  
  return Dif;
  
}

// [[Rcpp::export]]
arma::mat OmegaMat(const int& n, const arma::colvec& order, const arma::colvec& phi, const arma::colvec& theta){
  arma::mat out;
  arma::mat phiMat(n - order(1), n - order(1), arma::fill::eye);
  arma::mat thetaMat(n - order(1), n - order(1), arma::fill::eye);
  arma::mat Dif;
  //arma::mat LeftMat;
  //arma::mat RhoMat;
  arma::mat OmegaMat;
  //arma::vec gamma0;
  
  if (order(0) > 0) {
    phiMat = parsMat(n - order(1), -phi);
  } 
  
  if (order(2) > 0) {
    thetaMat = parsMat(n - order(1), theta);
  } 
  
  if (order(1) > 0) {
    
    Dif = DifMat(n, order);
    
    //LeftMat = phiMat * Dif;
    OmegaMat = arma::pinv(phiMat * Dif) * thetaMat;
    
  } else {
    
    //LeftMat = phiMat;
    OmegaMat = arma::inv(phiMat) * thetaMat;
    
  }
  
  //RhoMat = arma::inv(thetaMat) * LeftMat;
  //out = OmegaMat * OmegaMat.t();
  //gamma0 = diagvec(out);
  
  return OmegaMat;
  
}


// [[Rcpp::export]]
Rcpp::List SigmaMat(const int& n, const arma::colvec& order, const arma::colvec& phi, const arma::colvec& theta){
  arma::mat out;
  arma::mat phiMat(n - order(1), n - order(1), arma::fill::eye);
  arma::mat thetaMat(n - order(1), n - order(1), arma::fill::eye);
  arma::mat Dif;
  arma::mat LeftMat;
  arma::mat RhoMat;
  arma::mat OmegaMat;
  arma::vec gamma0;
  
  if (order(0) > 0) {
    phiMat = parsMat(n - order(1), -phi);
  } 
  
  if (order(2) > 0) {
    thetaMat = parsMat(n - order(1), theta);
  } 
  
  if (order(1) > 0) {
    
    Dif = DifMat(n, order);
    
    
    LeftMat = phiMat * Dif;
    OmegaMat = arma::pinv(phiMat * Dif) * thetaMat;
    
  } else {
    
    LeftMat = phiMat;
    OmegaMat = arma::inv(phiMat) * thetaMat;
    
  }
  
  RhoMat = arma::inv(thetaMat) * LeftMat;
  out = OmegaMat * OmegaMat.t();
  gamma0 = diagvec(out);
  
  return Rcpp::List::create(Rcpp::Named("SigmaMat") = out,
                            Rcpp::Named("gamma0")   = gamma0,
                            Rcpp::Named("OmegaMat") = OmegaMat,
                            Rcpp::Named("RhoMat") = RhoMat);
  
}


arma::colvec simInnov(const int& n, const Rcpp::String& XSim, const Rcpp::NumericVector& XPars) {
  
  arma::colvec out(n);
  
  if (XSim == "norm") {
    out.randn();
    out = out * std::sqrt(XPars[1]) + XPars[0];
  }
  
  return out;
  
}

// [[Rcpp::export]]
arma::colvec simARIMA(const double& mu, const double& sigma2, const Rcpp::String& innovDist, const arma::mat& OmegaMat, 
                      const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta) {
  
  Rcpp::NumericVector XPars = {0, 1};
  
  int n = OmegaMat.n_cols;
  
  arma::colvec out = OmegaMat * (simInnov(n, innovDist, XPars) * std::sqrt(sigma2)) + mu;

  if (BetaFlg == 1) {
    out = out + W * Beta;
  }
  
  return out;
  
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::colvec BoxCox(const arma::colvec& Y, const double& lambda) {
  
  arma::colvec out;
  
  if (lambda > 0) {
    
    out = (arma::pow(Y, lambda) - 1) / lambda;
    
  } else if (lambda == 0) {
    
    out = arma::log(Y);
    
  }
  
  return out;
  
} 

// [[Rcpp::export]]
arma::colvec invBoxCox(const arma::colvec& X, const double& lambda) {
  
  arma::colvec out;
  
  if (lambda > 0) {
    
    out = arma::exp(arma::log(lambda * X + 1) / lambda);
    
  } else if (lambda == 0) {
    
    out = arma::exp(X);
    
  }
  
  return out;
  
} 


// [[Rcpp::export]]
arma::colvec BoxCox1stDeriv(const arma::colvec& Y, const double& lambda) {
  
  arma::colvec out;
  
  if (lambda > 0) {
    
    out = arma::pow(Y, lambda - 1);
    
  } else if (lambda == 0) {
    
    out = arma::pow(Y, -1);
    
  }
  
  return out;
  
}

// [[Rcpp::export]]
arma::colvec invBoxCox2ndDeriv(const arma::colvec& X, const double& lambda) {
  
  arma::colvec out;
  
  if (lambda > 0) {

    out = std::pow(lambda, 3) * (lambda + 1) * arma::pow((lambda * X + 1), -lambda - 2);
    
  } else if (lambda == 0) {
    
    out = arma::exp(X);
    
  }
  
  return out;
  
}

arma::colvec biasajdBoxCox(const arma::colvec& pred, const double& lambda, const double& sigma2) {
  arma::colvec out = invBoxCox(pred, lambda) + invBoxCox2ndDeriv(pred, lambda) / 2 * sigma2;
  return out;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
Rcpp::List fastLm(const arma::mat& X, const arma::colvec& y) {
  int n = X.n_rows, k = X.n_cols;
  
  arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
  arma::colvec res  = y - X*coef;           // residuals
  
  // std.errors of coefficients
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);
  
  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));  
  
  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("stderr")       = std_err,
                            Rcpp::Named("df.residual")  = n - k);
}

//// [[Rcpp::export]]
//Rcpp::List estCSS(const arma::colvec& X, const arma::colvec& order, const int& muFlg) {
//  
//  int n = X.n_elem;
//  arma::colvec tmpX = X;
//  arma::colvec tmpX1;
//  arma::mat tmpA;
//  arma::mat tmpA1(n - order(1), order(0));
//  arma::mat tmpB;
//  arma::mat tmpB1(n - order(0) - order(1), order(0));
//  arma::colvec phi;
//  arma::colvec theta;
//  arma::colvec tmpEps;
//  double mu = 0;
//  arma::colvec sigma2;
//  
//  tmpEps = tmpX;
//  
//  if (muFlg == 1) {
//   mu = arma::mean(tmpX);
//   tmpX = tmpX - mu;
//   tmpEps = tmpX;
// }
//  
//  if (order(1) > 0) {
//    tmpX = DifMat(n, order) * tmpX;
//    tmpEps = tmpX;
//  } 
//  
//  if (order(0) > 0) {
//    for (int j = 0; j < order(0); j++) {
//      for (int i = 0; i < n - order(1); i++) {
//        if (i - order(0) >= 0) {
//          tmpA1(i, j) = tmpX(i - j - 1);
//        }
//      }
//    }
//    
//    tmpA = tmpA1.submat(order(0), 0, n - order(1) - 1, order(0) - 1);
//    tmpX1 = tmpX.subvec(order(0), n - order(1) - 1);
//    phi = arma::solve(tmpA, tmpX1);
//    tmpX = tmpX1 - tmpA * phi;
//    tmpEps = tmpX;
//  }
//  
//  if (order(2) > 0) {
//    for (int j = 0; j < order(2); j++) {
//      for (int i = 0; i < n - order(1) - order(0); i++) {
//        if (i - order(2) >= 0) {
//          tmpB1(i, j) = tmpX(i - j - 1);
//        }
//      }
//    }
//    
//    tmpB = tmpB1.submat(order(2), 0, n - order(1) - order(0) - 1, order(2) - 1);
//    tmpX1 = tmpX.subvec(order(2), n - order(1) - order(0) - 1);
//    theta = arma::solve(tmpB, tmpX1);
//    tmpEps = tmpX1 - tmpB * theta;
//  }
//  
//  //sigma2 = arma::mean(arma::pow(tmpEps, 2));
//  sigma2 = arma::sum(arma::pow(tmpEps, 2)) / tmpEps.n_elem;
//  
//  return Rcpp::List::create(Rcpp::Named("phi") = phi, 
//                            Rcpp::Named("theta") = theta,
//                            Rcpp::Named("mu") = mu, 
//                            Rcpp::Named("sigma2") = sigma2, 
//                            Rcpp::Named("eps") = tmpEps);
//  
//}

// [[Rcpp::export]]
Rcpp::List ArimaCpp(const arma::colvec& x, const arma::colvec& order, const int& include_mean) {
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("forecast"); 
  Rcpp::Function rfunction = pkg["Arima"];  
  Rcpp::List out = rfunction(Rcpp::Named("y", x), Rcpp::Named("order", order), Rcpp::Named("include.mean", include_mean));
  return out;
}

// [[Rcpp::export]]
double AIC(const double& loglik, const int& npars) {
  double out = loglik * (-2) + 2 * (npars);
  return out;
}

// [[Rcpp::export]]
double BIC(const double& loglik, const int& npars, const int& nobs) {
  double out = npars * std::log(nobs) - 2 * loglik;
  return out;
}

// [[Rcpp::export]]
double AICc(const double& loglik, const int& npars, const int& nobs) {
  double out = AIC(loglik, npars) + (2 * std::pow(npars, 2) + 2 * npars) / (nobs - npars - 1);
  return out;
}

// [[Rcpp::export]]
Rcpp::List loglik(const arma::colvec& Y, const arma::colvec& order, const int& include_mean, 
         const int& BoxCoxFlg, const double& lambda, const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta) {

  Rcpp::List out;
  double loglik;
  double loglik1;
  double loglik2;
  arma::colvec X;
  Rcpp::List model;
  
  int npars = order(0) + order(2) + 1;
  
  if (order(1) == 0) {
    npars = npars + include_mean;
  }
    
  X = Y;
    
  if (BoxCoxFlg == 0) {
    loglik2 = 0;
  } else {
    npars = npars + 1;
    X = BoxCox(X, lambda);
    loglik2 = arma::sum(arma::log(BoxCox1stDeriv(Y, lambda)));
  }
  
  if (BetaFlg == 1) {
    X = X - W * Beta;
  }  
  
  model = ArimaCpp(X, order, include_mean);
  loglik1 = model["loglik"];
  
  loglik = loglik1 + loglik2;

  out = model;
  out["loglik"] = loglik;
  out["aic"] = AIC(loglik, npars);
  out["bic"] = BIC(loglik, npars, out["nobs"]);
  out["aicc"] = AICc(loglik, npars, out["nobs"]);
  
  if (BoxCoxFlg == 1) {
    out["fitted"] = biasajdBoxCox(out["fitted"], lambda, out["sigma2"]);
  }
  
  return out;
  
}



// [[Rcpp::export]]
Rcpp::List loglikFoward(const arma::colvec& Y, const int& include_mean, 
                  const int& BoxCoxFlg, const double& lambda, const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, 
                  const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q) {
  
  Rcpp::List out;
  Rcpp::List tmp;
  Rcpp::List minModel;
  arma::colvec tmpOrder(5);
  arma::colvec minOrder(5);
  double mincrit = LONG_LONG_MAX;
  double tmpCrit;
  int errFlg;
  
  int n = Y.n_elem;
  int nobs;
  int npars;
  int i;
  int j;
  int k;
  int r = 0;
  
  for (i = 0; i < (max_p + 1); i++) {

    for (j = 0; j < (max_d + 1); j++) {

      for (k = 0; k < (max_q + 1); k++) {
        
        errFlg = 0;
        
        r = r + 1;
        tmpOrder = {i * 1.0, j * 1.0, k * 1.0, include_mean * 1.0, BoxCoxFlg * 1.0};
        
        nobs = n - j;
        npars = i + k + 1 + BoxCoxFlg;
        
        if (j == 0) {
          npars = npars + include_mean;
        }
        
        if (nobs > npars) {
          try{
            
            tmp = loglik(Y, tmpOrder.subvec(0, 2), include_mean, BoxCoxFlg, lambda, BetaFlg, W, Beta);
            
          } catch(...) {
            
            errFlg = 1;
            
          }
          
          if (errFlg == 0) {
            tmpCrit = tmp[crit];
            //Rcpp::Rcout << "Iter:" << r << ", minCrit:" << mincrit << ", crit:" << tmpCrit << "\n";
            
            if (tmpCrit < mincrit) {
              mincrit = tmpCrit;
              minOrder = tmpOrder;
              minModel = tmp;
            }
            
          }
          
        } else {
          Rcpp::warning("Sample size is not enough");
        }

      }
    }
  }
  
  out = Rcpp::List::create(Rcpp::Named("order") = minOrder, 
                           Rcpp::Named("coef") = minModel["coef"],
                           Rcpp::Named("sigma2") = minModel["sigma2"],
                           Rcpp::Named("lambda") = lambda,
                           Rcpp::Named("crit") = crit,
                           Rcpp::Named("critVal") = mincrit,
                           Rcpp::Named("loglik") = minModel["loglik"],
                           Rcpp::Named("fitted") = minModel["fitted"],
                           Rcpp::Named("model") = minModel);
  
  return out;
  
}

//// [[Rcpp::export]]
//Rcpp::List OptLambdaCritSecant(const arma::colvec& Y, const double& lambda0, const double& lambda1, const int& include_mean, 
//                     const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, const Rcpp::String& crit, 
//                     const int& max_p, const int& max_d, const int& max_q, const double& tol, const int& maxIter) {
//
//  double dif = LONG_LONG_MAX;
//  
//  double la0 = lambda0;
//  double la1 = lambda1;
//  double la2;
//  
//  Rcpp::List model0 = loglikFoward(Y, include_mean, 
//                                   1, la0, BetaFlg, W, Beta, 
//                                   crit, max_p, max_d, max_q);
//  
//  Rcpp::List model1 = loglikFoward(Y, include_mean, 
//                                   1, la1, BetaFlg, W, Beta, 
//                                   crit, max_p, max_d, max_q);
//  
//  double f0 = model0["critVal"];
//  double f1 = model1["critVal"];
//  double f2;
//  
//  //Rcpp::Rcout << "f0:" << f0 << ", f1:" << f1 << "\n";
//  
//  Rcpp::List model2;
//  
//  int iter = 0;
//  
//  while(dif >= tol && iter <= maxIter) {
//    
//    iter = iter + 1;
//    
//    la2 = la1 - f1 * (la1 - la0) / (f1 - f0);
//    
//    model2 = loglikFoward(Y, include_mean, 
//                      1, la2, BetaFlg, W, Beta, 
//                      crit, max_p, max_d, max_q);
//    
//    f2 = model2["critVal"];
//    
//    dif = std::abs(f2 - f1);
//    
//    Rcpp::Rcout << "Iter:" << iter << ", lambda:" << la2 << ", dif:" << dif << "\n";
//    
//    f0 = f1;
//    f1 = f2;
//    
//    model0 = model1;
//    model1 = model2;
//    
//  }
//  
//  return model0;
//  
//}


// [[Rcpp::export]]
Rcpp::List OptLambdaCritBisec(const arma::colvec& Y, const double& lowerLambda, const double& upperLambda, const int& breakPoint,
                              const int& include_mean, const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, 
                              const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q, const double& tol, 
                              const int& maxIter) {
  
  Rcpp::List model2 = 0;
  
  double la0 = lowerLambda;
  double la1 = upperLambda;
  double la2;
  double tmpla0;
  
  arma::colvec laVec(breakPoint + 1);
    
  Rcpp::List model0 = loglikFoward(Y, include_mean, 
                                   1, la0, BetaFlg, W, Beta, 
                                   crit, max_p, max_d, max_q);
  
  Rcpp::List model1 = loglikFoward(Y, include_mean, 
                                   1, la1, BetaFlg, W, Beta, 
                                   crit, max_p, max_d, max_q);
  
  double f0 = model0["critVal"];
  double f1 = model1["critVal"];
  double f2(breakPoint + 1);
  
  double dif = std::abs(f0 - f1);
  
  arma::colvec tmp(breakPoint + 1);
  
  double step;
  
  //Rcpp::Rcout << "f0:" << f0 << ", f1:" << f1 << "\n";
  
  
  
  int iter = 0;
  int i;
  while(dif >= tol && iter <= maxIter) {
    
    iter = iter + 1;
    
    //Rcpp::Rcout << "iter:" << iter << "\n";
    
    step = std::abs(la1 - la0) / breakPoint;
    
    if (la0 <= la1) {
      tmpla0 = la0;
    } else {
      tmpla0 = la1;
    }
    
    //Rcpp::Rcout << "step:" << step << "\n";
    for (i = 0; i < (breakPoint + 1); i++) {
      //Rcpp::Rcout << "i:" << i << ", step i:" << step * i << "\n";
      la2 = tmpla0 + step * i;
      model2 = loglikFoward(Y, include_mean, 
                            1, la2, BetaFlg, W, Beta, 
                            crit, max_p, max_d, max_q);
      f2 = model2["critVal"];
      
      if (f2 <= f1) {
        if (f2 <= f0) {
          f0 = f2;
          la0 = la2;
        } else {
          f1 = f2;
          la1 = la2;
        }
      }
      
      //Rcpp::Rcout << "la0:" << la0 << ", la1:" << la1 << ", la2:" << la2 << "\n";
      //Rcpp::Rcout << "f0:" << f0 << ", f1:" << f1 << ", f2:" << f2 << "\n";
      
    }

    dif = std::abs(f0 - f1);
    //Rcpp::Rcout << "dif:" << dif << "\n";
    
  }

  return model2;
  
}

// [[Rcpp::export]]
Rcpp::List optModel(const arma::colvec& Y, const double& lowerLambda, const double& upperLambda, const int& breakPoint,
                              const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, 
                              const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q, const double& tol, 
                              const int& maxIter) {
  
  
  double minCrit = LONG_LONG_MAX;
  double tmpCrit;
  Rcpp::List minModel;
  //int errCode1;
  //int errCode2;
  //int errCode3;
  //int errCode4;
  
  //// with mean
  try {
    Rcpp::List modelNoBoxCoxWithMean = loglikFoward(Y, 1, 0, 0, BetaFlg, W, Beta, 
                                                     crit, max_p, 0, max_q);
    
    if( modelNoBoxCoxWithMean.containsElementNamed("critVal") ){
      minCrit = modelNoBoxCoxWithMean["critVal"];
      minModel = modelNoBoxCoxWithMean;
    } 
    //int errCode1 = 0;
    //Rcpp::Rcout << "minCrit:" << minCrit << "\n";
  } catch(...) {
    //int errCode1 = 1;
  }
  
  try{
    Rcpp::List modelBoxCoxWithMean = OptLambdaCritBisec(Y, lowerLambda, upperLambda, breakPoint,
                                                         1, BetaFlg, W, Beta, 
                                                         crit, max_p, 0, max_q, tol, 
                                                         maxIter);
    
    if( modelBoxCoxWithMean.containsElementNamed("critVal") ){
      tmpCrit = modelBoxCoxWithMean["critVal"];
      //Rcpp::Rcout << "tmpCrit:" << tmpCrit << "\n";
      if (tmpCrit <= minCrit) {
        minCrit = tmpCrit;
        minModel = modelBoxCoxWithMean;
      }
    }
    //int errCode2 = 0;
  } catch(...) {
    //int errCode2 = 1;
  }
  
  
  
  //// without mean
  try{
    Rcpp::List modelNoBoxCoxNoMean = loglikFoward(Y, 0, 0, 0, BetaFlg, W, Beta, 
                                                   crit, max_p, max_d, max_q);
    
    if( modelNoBoxCoxNoMean.containsElementNamed("critVal") ){
      tmpCrit = modelNoBoxCoxNoMean["critVal"];
      //Rcpp::Rcout << "tmpCrit:" << tmpCrit << "\n";
      if (tmpCrit <= minCrit) {
        minCrit = tmpCrit;
        minModel = modelNoBoxCoxNoMean;
      }
    }
    //int errCode3 = 0;
  } catch(...) {
    //int errCode3 = 1;
  }
  
  try{
    Rcpp::List modelBoxCoxNoMean = OptLambdaCritBisec(Y, lowerLambda, upperLambda, breakPoint,
                                                       0, BetaFlg, W, Beta, 
                                                       crit, max_p, max_d, max_q, tol, 
                                                       maxIter);
    
    if( modelBoxCoxNoMean.containsElementNamed("critVal") ){
      tmpCrit = modelBoxCoxNoMean["critVal"];
      //Rcpp::Rcout << "tmpCrit:" << tmpCrit << "\n";
      if (tmpCrit <= minCrit) {
        minCrit = tmpCrit;
        minModel = modelBoxCoxNoMean;
      }
    }
    //int errCode4 = 0;
  } catch(...) {
    //int errCode4 = 1;
  }

  return minModel;
  
}

// [[Rcpp::export]]
double loglikRatio(const arma::colvec& Y1, const arma::colvec& Y2, const double& loglik0, 
                    const double& lowerLambda, const double& upperLambda, const int& breakPoint,
                    const int& BetaFlg, const arma::mat& W1, const arma::mat& W2, const arma::colvec& Beta, 
                    const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q, const double& tol, 
                    const int& maxIter) {
  
  //Rcpp::Rcout << "Y1:" << Y1 << "\n";
  //Rcpp::Rcout << "Y2:" << Y2 << "\n";
  //Rcpp::Rcout << "W1:" << W1 << "\n";
  //Rcpp::Rcout << "W2:" << W2 << "\n";
  
  Rcpp::List model1 = optModel(Y1, lowerLambda, upperLambda, breakPoint,
                    BetaFlg, W1, Beta, crit, max_p, max_d, max_q, tol, 
                    maxIter);
  
  double loglik1 = model1["loglik"];
  
  Rcpp::List model2 = optModel(Y2, lowerLambda, upperLambda, breakPoint,
                    BetaFlg, W2, Beta, crit, max_p, max_d, max_q, tol, 
                    maxIter);
  
  double loglik2 = model2["loglik"];
  
  double llr = -2 * (loglik0 - (loglik1 + loglik2));
  
  Rcpp::Rcout << "GLR:" << llr << ", LL0:" << loglik0 << ", LL1:" << loglik1 << ", LL2:" << loglik2 << "\n";
  
  return llr;
  
}

// [[Rcpp::export]]
Rcpp::List loglikRatioMax(const arma::colvec& Y, const int& minSize, 
                   const double& lowerLambda, const double& upperLambda, const int& breakPoint,
                   const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, 
                   const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q, const double& tol, 
                   const int& maxIter) {
  
  int t;
  int n = Y.n_elem;
  arma::colvec Y1;
  arma::colvec Y2;
  arma::mat W1;
  arma::mat W2;
  arma::colvec llrVec(n);
  llrVec.fill(LONG_LONG_MIN);
  
  int k = 0;
  
  if (BetaFlg == 1) {
    k = W.n_cols;
  }
  
  Rcpp::List model0 = optModel(Y, lowerLambda, upperLambda, breakPoint,
                               BetaFlg, W, Beta, crit, max_p, max_d, max_q, tol, 
                               maxIter);
  
  
  double loglik0 = model0["loglik"];
  double llr;
  double maxLlr = LONG_LONG_MIN;
  int maxT = 0;
  
  for (t = minSize - 1; t <= (n - minSize - 1); t++) {
    
    Rcpp::Rcout << "t:" << t + 1 << "\n";
    
    Y1 = Y.subvec(0, t);
    Y2 = Y.subvec(t + 1, n - 1);
    
    //Rcpp::Rcout << "Y1:" << Y1 << "\n";
    //Rcpp::Rcout << "Y2:" << Y2 << "\n";
    
    if (BetaFlg == 1) {
      W1 = W.submat(0, 0, t, k - 1);
      W2 = W.submat(t + 1, 0, n - 1, k - 1);
      //Rcpp::Rcout << "W1:" << W1 << "\n";
      //Rcpp::Rcout << "W2:" << W2 << "\n";
    }

    //Rcpp::Rcout << "W1:" << W1 << "\n";
    //Rcpp::Rcout << "W2:" << W2 << "\n";
    
    llr = loglikRatio(Y1, Y2, loglik0, lowerLambda, upperLambda, breakPoint,
                BetaFlg, W1, W2, Beta, crit, max_p, max_d, max_q, tol, maxIter);

    llrVec(t) = llr;
    
    if (llr >= maxLlr) {
      maxLlr = llr;
      maxT = t;
    }
    
    Rcpp::Rcout << "max LLR:" << maxLlr << ", maxT:" << maxT + 1 << "\n";
    
  }
  
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("llr") = maxLlr, 
                                      Rcpp::Named("t") = maxT,
                                      Rcpp::Named("llrVec") = llrVec,
                                      Rcpp::Named("order") = model0["order"], 
                                      Rcpp::Named("coef") = model0["coef"],
                                      Rcpp::Named("sigma2") = model0["sigma2"],
                                      Rcpp::Named("lambda") = model0["lambda"],
                                      Rcpp::Named("crit") = model0["crit"],
                                      Rcpp::Named("critVal") = model0["critVal"],
                                      Rcpp::Named("loglik") = model0["loglik"]);
  
  return out;
  
}

// [[Rcpp::export]]
Rcpp::List distPars(const int& n, const Rcpp::List& model0, const double& lowerLambda, const double& upperLambda, const int& breakPoint,
                 const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, 
                 const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q, const double& tol, 
                 const int& maxIter, const Rcpp::String& innovDist, const int& nsim, const int& maxErr) {
 
  Rcpp::List out, modelSim;
  int errCode;
  
  arma::colvec order = model0["order"];
  arma::colvec order0 = order.subvec(0, 2);
  arma::colvec phi;
  arma::colvec theta;
  arma::colvec coef = model0["coef"];
  double mu;
  arma::mat muBeta;
  double sigma2 = model0["sigma2"];
  double lambda = model0["lambda"];
  
  //Rcpp::List sigMat;
  arma::mat OMat;
  //Rcpp::NumericMatrix tmpOMat;
  
  arma::colvec X;
  arma::colvec Y;
  
  arma::mat orderMat(nsim, 5);
  arma::mat coefMat(nsim, max_p + max_q + 1);
  arma::colvec sigma2Vec(nsim);
  arma::colvec lambdaVec(nsim);                               
  arma::colvec loglikVec(nsim);
  
  Rcpp::NumericVector orderSim;
  Rcpp::NumericVector coefSim;
  double sigma2Sim;
  double lambdaSim;;
  double loglikSim;
  
  int i;
  int j;
  int iter;
  int k;
  
  int idx = 0;
  
  if (order(0) > 0) {
    phi = coef.subvec(idx, idx + order(0) - 1);
    idx = idx + order(0);
  }
  
  Rcpp::Rcout << "idx:" << idx << "\n";
  
  if (order(2) > 0) {
    theta = coef.subvec(idx, idx + order(2) - 1);
    idx = idx + order(2);
    //Rcpp::Rcout << "theta:" << theta << "\n";
  }
  
  //Rcpp::Rcout << "idx:" << idx << "\n";
  
  
  if (order(3) > 0) {
    mu = coef(idx);
    //Rcpp::Rcout << "mu:" << mu << "\n";
  }

  //if (BetaFlg == 1) {
  //  muBeta = W * Beta;
  //  mu = mu + muBeta(0, 0);
  //}


  OMat = OmegaMat(n, order0, phi, theta);

  for (i = 0; i < nsim; i++) {
   
    errCode = 1;
    iter = 0;
    
    Rcpp::Rcout << "Simulation of models and parameters is in progress:" << i + 1 << " / " << nsim << "\n";
    
    while (errCode == 1 && iter <= maxErr) {
      
      iter = iter + 1;  

      try{
        
        X = simARIMA(mu, sigma2, innovDist, OMat, BetaFlg, W, Beta);
        
        if (order(4) > 0) {
          Y = invBoxCox(X, lambda);
        } else {
          Y = X;
        }
        
        //Rcpp::Rcout << "X:" << X << "\n";
        //Rcpp::Rcout << "Y:" << Y << "\n";
        //Rcpp::Rcout << "i:" << i << ", iter:" << iter << ".has_nan():" << Y.has_nan() << ", Y:" << Y << "\n";
        
        //Rcpp::Rcout << "Y.has_nan:" << Y.has_nan() << "\n";
        
        if (Y.has_nan() == 0) {
          
          modelSim = optModel(Y, lowerLambda, upperLambda, breakPoint,
                              BetaFlg, W, Beta, 
                              crit, max_p, max_d, max_q, tol, 
                              maxIter);
          
          //Rcpp::Rcout << "modelSim.length:" << modelSim.length() << "\n";
          
          //Rcpp::Rcout << "i:" << i << ", iter:" << iter << "order:" << modelSim["order"] << "\n";
          
          //Rcpp::Rcout << "Y:" << Y << ", order:" << modelSim["order"] << "\n";
          
          orderSim = modelSim["order"];
          coefSim = modelSim["coef"];
          sigma2Sim = modelSim["sigma2"];
          lambdaSim = modelSim["lambda"];
          loglikSim = modelSim["loglik"];
          
          for (j = 0; j < 5; j++) {
            orderMat(i, j) = orderSim(j);
          }
          
          k = coefSim.length();
          for (j = 0; j < k; j++) {
            coefMat(i, j) = coefSim(j);
          }
          
          sigma2Vec(i) = sigma2Sim;
          
          lambdaVec(i) = lambdaSim;
          
          loglikVec(i) = loglikSim;
          
          errCode = 0;
          
        } else {
          
          errCode = 1;
          
        }

      } catch(...) {
        
        errCode = 1;
        
      }
      
      
    }
    
  }
  
  out = Rcpp::List::create(Rcpp::Named("orderMat") = orderMat, 
        Rcpp::Named("coefMat") = coefMat,
        Rcpp::Named("sigma2Vec") = sigma2Vec,
        Rcpp::Named("lambdaVec") = lambdaVec,
        Rcpp::Named("loglikVec") = loglikVec);
  
  return out;
  //return modelSim;
}

// [[Rcpp::export]]
arma::colvec distLoglikRatio(const int& n, const int& t, const Rcpp::List& distPars0, const double& lowerLambda, const double& upperLambda, 
                        const int& breakPoint, const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, 
                        const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q, const double& tol, 
                        const int& maxIter, const Rcpp::String& innovDist, const int& nsim, const int& maxErr) {
  
  Rcpp::NumericMatrix orderMat = distPars0["orderMat"];
  Rcpp::NumericMatrix coefMat = distPars0["coefMat"];
  Rcpp::NumericVector sigma2Vec = distPars0["sigma2Vec"];
  Rcpp::NumericVector lambdaVec = distPars0["lambdaVec"];
  
  arma::colvec order0(3);
  arma::colvec phi;
  arma::colvec theta;
  double mu;
  double sigma2;
  double lambda;
  arma::mat OMat;
  
  int sta;
  int ed;
  
  arma::colvec X;
  arma::colvec Y;
  arma::colvec Y1;
  arma::colvec Y2;
  arma::mat W1;
  arma::mat W2;
  
  double loglik0;
  //double llr0;
  double llr;
  
  Rcpp::List modelSim0;
    
  int ndim = orderMat.nrow();
  int i;
  int j;
  int k;
  int a;
  int kk = 0;
  
  if (BetaFlg == 1) {
    kk = W.n_cols;
  }
  
  arma::colvec out(ndim * nsim);
  
  int r = 0;
  int iter;
  int errCode;
  
  for (i = 0; i < ndim; i++) {

    a = 0;
    
    for (a = 0; a < 3; a++) {
      
      //Rcpp::Rcout << "i:" << i << ", a: " << a << "\n";
      
      order0(a) = orderMat(i, a);
    }
    
    int idx = 0;
     
    if (orderMat(i, 0) > 0) {
      sta = idx;
      ed = idx + orderMat(i, 0);
      
      a = 0;
      
      //Rcpp::Rcout << "sta:" << sta << ", ed: " << ed << "\n";
      //Rcpp::Rcout << "coef:" << coefMat(i, 0) << "\n";
      
      phi = arma::colvec(orderMat(i, 0));
      
      for (k = sta; k < ed; k++) {
        
        //Rcpp::Rcout << "a:" << a << "\n";
        
        phi(a) = coefMat(i, k);
        a = a + 1;
      }
      
      //Rcpp::Rcout << "phi:" << phi << "\n";
      
      idx = idx + orderMat(i, 0) + 1;
    }
     
    if (orderMat(i, 2) > 0) {
      sta = idx;
      ed = idx + orderMat(i, 2);
      
      a = 0;
      
      theta = arma::colvec(orderMat(i, 2));
      
      for (k = sta; k < ed; k++) {
        theta(a) = coefMat(i, k);
        a = a + 1;
      }
      
      //Rcpp::Rcout << "theta:" << theta << "\n";
      idx = idx + orderMat(i, 2) + 1;
    }
     
    if (orderMat(i, 3) > 0) {
      mu = coefMat(i, idx);
      //Rcpp::Rcout << "mu:" << mu << "\n";
    }
     
    OMat = OmegaMat(n, order0, phi, theta);
  
    sigma2 = sigma2Vec(i);
    lambda = lambdaVec(i);
 
    for (j = 0; j < nsim; j++) {
      
      errCode = 1;
      iter = 0;
      
      Rcpp::Rcout << "Simulation of GLRs is in progress:" << r + 1 << " / " << (ndim * nsim) << "\n";
      
      while (errCode == 1 && iter <= maxErr) {
        
        try{
  
          X = simARIMA(mu, sigma2, innovDist, OMat, BetaFlg, W, Beta);
          
          if (orderMat(i, 4) > 0) {
            Y = invBoxCox(X, lambda);
          } else {
            Y = X;
          }
          
          if (Y.has_nan() == 0) {
          
            modelSim0 = optModel(Y, lowerLambda, upperLambda, breakPoint,
                   BetaFlg, W, Beta, crit, max_p, max_d, max_q, tol, 
                   maxIter);
          
            loglik0 = modelSim0["loglik"];
          
            Y1 = Y.subvec(0, t);
            Y2 = Y.subvec(t + 1, n - 1);
            
            //Rcpp::Rcout << "Y1:" << Y1 << "\n";
            //Rcpp::Rcout << "Y2:" << Y2 << "\n";
            
            if (BetaFlg == 1) {
              W1 = W.submat(0, 0, t, kk - 1);
              W2 = W.submat(t + 1, 0, n - 1, kk - 1);
              //Rcpp::Rcout << "W1:" << W1  << "\n";
              //Rcpp::Rcout << "W2:" << W2  << "\n";
            }
            
            llr = loglikRatio(Y1, Y2, loglik0, lowerLambda, upperLambda, breakPoint,
                              BetaFlg, W1, W2, Beta, crit, max_p, max_d, max_q, tol, maxIter);
            
            out(r) = llr;
            
            errCode = 0;
            
          } else {
            
            errCode = 1;
            
          }
          
        } catch(...) {
          
          errCode = 1;
          
        }
        
      }
      
      r = r + 1; 
      
    }
  }
  
  //return phi;
  return out;
  
}


// [[Rcpp::export]]
Rcpp::NumericVector returnCritValKDECpp(const arma::colvec& x, const double& alpha) {
  Rcpp::Environment rEnv = Rcpp::Environment::global_env();
  Rcpp::Function rfun = rEnv["returnCritValKDE"];
  Rcpp::NumericVector out = rfun(Rcpp::Named("x", x), Rcpp::Named("alpha", alpha));
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector quantileCpp(const arma::colvec& x, const double& alpha) {
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("stats"); 
  Rcpp::Function rfunction = pkg["quantile"];  
  Rcpp::NumericVector out = rfunction(Rcpp::Named("x", x), Rcpp::Named("probs", 1 - alpha));
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector idxFinder(const Rcpp::LogicalVector& x){
  
  int n = x.length();
  Rcpp::NumericVector out(n);
  double tmp = 0;
  
  for (int i = 0; i < n; i++) {
    out(i) = tmp;
    tmp = tmp + 1;
  }
  
  out = out[x];
  
  return out;
  
}



// [[Rcpp::export]]
Rcpp::List binSeg(const arma::colvec& Y, const double & alpha, const int& GLRSApprox, const int& minSize,
                  const double& lowerLambda, const double& upperLambda, const int& breakPoint,
                  const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, 
                  const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q, const double& tol, 
                  const int& maxIter, const Rcpp::String& innovDist, const int& nsim1, const int& nsim2, const int& maxErr) {
  
  Rcpp::List out;
  
  int n = Y.n_elem;
  
  Rcpp::NumericVector tmpY;
  tmpY = Rcpp::rep(0, n);
  int i = 0;
  
  for (i = 0; i < n; i++) {
    tmpY[i] = Y(i);
  }

  Rcpp::NumericVector AvailCheckVec;
  AvailCheckVec = Rcpp::rep(1, n);
  
  Rcpp::NumericVector GroupVec;
  GroupVec = Rcpp::rep(0, n);
  
  Rcpp::NumericVector checkGroupVec;
  Rcpp::NumericVector checkGroupCntVec;
  
  int nCheckGroupVec;
  int workGroup;
  Rcpp::NumericVector workY;
  int nWorkY;

  Rcpp::LogicalVector workIdx;
  Rcpp::NumericVector workIdxFinder;
  Rcpp::NumericVector llrVec;
    
  int k = W.n_cols;
  Rcpp::NumericMatrix tmpW(n, k);
  //Rcpp::NumericMatrix workW(n, k);
  
  
  //Rcpp::Rcout << "k:" << k << "\n";  
  //Rcpp::Rcout << "is_vec:" << W.is_vec() << "\n";  
  int j = 0;
  
  if (BetaFlg == 1) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < k; j++) {
        
        tmpW(i, j) = W(i, j);
        
      }
    }
  }
  
  //Rcpp::Rcout << "tmpW:" << tmpW << "\n";

  Rcpp::List step1;
  Rcpp::List step2;
  Rcpp::NumericVector step3;
  Rcpp::NumericVector step4;
  double criticalVal = 0;
  double maxllr = 0;
  
  int minUpdateIdx, maxUpdateIdx;
  
  int level = 0;
  
  int flg = 0;
  int g = 0;
  //double tmpT = 0;
  int t = 0;
  int Finder = 0;
  
  int groupIdx = 0;
  
  while(flg == 0) {
      
    level = level + 1;
    
    checkGroupVec = GroupVec[AvailCheckVec == 1];
    checkGroupCntVec = Rcpp::table(checkGroupVec);
    checkGroupVec = Rcpp::sort_unique(checkGroupVec);
    checkGroupVec = checkGroupVec[checkGroupCntVec >= 2 * minSize];

    nCheckGroupVec = checkGroupVec.length();
    
    Rcpp::Rcout << "checkGroupCntVec:" << checkGroupCntVec << "\n";
    Rcpp::Rcout << "checkGroupVec:" << checkGroupVec << "\n";
    Rcpp::Rcout << "nCheckGroupVec:" << nCheckGroupVec << "\n";
    
    if (nCheckGroupVec > 0) {
    
      for (g = 0; g < nCheckGroupVec; g++) {
        
        
        
        workGroup = checkGroupVec[g];
        //Rcpp::Rcout << "workGroup:" << workGroup << "\n";
        
        workIdx = GroupVec == workGroup;
        workIdxFinder = idxFinder(workIdx);
        //Rcpp::Rcout << "workIdx:" << workIdx << "\n";
        Rcpp::Rcout << "workIdxFinder:" << workIdxFinder << "\n";
        
        workY = tmpY[workIdx];
        //Rcpp::Rcout << "workY:" << workY << "\n";
        nWorkY = workY.length();
        //Rcpp::Rcout << "nWorkY:" << nWorkY << "\n";
        
        arma::colvec tmpworkY(nWorkY, arma::fill::zeros);
        //Rcpp::Rcout << "tmpworkY:" << tmpworkY << "\n";
        for (i = 0; i < nWorkY; i++) {
          tmpworkY(i) = workY(i);
        }
        //Rcpp::Rcout << "tmpworkY:" << tmpworkY << "\n";
        
        arma::mat tmpworkW(nWorkY, k, arma::fill::zeros);
        if (BetaFlg == 1) {
          for (i = 0; i < nWorkY; i++) {
            Finder = workIdxFinder(i);
            for (j = 0; j < k; j++) {
              tmpworkW(i, j) = tmpW(Finder, j);
            }
          }
        }

        //Rcpp::Rcout << "tmpworkW:" << tmpworkW << "\n";
        
        ////Step1: check the max of GLRs
        step1 = loglikRatioMax(tmpworkY, minSize, lowerLambda, upperLambda, breakPoint, 
                                BetaFlg, tmpworkW, Beta, crit, max_p, max_d, max_q, tol, maxIter);
        //tmpT = step1["t"];
        t = Rcpp::as<int>(step1["t"]);
        maxllr = step1["llr"];
        llrVec = step1["llrVec"];
        
        //Rcpp::Rcout << "t:" << t << ", maxllr:" << maxllr << "\n";
        //Rcpp::Rcout << "step1:" << step1 << "\n";
        
        ////Step2: simulate the distribution of models and parameters
        step2 = distPars(nWorkY, step1, lowerLambda, upperLambda, breakPoint, 
                          BetaFlg, tmpworkW, Beta, crit, max_p, max_d, max_q, tol, maxIter, "norm", nsim1, maxErr);
        //Rcpp::Rcout << "step2:" << step2 << "\n";
        
        ////Step3: simulate the distribution of GLRs
        step3 = distLoglikRatio(nWorkY, t, step2, lowerLambda, upperLambda, breakPoint, 
                                 BetaFlg, tmpworkW, Beta, crit, max_p, max_d, max_q, tol, maxIter, "norm", nsim2, maxErr);
        //Rcpp::Rcout << "step3:" << step3 << "\n";
        
        ////Step4: check the max of GLRs
        if (GLRSApprox == 1) {
          step4 = returnCritValKDECpp(step3, alpha);
          //Rcpp::Rcout << "step4:" << step4 << "\n";
        } else {
          step4 = quantileCpp(step3, alpha);
          //Rcpp::Rcout << "step4:" << step4 << "\n";
        }
        
        
        //Rcpp::Rcout << "criticalVal:" << criticalVal << "\n";
        
        // update group vector and availability vector
        criticalVal = step4(0);
        minUpdateIdx = Rcpp::min(workIdxFinder);
        maxUpdateIdx = minUpdateIdx + t;
        
        if (maxllr >= criticalVal) {
          
          //GroupVec[workIdx] = groupIdx;
          
          groupIdx = groupIdx + 1;
          for (i = minUpdateIdx; i <= maxUpdateIdx; i++) {
            GroupVec(i) = groupIdx;
          }
          
          //Rcpp::Rcout << "GroupVec:" << GroupVec << "\n";
          
        } else {
          
          for (i = minUpdateIdx; i <= Rcpp::max(workIdxFinder); i++) {
            AvailCheckVec(i) = 0;
          }
          
        }
        
        Rcpp::Rcout << "AvailCheckVec:" << AvailCheckVec << "\n";
        Rcpp::Rcout << "GroupVec:" << GroupVec << "\n";

      }
      
    } else {
      
      AvailCheckVec = Rcpp::rep(0, n);
      flg = 1;
      
    }
      
  }
    
  return out = Rcpp::List::create(Rcpp::Named("GroupVec") = GroupVec);

}


//double loglik(const arma::colvec& Y, const arma::mat& SigMat, const double& mu, const double& sigma2, 
//            const int& BoxCoxFlg, const double& lambda, const int& betaFlg, const arma::mat& W, const arma::colvec& Beta) {
//  
//  double pi = 3.141592653589;
//  
//  double out;
//  arma::mat out1;
//  //double out;
//  arma::mat loglik1;
//  double loglik2;
//  
//  arma::colvec X;
//  int n = Y.n_elem;
//    
//  arma::mat SigMat1 = SigMat * sigma2;
//    
//  X = Y;
//    
//  if (BoxCoxFlg == 0) {
//    X = X - mu;
//    loglik2 = 0;
//  } else {
//    X = BoxCox(X, lambda) - mu;
//    loglik2 = arma::prod(BoxCox1stDeriv(Y, lambda));
//  }
//  
//  if (betaFlg == 1) {
//    X = X - W * Beta;
//  }
//  
//  loglik1 = (std::log(2 * pi) * n / 2 + 0.5 * std::log(arma::det(SigMat1)) + 0.5 * X.t() * arma::inv(SigMat1) * X) * (-1);
//  
//  out1 = loglik1 + loglik2;
//  
//  return out = out1(0, 0);
//  
//}









/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


