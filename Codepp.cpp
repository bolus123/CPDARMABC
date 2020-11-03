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
    Rcpp::stop("The parameter vector is not invertible/stationary");
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


arma::colvec simInnov(const int& n, Rcpp::String& XSim, Rcpp::NumericVector& XPars) {
  
  arma::colvec out(n);
  
  if (XSim == "norm") {
    out.randn();
    out = out * std::sqrt(XPars[1]) + XPars[0];
  }
  
  return out;
  
}

// [[Rcpp::export]]
arma::colvec simARIMA(const int& n, const double& mu, const double& sigma2, Rcpp::String& innovDist, arma::mat& OmegaMat) {
  
  Rcpp::NumericVector XPars = {0, 1};
  
  arma::colvec out = OmegaMat * (simInnov(n, innovDist, XPars) * std::sqrt(sigma2) + mu);

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

// [[Rcpp::export]]
Rcpp::List estCSS(const arma::colvec& X, const arma::colvec& order, const int& muFlg) {
  
  int n = X.n_elem;
  arma::colvec tmpX = X;
  arma::colvec tmpX1;
  arma::mat tmpA;
  arma::mat tmpA1(n - order(1), order(0));
  arma::mat tmpB;
  arma::mat tmpB1(n - order(0) - order(1), order(0));
  arma::colvec phi;
  arma::colvec theta;
  arma::colvec tmpEps;
  double mu = 0;
  arma::colvec sigma2;
  
  tmpEps = tmpX;
  
  if (muFlg == 1) {
    mu = arma::mean(tmpX);
    tmpX = tmpX - mu;
    tmpEps = tmpX;
  }
  
  if (order(1) > 0) {
    tmpX = DifMat(n, order) * tmpX;
    tmpEps = tmpX;
  } 
  
  if (order(0) > 0) {
    for (int j = 0; j < order(0); j++) {
      for (int i = 0; i < n - order(1); i++) {
        if (i - order(0) >= 0) {
          tmpA1(i, j) = tmpX(i - j - 1);
        }
      }
    }
    
    tmpA = tmpA1.submat(order(0), 0, n - order(1) - 1, order(0) - 1);
    tmpX1 = tmpX.subvec(order(0), n - order(1) - 1);
    phi = arma::solve(tmpA, tmpX1);
    tmpX = tmpX1 - tmpA * phi;
    tmpEps = tmpX;
  }
  
  if (order(2) > 0) {
    for (int j = 0; j < order(2); j++) {
      for (int i = 0; i < n - order(1) - order(0); i++) {
        if (i - order(2) >= 0) {
          tmpB1(i, j) = tmpX(i - j - 1);
        }
      }
    }
    
    tmpB = tmpB1.submat(order(2), 0, n - order(1) - order(0) - 1, order(2) - 1);
    tmpX1 = tmpX.subvec(order(2), n - order(1) - order(0) - 1);
    theta = arma::solve(tmpB, tmpX1);
    tmpEps = tmpX1 - tmpB * theta;
  }
  
  //sigma2 = arma::mean(arma::pow(tmpEps, 2));
  sigma2 = arma::sum(arma::pow(tmpEps, 2)) / tmpEps.n_elem;
  
  return Rcpp::List::create(Rcpp::Named("phi") = phi, 
                            Rcpp::Named("theta") = theta,
                            Rcpp::Named("mu") = mu, 
                            Rcpp::Named("sigma2") = sigma2, 
                            Rcpp::Named("eps") = tmpEps);
  
}

// [[Rcpp::export]]
Rcpp::List arimaCpp(const arma::colvec& x, const arma::colvec& order, const Rcpp::LogicalVector& include_mean) {
  Rcpp::Function arima_("arima");
  Rcpp::List out = arima_(Rcpp::Named("x", x),
                          Rcpp::Named("order", order),
                          Rcpp::Named("include.mean", include_mean));
  return(out);
}

// [[Rcpp::export]]
double loglik(const arma::colvec& Y, const arma::colvec& order, const Rcpp::LogicalVector& include_mean) {

  double out;
  arma::mat out1;
  //double out;
  arma::mat loglik1;
  double loglik2;
  
  arma::colvec X;
  int n = Y.n_elem;
    
  X = Y;
    
  if (BoxCoxFlg == 0) {
    X = X - mu;
    loglik2 = 0;
  } else {
    X = BoxCox(X, lambda) - mu;
    loglik2 = arma::prod(BoxCox1stDeriv(Y, lambda));
  }
  
  if (betaFlg == 1) {
    X = X - W * Beta;
  }  
  
  
  
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



// [[Rcpp::export]]
double AIC(const double& loglik, const int& npars) {
  double out = loglik * (-2) + 2 * (npars);
  return out;
}





/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


