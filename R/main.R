YJTrans <- function(Y, lambda) {
  
  fun1 <- function(Y, lambda) {
    sign(Y) * ((abs(Y) + 1) ^ (sign(Y) * (lambda - 1) + 1) - 1) /
      (sign(Y) * (lambda - 1) + 1)
  }
  
  fun2 <- function(Y) {
    sign(Y) * log(abs(Y) + 1)
  }
  
  if (lambda != 0 & lambda != 2) {
    out <- fun1(Y, lambda)
  } else {
    if (lambda == 0) {
      out <- ifelse(Y >= 0, fun2(Y), fun1(Y, lambda))
    } else if (lambda == 2) {
      out <- ifelse(Y < 0, fun2(Y), fun1(Y, lambda))
    }
    
  }
  
  out
  
}


loglikFunc <- function(Y, order, pars, meanFlg = TRUE, YJFlg = TRUE) {
  
  X <- Y
  
  mu <- 0
  pos <- 1
  if (YJFlg == TRUE) {
    lambda <- pars[pos]
    X <- YJTrans(X, lambda)
    pos <- pos + 1
  }
  
  if (meanFlg == TRUE) {
    mu <- pars[pos]
    X <- X - mu
    pos <- pos + 1
  }
  
  sigma2 <- pars[pos]
  pos <- pos + 1
  
  if (order[1] > 0) {
    psiVec <- pars[pos:(pos + order[1] - 1)]
    pos <- pos + order[1]
  }
  
  if (order[3] > 0) {
    thetaVec <- pars[pos:(pos + order[3] - 1)]
    pos <- pos + order[3]
  }
  
  
  if (order[2] > 0) {
    for (i in 1:order[2]) {
      X <- diff(X, 1)
    }
  }
  
  n <- length(X)
  eps <- rep(NA, n)
  
  l1 <- 0
  if (YJFlg == TRUE) {
    l1 <- (lambda - 1) * sum(sign(Y) * log(abs(Y) + 1))
  }
  
  if (order[1] == 0 & order[3] == 0) {
    eps <- X
  } else {
    for (i in 1:n) {
      if (order[1] > 0) {
        tmpX <- rep(0, order[1])
        selectX <- (i - 1):(i - order[1])
        whichX <- which(selectX > 0)
        lengthX <- length(whichX)
        if (lengthX > 0) {
          tmpX[1:lengthX] <- X[(i - 1):(i - lengthX)]
        }
        eps[i] <- X[i] - psiVec %*% tmpX
      } else {
        eps <- X
      }
      
      if (order[3] > 0) {
        tmpEps <- rep(0, order[3])
        selectEps <- (i - 1):(i - order[3])
        whichEps <- which(selectEps > 0)
        lengthEps <- length(whichEps)
        if (length(whichEps) > 0) {
          tmpEps[1:lengthEps] <- eps[(i - 1):(i - lengthEps)]
        }
        eps[i] <- eps[i] - thetaVec %*% tmpEps
      }
      
    }
  }
  
  
  
  l2 <- sum(log(dnorm(eps, sd = sqrt(sigma2))))
  
  loglik <- l1 + l2
  
  list(loglik = loglik, residuals = eps)
  
  #loglik
  
}


getEsts <- function(Y, order, pars0, meanFlg = TRUE, YJFlg = TRUE) {
  
  fnGradDesc <- function(Y, order, pars, meanFlg, YJFlg) {
    -loglikFunc(Y, order, pars, meanFlg, YJFlg)$loglik
  }
  
  if (length(pars0) > 1) {
    out <- optim(pars0, fn = fnGradDesc, Y = Y, order = order, 
                 meanFlg = meanFlg, YJFlg = YJFlg, method = "Nelder-Mead")
  } else {
    out <- optim(pars0, fn = fnGradDesc, Y = Y, order = order, 
                 meanFlg = meanFlg, YJFlg = YJFlg, method = "Brent")
  }
  
  
  
  list(Ests = out$par, loglik = -out$value)
  
}

AIC <- function(loglik, npars) {
  ifelse(is.finite(loglik), loglik * (-2) + 2 * (npars), NA)
}

BIC <- function(loglik, npars, nobs) {
  ifelse(is.finite(loglik), npars * log(nobs) - 2 * loglik, NA)
}

AICc <- function(loglik, npars, nobs) {
  ifelse(is.finite(loglik), AIC(loglik, npars) + (2 * npars ^ 2 + 2 * npars) / (nobs - npars - 1), NA)
}

ARMACheck <- function (pars = NULL) {
  if (length(pars) > 0) {
    phi = c(1, pars)
    out <- ifelse(all(abs(polyroot(phi)) > 1), 1, 0)
  } else {
    out <- 1
  }
  out
}


getBestModel <- function(Y, halfWinLength, maxOrder, method = 'bic', searchYJ = TRUE, lambda0 = 1, searchM1M2 = TRUE) {
  
  if (searchYJ == TRUE) {
    orderMat <- expand.grid(0:maxOrder[1], 0:maxOrder[2], 0:maxOrder[3], 
                            c(FALSE, TRUE), c(FALSE, TRUE))
  } else {
    orderMat <- expand.grid(0:maxOrder[1], 0:maxOrder[2], 0:maxOrder[3], 
                            c(FALSE, TRUE), c(FALSE))
  }
  
  orderMat <- orderMat[-which(orderMat[, 2] > 0 & orderMat[, 4] == TRUE), ]
  
  k <- dim(orderMat)[1]
  
  n <- length(Y)
  
  llVec <- rep(-Inf, k)
  
  estsList <- list()
  
  OrderMat <- matrix(0, nrow = k, ncol = 3)
  
  meanFlgVec <- rep(NA, k)
  YJFlgVec <- rep(NA, k)
  
  critVec <- rep(Inf, k)
  minCrit <- Inf
  minR <- 0
  
  ARCheck <- 1
  MACheck <- 1
  
  for (r in 1:k) {
    #cat('r:', r, '\n')
    tmpOrder <- as.numeric(orderMat[r, 1:3])
    tmpMeanFlg <- orderMat[r, 4]
    tmpYJFlg <- orderMat[r, 5]
    
    npars <- tmpOrder[1] + tmpOrder[3] + tmpMeanFlg + tmpYJFlg + 1
    tmpPars <- rep(NA, npars)
    
    #nobs <- length(Y) - tmpOrder[2]
    
    tmpObs <- length(Y) - tmpOrder[2]
    nobs <- ifelse(searchM1M2 == TRUE, length(Y) - tmpOrder[2], halfWinLength - tmpOrder[2])
    
    if (nobs > npars) {
      #if (halfWinLength > npars) {
      
      if (tmpYJFlg == TRUE) {
        tmpModel <- try(arima(Y, order = tmpOrder, 
                              include.mean = tmpMeanFlg, method = "CSS"),
                        silent = TRUE)
      } else {
        tmpModel <- try(arima(Y, order = tmpOrder, 
                              include.mean = tmpMeanFlg, method = "CSS-ML"),
                        silent = TRUE)
        
      }
      
      if (class(tmpModel) != 'try-error') {
        XResi <- tmpModel$residuals
        pos <- 1
        if (tmpYJFlg == TRUE) {
          tmpPars[pos] <- lambda0
          pos <- pos + 1
        }
        
        if (tmpMeanFlg == TRUE) {
          tmpPars[pos] <- tmpModel$coef[tmpOrder[1] + tmpOrder[3] + 1]
          pos <- pos + 1
        }
        
        tmpPars[pos] <- tmpModel$sigma2
        pos <- pos + 1
        
        if (tmpOrder[1] > 0) {
          tmpPars[pos:(pos + tmpOrder[1] - 1)] <- tmpModel$coef[1:tmpOrder[1]]
          ARCheck <- ARMACheck(pars = -tmpModel$coef[1:tmpOrder[1]])
          pos <- pos + tmpOrder[1]
        }
        
        if (tmpOrder[3] > 0) {
          tmpPars[pos:(pos + tmpOrder[3] - 1)] <- tmpModel$coef[(tmpOrder[1] + 1):(tmpOrder[1] + tmpOrder[3])]
          MACheck <- ARMACheck(pars = tmpModel$coef[(tmpOrder[1] + 1):(tmpOrder[1] + tmpOrder[3])])
          pos <- pos + tmpOrder[3]
        }
        
        if (tmpYJFlg == TRUE) {
          model <- try(getEsts(Y = Y, order = tmpOrder, pars0 = tmpPars, 
                               meanFlg = tmpMeanFlg, YJFlg = tmpYJFlg), silent = TRUE)
          
          if (class(model) != 'try-error') {
            
            pos <- 1
            if (tmpYJFlg == TRUE) {
              pos <- pos + 1
            }
            
            if (tmpMeanFlg == TRUE) {
              pos <- pos + 1
            }
            
            pos <- pos + 1
            
            if (tmpOrder[1] > 0) {
              tmpAR <- model$Ests[pos:(pos + tmpOrder[1] - 1)]
              #tmpPars[pos:(pos + tmpOrder[1] - 1)] <- tmpModel$coef[1:tmpOrder[1]]
              ARCheck <- ARMACheck(pars = -tmpAR)
              pos <- pos + tmpOrder[1]
            }
            
            if (tmpOrder[3] > 0) {
              tmpMA <- model$Ests[pos:(pos + tmpOrder[3] - 1)]
              #tmpPars[pos:(pos + tmpOrder[3] - 1)] <- tmpModel$coef[(tmpOrder[1] + 1):(tmpOrder[1] + tmpOrder[3])]
              MACheck <- ARMACheck(pars = tmpMA)
              pos <- pos + tmpOrder[3]
            }
            
            
            
            estsList[[r]] <- model$Ests
            
            ##
            ## control YJ parameter must between 0 and 2
            ##
            
            if (0 <= model$Ests[1] & model$Ests[1] <= 2) {
              
              ##
              ## Both casuality and invertibility need to be satisfied
              ##
              
              if (ARCheck == 1 & MACheck == 1) {
                llVec[r] <- model$loglik
                XResi <- loglikFunc(Y, tmpOrder, model$Ests, meanFlg = tmpMeanFlg, YJFlg = tmpYJFlg)$residuals
              }
            } 
            
          }
          
          
          
        } else {
          ##
          ## Both casuality and invertibility need to be satisfied
          ##
          
          estsList[[r]] <- tmpPars
          
          if (ARCheck == 1 & MACheck == 1) {
            llVec[r] <- tmpModel$loglik
          }
          
          
        }
        
        OrderMat[r, ] <- tmpOrder
        meanFlgVec[r] <- tmpMeanFlg
        YJFlgVec[r] <- tmpYJFlg
        
        if (method == 'aic') {
          critVec[r] <- AIC(llVec[r], npars)
        } else if (method == 'bic') {
          critVec[r] <- BIC(llVec[r], npars, tmpObs)
        } else if (method == 'aicc') {
          critVec[r] <- AICc(llVec[r], npars, tmpObs)
        }
        
        if (!is.na(critVec[r])) {
          if (critVec[r] < minCrit) {
            minCrit <- critVec[r]
            minR <- r
            minXResi <- XResi
            
          }
        }
        
      }
    }
    
    
  }
  
  ##if (YJFlgVec[minR] == TRUE) {
  ##  
  ##  minlambda <- estsList[[minR]][1]
  ##  X <- YJTrans(Y, minlambda)
  ##  
  ##  if ()
  ##  
  ##  Xfitted <- X - minXResi
  ##  
  ##  if (meanFlgVec[minR] == TRUE) {
  ##    minSigma2 <- estsList[[minR]][3]
  ##  } else {
  ##    minSigma2 <- estsList[[minR]][2]
  ##  }
  ##
  ##  fitted <- InvYJTrans(Xfitted, minlambda) + InvYJTrans2ndDerv(Xfitted, minlambda) / 2 * minSigma2
  ##  
  ##} else {
  ##  
  ##  fitted <- Y - minXResi
  ##  
  ##}
  
  
  list(loglik = llVec[minR], order = OrderMat[minR, ], 
       ests = estsList[[minR]], meanFlg = meanFlgVec[minR], 
       YJFlg = YJFlgVec[minR], XResi = minXResi,
       critVec = critVec, llVec = llVec, 
       OrderMat = OrderMat, estsList = estsList,
       meanFlgVec = meanFlgVec, YJFlgVec = YJFlgVec)
  
  #list(llVec = llVec, estsList = estsList, meanFlgVec = meanFlgVec, YJFlgVec = YJFlgVec)
  
}

getFitted <- function(Y, XResi, order, ests, YJFlg = TRUE, meanFlg = TRUE) {
  
  n <- length(Y)
  
  if (YJFlg == TRUE) {
    lambda <- ests[1]
    X <- YJTrans(Y, lambda)
    if (order[2] > 0) {
      D <- matrix(0, nrow = order[2] + 1, ncol = n)
      D[1, ] <- X
      for (i in 1:order[2]) {
        D[i + 1, (1): (n - i)] <- D[i, 2:(n - i + 1)] - D[i, 1:(n - i)]
      }
      X <- D[i + 1, -c((n - i + 1):n)]
    }
    Xfitted <- X - XResi
    
    if (order[2] > 0) {
      D[order[2] + 1, 1:(n - order[2])] <- Xfitted
      for (i in order[2]:1) {
        D[i, 1:(n - i + 1)] <- colSums(D[(i + 1):i, 1:(n - i + 1)])
        #D[order[2] - i + 1, (i + 1):n] <- D[order[2] - i + 1, c((order[2] - i + 1):(n - i))] + D[order[2] - i + 2, c(order[2] - i + 2)]
      }
      Xfitted <- D[1, ]
    }
    if (meanFlg == TRUE) {
      sigma2 <- ests[3]
    } else {
      sigma2 <- ests[2]
    }
    fitted <-  biasCorrect(Xfitted, lambda, sigma2)#InvYJTrans(Xfitted, lambda) + InvYJTrans2ndDerv(Xfitted, lambda) / 2 * sigma2
  } else {
    X <- Y
    fitted <- X - XResi
  }
  
  
  fitted
  
}

getCI <- function(alpha, Y, XResi, order, ests, YJFlg = TRUE, meanFlg = TRUE) {
  
  n <- length(Y)
  
  if (YJFlg == TRUE) {
    lambda <- ests[1]
    X <- YJTrans(Y, lambda)
    if (order[2] > 0) {
      D <- matrix(0, nrow = order[2] + 1, ncol = n)
      D[1, ] <- X
      for (i in 1:order[2]) {
        D[i + 1, (1): (n - i)] <- D[i, 2:(n - i + 1)] - D[i, 1:(n - i)]
      }
      X <- D[i + 1, -c((n - i + 1):n)]
    }
    Xfitted <- X - XResi
    
    if (order[2] > 0) {
      D[order[2] + 1, 1:(n - order[2])] <- Xfitted
      for (i in order[2]:1) {
        D[i, 1:(n - i + 1)] <- colSums(D[(i + 1):i, 1:(n - i + 1)])
        #D[order[2] - i + 1, (i + 1):n] <- D[order[2] - i + 1, c((order[2] - i + 1):(n - i))] + D[order[2] - i + 2, c(order[2] - i + 2)]
      }
      Xfitted <- D[1, ]
    }
    if (meanFlg == TRUE) {
      sigma2 <- ests[3]
    } else {
      sigma2 <- ests[2]
    }
    
    lower <- Xfitted - qnorm(1 - alpha / 2) * sqrt(sigma2)
    upper <- Xfitted + qnorm(1 - alpha / 2) * sqrt(sigma2)
    lower <-  biasCorrect(lower, lambda, sigma2)#InvYJTrans(Xfitted, lambda) + InvYJTrans2ndDerv(Xfitted, lambda) / 2 * sigma2
    upper <-  biasCorrect(upper, lambda, sigma2)
  } else {
    X <- Y
    fitted <- X - XResi
    
    if (meanFlg == TRUE) {
      sigma2 <- ests[2]
    } else {
      sigma2 <- ests[1]
    }
    
    lower <- Xfitted - qnorm(1 - alpha / 2) * sqrt(sigma2)
    upper <- Xfitted + qnorm(1 - alpha / 2) * sqrt(sigma2)
  }
  
  
  list(lower = lower, upper = upper)
  
}

biasCorrect <- function(mu, lambda, sigma2) {
  
  InvYJTrans(mu, lambda) + InvYJTrans2ndDerv(mu, lambda) / 2 * sigma2
  
}


generalizedloglikRatio <- function(Y1, Y2, halfWinLength, maxOrder, method = 'bic', searchYJ = TRUE, lambda0 = 1, searchM1M2 = TRUE) {
  Y0 <- c(Y1, Y2)
  
  model0 <- getBestModel(Y0, halfWinLength, maxOrder, method, searchYJ, lambda0, searchM1M2)
  
  loglik0 <- model0$loglik
  
  loglik1 <- NA
  loglik2 <- NA
  
  if (searchM1M2 == TRUE) {
    model1 <- getBestModel(Y1, halfWinLength, maxOrder, method, searchYJ, lambda0, searchM1M2)
    model2 <- getBestModel(Y2, halfWinLength, maxOrder, method, searchYJ, lambda0, searchM1M2)
  } else {
    
    tmpPars <- vector()
    tmpOrder <- model0$order
    tmpYJFlg <- model0$YJFlg
    tmpMeanFlg <- model0$meanFlg
    
    if (tmpYJFlg == TRUE) {
      tmpModel1 <- try(arima(Y1, order = tmpOrder, 
                             include.mean = tmpMeanFlg, method = "CSS"),
                       silent = TRUE)
      
      if (class(tmpModel1) != 'try-error') {
        pos <- 1
        if (tmpYJFlg == TRUE) {
          tmpPars[pos] <- lambda0
          pos <- pos + 1
        }
        
        if (tmpMeanFlg == TRUE) {
          tmpPars[pos] <- tmpModel1$coef[tmpOrder[1] + tmpOrder[3] + 1]
          pos <- pos + 1
        }
        
        tmpPars[pos] <- tmpModel1$sigma2
        pos <- pos + 1
        
        if (tmpOrder[1] > 0) {
          tmpPars[pos:(pos + tmpOrder[1] - 1)] <- tmpModel1$coef[1:tmpOrder[1]]
          ARCheck <- ARMACheck(pars = -tmpModel1$coef[1:tmpOrder[1]])
          pos <- pos + tmpOrder[1]
        }
        
        if (tmpOrder[3] > 0) {
          tmpPars[pos:(pos + tmpOrder[3] - 1)] <- tmpModel1$coef[(tmpOrder[1] + 1):(tmpOrder[1] + tmpOrder[3])]
          MACheck <- ARMACheck(pars = tmpModel1$coef[(tmpOrder[1] + 1):(tmpOrder[1] + tmpOrder[3])])
          pos <- pos + tmpOrder[3]
        }
        
        
        model1 <- try(getEsts(Y = Y1, order = tmpOrder, pars0 = tmpPars, 
                              meanFlg = tmpMeanFlg, YJFlg = tmpYJFlg), silent = TRUE)
        
        
      } else {
        model1 <- tmpModel1
      }
      
    } else {
      tmpModel1 <- try(arima(Y1, order = tmpOrder, 
                             include.mean = tmpMeanFlg, method = "CSS-ML"),
                       silent = TRUE)
      
      model1 <- tmpModel1
      
    }
    
    
    #######################################################
    
    tmpPars <- vector()
    
    if (tmpYJFlg == TRUE) {
      tmpModel2 <- try(arima(Y2, order = tmpOrder, 
                             include.mean = tmpMeanFlg, method = "CSS"),
                       silent = TRUE)
      
      if (class(tmpModel2) != 'try-error') {
        pos <- 1
        if (tmpYJFlg == TRUE) {
          tmpPars[pos] <- lambda0
          pos <- pos + 1
        }
        
        if (tmpMeanFlg == TRUE) {
          tmpPars[pos] <- tmpModel2$coef[tmpOrder[1] + tmpOrder[3] + 1]
          pos <- pos + 1
        }
        
        tmpPars[pos] <- tmpModel2$sigma2
        pos <- pos + 1
        
        if (tmpOrder[1] > 0) {
          tmpPars[pos:(pos + tmpOrder[1] - 1)] <- tmpModel2$coef[1:tmpOrder[1]]
          ARCheck <- ARMACheck(pars = -tmpModel2$coef[1:tmpOrder[1]])
          pos <- pos + tmpOrder[1]
        }
        
        if (tmpOrder[3] > 0) {
          tmpPars[pos:(pos + tmpOrder[3] - 1)] <- tmpModel2$coef[(tmpOrder[1] + 1):(tmpOrder[1] + tmpOrder[3])]
          MACheck <- ARMACheck(pars = tmpModel2$coef[(tmpOrder[1] + 1):(tmpOrder[1] + tmpOrder[3])])
          pos <- pos + tmpOrder[3]
        }
        
        
        model2 <- try(getEsts(Y = Y2, order = tmpOrder, pars0 = tmpPars, 
                              meanFlg = tmpMeanFlg, YJFlg = tmpYJFlg), silent = TRUE)
        
      } else {
        model2 <- tmpModel2
      }
      
    } else {
      tmpModel2 <- try(arima(Y2, order = tmpOrder, 
                             include.mean = tmpMeanFlg, method = "CSS-ML"),
                       silent = TRUE)
      model2 <- tmpModel2
    }
    
    
    
  }
  
  loglik1 <- ifelse(class(model1) != 'try-error', model1$loglik, NA)
  loglik2 <- ifelse(class(model2) != 'try-error', model2$loglik, NA)
  
  
  glr <- -2 * (loglik0 - (loglik1 + loglik2))
  
  cat('GLR:', glr, ', ll0:', loglik0, ', ll1:', loglik1, ', ll2:', loglik2, '\n')
  
  list(GLR = glr, model0 = model0, model1 = model1, model2 = model2)
}


SldWin <- function(Y, halfWinLength, maxOrder, method = 'bic', searchYJ = TRUE, lambda0 = 1, searchM1M2 = TRUE) {
  n <- length(Y)
  glrVec <- rep(NA, n - halfWinLength * 2 + 1)
  model0 <- list()
  model1 <- list()
  model2 <- list()
  
  t <- 0
  for (j in halfWinLength:(n - halfWinLength)) {
    t <- t + 1
    cat('t:', j + 1, '\n')
    Y1 <- Y[t:j]
    Y2 <- Y[(j + 1):(j + halfWinLength)]
    models <- suppressWarnings(
      generalizedloglikRatio(Y1, Y2, halfWinLength, maxOrder, method, searchYJ, lambda0, searchM1M2))
    glrVec[t] <- models$GLR
    model0[[t]] <- models$model0
    model1[[t]] <- models$model1
    model2[[t]] <- models$model2
  }
  list(GLRVec = glrVec, model0 = model0, model1 = model1, model2 = model2)
}

#ScanGLRVec.old <- function(SldWin, n, halfWinLength) {
#  tmpGLRVec <- rep(-Inf, n)
#  tmpGLRVec[(halfWinLength + 1):(n - halfWinLength + 1)] <- SldWin$GLRVec
#  GLRVecDiff1 <- sign(diff(tmpGLRVec))
#  GLRVecDiff2 <- diff(GLRVecDiff1)
#  Idx <- (which(GLRVecDiff2 == -2) + 1) - halfWinLength
#  tmp <- SldWin$GLRVec[Idx]
#  ord <- order(tmp, decreasing = TRUE)
#  Idx[ord]
#}

ScanGLRVec <- function(SldWin, n, halfWinLength, minpeakdistance = halfWinLength, minpeakheight = -.Machine$double.xmax) {
  tmpGLRVec <- SldWin$GLRVec
  tmpGLRVec[is.na(SldWin$GLRVec)] <- minpeakheight
  tmpGLRVec[is.infinite(SldWin$GLRVec)] <- minpeakheight
  tmpGLRVec <- c(minpeakheight, tmpGLRVec, minpeakheight)
  pracma::findpeaks(tmpGLRVec, sortstr = TRUE, minpeakdistance = minpeakdistance, minpeakheight = minpeakheight)[, 2] - 1
}



BFunc <- function(X, lambda) {
  (sign(X) * (lambda - 1) + 1) * abs(X) + 1
}

AFunc <- function(X, lambda) {
  exp(1 / (sign(X) * (lambda - 1) + 1) * log(BFunc(X, lambda)))
}

InvYJTrans <- function(X, lambda) {
  if (0 < lambda & lambda < 2) {
    out <- sign(X) * (AFunc(X, lambda) - 1)
  } else if (lambda == 0 & lambda == 2) {
    out <- sign(X) * (exp(abs(X)) - 1)
  } else if (lambda < 0) {
    out <- ifelse(0 <= X & X <= - 1 / lambda, 
                  sign(X) * (AFunc(X, lambda) - 1), NA)
  } else if (lambda > 2) {
    out <- ifelse(1 / (2 - lambda) < X & X < 0, 
                  sign(X) * (AFunc(X, lambda) - 1), NA)
  }
  out
}

#doubleFactorial <- function(j) {
#  prod(seq(j, 0, -2))
#}
#
#normalMoment <- function(j, sigma2) {
#  
#  sigma2 ^ (j / 2) * doubleFactorial(j - 1)
#  
#}
#
#normalMoment <- Vectorize(normalMoment, vectorize.args = c('j', 'sigma2'))
#
#InvYJTransJthDerv <- function(X, lambda, j) {
#  
#  BpartFunc <- function(jVec, lambda, X) {
#    prod(1 - jVec * sign(X) * (sign(X) * (lambda - 1) + 1))
#  }
#  
#  jVec <- seq(0, j - 1, 1)
#  
#  Apart <- AFunc(X, lambda) / BFunc(X, lambda) ^ j
#  Bpart <- BpartFunc(jVec, lambda, X)
#  
#  if (0 < lambda & lambda < 2) {
#    out <- Apart * Bpart
#  } else if (lambda == 0 & lambda == 2) {
#    out <- (-1) ^ (j - 1) * exp(abs(X))
#  } else if (lambda < 0) {
#    out <- ifelse(0 <= X & X <= - 1 / lambda, 
#                  Apart * Bpart, NA)
#  } else if (lambda > 2) {
#    out <- ifelse(1 / (2 - lambda) < X & X < 0, 
#                  Apart * Bpart, NA)
#  }
#  out
#}
#
#InvYJTransJthDerv <- Vectorize(InvYJTransJthDerv, vectorize.args = c('X', 'lambda', 'j'))
#
InvYJTrans1stDerv <- function(X, lambda){
  if (0 < lambda & lambda < 2) {
    out <- AFunc(X, lambda) / BFunc(X, lambda)
  } else if (lambda == 0 & lambda == 2) {
    out <- exp(sign(X))
  } else if (lambda < 0) {
    out <- ifelse(0 <= X & X <= - 1 / lambda, 
                  AFunc(X, lambda) / BFunc(X, lambda), NA)
  } else if (lambda > 2) {
    out <- ifelse(1 / (2 - lambda) < X & X < 0, 
                  AFunc(X, lambda) / BFunc(X, lambda), NA)
  }
  out
}

InvYJTrans2ndDerv <- function(X, lambda){
  if (0 < lambda & lambda < 2) {
    out <- AFunc(X, lambda) / BFunc(X, lambda) ^ 2 * (1 - sign(X) * (sign(X) * (lambda - 1) + 1))
  } else if (lambda == 0 & lambda == 2) {
    out <- sign(X) * exp(sign(X))
  } else if (lambda < 0) {
    out <- ifelse(0 <= X & X <= - 1 / lambda, 
                  AFunc(X, lambda) / BFunc(X, lambda) ^ 2 * (1 - sign(X) * (sign(X) * (lambda - 1) + 1)), NA)
  } else if (lambda > 2) {
    out <- ifelse(1 / (2 - lambda) < X & X < 0, 
                  AFunc(X, lambda) / BFunc(X, lambda) ^ 2 * (1 - sign(X) * (sign(X) * (lambda - 1) + 1)), NA)
  }
  out
}


varXARMA11 <- function(sigma2, ar, ma) {
  (1 + 2 * ar * ma + ma ^ 2) * sigma2 / (1 - ar ^ 2)
}


MuYARMA11 <- function(lambda, muX, sigma2, ar, ma) {
  
  InvYJTrans(muX, lambda) + InvYJTrans2ndDerv(muX, lambda) / 2 * 
    varXARMA11(sigma2, ar, ma)
  
}


VarYARMA11 <- function(lambda, muX, sigma2, ar, ma) {
 
  InvYJTrans1stDerv(muX, lambda) ^ 2 * varXARMA11(sigma2, ar, ma) + 
    1 / 2 * (InvYJTrans2ndDerv(muX, lambda) * varXARMA11(sigma2, ar, ma)) ^ 2
   
}

deltaARMA11Func <- function(muX, sigma2, ar, ma) {
  muX / sqrt(varXARMA11(sigma2, ar, ma))
}

muXARMA11Func <- function(delta, sigma2, ar, ma) {
  delta * sqrt(varXARMA11(sigma2, ar, ma))
}

muXARMA11Func <- Vectorize(muXARMA11Func, vectorize.args = c('delta', 'sigma2', 'ar', 'ma'))

arimaYJSim <- function(n, model, initial = 0) {
  
  mo <- list(order = model$order)
  mo$order[2] <- 0
  
  pos <- 1
  
  if (model$YJFlg == TRUE) {
    lambda <- as.numeric(model$ests[pos])
    pos <- pos + 1
  } 
  
  if (model$meanFlg == TRUE) {
    mu <- as.numeric(model$ests[pos])
    pos <- pos + 1
  } else {
    mu <- 0
  }
  
  sigma2 <- as.numeric(model$ests[pos])
  pos <- pos + 1
  
  if (model$order[1] > 0) {
    arVec <- as.numeric(model$ests[pos:(pos + model$order[1] - 1)])
    mo[["ar"]] <- arVec
    if (sum(arVec) == 0) {
      mo$order[1] <- 0
      mo$ar <- NULL
    } 
    pos <- pos + model$order[1]
  } 
  
  if (model$order[3] > 0) {
    maVec <- as.numeric(model$ests[pos:(pos + model$order[3] - 1)])
    mo[["ma"]] <- maVec
    if (sum(maVec) == 0) {
      mo$order[3] <- 0
      mo$ma <- NULL
    }
    
  } 
  
  
  if (mo$order[1] + mo$order[3] > 0) {
    start.innov <- rep(0, model$order[1] + model$order[3])
  } else {
    start.innov <- 0
  }

  
  tmpSim <- arima.sim(
    mo, n = n - model$order[2], n.start = mo$order[1] + mo$order[3],
    start.innov = start.innov, 
    mean = 0, sd = sqrt(sigma2))

  if (model$order[2] > 0) {
    nn <- length(tmpSim)
    
    tmpSim <- tmpSim[(nn - (n - model$order[2]) + 1):nn]
    for (j in 1:model$order[2]) {
      
      if (j != model$order[2]) {
        init <- 0
      } else {
        init <- initial
      }
      
      tmpSim <- filter(c(init, tmpSim), 1, method = 'recursive')
      
    }
  }
  
  tmpSim <- tmpSim + mu

  if (model$YJFlg == TRUE) {
    out <- InvYJTrans(tmpSim, lambda) 
  } else {
    out <- tmpSim
  }
  
  out

}

distGLR <- function(halfWinLength, model0, maxOrder, method = 'bic', 
                    searchYJ = TRUE, lambda0 = 1, searchM1M2 = TRUE, 
                    nsim1 = 10, nsim2 = 10, ntry = 10) {
  
  GLRVec <- rep(NA, nsim1 * nsim2)
  
  nt1 <- 0
  w <- 0
  h1 <- 0
  
  while(nsim1 > h1 & ntry > nt1) {
    nt2 <- 0
    h2 <- 0
    tmp <- try(arimaYJSim(2 * halfWinLength, model0), silent = TRUE)
    nt1 <- nt1 + 1
    if (class(tmp) != 'try-error' & sum(!is.na(tmp)) == 2 * halfWinLength) {
      nt1 <- nt1 - 1
      model1 <- try(getBestModel(tmp, halfWinLength, maxOrder, method, searchYJ, lambda0, searchM1M2), silent = TRUE)
      nt1 <- nt1 + 1
      if (class(model1) != 'try-error') {
        nt1 <- nt1 - 1
        while(nsim2 > h2 & ntry > nt2) {
          tmp1 <- try(arimaYJSim(2 * halfWinLength, model1), silent = TRUE)
          
          nt2 <- nt2 + 1
          if (class(tmp1) != 'try-error' & sum(!is.na(tmp1)) == 2 * halfWinLength) {
            
            nt2 <- nt2 - 1
            Y1 <- tmp1[1:halfWinLength]
            Y2 <- tmp1[-c(1:halfWinLength)]
            GLR <- try(generalizedloglikRatio(Y1, Y2, halfWinLength, maxOrder, method, searchYJ, lambda0, searchM1M2), 
                       silent = TRUE)
            
            nt2 <- nt2 + 1
            if (class(GLR) != 'try-error') {
              nt2 <- nt2 - 1
              w <- w + 1
              h2 <- h2 + 1
              if (h2 == 1) {h1 <- h1 + 1} 
              GLRVec[w] <- GLR$GLR
              
            }
            
          }
          
        }
        
        
      }
      
    }
    
  }
  
  
  GLRVec
  
}

detectCP <- function(alpha, SldWin, ScanIdx, halfWinLength, maxOrder, method = 'bic',
                     searchYJ = TRUE, lambda0 = 1, searchM1M2 = TRUE, asymptotic = FALSE,
                     nsim1 = 10, nsim2 = 10, ntry = 10, tolProp = 0.5, earlyStop = TRUE) {
  
  rootFinding <- function(p0, GLR, dGLR) {
    
    GLR - quantile(dGLR, p0, na.rm = TRUE)
    
  }
  
  n <- length(ScanIdx)
  out <- matrix(NA, nrow = n, ncol = 6)
  
  for (j in 1:n) {
    tmpIdx <- ScanIdx[j]
    tt <- tmpIdx + halfWinLength
    out[j, 1] <- tmpIdx
    out[j, 2] <- tt
    cat('Check the measurement at t:', tt, '\n')
    tmpModel <- SldWin$model0[[tmpIdx]]
    tmpGLR <- SldWin$GLRVec[tmpIdx]
    
    
    if (asymptotic == FALSE) {
      cat('Simulating the reference distribution', '\n')
      dGLR <- suppressWarnings(distGLR(halfWinLength, tmpModel, maxOrder, method, 
                                       searchYJ, lambda0, searchM1M2, nsim1, nsim2, ntry))
      if (sum(!is.na(dGLR)) > tolProp * nsim1 * nsim2) {
        #if (approx == TRUE) {
        kdeGLR <- kde1d::kde1d(dGLR)
        pValue1 <- 1 - kde1d::pkde1d(tmpGLR, kdeGLR)
        out[j, 3] <- pValue1
        cat('approximate pValue:', pValue1, '\n')
        #} else {
        #critValue <- quantile(dGLR, 1 - alpha)
        if (tmpGLR > max(dGLR)) {
          pValue2 <- 0
        } else if (tmpGLR < min(dGLR)) {
          pValue2 <- 1
        } else {
          pValue2 <- 1 - uniroot(rootFinding, interval = c(0, 1), GLR = tmpGLR, dGLR = dGLR)$root
        }
        out[j, 4] <- pValue2
        cat('empirical pValue:', pValue2, '\n')
        #}
        
        if (pValue1 >= alpha & pValue2 >= alpha) {
          #
          out[j, c(5, 6)] <- 0
          if (earlyStop == TRUE) {
            cat('The change point detection stops at t:', tt, '\n')
            break
          }
          #
        } else {
          if (pValue1 < alpha) {
            out[j, 5] <- 1
            cat('approximate method detects a change point at t:', tt, '\n')
          } else {
            out[j, 5] <- 0
          }
          
          if (pValue2 < alpha) {
            out[j, 6] <- 1
            cat('empirical method detects a change point at t:', tt, '\n')
          } else {
            out[j, 6] <- 0
          }
        }
        
      } else {
        cat('fail to get the reference distribution at t:', tt, '\n')
      }
    } else {
      
      pValue1 <- 1 - pchisq(tmpGLR, length(tmpModel$ests))
      out[j, 3] <- pValue1
      
      if (pValue1 < alpha) {
        out[j, 5] <- 1
        cat('asymptotic method detects a change point at t:', tt, '\n')
      } else {
        out[j, 5] <- 0
        cat('The change point detection stops at t:', tt, '\n')
      }
    }
    
  }
  out
}


detectCPPar <- function(alpha, SldWin, ScanIdx, halfWinLength, maxOrder, method = 'bic',
                        searchYJ = TRUE, lambda0 = 1, searchM1M2 = TRUE, asymptotic = FALSE,
                        nsim1 = 10, nsim2 = 10, ntry = 10, tolProp = 0.5, cluster = parallel::makeCluster(2)) {
  
  rootFinding <- function(p0, GLR, dGLR) {
    
    GLR - quantile(dGLR, p0, na.rm = TRUE)
    
  }
  
  n <- length(ScanIdx)
  #out <- matrix(NA, nrow = n, ncol = 6)
  alpha <- alpha
  SldWin <- SldWin
  ScanIdx <- ScanIdx
  halfWinLength <- halfWinLength
  maxOrder <- maxOrder
  method <- method
  searchYJ <- searchYJ
  lambda0 <- lambda0
  searchM1M2 <- searchM1M2
  asymptotic <- asymptotic
  nsim1 <- nsim1
  nsim2 <- nsim2
  ntry <- ntry
  tolProp <- tolProp
  
  clusterExport(cl = cluster, c('distGLR', 'arimaYJSim', 'muXARMA11Func',
                                'deltaARMA11Func', 'InvYJTrans2ndDerv', 'InvYJTrans',
                                'AFunc', 'BFunc', 'ScanGLRVec',
                                'SldWin', 'generalizedloglikRatio', 'getFitted',
                                'getBestModel', 'ARMACheck', 'AICc',
                                'BIC', 'AIC', 'getEsts',
                                'loglikFunc', 'YJTrans'), envir = .GlobalEnv)
  clusterExport(cl = cluster, c('rootFinding', 'alpha', 'SldWin',
                                'ScanIdx', 'halfWinLength', 'maxOrder',
                                'method', 'searchYJ', 'lambda0',
                                'searchM1M2', 'asymptotic', 'nsim1',
                                'nsim2', 'ntry', 'tolProp'), envir = environment())
  
  tmp <- parLapplyLB(cl = cluster, X = 1:n, fun = function(X) {
    tmpIdx <- ScanIdx[X]
    tt <- tmpIdx + halfWinLength
    out[X, 1] <- tmpIdx
    out[X, 2] <- tt
    cat('Check the measurement at t:', tt, '\n')
    tmpModel <- SldWin$model0[[tmpIdx]]
    tmpGLR <- SldWin$GLRVec[tmpIdx]
    
    
    if (asymptotic == FALSE) {
      cat('Simulating the reference distribution', '\n')
      dGLR <- suppressWarnings(distGLR(halfWinLength, tmpModel, maxOrder, method, 
                                       searchYJ, lambda0, searchM1M2, nsim1, nsim2, ntry))
      if (sum(!is.na(dGLR)) > tolProp * nsim1 * nsim2) {
        #if (approx == TRUE) {
        kdeGLR <- kde1d::kde1d(dGLR)
        pValue1 <- 1 - kde1d::pkde1d(tmpGLR, kdeGLR)
        out[X, 3] <- pValue1
        cat('approximate pValue:', pValue1, '\n')
        #} else {
        #critValue <- quantile(dGLR, 1 - alpha)
        if (tmpGLR > max(dGLR)) {
          pValue2 <- 0
        } else if (tmpGLR < min(dGLR)) {
          pValue2 <- 1
        } else {
          pValue2 <- 1 - uniroot(rootFinding, interval = c(0, 1), GLR = tmpGLR, dGLR = dGLR)$root
        }
        out[X, 4] <- pValue2
        cat('empirical pValue:', pValue2, '\n')
        #}
        
        if (pValue1 >= alpha & pValue2 >= alpha) {
          #cat('The change point detection stops at t:', tt, '\n')
          out[j, c(5, 6)] <- 0
          #break
        } else {
          if (pValue1 < alpha) {
            out[X, 5] <- 1
            cat('approximate method detects a change point at t:', tt, '\n')
          } else {
            out[X, 5] <- 0
          }
          
          if (pValue2 < alpha) {
            out[X, 6] <- 1
            cat('empirical method detects a change point at t:', tt, '\n')
          } else {
            out[X, 6] <- 0
          }
        }
        
      } else {
        cat('fail to get the reference distribution at t:', tt, '\n')
      }
    } else {
      
      pValue1 <- 1 - pchisq(tmpGLR, length(tmpModel$ests))
      out[X, 3] <- pValue1
      
      if (pValue1 < alpha) {
        out[X, 5] <- 1
        cat('asymptotic method detects a change point at t:', tt, '\n')
      } else {
        out[X, 5] <- 0
        cat('The change point detection stops at t:', tt, '\n')
      }
    }
  })
  
  tmp
}


#load(file = 'C:/Users/yyao17/Desktop/example.Rdata')

#debug(SldWin)
#debug(generalizedloglikRatio)
#debug(getBestModel)

CPD <- function(Y, alpha, halfWinLength, maxOrder = c(2, 2, 2), method = 'bic', searchYJ = TRUE, lambda0 = 1, 
                searchM1M2 = TRUE, asymptotic = FALSE, nsim1 = 100, nsim2 = 1, ntry = 10, tolProp = 0.5, earlyStop = TRUE,
                minpeakdistance = halfWinLength / 2, minpeakheight = 0, kdeFlg = TRUE, 
                graph = TRUE, ylab = 'Observation', xlab = NULL) {
  
  n <- length(Y)
  
  step1 <- SldWin(Y, halfWinLength = halfWinLength, maxOrder = maxOrder, 
                  method = method, searchYJ = searchYJ, lambda0 = lambda0, searchM1M2 = searchM1M2)
  step2 <- ScanGLRVec(step1, n, halfWinLength, minpeakdistance = minpeakdistance, minpeakheight = minpeakheight)
  
  if (length(step2) > 0) {
    step3 <- detectCP(alpha = alpha, SldWin = step1, ScanIdx = step2, halfWinLength = halfWinLength, 
                      maxOrder = maxOrder, method = method,
                      searchYJ = searchYJ, lambda0 = lambda0, searchM1M2 = searchM1M2, asymptotic = asymptotic,
                      nsim1 = nsim1, nsim2 = nsim2, ntry = ntry, tolProp = tolProp, earlyStop = earlyStop)
    
    step3 <- na.omit(step3)
    
    if (kdeFlg == TRUE) {
      idx1 <- step3[which(step3[, 5] > 0), 1]
      idx2 <- step3[which(step3[, 5] > 0), 2]
    } else {
      idx1 <- step3[which(step3[, 6] > 0), 1]
      idx2 <- step3[which(step3[, 6] > 0), 2]
    }
    
    
    if (length(idx2) > 0) {
      idx <- sort(idx2)
      nidx <- length(idx)
      models <- list()
      Pred <- list()
      fitted <- list()
      
      for (i in 1:(nidx + 1)) {
        
        
        if (1 < i & i <= nidx) {
          tmpY <- Y[idx[i - 1]:(idx[i] - 1)]
        } else if (i == 1) {
          tmpY <- Y[1:(idx[i] - 1)]
        } else if (i > nidx ) {
          tmpY <- Y[idx[i - 1]:n]
        }
        
        models[[i]] <- suppressWarnings(
          try(getBestModel(tmpY, halfWinLength = halfWinLength, maxOrder = maxOrder), 
              silent = TRUE))
        fitted[[i]] <- getFitted(tmpY, models[[i]]$XResi, models[[i]]$order, 
                                 models[[i]]$ests, models[[i]]$YJFlg, models[[i]]$meanFlg)
      }
      
      means <- vector()
      updowns <- vector()
      
      for (i in 1:(nidx + 1)) {
        means[i] <- mean(fitted[[i]], na.rm = TRUE)
      }
      
      fitted <- unlist(fitted)
      
      out <- list(cp = idx, means = means, models = models, fitted = fitted)
    }
    
    if (graph == TRUE) {
      par(mfrow = c(2, 1))
      
      plot(Y, type = 'l', ylab = ylab, xlab = xlab, main = 'Time series plot with change points and fitted values')
      
      if (length(idx2) > 0) {
        #points(fitted)
        abline(v = idx2, lty = 2)
      }
      
      
      plot(c(1,n), c(min(step1$GLRVec), max(step1$GLRVec)), type = 'n', xlab = xlab, ylab = 'GLR', main = 'GLR plot')
      points((halfWinLength + 1):(n - halfWinLength + 1), step1$GLRVec, type = 'l', lty = 1)
      if (length(idx2) > 0) {
        points(idx2, step1$GLRVec[idx1])
      }
      par(mfrow = c(1, 1))
    }
  }
  
  return(out)
  
}


parCPD <- function(Y, alpha, halfWinLength, maxOrder = c(2, 2, 2), 
                   method = 'bic', searchYJ = TRUE, lambda0 = 1, 
                   searchM1M2 = TRUE, asymptotic = FALSE, 
                   nsim1 = 100, nsim2 = 1, ntry = 10, tolProp = 0.5, #earlyStop = TRUE,
                   minpeakdistance = halfWinLength / 2, minpeakheight = 0, kdeFlg = TRUE, 
                   graph = TRUE, ylab = 'Observation', xlab = "Index", cluster = parallel::makeCluster(2)) {
  
  
  Y <- Y
  alpha <- alpha
  halfWinLength <- halfWinLength
  maxOrder <- maxOrder
  method <- method
  searchYJ <- searchYJ
  lambda0 <- lambda0
  searchM1M2 <- searchM1M2
  asymptotic <- asymptotic
  nsim1 <- nsim1
  nsim2 <- nsim2
  ntry <- ntry
  tolProp <- tolProp
  #earlyStop <- earlyStop
  minpeakdistance <- minpeakdistance
  minpeakheight <- minpeakheight
  
  cl <- cluster
  
  parallel::clusterExport(cl = cl, c('distGLR', 'arimaYJSim', 'muXARMA11Func',
                                     'deltaARMA11Func', 'InvYJTrans2ndDerv', 'InvYJTrans',
                                     'AFunc', 'BFunc', 'ScanGLRVec',
                                     'SldWin', 'generalizedloglikRatio', 'getFitted',
                                     'getBestModel', 'ARMACheck', 'AICc',
                                     'BIC', 'AIC', 'getEsts',
                                     'loglikFunc', 'YJTrans', 'detectCP'), envir = .GlobalEnv)
  parallel::clusterExport(cl = cl, c('Y', 'alpha', 'halfWinLength', 'maxOrder',
                                     'method', 'searchYJ', 'lambda0',
                                     'searchM1M2', 'asymptotic', 'nsim1',
                                     'nsim2', 'ntry', 'tolProp', ##'earlyStop',
                                     'minpeakdistance', 'minpeakheight'), envir = environment())
  
  n <- length(Y)
  
  v <- n - 2 * halfWinLength + 1
  
  out <- NA
  
  step1 <- parallel::parLapplyLB(cl = cluster, 
                                 X = (halfWinLength + 1):(n - halfWinLength + 1), fun = function(X) {
                                   tmp <- Y[(X - halfWinLength):(X + halfWinLength - 1)]
                                   try(SldWin(tmp, halfWinLength = halfWinLength, maxOrder = maxOrder, 
                                              method = method, searchYJ = searchYJ, lambda0 = lambda0, 
                                              searchM1M2 = searchM1M2), silent = TRUE)
                                 })
  
  #step1 <- SldWin(Y, halfWinLength = halfWinLength, maxOrder = maxOrder, 
  #                method = method, searchYJ = searchYJ, lambda0 = lambda0, 
  #                searchM1M2 = searchM1M2)
  
  #glrvec <- rep(NA, v)
  
  models0 <- list()
  tmpGLR <- rep(NA, v)
  models0[["model0"]] <- list()
  for (i in 1:v) {
    tmpGLR[i] <- step1[[i]]$GLRVec
    models0$model0[i] <- step1[[i]]$model0
    #glrvec[i] <- step1[[i]]$GLRVec
  }
  
  models0[["GLRVec"]] <- tmpGLR
  step2 <- ScanGLRVec(models0, n, halfWinLength, minpeakdistance = minpeakdistance, minpeakheight = minpeakheight)
  
  vv <- length(step2)
  
  if (length(step2) > 0) {
    
    parallel::clusterExport(cl = cl, c('step1'), envir = environment())
    
    #step3 <- detectCP(alpha = alpha, SldWin = models, ScanIdx = step2[1], halfWinLength = halfWinLength, 
    #      maxOrder = maxOrder, method = method,
    #      searchYJ = searchYJ, lambda0 = lambda0, searchM1M2 = searchM1M2, asymptotic = asymptotic,
    #      nsim1 = nsim1, nsim2 = nsim2, ntry = ntry, tolProp = tolProp, earlyStop = earlyStop)
    
    tmpstep3 <- parallel::parLapplyLB(cl = cluster, X = step2, fun = function(X) {
      try(detectCP(alpha = alpha, SldWin = models0, ScanIdx = X, halfWinLength = halfWinLength, 
                   maxOrder = maxOrder, method = method,
                   searchYJ = searchYJ, lambda0 = lambda0, searchM1M2 = searchM1M2, asymptotic = asymptotic,
                   nsim1 = nsim1, nsim2 = nsim2, ntry = ntry, tolProp = tolProp, earlyStop = FALSE), silent = TRUE)
    })
    
    #step3 <- detectCP(alpha = alpha, SldWin = step1, ScanIdx = step2, halfWinLength = halfWinLength, 
    #                  maxOrder = maxOrder, method = method,
    #                  searchYJ = searchYJ, lambda0 = lambda0, searchM1M2 = searchM1M2, asymptotic = asymptotic,
    #                  nsim1 = nsim1, nsim2 = nsim2, ntry = ntry, tolProp = tolProp, earlyStop = earlyStop)
    
    step3 <- matrix(NA, nrow = vv, ncol = 6)
    
    for (i in 1:vv) {
      
      step3[i, ] <- tmpstep3[[i]]
      
    }
    
    step3 <- na.omit(step3)
    
    if (kdeFlg == TRUE) {
      idx1 <- step3[which(step3[, 5] > 0), 1]
      idx2 <- step3[which(step3[, 5] > 0), 2]
    } else {
      idx1 <- step3[which(step3[, 6] > 0), 1]
      idx2 <- step3[which(step3[, 6] > 0), 2]
    }
    
    
    models <- list()
    fitted <- list()
    
    if (length(idx2) > 0) {
      idx <- sort(idx2)
      nidx <- length(idx)
      
      Pred <- list()
      
      
      for (i in 1:(nidx + 1)) {
        
        
        if (1 < i & i <= nidx) {
          tmpY <- Y[idx[i - 1]:(idx[i] - 1)]
        } else if (i == 1) {
          tmpY <- Y[1:(idx[i] - 1)]
        } else if (i > nidx ) {
          tmpY <- Y[idx[i - 1]:n]
        }
        
        models[[i]] <- suppressWarnings(
          try(getBestModel(tmpY, halfWinLength = halfWinLength, maxOrder = maxOrder), 
              silent = TRUE))
        #fitted[[i]] <- getFitted(tmpY, models[[i]]$XResi, models[[i]]$order, 
        #                         models[[i]]$ests, models[[i]]$YJFlg, models[[i]]$meanFlg)
      }
      
      #means <- vector()
      #updowns <- vector()
    #  
      ##for (i in 1:(nidx + 1)) {
      #  means[i] <- mean(fitted[[i]], na.rm = TRUE)
      #}
      
      #fitted <- unlist(fitted)
      
      #out <- list(cp = idx, GLR = tmpGLR, means = means, models = models, fitted = fitted)
      out <- list(cp = idx, GLR = tmpGLR, models = models)
    } else {
      
      models[[1]] <- suppressWarnings(
        try(getBestModel(Y, halfWinLength = halfWinLength, maxOrder = maxOrder), 
            silent = TRUE))
      #fitted[[1]] <- getFitted(Y, models[[1]]$XResi, models[[1]]$order, 
      #                         models[[1]]$ests, models[[1]]$YJFlg, models[[1]]$meanFlg)
      
      #out <- list(cp = NA, GLR = tmpGLR, means = means, models = models, fitted = fitted)
      out <- list(cp = NA, GLR = tmpGLR, models = models)
      
    }
    
    if (graph == TRUE) {
      par(mfrow = c(2, 1))
      
      plot(Y, type = 'l', ylab = ylab, xlab = xlab, main = 'Time series plot with change points and fitted values')
      
      if (length(idx2) > 0) {
        #points(fitted)
        abline(v = idx2, lty = 2)
      }
      
      
      plot(c(1,n), c(min(tmpGLR), max(tmpGLR)), type = 'n', xlab = xlab, ylab = 'GLR', main = 'GLR plot')
      points((halfWinLength + 1):(n - halfWinLength + 1), tmpGLR, type = 'l', lty = 1)
      if (length(idx2) > 0) {
        points(idx2, tmpGLR[idx1])
      }
      par(mfrow = c(1, 1))
    }
  }
  
  parallel::stopCluster(cl)
  
  return(out)
  
}




getRef <- function(h, Y, order, ests, 
                YJFlg = TRUE, meanFlg = TRUE, nsim = 10000) {
  
  pred <- rep(NA, h)
  
  arVec <- numeric()
  maVec <- numeric()
  
  n <- length(Y)
  
  pos <- 1
  lambda <- 1
  mu <- 0
  
  if (YJFlg == TRUE) {
    lambda <- ests[pos]
    X <- YJTrans(Y, lambda)
    pos <- pos + 1
  } else {
    X <- Y
  }
  
  if (meanFlg == TRUE) {
    mu <- ests[pos]
    pos <- pos + 1
  }
  
  sigma2 <- ests[pos]
  pos <- pos + 1
  
  if (order[1] > 0) {
    arVec <- ests[(pos):(pos + order[1] - 1)]
    pos <- pos + order[1]
  }
  
  if (order[3] > 0) {
    maVec <- ests[(pos):(pos + order[3] - 1)]
    pos <- pos + order[3]
  }
  
  
  maMat <- matrix(0, nrow = n + h - order[2], ncol = n + h - order[2])
    
  for (i in 1:(n + h - order[2])) {
    if (i > 1) {
      maEst <- c(1, ARMAtoMA(ar = arVec, ma = maVec, lag.max = i - 1))
    } else {
      maEst <- 1
    }
    maMat[i, 1:i] <- maEst
  }
  
  sigma2Vec <- rep(NA, n + h - order[2])
  
  j <- 0
  for (i in 1:(n + h - order[2])) {
    if (i <= (n - order[2])) {
      sigma2Vec[i] <- sigma2
    } else {
      j <- j + 1
      sigma2Vec[i] <- sigma2 * sum(maMat[i, 1:j] ^ 2)
    }
  }
  
  ref <- matrix(NA, nrow = n + h, ncol = nsim)
  
  for (sim in 1:nsim) {
    fitted <- rep(NA, n + h - order[2])
    Eps <- rnorm(n + h - order[2], 0, sqrt(sigma2Vec))
    for (i in 1:(n + h - order[2])) {
      #if (i > 1) {
      #  maEst <- c(1, ARMAtoMA(ar = arVec, ma = maVec, lag.max = i - 1))
      #} else {
      #  maEst <- 1
      #}
      
      fitted[i] <- sum(Eps[1:i] * maMat[i, i:1])
      
      #j <- i + order[2]
      #
      #if (j > n) {
      #  ME[j] <- qnorm(1 - alpha / 2) * sqrt(sigma2) * sqrt(sum(maEst ^ 2))
      #}
      
    }
    
    fitted <- fitted + mu
    #lower <- fitted - ME
    #upper <- fitted + ME
    
    
    if (order[2] > 0) {
      fitted <- c(fitted, rep(0, order[2]))
      D <- matrix(0, nrow = order[2] + 1, ncol = n + h)
      D[1, 1:n] <- X
      for (i in 2:(order[2] + 1)) {
        for (j in 1:(n - i + 1)) {
          D[i, j] <- D[i - 1, j + 1] - D[i - 1, j]
        }
      }
      
      D[order[2] + 1, ] <- fitted
      
      for (i in order[2]:1) {
        for (j in 1:(n + h - i + 1)) {
          if (j > 1) {
            D[i, j] <- D[i + 1, j - 1] + D[i, j - 1]
          }
        }
      }
      
      fitted <- D[1, ]
      
    }
    
    if (YJFlg == TRUE) {
      fitted <- InvYJTrans(fitted, lambda)#biasCorrect(fitted, lambda, sigma2)
    }
    
    ref[, sim] <- fitted
    
  }
  

  ref
  
  
} 

getInterval <- function(Ref, alpha) {
  
  optFunc <- function(a, Ref, alpha) {
    
    q1 <- quantile(Ref, a)
    q2 <- quantile(Ref, 1 - (alpha - a))
    
    q2 - q1
    
  }
  
  n <- dim(Ref)[1]
  
  p1Vec <- rep(NA, n)
  
  for (i in 1:n) {
    p1Vec[i] <- optimize(f = optFunc, interval = c(0, alpha), Ref = Ref[i, ], alpha = alpha)$minimum
  }
  
  p2Vec <- 1 - (alpha - p1Vec)
  
  med <- rep(NA, n)
  lower <- rep(NA, n)
  upper <- rep(NA, n)
  
  for (i in 1:n) {
    med[i] <- quantile(Ref[i, ], 0.5)
    lower[i] <- quantile(Ref[i, ], p1Vec[i])
    upper[i] <- quantile(Ref[i, ], p2Vec[i])
  }
  list(fitted = rowMeans(Ref), med = med, lower = lower, upper = upper)
}

