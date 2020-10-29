library(forecast)
library(pracma)

##############################################
#Matrix form
##############################################
InvertQ <- function(coef){
  stopifnot(class(coef)=="numeric"||class(coef)=="matrix"||(class(coef)=="array" && (dim(coef)[1]==dim(coef)[2])))
  if (class(coef) == "numeric")
    coef <- array(coef,dim=c(1,1,length(coef)))
  if (class(coef) == "matrix")
    coef <- array(coef,dim=c(NROW(coef),NROW(coef),1))
  k <- dim(coef)[1]
  order <- dim(coef)[3]
  if (order==1)
    ans <- eigen(coef[,,1], symmetric=FALSE, only.values =TRUE)$value
  else{
    blockMat <- matrix(numeric((k*order)^2),k*order,k*order)
    blockMat[1:k,] <- coef
    Imat <- diag((order-1)*k)
    blockMat[((k+1):(k*order)),1:((order-1)*k)] <- Imat
    ans <- eigen(blockMat, symmetric=FALSE, only.values =TRUE)$value
  }
  MaxEigenvalue <- max(Mod(ans))
  if (MaxEigenvalue >= 1)
    return( warning("check stationary/invertibility condition !"))
}

parsMat <- function(n, parsVec, norder = 1) {
  Check <- InvertQ(parsVec)
  if (class(Check) == 'character') {
    NULL
  } else if (is.null(Check)) {
    Mat <- diag(n)
    for (i in 1:norder) {
      Mat <- Mat + Diag(rep(parsVec[i], n - i), k = -i)
    }
    Mat
  }
}

SigmaMat <- function(n, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL) {
  if (order[1] == 0) {
    phiMat <- diag(n - order[2])
  } else {
    phiMat <- parsMat(n - order[2], -phiVec, norder = order[1])
  }
  
  if (order[3] == 0) {
    thetaMat <- diag(n - order[2])
  } else {
    thetaMat <- parsMat(n - order[2], thetaVec, norder = order[3])
  }
  
  if (order[2] > 0) {
    Dif <- diag(n)
    for (k in 1:order[2]) {
      Dif1 <- Diag(rep(1, n - k + 1)) + Diag(rep(-1, n - k), -1)
      Dif1 <- Dif1[-1, ]
      Dif <- Dif1 %*% Dif
    }
  }

  if (order[2] == 0) {
    out1 <- solve(phiMat) %*% thetaMat
  } else {
    out1 <- pinv(phiMat %*% Dif) %*% thetaMat
  }
  
  out <- out1 %*% t(out1)

  gamma0 <- diag(out)

  list(SigmaMat = out, gamma0 = gamma0, OmegaMat = out1)
}

##############################################
  #process simulation
##############################################

# Simulate innovations
simInnov <- function(n, XSim = 'norm', XPars = c(0, 1)) {

  if (XSim == "norm") {

    me <- XPars[1];
    std <- sqrt(XPars[2]);

    pars <- c(me, std)

    rgen <- rnorm

  } else if (XSim == "t") {

    me <- 0;
    std <- sqrt(XPars[1] / (XPars[1] - 2));

    pars <- c(XPars[1])

    rgen <- rt

  } else if (XSim == "gamma") {

    me <- XPars[1] * XPars[2];
    std <- sqrt(XPars[1] * XPars[2] ^ 2);

    pars <- c(XPars[1], 1 / XPars[2])

    rgen <- rgamma


  } else if (XSim == "beta") {

    me <- XPars[1] / (XPars[1] + XPars[2]);
    std <- sqrt(XPars[1] * XPars[2] / ((XPars[1] + XPars[2]) ^ 2) / (XPars[1] + XPars[2] + 1));

    pars <- c(XPars[1], XPars[2])

    rgen <- rbeta

  }

  if (XSim != 't') {

    out <- rgen(n, pars[1], pars[2])

  } else {

    out <- rgen(n, pars[1])

  }

  (out - me) / std

}


simARMAProcess <- function(n, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, mu = 0, sigma2 = 1, innovDist = 'norm', innovPars = c(0, 1),
                         SigMat = SigmaMat(n = n, order = order, phiVec = phiVec, thetaVec = thetaVec)) {
  
    as.vector(SigMat$OmegaMat %*% (matrix(simInnov(n - order[2], XSim = innovDist, XPars = innovPars), ncol = 1) * sqrt(sigma2) - mu))

}


##############################################

boxcoxDeriv <- function(y, lambda) {
  
  n <- length(y)
  
  if (lambda >= 0) {
    
    if (lambda > 0) {
      
      res <- y ^ (lambda - 1)
      
    } else {
      
      res <- y ^ (-1)
      
    }
    
  } else {
    
    res <- y
    
  }
  
  return(res)
  
}

loglikForward <- function(y, max.p = 2, max.d = 2, max.q = 2, boxcox = TRUE, lambda = 0, nred = 6) {
  
  pars <- expand.grid(c(0:max.p), c(0:max.d), c(0:max.q), c(FALSE, TRUE))
  
  npars <- rowSums(pars[, c(-2)]) + 1 + boxcox
  
  pars <- cbind(pars, npars)#[10, ]
  
  nsample <- length(y) - pars[, 2]
  
  pars <- pars[which((nsample - pars[, 5] - nred) >= 0), ]
  
  kk <- dim(pars)[1]
  
  phi <- matrix(NA, nrow = kk, ncol = max.p)
  
  theta <- matrix(NA, nrow = kk, ncol = max.q)
  
  mu <- rep(0, kk)
  sigma2 <- rep(NA, kk)
  
  fitted <- matrix(NA, nrow = kk, ncol = length(y))
  
  loglik <- base::rep(-Inf, kk)

  for (jj in 1:kk) {
  
    #if (pars[jj, 5] < nsample[jj]) {
    
      if (boxcox == TRUE) {
        
        y.boxcox <- BoxCox(y, lambda = lambda)
        
        model <- try(Arima(y.boxcox, order = c(pars[jj, 1], pars[jj, 2], pars[jj, 3]), include.mean = pars[jj, 4]), silent = TRUE)
        
        if (class(model)[1] != 'try-error') {
        
          loglik[jj] <- sum(log(boxcoxDeriv(y, lambda = lambda))) + model$loglik
          
        } 
        
      } else {
        
        model <- try(Arima(y, order = c(pars[jj, 1], pars[jj, 2], pars[jj, 3]), include.mean = pars[jj, 4]), silent = TRUE)
        
        if (class(model)[1] != 'try-error') {
          
          loglik[jj] <- model$loglik
          
        } 
        
      }
      
      if (class(model)[1] != 'try-error') {
        
        #cat(pars[jj, 1], pars[jj, 3], pars[jj, 4], '\n')
        
        if (pars[jj, 1] > 0) {
          phi[jj, 1:pars[jj, 1]] <- model$model$phi
        }
        
        if (pars[jj, 3] > 0) {
          theta[jj, 1:pars[jj, 3]] <- model$model$theta[which(abs(model$model$theta) != 0)]
        }
        
        if (pars[jj, 4] == TRUE) {
          mu[jj] <- model$coef[pars[jj, 1] + pars[jj, 3] + 1]
        }
        
        sigma2[jj] <- model$sigma2

		fitted[jj, ] <- model$fitted

      }
      
    #} 
    
  }
  
  nfinit <- which(!is.infinite(loglik))
  
  out <- list(pars = pars[nfinit, ], loglik = loglik[nfinit], phi = phi[nfinit, ], theta = theta[nfinit, ], 
              mu = mu[nfinit], sigma2 = sigma2[nfinit], y.length = nsample[nfinit], fitted = fitted[nfinit, ])
  
  return(out)
  
}

BIC <- function(model) {
  
  #loglik <- loglikForward(y, max.p = max.p, max.d = max.d, max.q = max.q, boxcox = boxcox, lambda = lambda)
  
  
  
  npars <- model$pars[, 5]
  nsample <- length(model$y.length) - model$pars[, 2]
  loglik <- model$loglik
  
  npars * log(nsample) - 2 * loglik
  
}

optBIC <- function(y, initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = 6, lambda.round.digit = 8) {
  
  optimizerBoxCox <- function(lambda, y, max.p = 2, max.d = 2, max.q = 2, nred = 6) {
    
    #y.length = length(y)
    
    
    
    model <- loglikForward(y = y, 
                           max.p = max.p, max.d = max.d, max.q = max.q, 
                           boxcox = TRUE, lambda = lambda, nred = nred)
    
    out <- min(BIC(model = model))
    
    #cat('lambda:', lambda, 'and min BIC:', out, '\n')
    
    return(out)
    
  }
  
  #######################################################
  minBIC1 <- Inf
  minBIC2 <- Inf
  
  #######################################################
  #Without Box Cox
  #######################################################
  
  model1 <- loglikForward(y = y, 
                          max.p = max.p, max.d = max.d, max.q = max.q, boxcox = FALSE, lambda = 0, nred = nred)
  BIC1 <- BIC(model = model1)
  numMinBIC1 <- which.min(BIC1)[1]
  minBIC1 <- BIC1[numMinBIC1]
  
  #######################################################
  #With Box Cox
  #######################################################
  
  lambda <- try(optim(par = initLambda, fn = optimizerBoxCox, method = 'Brent', lower = lower, upper = upper, 
                      y = y, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred), silent = TRUE)#, hessian = TRUE, control = list(parscale = 1e-2))
  
  #lambda <- optim(par = initLambda, fn = optimizerBoxCox, method = 'L-BFGS-B', lower = lower, upper = upper, 
  #                y = y, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred, hessian = TRUE)#, control = list(parscale = 1e-2))
  
  if (class(lambda) != 'try-error' & !is.na(lambda$value)) {
    model2 <- loglikForward(y = y, 
                            max.p = max.p, max.d = max.d, max.q = max.q, boxcox = TRUE, lambda = round(lambda$par, lambda.round.digit), nred = nred)
    
    BIC2 <- BIC(model = model2)
    numMinBIC2 <- which.min(BIC2)[1]
    minBIC2 <- BIC2[numMinBIC2]  
  }
  
  #######################################################
  
  if (minBIC1 <= minBIC2) {
    
	outmodel <- model1
	numMinBIC <- numMinBIC1
	minBIC <- minBIC1
	boxcox = FALSE
	lambda = NA
	
  } else {
    
    outmodel <- model2
	numMinBIC <- numMinBIC2
	minBIC <- minBIC2
	boxcox = TRUE
    lambda = lambda$par
	
  }
  
  if (dim(outmodel$pars)[1] == 1) {
	
	pars <- outmodel$pars
	phi <- outmodel$phi
	theta <- outmodel$theta
	loglik <- outmodel$loglik
	mu <- outmodel$mu
	sigma2 <- outmodel$sigma2
	fitted <- outmodel$fitted
	
  } else {
	if (max.p <= 1) {
      
      phi <- outmodel$phi[numMinBIC]
      
    } else {
      
      phi <- outmodel$phi[numMinBIC, ]
      
    }
    
    if (max.q <= 1) {
      
      theta <- outmodel$theta[numMinBIC]
      
    } else {
      
      theta <- outmodel$theta[numMinBIC, ]
      
    }
    
	pars <- outmodel$pars[numMinBIC, ]
	loglik <- outmodel$loglik[numMinBIC]
	mu <- outmodel$mu[numMinBIC]
	sigma2 <- outmodel$sigma2[numMinBIC]
	#y.length <- outmodel$y.length[numMinBIC]
	fitted <- outmodel$fitted[numMinBIC, ]
	
  }
  
  out <- list(minBIC = minBIC, pars = pars, loglik = loglik, 
                phi = phi, theta = theta, 
                mu = mu, sigma2 = sigma2, #y.length = , 
                boxcox = boxcox, lambda = lambda, fitted = fitted)
  
  out
  
}


loglikRatio <- function(y, t, loglik0, initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = 6) {
  
  n <- length(y)
  
  y1 <- y[1:t]
  model1 <- optBIC(y1, initLambda = initLambda, lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred)
  
  y2 <- y[(t + 1):n]
  model2 <- optBIC(y2, initLambda = initLambda, lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred)
  
  loglikRatio <- -2 * (loglik0 - (model1$loglik + model2$loglik))
  
  out <- list(loglikRatio = loglikRatio, model1 = model1, model2 = model2)
  
  cat('t:', t, 'model0:', loglik0, 'model1:', model1$loglik, 'model2:', model2$loglik, 'loglikRatio:', loglikRatio, '\n')
  
  return(out)
  
}

loglikRatioMax <- function(y, minsamp = 9, initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = minsamp) {
  
  n <- length(y)
  
  start <- minsamp
  end <- n - minsamp + 1
  
  loglikRatio <- rep(NA, end - start + 1)
  
  model0 <- optBIC(y, initLambda = initLambda, lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred)
  
  tvec <- start:end
  
  maxloglikRatio <- -Inf
  
  r <- 0
  
  for (t in tvec){
    
    r <- r + 1
    
    loglikRatioModel <- loglikRatio(y, t, loglik0 = model0$loglik, initLambda = initLambda, 
                                       lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q)
    
    loglikRatio[r] <- loglikRatioModel$loglikRatio
    
    if (loglikRatioModel$loglikRatio > maxloglikRatio) {
      
      maxloglikRatio <- loglikRatioModel$loglikRatio
      model1 <- loglikRatioModel$model1
      model2 <- loglikRatioModel$model2
      
    }
    
  }
  
  numMaxlikRatio <- which.max(loglikRatio)
  
  out <- list(
    t = tvec[numMaxlikRatio], 
    maxloglikRatio = maxloglikRatio, 
    model0 = model0, 
    model1 = model1, 
    model2 = model2,
    n = n)
  
  out
  
}

parsVecs <- function(model) {
  
  if (!is.na(unlist(model$lambda[1]))) {
    lambda <- unlist(model$lambda[1])
  } else {
    lambda <- NULL
  }
  
  if (unlist(model$pars[1]) > 0) {
    phi <- unlist(model$phi[1:unlist(model$pars[1])])
  } else {
    phi <- NULL
  }
  
  if (unlist(model$pars[3]) > 0) {
    theta <- unlist(model$theta[1:unlist(model$pars[3])])
  } else {
    theta <- NULL
  }
  
  #phi <- ifelse(unlist(model$pars[1]) > 0, unlist(model$phi[1:unlist(model$pars[1])]), 0)
  #theta <- ifelse(unlist(model$pars[3]) > 0, unlist(model$theta[1:unlist(model$pars[3])]), 0)
  order <- unlist(model$pars[1:3])
  mu <- unlist(model$mu)
  sigma2 <- unlist(model$sigma2)
  
  list(lambda = lambda, phi = phi, theta = theta, order = order, mu = mu, sigma2 = sigma2)
  
}

distPars0 <- function(n, model = list(pars = c(0, 0, 0, FALSE, 2), phi = NULL, theta = NULL, mu = 0, sigma2 = 1, boxcox = TRUE, lambda = 1/3), 
                  initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = 6, nsim = 100, seed = 12345) {
  
  set.seed(seed)
  
  pars0 <- parsVecs(model)
  
  lambda0 <- pars0$lambda
  phi0 <- pars0$phi
  theta0 <- pars0$theta
  order0 <- pars0$order
  mu0 <- pars0$mu
  sigma20 <- pars0$sigma2
  
  SigMat0 <- SigmaMat(n = n, order = order0, phiVec = phi0, thetaVec = theta0)
  
  lambda01 <- matrix(NA, nrow = nsim, ncol = 1)
  phi01 <- matrix(NA, nrow = nsim, ncol = max.p)
  theta01 <- matrix(NA, nrow = nsim, ncol = max.q)
  order01 <- matrix(0, nrow = nsim, ncol = 3)
  mu01 <- matrix(0, nrow = nsim, ncol = 1)
  sigma201 <- matrix(1, nrow = nsim, ncol = 1)
  loglik0 <- matrix(1, nrow = nsim, ncol = 1)
    
  testPhi <- NULL
  testTheta <- NULL
	
  for (ii in 1:nsim) {
    
    cat('process:', ii, '/', nsim, '\n')
    
    flg <- 0
    
    while(flg == 0) {
      
      #simVec <- arima.sim(
      #  model = list(
      #    ar = phi0, 
      #    ma = theta0,
      #    order = order0
      #  ),
      #  mean = mu0,
      #  sd = sqrt(sigma20),
      #  n = n
      #)
	  
	  simVec <- simARMAProcess(n = n, order = order0, phiVec = phi0, thetaVec = theta0, sigma2 = sigma20, SigMat = SigMat0)
      
      if (model$boxcox == TRUE) {
        simVec <- InvBoxCox(simVec, lambda0)
      }
      
      mo <- try(optBIC(simVec, initLambda = initLambda, lower = lower, upper = upper, 
                        max.p = max.p, max.d = max.d, max.q = max.q, nred = nred), silent = TRUE)
      
      if (class(mo)[1] != 'try-error') {
        
        pars1 <- parsVecs(mo)
        
		lambda1 <- pars1$lambda
		phi1 <- pars1$phi
		theta1 <- pars1$theta
		order1 <- pars1$order
		mu1 <- pars1$mu
		sigma21 <- pars1$sigma2
		
		if (order1[1] > 0) {
			testPhi <- InvertQ(-phi1)
		}
		
		if (order1[3] > 0) {
			testTheta <- InvertQ(theta1)
		}
		
		if (is.null(testPhi) & is.null(testTheta)) {
		
			flg <- 1
			
			if (order1[1] > 0) {
				phi01[ii, 1:order1[1]] <- phi1
			}
			
			if (order1[3] > 0) {
				theta01[ii, 1:order1[3]] <- theta1
			}
			
			lambda01[ii] <- ifelse(is.null(pars1$lambda[1]), NA, pars1$lambda[1])
			order01[ii, ] <- pars1$order
			mu01[ii] <- pars1$mu
			sigma201[ii] <- pars1$sigma2
			
		}
        
      }
      
    }
    
  }
  
  out <- list(order = order01, phi = phi01, theta = theta01, mu = mu01, sigma2 = sigma201, lambda = lambda01)
  out
  
}


#debug(distPars0)


distLoglikRatio <- function(distPars0, t, n, initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = 6, 
                            nsim = 100, seed = 12345) {
  
  set.seed(seed)
  
  npars <- dim(distPars0$order)[1]
  #Lambda <- matrix(NA, nrow = nsim, ncol = npars)
  
  totsim <- npars * nsim
  out <- rep(NA, totsim)
  
  a <- 0
  
  for (r in 1:npars) {
    
    if (dim(distPars0$phi)[2] > 1) {
      
      if (distPars0$order[r, 1] >= 1) {
        phi1 <- distPars0$phi[r, 1:distPars0$order[r, 1]]
      } else {
        phi1 <- NULL
      }
      
    } else {
      phi1 <- ifelse(is.na(distPars0$phi[r]), NULL, distPars0$phi[r])
    }
    
    if (dim(distPars0$theta)[2] > 1) {
      
      if (distPars0$order[r, 3] >= 1) {
        theta1 <- distPars0$theta[r, 1:distPars0$order[r, 3]]
      } else {
        theta1 <- NULL
      }
      
    } else {
      theta1 <- ifelse(is.na(distPars0$theta[r]), NULL, distPars0$theta[r])
    }
    
    lambda1 <- distPars0$lambda[r]
    #phi1 <- distPars0$phi
    #theta1 <- distPars0$theta
    order1 <- distPars0$order[r, ]
    mu1 <- distPars0$mu[r]
    sigma21 <- distPars0$sigma2[r]
    
	SigMat1 <- SigmaMat(n = n, order = order1, phiVec = phi1, thetaVec = theta1)
	
	for (v in 1:nsim) {
	  
	  a <- a + 1
    
	  cat('process:', a, '/', totsim, '\n')
	  
	  flg <- 0
	  flgPhi <- 0
	  flgTheta <- 0
	  
      while(flg == 0) {
       
        #simVec2 <- arima.sim(
        #  model = list(
        #    ar = phi1, 
        #    ma = theta1,
        #    order = order1
        #  ),
        #  mean = mu1,
        #  sd = sqrt(sigma21),
        #  n = n
        #)
        #
        #if (!is.na(lambda1)) {
        #  simVec2 <- InvBoxCox(simVec2, lambda1)
        #}
         
		simVec2 <- simARMAProcess(n = n, order = order1, phiVec = phi1, thetaVec = theta1, sigma2 = sigma21, SigMat = SigMat1)
      
		#if (model$boxcox == TRUE) {
		if (!is.na(lambda1)) {
			simVec2 <- InvBoxCox(simVec2, lambda1)
		} 
		 
        model0 <- try(optBIC(simVec2, initLambda = initLambda, lower = lower, upper = upper, 
                    max.p = max.p, max.d = max.d, max.q = max.q, nred = nred), silent = TRUE)
        
        if (class(model0)[1] != 'try-error') {
          
          if (model0$pars[1] > 0) {
            
            testPhi <- InvertQ(model0$phi[1:unlist(model0$pars[1])])
            
            if (class(testPhi) == 'NULL') {
              
              flgPhi <- 1
              
            }
            
          } else {
            
            flgPhi <- 1
            
          }
          
          if (model0$pars[3] > 0) {
            
            testTheta <- InvertQ(model0$theta[1:unlist(model0$pars[3])])
            
            if (class(testTheta) == 'NULL') {
              
              flgTheta <- 1
              
            }
            
          } else {
            
            flgTheta <- 1
            
          }
          
          if (flgPhi * flgTheta == 1) {
            
            llRatio <- try(loglikRatio(simVec2, t, model0$loglik, initLambda = initLambda, 
                                       lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred)$loglikRatio, silent = TRUE)
            
            if (class(llRatio)[1] != 'try-error') {
              #Lambda[v, r] <- llRatio
              out[a] <- llRatio
              flg <- 1
            }
            
          }
          
        }
        
        
      }
      
    }
    
  }
  
  #out <- as.vector(Lambda)
  out
  
}






#qq <- distPars0(36, max.p = 2, max.q = 2, nsim = 10)

#debug(distLoglikRatio)
#qq1 <- distLoglikRatio(qq, t = 18, n = 36, max.p = 2, max.q = 2, nsim = 10)

critLikelihoodRatio <- function(distLoglikRatio, alpha = 0.05) {
  
  quantile(distLoglikRatio, 1 - alpha)
  
}


binarySegmentation <- function(y, alpha = 0.05, minsamp = 9, initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = minsamp, nsim1 = 100, nsim2 = 10) {

	moVec <- rep(0, length(y))
	
	flg <- 0
	
	moInd <- 0
	
	while(flg == 0) {
	
		sumMoVec <- table(moVec)
		checkMoVec <- sumMoVec[sumMoVec >= 2 * minsamp]
		workMoVec <- as.numeric(names(checkMoVec]))
	
		if (sum(checkMoVec) > 0) {
		
			vv <- 0
			
			for (workMo in workMoVec) {
			
				vv <- vv + 1
				
				workMo <- moVec == workMoVec[vv]
				workY <- y[workmoVec]
				
				nn <- length(workY)
				
				step1 <- loglikRatioMax(y = workY, minsamp = minsamp)
				
				step2 <- distPars0(n = nn, 
					model = list(pars = step1$model0$pars, phi = step1$model0$phi, theta = step1$model0$theta, 
                              mu = step1$model0$mu, sigma2 = step1$model0$sigma2, boxcox = step1$model0$boxcox, lambda = step1$model0$lambda), 
					nsim = nsim1)
				
				step3 <- distLoglikRatio(step2, t = step1$t, n = nn, nsim = nsim2)
				
				step4 <- critLikelihoodRatio(step3, alpha = alpha)
				
				if (step1$maxloglikRatio > step4) {
					moInd <- moInd + 1
					ind.begin <- min(which(workmoVec))
					ind.end <- max(which(workmoVec))
					moVec[ind.begin:(ind.begin + step1$t)] <- moInd
					
					moInd <- moInd + 1
					moVec[(ind.begin + step1$t + 1):(ind.end)] <- moInd
				}

			}
	
		} else {
	
			flg <- 1
		}
		
	}

}


##############################################################

invBoxcoxDeriv2nd <- function(x, lambda) {
  
  n <- length(x)
  
  if (lambda >= 0) {
    
    if (lambda > 0) {
      
      res <- (lambda ^ 3) * (lambda + 1) * (lambda * x + 1) ^ (- lambda - 2)
      
    } else {
      
      res <- exp(lambda)
      
    }
    
  } else {
    
    res <- rep(0, n)
    
  }
  
  return(res)
  
}

predAdjTaylor <- function(predx, lambda, sigma2x) {
 
  InvBoxCox(predx, lambda) + invBoxcoxDeriv2nd(predx, lambda) / 2 * sigma2x
   
}


optBICEM <- function(y, initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = 6, tol = 1e-4, maxiter = 100) {
  
  MissVec <- is.na(y)
  y1 <- y
  
  flg <- 0
  
  predy <- mean(y, na.rm = TRUE)
  oldLoglik <- -Inf
  
  iter <- 0
  
  while(flg == 0) {
    
    iter <- iter + 1
    
    if (iter > maxiter) {
      flg <- 1
      stop('reach maximum of iteration')
    }
    
    y1[MissVec] <- predy
    
    # M step
    
    newmo <- optBIC(y1, initLambda = 1/3, lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred)
    newLoglik <- newmo$loglik
    
    dif <- newLoglik - oldLoglik
    
    cat('iter:', iter, 'dif:', dif, 'oldLoglik:', oldLoglik, 'newLoglik:', newLoglik, '\n')
    
    if (dif < tol) {
      
      flg <- 1
      
      if (dif < 0) {
        
        out <- oldmo
        
      } else {
        
        out <- newmo
        
      }
      
    } else {
      
      oldmo <- newmo
      
      oldLoglik <- newLoglik
      
      if (newmo$boxcox == TRUE) {
        
        # E step
        
        predy <- predAdjTaylor(newmo$mu, newmo$lambda, newmo$sigma2)
        
      } else {
        
        predy <- newmo$mu
        
      }
      
    }
      
  }
  
  return(out)
  
}


loglikRatioEM <- function(y, t, loglik0, initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = 6) {
  
  n <- length(y)
  
  y1 <- y[1:t]
  model1 <- optBICEM(y1, initLambda = initLambda, lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred)
  
  y2 <- y[(t + 1):n]
  model2 <- optBICEM(y2, initLambda = initLambda, lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred)
  
  loglikRatio <- -2 * (loglik0 - (model1$loglik + model2$loglik))
  
  out <- list(loglikRatio = loglikRatio, model1 = model1, model2 = model2)
  
  cat('t:', t, 'model0:', loglik0, 'model1:', model1$loglik, 'model2:', model2$loglik, 'loglikRatio:', loglikRatio, '\n')
  
  return(out)
  
}

loglikRatioMaxEM <- function(y, minsamp = 9, initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = minsamp) {
  
  n <- length(y)
  
  start <- minsamp
  end <- n - minsamp + 1
  
  loglikRatio <- rep(NA, end - start + 1)
  
  model0 <- optBICEM(y, initLambda = initLambda, lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred)
  
  tvec <- start:end
  
  maxloglikRatio <- -Inf
  
  r <- 0
  
  for (t in tvec){
    
    r <- r + 1
    
    loglikRatioModel <- loglikRatioEM(y, t, loglik0 = model0$loglik, initLambda = initLambda, 
                                      lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q)
    
    loglikRatio[r] <- loglikRatioModel$loglikRatio
    
    if (loglikRatioModel$loglikRatio > maxloglikRatio) {
      
      maxloglikRatio <- loglikRatioModel$loglikRatio
      model1 <- loglikRatioModel$model1
      model2 <- loglikRatioModel$model2
      
    }
    
  }
  
  numMaxlikRatio <- which.max(loglikRatio)
  
  out <- list(
    t = tvec[numMaxlikRatio], 
    maxloglikRatio = maxloglikRatio, 
    model0 = model0, 
    model1 = model1, 
    model2 = model2,
    n = n)
  
  out
  
}

distPars0EM <- function(n, MissVec, model = list(pars = c(0, 0, 0, FALSE, 2), phi = NULL, theta = NULL, mu = 0, sigma2 = 1, boxcox = TRUE, lambda = 1/3), 
                        initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = 6, nsim = 100, seed = 12345) {
  
  set.seed(seed)
  
  pars0 <- parsVecs(model)
  
  lambda0 <- pars0$lambda
  phi0 <- pars0$phi
  theta0 <- pars0$theta
  order0 <- pars0$order
  mu0 <- pars0$mu
  sigma20 <- pars0$sigma2
  
  lambda01 <- matrix(NA, nrow = nsim, ncol = 1)
  phi01 <- matrix(NA, nrow = nsim, ncol = max.p)
  theta01 <- matrix(NA, nrow = nsim, ncol = max.q)
  order01 <- matrix(0, nrow = nsim, ncol = 3)
  mu01 <- matrix(0, nrow = nsim, ncol = 1)
  sigma201 <- matrix(1, nrow = nsim, ncol = 1)
  loglik0 <- matrix(1, nrow = nsim, ncol = 1)
  
  for (ii in 1:nsim) {
    
    flg <- 0
    flgPhi <- 0
    flgTheta <- 0
    
    while(flg == 0) {
      
      simVec <- arima.sim(
        model = list(
          ar = phi0, 
          ma = theta0,
          order = order0
        ),
        mean = mu0,
        sd = sqrt(sigma20),
        n = n
      )
      
      if (model$boxcox == TRUE) {
        simVec <- InvBoxCox(simVec, lambda0)
      }
      
      simVec[MissVec] <- NA
      
      mo <- try(optBICEM(simVec, initLambda = initLambda, lower = lower, upper = upper, 
                       max.p = max.p, max.d = max.d, max.q = max.q, nred = nred), silent = TRUE)
      
      if (class(mo)[1] != 'try-error') {
        
        pars1 <- parsVecs(mo)
        
        if (max.p <= 1) {
          
          phi01[ii] <- ifelse(is.null(pars1$phi), NA, pars1$phi)
          
          # Stationary test
          testPhi <- InvertQ(phi01[ii])
          
        } else {
          
          if (pars1$order[1] > 0) {
            phi01[ii, 1:pars1$order[1]] <- pars1$phi[1:pars1$order[1]]
            
            # Stationary test
            testPhi <- InvertQ(phi01[ii, 1:pars1$order[1]])
          }
          
        }
        
        if (max.q <= 1) {
          
          theta01[ii] <- ifelse(is.null(pars1$theta), NA, pars1$theta)
          
          # Stationary test
          testTheta <- InvertQ(theta01[ii])
          
        } else {
          
          if (pars1$order[3] > 0) {
            theta01[ii, 1:pars1$order[3]] <- pars1$theta[1:pars1$order[3]]
            
            # Stationary test
            testTheta <- InvertQ(theta01[ii, 1:pars1$order[3]])
          }
          
        }
        
        lambda01[ii] <- pars1$lambda[1]
        #phi01[ii, ] <- ifelse(is.null(pars1$phi), NA, pars1$phi)
        #theta01[ii, ] <- ifelse(is.null(pars1$theta), NA, pars1$theta)
        order01[ii, ] <- pars1$order
        mu01[ii] <- pars1$mu
        sigma201[ii] <- pars1$sigma2
        
        if (exists('testPhi')) {
          
          if (class(testPhi) == 'NULL') {
            flgPhi <- 1
          }
          
        } else {
          flgPhi <- 1
        }
        
        if (exists('testTheta')) {
          
          if (class(testTheta) == 'NULL') {
            flgTheta <- 1
          }
          
        } else {
          flgTheta <- 1
        }
        
        flg <- flgPhi * flgTheta
        
      }
      
    }
    
  }
  
  out <- list(order = order01, phi = phi01, theta = theta01, mu = mu01, sigma2 = sigma201, lambda = lambda01)
  out
  
}


#abc <- distPars0EM(n = 60, MissVec = is.na(simEM), model = list(pars = c(0, 0, 0, FALSE, 2), phi = NULL, theta = NULL, mu = 0, sigma2 = 1, boxcox = TRUE, lambda = 1/3), 
#                        initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = 6, nsim = 10, seed = 12345)


distLoglikRatioEM <- function(distPars0EM, t, n, MissVec, initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = 6, 
                              nsim = 100, seed = 12345) {
  
  set.seed(seed)
  
  npars <- dim(distPars0EM$order)[1]
  totsim <- npars * nsim
  out <- rep(NA, totsim)
  
  a <- 0
  
  for (v in 1:nsim) {
    
    for (r in 1:npars) {
      
      a <- a + 1
      
      cat('process:', a, '/', totsim, '\n')
      
      flg <- 0
      flgPhi <- 0
      flgTheta <- 0
      
      if (dim(distPars0EM$phi)[2] > 1) {
        
        if (distPars0EM$order[r, 1] >= 1) {
          phi1 <- distPars0EM$phi[r, 1:distPars0EM$order[r, 1]]
        } else {
          phi1 <- NULL
        }
        
      } else {
        phi1 <- ifelse(is.na(distPars0EM$phi[r]), NULL, distPars0EM$phi[r])
      }
      
      if (dim(distPars0EM$theta)[2] > 1) {
        
        if (distPars0EM$order[r, 3] >= 1) {
          theta1 <- distPars0EM$theta[r, 1:distPars0EM$order[r, 3]]
        } else {
          theta1 <- NULL
        }
        
      } else {
        theta1 <- ifelse(is.na(distPars0EM$theta[r]), NULL, distPars0EM$theta[r])
      }
      
      lambda1 <- distPars0EM$lambda[r]
      #phi1 <- distPars0EM$phi
      #theta1 <- distPars0EM$theta
      order1 <- distPars0EM$order[r, ]
      mu1 <- distPars0EM$mu[r]
      sigma21 <- distPars0EM$sigma2[r]
      
      while(flg == 0) {
        
        simVec2 <- arima.sim(
          model = list(
            ar = phi1, 
            ma = theta1,
            order = order1
          ),
          mean = mu1,
          sd = sqrt(sigma21),
          n = n
        )
        
        if (!is.na(lambda1)) {
          simVec2 <- InvBoxCox(simVec2, lambda1)
        }
        
        simVec2[MissVec] <- NA
        
        model0 <- try(optBICEM(simVec2, initLambda = initLambda, lower = lower, upper = upper, 
                             max.p = max.p, max.d = max.d, max.q = max.q, nred = nred), silent = TRUE)
        
        if (class(model0)[1] != 'try-error') {
          
          if (model0$pars[1] > 0) {
            
            testPhi <- InvertQ(model0$phi[1:unlist(model0$pars[1])])
            
            if (class(testPhi) == 'NULL') {
              
              flgPhi <- 1
              
            }
            
          } else {
            
            flgPhi <- 1
            
          }
          
          if (model0$pars[3] > 0) {
            
            testTheta <- InvertQ(model0$theta[1:unlist(model0$pars[3])])
            
            if (class(testTheta) == 'NULL') {
              
              flgTheta <- 1
              
            }
            
          } else {
            
            flgTheta <- 1
            
          }
          
          if (flgPhi * flgTheta == 1) {
            
            llRatio <- try(loglikRatioEM(simVec2, t, model0$loglik, initLambda = initLambda, 
                                       lower = lower, upper = upper, max.p = max.p, max.d = max.d, max.q = max.q, nred = nred)$loglikRatio, silent = TRUE)
            
            if (class(llRatio)[1] != 'try-error') {
              out[a] <- llRatio
              flg <- 1
            }
            
          }
          
        }
        
        
      }
      
    }
    
  }
  
  #out <- as.vector(Lambda)
  out
  
}

#debug(distLoglikRatioEM)

#abc11 <- distLoglikRatioEM(abc, t = 37, n = 60, MissVec = is.na(simEM), initLambda = 1/3, lower = 0, upper = 1, max.p = 2, max.d = 2, max.q = 2, nred = 6, 
#                              nsim = 10, seed = 12345)



