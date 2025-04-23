test_that("multiplication works", {
  
  set.seed(59542234)
  n <- 10
  J <- 50
  X <- cbind(1,rep(c(0,1),each = n/2))
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b1 <- b1 - mean(b1)
  b0 <- b0 - mean(b0)
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = b0,
                              b1 = b1,
                              distn = "Poisson",
                              mean_z = 8)
  
  k_constr <- 2
  j_constr <- 1
  p <- 2
  
  # constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  constraint_fn <- rep(list(function(x){mean(x)}), 2)
  
  ##### Arguments to fix:
  
  # constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)
  constraint_grad_fn <- rep(list(function(x){ rep(1/length(x), length(x))}), 2)
  rho_init = 1
  tau = 1.2
  kappa = 0.8
  obj_tol = 100
  B_tol <- 1e-3
  constraint_tol = 1e-5
  init_tol = 1e6
  c1 = 1e-4
  maxit = 1000
  inner_maxit = 25
  
  full_fit <- #suppressMessages(
    emuFit_micro_penalized(X = X,
                           Y = Y,
                           B = NULL,
                           constraint_fn = rep(list(mean), 2), 
                           tolerance = 1e-5,
                           verbose = FALSE)
  
  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented
  
  X_cup <- X_cup_from_X(X,J)
  
  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  constraint_grad_fn <- function(x){ dpseudohuber_center_dx(x,0.1)}
  
  
  j_ref <- 5
  null_time <- proc.time()[3]
  null_fit <- fit_null(B = B,
                       Y = Y_aug,
                       X = X ,
                       X_cup = X_cup,
                       k_constr = k_constr,
                       j_constr = j_constr,
                       j_ref = j_ref,
                       constraint_fn = list(constraint_fn,constraint_fn),
                       constraint_tol = 1e-5,
                       B_tol = 1e-5,
                       constraint_grad_fn = list(constraint_grad_fn,
                                                 constraint_grad_fn),
                       verbose = TRUE,
                       trackB = (j_ref == 4)) ## just track for one j
  null_time <- proc.time()[3] - null_time
  
  null_repar_time <- proc.time()[3]
  null_repar_fit <- null_repar_cd(Y = Y_aug,
                                  X = X,
                                  B = B,
                                  j_constr = j_constr,
                                  k_constr = k_constr,
                                  j_ref = j_ref,
                                  constraint_fn = constraint_fn,
                                  constraint_grad_fn = constraint_grad_fn,
                                  tolerance = 3e-5,
                                  maxit = 1000)
  null_repar_time <- proc.time()[3] - null_repar_time
  
  min_mse_lag(X = X,
              Y = Y_aug,
              B = null_fit$B,
              constraint_grad_fn = constraint_grad_fn,
              k_constr = k_constr,
              j_constr = j_constr,
              j_ref = j_ref)
  
  min_mse_lag(X = X,
              Y = Y_aug,
              B = null_repar_fit$B,
              constraint_grad_fn = constraint_grad_fn,
              k_constr = k_constr,
              j_constr = j_constr,
              j_ref = j_ref)
  plot(null_repar_fit$B[2,])  
  points(null_fit$B[2,],col= "red",pch = 20)
  
  set.seed(8843)
  js <- replicate(1000,sample(setdiff(1:J,c(j_ref,j_constr)),2))
  diff_mse <- numeric(1000)
  for(jind in 1:1000){
    print(jind)
  an_deriv <- null_repar_ll_gr(x = rep(0,5),js = js[,jind],
                   B = null_fit$B,
                   z = update_z(Y_aug,X,null_fit$B),
                   p = 2,
                   Y = Y_aug,
                   X = X,
                   j_constr = j_constr,
                   k_constr = k_constr,
                   constraint_fn = constraint_fn,
                   constraint_grad_fn = constraint_grad_fn,
                               return_hess = FALSE)
  
  
  n_deriv <- numDeriv::grad(function(x) null_repar_ll(x,js = js[,jind],
                   B = null_fit$B,
                   z = update_z(Y_aug,X,null_fit$B),
                   p = 2,
                   Y = Y_aug,
                   X = X,
                   j_constr = j_constr,
                   k_constr = k_constr,
                   constraint_fn = constraint_fn),
                 x = rep(0,5))
  diff_mse[jind] <- sum((an_deriv - n_deriv)^2)
  }
  
  plot(1:1000,log(diff_mse))
  
  an_deriv <- null_repar_ll_gr(x = rep(0,5),js = js[,jind],
                               B = null_fit$B,
                               z = update_z(Y_aug,X,null_fit$B),
                               p = 2,
                               Y = Y_aug,
                               X = X,
                               j_constr = j_constr,
                               k_constr = k_constr,
                               constraint_fn = constraint_fn,
                               constraint_grad_fn = constraint_grad_fn,
                               return_hess = TRUE)
  

  n_hess <- numDeriv::hessian(function(x) null_repar_ll(x,js = js[,jind],
                                                      B = null_fit$B,
                                                      z = update_z(Y_aug,X,null_fit$B),
                                                      p = 2,
                                                      Y = Y_aug,
                                                      X = X,
                                                      j_constr = j_constr,
                                                      k_constr = k_constr,
                                                      constraint_fn = constraint_fn),
                            x = rep(0,5))
  
  an_deriv$info_inv/solve(n_hess)
  
})
