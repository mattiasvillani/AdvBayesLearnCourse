library(RTMB)

f <- function(pars) {
  x <- pars$x          
  sum(x^2) + (x[1]-1)*x[2]
}

objfun <- MakeADFun(f, list(x = c(1, 2, 3)), silent = TRUE)
objfun$fn(c(1, 2, 3))   # f(c(1,2,3))
objfun$gr(c(1, 2, 3))   # gradient vector
objfun$he(c(1, 2, 3))   # 3x3 Hessian matrix

initVal = c(1,2,0)
optRes <- optim(initVal, objfun$fn, gr = objfun$gr, method = "BFGS",
                hessian = TRUE)
optRes$par
optRes$hessian
