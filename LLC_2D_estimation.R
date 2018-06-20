#--- Code to compute LLC density estimation in dimension 2.
rm(list = ls())
setwd("~/gitfolders/MTP_2_density_estimation")
library(MASS)
library(Rglpk)
library(lpSolve)
library(geometry)
library(scatterplot3d)
library(plot3D)
library(fields)
library(sos)
library(utils)
library(rgl)
library(matlib)
library(rmutil)
library(fastmatch)
#--- load code for functions

source("utilities.R")

#=================
#----------------
#---- Code for 2D ell-convex
#---- Assume points in X are integral. If not: will do rational approximations and scale
#example 1: X ~ points from a binomial, with correct density weights
if(FALSE){
  n = 11
  X = matrix(NA, nrow = 1, ncol = 2)
  w = c()
  for(i in 0:10){
    for(j in 0:10){
      X = rbind(X, c(i,j))
  w = append(w,choose(n,i)*choose(n,j))
    }
  }
  X = X[-1,]
  w = w/sum(w)
  #-- do the LLC estimation
  C = polytrope(X)
  V = polVertex.2d(C)
  info = get.ell.ineq(X,w,C,V)
  #start with an initial guess = log of independent Gaussian evaluated at these points
  initialGuess = -1/100*apply(info$xlist, 1, function(y) sum(y**2))
  #check that initial guess satisfies the constraints
  mult = info$ineq.mat%*%initialGuess
  eps = 0.001
  S = cgm(initialGuess, sigmafuncS(info$xlist,info$w.vec, eps), gradsigmafuncS(info$xlist,info$w.vec, eps), info$ineq.mat, info$w.vec,tol=1e-12)
  #double-check solution sums up to 1. This indicates quality of the convergence as well.
  S$value + S$par %*% c(w, rep(0, length(S$par) - length(w)))
#rescale should the solution not sum to 1 already
  S2 <- S
  S2$par = S$par - log(S$par %*% info$w.vec + sigmafunc(info$xlist,S$par,info$w.vec,eps))
  #-- print out the integral of the density
  S2$value + S2$par %*% c(w, rep(0, length(S2$par) - length(w)))
  #plot the results
  scatterplot3d(info$xlist[,1],info$xlist[,2],exp(S$par), highlight.3d=TRUE)
  #plot initial guess for comparison
  scatterplot3d(info$xlist[,1],info$xlist[,2],exp(initialGuess), highlight.3d=TRUE)
}

#---- Example 2: reproduce the plot in the paper
#load(file = "npoints100date2017-11-09.Rdata")
load(file="50PointsJun19.Rdata")
X = LS[[1]]
npoints = dim(X)[1]
X = round(X,1)*10
#weight vectors
w = rep(1/npoints,npoints)
C = polytrope(X)
V = polVertex.2d(C)
info = get.ell.ineq(X,w,C,V)
#start with an initial guess = log of independent Gaussian evaluated at these points
initialGuess = -1/100*apply(info$xlist, 1, function(y) sum(y**2))
eps = 0.001
#run the optimization
S = cgm(initialGuess, sigmafuncS(info$xlist,info$w.vec, eps), gradsigmafuncS(info$xlist,info$w.vec, eps), info$ineq.mat, info$w.vec,tol=1e-12)
#plot the results
scatterplot3d(info$xlist[,1],info$xlist[,2],exp(S$par), highlight.3d=TRUE)
#save
save(list = ls(), file="llc-output-June20.Rdata")
#--- fancier plot code
#pdf(file = "mpt2-vs-ell.pdf")
pdf(file = "mpt2-vs-ell-June20.pdf")
par(mfrow = c(2,2))
g = LS[[4]]
M = mesh(g$x,g$y)
scatterplot3d(c(M$x), c(M$y), c(g$z), highlight.3d=TRUE)
scatterplot3d(info$xlist[,1]/10,info$xlist[,2]/10,S$par, highlight.3d=TRUE)
scatterplot3d(c(M$x), c(M$y), exp(c(g$z)), highlight.3d=TRUE)
scatterplot3d(info$xlist[,1]/10,info$xlist[,2]/10,exp(S$par), highlight.3d=TRUE)
dev.off()




