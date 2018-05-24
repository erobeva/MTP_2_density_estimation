#--- utility functions used by LLC_2D
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
#--- additional library
library(fastmatch)

#-- functions share with MTP2
integral_appr <- function(y) {
  d = length(y)
  ymean = mean(y)
  z = y - ymean
  tmp = 1
  for (i in 1:d){
    tmp = tmp + (z[i]*z[i])/(2*(d)*(d+1)) + (z[i]*z[i]*z[i])/
      (3*(d)*(d+1)*(d+2));
  }
  for (i in 1:(d-1)) {
    tmp = tmp / i;
  }
  tmp = tmp*exp(ymean);
  return(tmp);
}

integral <- function(y, eps) {
  y = sort(y);
  d = length(y);
  tmp = 0;
  if(y[d] - y[1] < eps) {
    return(integral_appr(y));
  } else {
    if(d == 2) {
      tmp = (exp(y[1]) - exp(y[2]))/(y[1] - y[2]);
      return(tmp);
    } else {
      e1 = integral(y[2:d], eps);
      e2 = integral(y[1:(d-1)], eps);
      return((e1-e2)/(y[d] - y[1]));
    }
  }
}

integral_deriv <- function(y, i, eps) {
  z = c(y, y[i])
  return(integral(z, eps));
}

upper_chull <- function(x,y) {
  Q = cbind(x,y)
  if (ncol(Null(t(Q) - Q[1,])) > 0) {
    Q = cbind(x,y+rnorm(length(y), mean=0, sd=1))
  }
  C = convhulln(Q)
  Cupper = 0
  Q1 = cbind(Q, rep(1,nrow(x)))
  for (i in 1:nrow(C)) {
    N = Null(t(cbind(x,rep(1,nrow(x)))[C[i,],]));
    if(ncol(N) > 0) { next; }
    K = Null(t(Q1[C[i,],]))
    if (K[ncol(Q)] < 1e-08 && K[ncol(Q)] > -1e-08) {
      next;
    }
    if (K[ncol(Q)] < 0) {
      K = -K
    }
    L = Q1 %*% K;
    if (sum(L) < 0) {
      Cupper = rbind(Cupper, C[i,])
    }
  }
  return(Cupper[2:nrow(Cupper),]);
}

sigmafunc <- function(x, y, w, eps) {
  convh = upper_chull(x,y)
  fval = -y%*% w;
  #computing the integral
  for (i in 1:nrow(convh)) {
    integ = integral(y[convh[i,]], eps);
    A = cbind(x[convh[i,],], rep(1,ncol(x)+1))
    fval = fval + integ*abs(det(A));
  }
  return(fval);
}

sigmafuncS <- function(x,w,eps) {
  s <- function(y){
    return(sigmafunc(x,y,w,eps));
  }
  return(s);
}

gradsigmafunc <- function(x,y,w, eps) {
  convh = upper_chull(x,y)
  gradval = -w;
  for (j in 1:length(y)) {
    for (i in 1:nrow(convh)) {
      a = -1;
      for (k in 1:ncol(convh)) {
        if (convh[i,k] == j) {
          a = k;
        }
      }
      if (a >= 1) {
        integ = integral_deriv(y[convh[i,]], k, eps);
        A = cbind(x[convh[i,],], rep(1,ncol(x)+1))
        gradval[j] = gradval[j] + integ*abs(det(A));
      }
    }
  }
  return(gradval);
}

gradsigmafuncS <- function(x,w,eps){
  s <- function(y){
    return(gradsigmafunc(x,y,w,eps));
  }
  return(s);
}

mini <- function(x,y) {
  res = x;
  for (i in 1:length(x)) {
    res[i] = min(x[i], y[i]);
  }
  return(res);
}

maxi <- function(x,y) {
  res = x;
  for (i in 1:length(x)) {
    res[i] = max(x[i], y[i]);
  }
  return(res);
}

minmaxx <- function(x) {
  flag = 1;
  mm = x;
  sortx = order(sort(x[,1], index.return = TRUE)$ix);
  sorty = order(sort(x[,2], index.return = TRUE)$ix);
  index_table = matrix(0, nrow(x), nrow(x));
  for (i in 1:nrow(x)) {
    index_table[sortx[i], sorty[i]] = i;
  }
  while(flag == 1) {
    flag = 0;
    mmnew = 0;
    old = 0;
    #print(mm);
    for (i in 1:(nrow(mm)-1)){
      for (j in max(old+1,i+1):nrow(mm)) {
        mi_ind_x = min(sortx[i], sortx[j]);
        ma_ind_x = max(sortx[i], sortx[j]);
        mi_ind_y = min(sorty[i], sorty[j]);
        ma_ind_y = max(sorty[i], sorty[j]);
        if(index_table[mi_ind_x,mi_ind_y] != 0 && index_table[ma_ind_x,ma_ind_y] != 0) {
          next;
        }
        if (sum(mm[i,] <= mm[j,]) == ncol(mm) ) {
          next;
        }
        if (sum(mm[j,] <= mm[i,]) == ncol(mm)) {
          next;
        }
        
        mi = mini(mm[i,], mm[j,]);
        ma = maxi(mm[i,], mm[j,]);
        foundmi = 0;
        foundma = 0;
        if (index_table[mi_ind_x, mi_ind_y] != 0) {
          foundmi = 1;
        }
        if (index_table[ma_ind_x, ma_ind_y] != 0) {
          foundma = 1;
        }
        if(foundmi == 0) {
          flag = 1;
          mmnew = rbind(mmnew, mi);
          sortx = c(sortx, mi_ind_x);
          sorty = c(sorty, mi_ind_y);
          index_table[mi_ind_x, mi_ind_y] = nrow(mm) + nrow(mmnew) - 1;
        }
        if(foundma == 0) {
          flag = 1;
          mmnew = rbind(mmnew, ma);
          sortx = c(sortx, ma_ind_x);
          sorty = c(sorty, ma_ind_y);
          index_table[ma_ind_x, ma_ind_y] = nrow(mm) + nrow(mmnew) - 1;
        }
      }
    }
    if (flag == 1) {
      old = nrow(mm);
      mm = rbind(mm, mmnew[2:nrow(mmnew),]);
    }
  }
  return(list(mm = mm, sortx = sortx, sorty = sorty, index_table = index_table));
}

cgm <- function(initialGuess, obj_fun, gr_obj_fun, Amat, w, maxSteps = 500,tol=1e-10) {
  guess = initialGuess/sum(initialGuess);
  #print("Adjusting so that integral equals 1.");
  integral_value = obj_fun(guess) + guess %*% w;
  #cat("obj_fun(guess) = ", obj_fun(guess), " and guess %*% w = ", guess %*% w);
  #cat("Found the integral value. It equals ", integral_value);
  guess = guess - log(integral_value);
  #print("Adjusted the new guess.");
  val = obj_fun(guess);
  #print("Found initial value.");
  oldval = 1e5
  improve = 1e4  
  i = 1
  while(improve > tol & i < maxSteps){
    #cat("Step ", i, " of CGM.\n");
    gr = gr_obj_fun(guess);
    #print("Gradient found. Starting LP.");
    furthest_point = Rglpk_solve_LP(obj = gr, mat = Amat, dir = rep(">=", nrow(Amat)), rhs = rep(0, nrow(Amat)),
                                    bounds=list(lower=list(ind = 1:ncol(Amat), val = rep(-30, ncol(Amat))), 
                                                upper=list(ind = 1:ncol(Amat), val = rep(6, ncol(Amat)))));
    #print("Furthest point found. Starting line search.");
    gamma = 2 / (2 + i);
    new_guess = gamma * furthest_point$solution + (1 - gamma) * guess;
    new_integral_value = obj_fun(new_guess) + new_guess %*% w;
    new_guess = new_guess - log(new_integral_value);
    new_val = obj_fun(new_guess);
    while (new_val > val) {
      gamma = 0.5 * gamma;
      new_guess = gamma * furthest_point$solution + (1 - gamma) * guess;
      new_integral_value = obj_fun(new_guess) + new_guess %*% w;
      new_guess = new_guess - log(new_integral_value);
      new_val = obj_fun(new_guess);
    }
    guess = new_guess;
    val = new_val;
    #if (i %% 10 == 0) {
    cat("iter: ", i, "objective value:", val, ".\n");
    #}
    improve = abs(val-oldval)
    oldval = val
    i = i+1
  }
  return(list(par = guess, value = val));
}

#--- functions specific to LLC
#-- find smallest polytrope that contains X
#input: X = matrix, row vectors = points
#output: matrix C, where x_i - x_j >= c_{ij}
polytrope <- function(X){
  points = cbind(0,X)
  m = dim(points)[2]
  C = matrix(0, nrow = m, ncol = m)
  for(i in 1:m){
    for(j in 1:m){
      C[i,j] = min(points[,i] - points[,j])
    }
  }
  return(C)
}

#-- find: ordinary vertices of the polytrope C
polVertex.2d <- function(C){
  mat = (C[,1] - C)[,c(2,3)]
  #add points there are new vertices
  x = mat[,1]
  y = mat[,2]
  if(y[1] > y[2]){
    mat = rbind(mat,c(x[2]+y[1]-y[2],y[1]))
  }
  if(x[1] > x[3]){
    mat = rbind(mat,c(x[1],y[3]+x[1]-x[3]))
  }
  if(x[2] < x[3] & y[2] > y[3]){
    mat = rbind(mat,c(x[2],y[3]))
  }
  return(mat)
}
#function: check if a point x is in the polytrope(C). (Append first index by 0)
isInPolytrope <- function(x,C){
  for(i in 1:2){
    for(j in i:3){
      if(x[i]-x[j] < C[i,j]){
        return(FALSE)
      }
      if(x[j]-x[i] < C[j,i]){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

#function: returns the weight of a point in list. 
#assume: points are NOT appended by 0
getWeight <- function(x,X,w){
id = fmatch(2,apply(apply(X,1,function(y) y == x),2,sum))
if(!is.na(id)){
  return(w[id])
}else{
  return(0)
}
}


#function to convert a point in the index vector to the point on the grid.
idx <- function(i,j,n1,n2){
  return(i+(j-1)*n1)
}

#--- function to generate inequalities, weights and added points
#COMMENT: produce inequalities of the kind
#blah >= 0.
get.ell.ineq <- function(X,w,C,V){
  #compute bounding box
  n1 = max(V[,1]) - min(V[,1])+1
  n2 = max(V[,2]) - min(V[,2])+1  
  #subtract off (1,1) since R indices start at (1,1)
  corner = c(0,apply(V,2,min))-c(0,1,1)
  #create a big matrix of n1xn2 points, and corresponding weight vectors
  ineq.mat = matrix(NA, nrow = 1, ncol = n1*n2)
  w.vec = rep(0,n1*n2)
  #function tested ok. 

  #--- list L-convex inequalities:
  #for each of the three triangle type: anchor point = minimum point. 
  #so for each point: generate at most 3 unique inequalities. 
  xlist = matrix(NA,nrow=1,ncol=2)
  relevant.idx = c()
  for(j in 1:n2){
    print(paste("j = ", j, "out of ", n2))
    for(i in 1:n1){
      x = corner + c(0,i,j)
      if(isInPolytrope(x,C)){
        xlist <- rbind(xlist, x[-1])
        relevant.idx <- append(relevant.idx,idx(i,j,n1,n2))
        w.vec[idx(i,j,n1,n2)] = getWeight(x[-1],X,w)
        x10 = x + c(0,1,0)
        x01 = x + c(0,0,1)
        x11 = x + c(0,1,1)
        x12 = x + c(0,1,2)
        x21 = x + c(0,2,1)
        xneighbor = matrix(c(x10,x01,x11,x12,x21),ncol =3,byrow = T)
        inC = apply(xneighbor,1,isInPolytrope,C=C)
        if(inC[1] && inC[2] && inC[3]){
          ineq = rep(0, n1*n2)
          ineq[idx(i,j,n1,n2)] = 1
          ineq[idx(i+1,j+1,n1,n2)] = 1
          ineq[idx(i+1,j,n1,n2)] = -1
          ineq[idx(i,j+1,n1,n2)] = -1
          ineq.mat <- rbind(ineq.mat,ineq)
        }      
        if(inC[2] && inC[3] && inC[4]){
          ineq = rep(0, n1*n2)
          ineq[idx(i,j,n1,n2)] = -1
          ineq[idx(i+1,j+2,n1,n2)] = -1
          ineq[idx(i,j+1,n1,n2)] = 1
          ineq[idx(i+1,j+1,n1,n2)] = 1
          ineq.mat <- rbind(ineq.mat,ineq)        
        }
        if(inC[1] && inC[3] && inC[5]){
          ineq = rep(0, n1*n2)
          ineq[idx(i,j,n1,n2)] = -1
          ineq[idx(i+2,j+1,n1,n2)] = -1
          ineq[idx(i+1,j,n1,n2)] = 1
          ineq[idx(i+1,j+1,n1,n2)] = 1
          ineq.mat <- rbind(ineq.mat,ineq)        
        }
      }
    }
  }
  xlist <- xlist[-1,]
  ineq.mat <- ineq.mat[-1,]
  #only keep the relevant variables
  ineq.mat <- ineq.mat[,relevant.idx]
  w.vec <- w.vec[relevant.idx]
  m <- list()
  m$ineq.mat <- ineq.mat
  m$xlist <- xlist
  m$w.vec <- w.vec
  return(m)
}

#--- interpret output
interp_ell_mle <- function(info, S) {
  mm = info$xlist
  min_coord = rep(0, ncol(mm));
  max_coord = rep(0, ncol(mm));
  for(i in 1:ncol(mm)){
    min_coord[i] = min(mm[,i]);
    max_coord[i] = max(mm[,i]);
  }
  x = seq(min_coord[1]-0.1, max_coord[1]+0.1, length.out = 100);
  y = seq(min_coord[2]-0.1, max_coord[2]+0.1, length.out = 100);
  z = matrix(0, ncol=100, nrow=100);
  triangs = upper_chull(mm, S$par);
  for (i in 1:length(x)) {
    for (j in 1:length(y)) {
      triang_num = 0;
      ker = c(0,0,0);
      flag = 0;
      for (k in 1:nrow(triangs)){
        ker = Null(rbind(mm[triangs[k,1],] - c(x[i], y[j]), mm[triangs[k,2],]- c(x[i], y[j]), mm[triangs[k,3],]- c(x[i], y[j])))[,1];
        ker = ker/sum(ker);
        if(ker[1] >=-1e-06 && ker[2] >=-1e-06 && ker[3] >=-1e-06) {
          triang_num = k;
          flag = 1;
          break;
        }
      }
      if (flag == 0) {
        z[i,j] = -Inf;
      } else {
        z[i,j] = ker %*% c(S$par[triangs[triang_num,1]], S$par[triangs[triang_num,2]], S$par[triangs[triang_num,3]]);
      }
    }
  }
  return(list(x = x, y = y, z = z));
}

