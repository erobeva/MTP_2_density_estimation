install.packages("MASS")
install.packages("geometry")
install.packages("scatterplot3d")
install.packages("plot3D")
install.packages("fields")
install.packages("sos")
install.packages("utils")
install.packages("rgl")
install.packages("lpSolve")
install.packages("Rglpk")
install.packages("matlib")
install.packages("rmutil")
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

# Computes (an approximation of) the integral of 
# exp(tent function) above a simplex where the tent
# function is linear
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

# Computes (an approximation of) the integral of
# exp(tent function)
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

# Computes the derivative of the integral of
# exp(tent function)
integral_deriv <- function(y, i, eps) {
  z = c(y, y[i])
  return(integral(z, eps));
}

# Computes the upper convex hull of points (x_i,y_i)
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

# Computes the objective function to be minimized
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

# Returns the objective function to be minimized
# for specified values of x, w, epsilon
sigmafuncS <- function(x,w,eps) {
  s <- function(y){
    return(sigmafunc(x,y,w,eps));
  }
  return(s);
}

# Computes the derivative of the objective funcion
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

# Returns the integral of the objective function
# for specified values of x, w, epsilon
gradsigmafuncS <- function(x,w,eps){
  s <- function(y){
    return(gradsigmafunc(x,y,w,eps));
  }
  return(s);
}

# Computes the coordinate-wise minimum of x and y
mini <- function(x,y) {
  res = x;
  for (i in 1:length(x)) {
    res[i] = min(x[i], y[i]);
  }
  return(res);
}

# Computes the coordinate-wise maximum of x and y
maxi <- function(x,y) {
  res = x;
  for (i in 1:length(x)) {
    res[i] = max(x[i], y[i]);
  }
  return(res);
}

# Computes the min-max convex hull of a set of points x
# in R^2
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

# Creates the matrix of supermodular inequalities
create_Amat <- function(m1) {
  mm = m1$mm;
  index_table = m1$index_table;
  sortx = m1$sortx;
  sorty = m1$sorty;
  n = nrow(index_table);
  m = nrow(mm);
  Amat = 0;
  flag = 0;
  for (i in 1:(n-1)) {
    #print(i);
    for (j in 1:(n-1)) {
      if(index_table[i, j] == 0 || index_table[i+1,j] ==0 || index_table[i,j+1] == 0 || index_table[i+1,j+1] == 0) {
        next;
      }
      # print(j);
      v = rep(0,m)
      v[index_table[i+1,j]] = -1;
      v[index_table[i,j+1]] = -1;
      v[index_table[i,j]] = 1;
      v[index_table[i+1,j+1]] = 1;
      Amat = rbind(Amat, v);
      flag = 1;
    }
  }
  if(flag == 1) {
    return(Amat[2:nrow(Amat),]);
  } else {
    return(0);
  }
}

# Computes the MTP_2 and log-concave MLE using the spg
# optimization algorithm
mtp2_opt <- function(x, m1,w) {
  eps=0.001;
  mm = m1$mm;
  w = c(w, rep(0, nrow(mm)-nrow(x)))
  Amat = create_Amat(m1);
  initialGuess = c(rep(0, nrow(x)), rep(1, nrow(mm)-nrow(x)));
  if(!is.vector(Amat) && !is.matrix(Amat)) {
    S = spg(par = initialGuess, fn = sigmafuncS(mm,w, eps), gr = gradsigmafuncS(mm,w, eps), control = list(checkGrad.tol = 100));
    return(S);
  } else {
    if(is.vector(Amat)){
      b=0
      meq = 0
    }
    if(is.matrix(Amat)){
      b = rep(0, nrow(Amat));
      meq = 0;
    }
    for (i in 1:3) {
      S = spg(par = initialGuess, fn = sigmafuncS(mm,w, eps), gr = gradsigmafuncS(mm,w, eps), project = "projectLinear", 
              projectArgs=list(A=Amat, b=b, meq=meq), control = list(checkGrad.tol = 1000, maxit=2000));
      integral_value = S$value + S$par %*% c(w, rep(0, length(S$par) - length(w)));
      S$par = S$par - log(integral_value);
      S$value = sigmafunc(mm, S$par, w, eps);
      initialGuess = S$par;
    }
    #S = spg(par = initialGuess, fn = sigmafuncS(mm,w, eps), project = "projectLinear", 
    #        projectArgs=list(A=Amat, b=b, meq=meq), control = list(checkGrad.tol = 1000, maxit=3000));
    return(S);
  }
}

# Prepares the computed MLE for plotting
interp_mtp2_mle <- function(m1, S) {
  mm = m1$mm;
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

# The conditional gradient method
cgm <- function(initialGuess, obj_fun, gr_obj_fun, Amat, w, maxSteps = 500) {
  guess = initialGuess/sum(initialGuess);
  integral_value = obj_fun(guess) + guess %*% w;
  guess = guess - log(integral_value);
  val = obj_fun(guess);
  for (i in 1:maxSteps) {
    gr = gr_obj_fun(guess);
    furthest_point = Rglpk_solve_LP(obj = gr, mat = Amat, dir = rep(">=", nrow(Amat)), rhs = rep(0, nrow(Amat)),
                                    bounds=list(lower=list(ind = 1:ncol(Amat), val = rep(-30, ncol(Amat))), 
                                                upper=list(ind = 1:ncol(Amat), val = rep(6, ncol(Amat)))));
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
    cat("iter: ", i, "objective value:", val, ".\n");
  }
  return(list(par = guess, value = val));
}

# Finds an interior point in the convex body specified by the
# inequalities Amat.x >= 0
find_interior_point <- function(m1, Amat) {
  eps = 0.5;
  index_table = m1$index_table;
  n = nrow(index_table);
  guess = rep(10, ncol(Amat));
  for(i in 1:(n-1)) {
    for (j in 1:(n-1)) {
      if(index_table[i,j] == 0 || index_table[i+1,j] == 0 || index_table[i+1,j+1] == 0 || index_table[i, j+1] == 0) {
        next;
      }
      alpha = guess[index_table[i,j]] + guess[index_table[i+1,j+1]] - guess[index_table[i, j+1]]-guess[index_table[i+1,j]];
      if (alpha <= 0){
        guess[index_table[i+1,j+1]] = guess[index_table[i+1,j+1]] - alpha + eps;
      }
    }
  }
  if(sum(Amat %*% guess > 0) < nrow(Amat)){
    return(0);
  }
  return(guess);
}

# Computes the MTP_2 and log-concave MLE using the conditional
# gradient method
mtp2_cgm <- function(x, m1, w) {
  eps = 0.001;
  Amat = create_Amat(m1);
  mm = m1$mm;
  w = c(w, rep(0, nrow(mm)-nrow(x)));
  initialGuess = find_interior_point(m1,Amat);
  S = cgm(initialGuess, sigmafuncS(mm,w, eps), gradsigmafuncS(mm,w, eps), Amat, w);
  return(S);
}


###################################
# To use the functions:
###################################
npoints = 5;
x = matrix(rnorm(2*npoints),ncol=2)
par(mar = c(1,1,1,1))
plot(x, xlim = c(-6, 6), ylim = c(-6,6))
m1 = minmaxx(x)
points(m1$mm)
w = 1/npoints*rep(1, npoints)
#######
# Running the optimization problem:
#######
lcd = mlelcd(x,w)
S = mtp2_cgm(x, m1, w)
S$value + S$par %*% c(w, rep(0, length(S$par) - length(w)))
g <- interp_mtp2_mle(m1, S);

##################################
# To plot the solution:
##################################
M = mesh(g$x, g$y)
scatterplot3d(c(M$x), c(M$y), c(g$z), highlight.3d=TRUE)
persp3d(g$x,g$y,exp(g$z), zlim=c(-0.0001,0.001), col=rainbow(40000))
points3d(cbind(L1$X/10, exp(S$par[1:nrow(L1$X/10)])))


##################################
# To save the experiment:
##################################
LS = list(x, m1, S, g);
filename=paste(dir_name,"npoints",npoints,"date",toString(Sys.time()), ".Rdata", sep="", collapse=NULL);
save(LS, file=filename);

load(file="npoints100date2017-11-09.Rdata")


##################################
# Finding the squared L_2 distance between densities:
##################################
squared_L2 <- function(p1, p2, L, R, N) {
  xcoord = seq(L, R, length.out = N);
  ycoord = seq(L, R, length.out = N);
  distance = 0;
  for (i in 1:N) {
    for (j in 1:N) {
      difference = (p1(c(xcoord[i], ycoord[j])) - p2(c(xcoord[i], ycoord[j])));
      distance = distance + difference*difference;
    }
  }
  return(distance*(R-L)*(R-L)/(N*N));
}

hellinger <- function(p1, p2, L, R, N) {
  xcoord = seq(L, R, length.out = N);
  ycoord = seq(L, R, length.out = N);
  distance = 0;
  for (i in 1:N) {
    for (j in 1:N) {
      difference = (sqrt(p1(c(xcoord[i], ycoord[j]))) - sqrt(p2(c(xcoord[i], ycoord[j]))));
      distance = distance + difference*difference;
    }
  }
  return(distance*(R-L)*(R-L)/(N*N));
}

standard2dGaussian <- function(x) {
  return(1/(2*pi)*exp(-(x[1]*x[1]+x[2]*x[2])/2));
}

standard2dGaussian2arg <- function(x1, x2) {
  return(standard2dGaussian(c(x1, x2)));
}

evaluateTentFunction <- function(x, y, z) {
  val = 0;
  z1 = c(z, 1);
  C = upper_chull(x,y);
  for (i in 1:nrow(C)) {
    M = t(x[C[i,],]);
    M = rbind(M, rep(1,ncol(C)));
    coords = solve(M) %*% z1;
    flag = TRUE;
    for (j in 1:ncol(C)) {
      if (coords[j] < 0) {
        flag = FALSE;
        break;
      }
    }
    if (flag == TRUE) {
      return(exp(t(coords[1:ncol(C)]) %*% y[C[i,]]));
    }
  }
  return(0);
}

tentFunction <- function(x, y) {
  s <- function(z) {
    return(evaluateTentFunction(x, y, z));
  }
  return(s);
}

L2differenceFunction <- function(x,y) {
  s <- function(a,b) {
    difference = 1/(2*pi)*exp(-(a*a+b*b)/2) - tentFunction(x, y)(c(a,b));
    return(difference*difference);
  }
  return(s);
}

squared_L2(standard2dGaussian, tentFunction(x, lcd$logMLE), -10, 10, 50)
squared_L2(standard2dGaussian, tentFunction(m1$mm, S$par), -10, 10, 50)

int2(L2differenceFunction(m1$mm, S$par)) # too slow
int2(L2differenceFunction(x,lcd$logMLE))


