library(mvtnorm)
# True paameter
alpha = 0.;
gamma = 0.;
beta = -1.;
delta = 0.9;
r = -0.7;
N = 1000; # data number
# Generate data
s = 500; # simulation number

estimate <- function(gamma, delta, alpha, beta, s, N, r){
  all_beta_ols = array();
  all_beta_2sls = array();
  all_alpha_ols = array();
  all_alpha_2sls = array();
  for (i in 1:s){
    # generate data
    uv = rmvnorm(N, mean=c(0,0), sigma=matrix(c(1,r,r,1), 2));
    u = uv[,1];
    v = uv[,2];
    p1 = runif(1);
    p2 = runif(1,min=0,max = 1 - p1);
    p3 = 1 - p1 - p2;
    z = sample(c(-1,0,1), N, prob = c(p1, p2, p3),replace = TRUE);
    y2 = (-gamma - delta * z) < v;
    y2 = as.integer(y2);
    y1 = (-alpha - beta * y2) < u;
    y1 = as.integer(y1);
    # estimate ols and 2sls
    beta_ols = cov(y1,y2) / var(y2);
    alpha_ols = mean(y1) - beta_ols * mean(y2);
    beta_2sls = cov(y1,z) / cov(y2,z);
    alpha_2sls = mean(y1) - beta_2sls * mean(y2);
    # append
    all_beta_ols[i] = beta_ols;
    all_alpha_ols[i] = alpha_ols;
    all_beta_2sls[i] = beta_2sls;
    all_alpha_2sls[i] = alpha_2sls;
  }
  return(list(all_beta_ols, all_alpha_ols, all_beta_2sls, all_alpha_2sls, y1, y2, z))
}

estimators = estimate(gamma, delta, alpha, beta, s, N, r);
all_beta_ols = unlist(estimators[1]);
all_alpha_ols = unlist(estimators[2]);
all_beta_2sls = unlist(estimators[3]);
all_alpha_2sls = unlist(estimators[4]);
y1 = unlist(estimators[5]);
y2 = unlist(estimators[6]);
z = unlist(estimators[7]);

# identify set
identify_set <- function(z, y1, y2){
  zs = c(-1,0,1);
  lb0s = array();
  ub0s = array();
  lb10s = array();
  ub10s = array();
  for (i in 1:3){
    Z = zs[i]
    ind = which(z==Z);
    n = length(ind);
    y1_z = y1[ind];
    y2_z = y2[ind];
    lb0 = sum(y1_z) / n;#lower bound of rho_0
    ub0 = 1 - sum((y1_z+y2_z) == 0) / n;#upper bound of rho_0
    lb10 = sum((y1_z + y2_z) == 2)/n;# lower bound of rho_1 + rho_0
    ub10 = sum(y1_z)/n;# upper bound of rho_0 + rho_1
    lb0s[i] = lb0;
    ub0s[i] = ub0;
    lb10s[i] = lb10;
    ub10s[i] = ub10;
  }
  lb0 = max(lb0s);
  ub0 = min(ub0s);
  lb10 = max(lb10s);
  ub10 = min(ub10s);
  ub1 = ub10 - lb0;
  lb1 = lb10 - ub0;
  
  return(list(lb0, ub0, lb1, ub1))
}

identify_set_pos <- function(z, y1, y2){
  zs = c(-1,0,1);
  lb0s = array();
  ub0s = array();
  lb10s = array();
  ub10s = array();
  for (i in 1:3){
    Z = zs[i]
    ind = which(z==Z);
    n = length(ind);
    y1_z = y1[ind];
    y2_z = y2[ind];
    lb0 = sum(y1_z[which(y2_z == 0)]) / n;#lower bound of rho_0
    ub0 =  sum(y1_z) / n;#upper bound of rho_0
    lb10 =  sum(y1_z) / n;# lower bound of rho_1 + rho_0
    ub10 = sum(y2_z[which(y1_z == 0)])/n;# upper bound of rho_0 + rho_1
    lb0s[i] = lb0;
    ub0s[i] = ub0;
    lb10s[i] = lb10;
    ub10s[i] = ub10;
  }
  lb0 = max(lb0s);
  ub0 = min(ub0s);
  lb10 = max(lb10s);
  ub10 = min(ub10s);
  ub1 = ub10 - lb0;
  lb1 = lb10 - ub0;
  
  return(list(lb0, ub0, lb1, ub1))
}



identify_bound = identify_set(z,y1,y2);
lb0 = as.numeric(identify_bound[1]);
ub0 = as.numeric(identify_bound[2]);
lb1 = as.numeric(identify_bound[3]);
ub1 = as.numeric(identify_bound[4]);

x_ = c(lb0,lb0,ub0,ub0);
y_ = c(lb1,ub1,ub1 - (ub0 - lb0),lb1 - (ub0-lb0));
plot(x_, y_, type = "n",xlim=c(0,1),ylim=c(-1,1),xlab = 'rho_0',ylab = 'rho_1');
polygon(x_ ,y_);


if (r > 0){
  identify_bound = identify_set_pos(z,y1,y2);
  lb0 = as.numeric(identify_bound[1]);
  ub0 = as.numeric(identify_bound[2]);
  lb1 = as.numeric(identify_bound[3]);
  ub1 = as.numeric(identify_bound[4]);
  x_ = c(lb0,lb0,ub0,ub0);
  y_ = c(lb1,ub1,ub1 - (ub0 - lb0),lb1 - (ub0-lb0));
  par(new = TRUE);
  plot(x_, y_, type = "n",xlim=c(0,1),ylim=c(-1,1),xlab = 'rho_0',ylab = 'rho_1');
  polygon(x_ ,y_)
}

par(new = TRUE);
plot(all_alpha_ols, all_beta_ols,col="blue",xlim=c(0,1),ylim=c(-1,1),xlab = 'rho_0',ylab = 'rho_1');
par(new = TRUE);
plot(all_alpha_2sls, all_beta_2sls,col="red",xlim=c(0,1),ylim=c(-1,1),xlab = 'rho_0',ylab = 'rho_1');
legend(0.8,0.8,legend = c('ols','2sls'),fill=c("blue","red"))

