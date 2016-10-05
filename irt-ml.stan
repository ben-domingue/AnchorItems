data {
  int<lower=1> J;              // number of students
  int<lower=1> K;              // number of questions
  int<lower=1> N;              // number of observations
  int<lower=1,upper=J> jj[N];  // student for observation n
  int<lower=1,upper=K> kk[N];  // question for observation n
  int<lower=0,upper=1> y[N];   // correctness for observation n
  //
  int<lower=1,upper=2> g[J];   // group for observation J
}
transformed data {
  vector[2] zeros;
  zeros = rep_vector(0,2);
}
parameters {    
  real mu_alpha[2];                  //mean of group specific means
  real sigma_alpha[2];                  //sd of group specific means
  real alpha[J];               // ability of student j - mean ability
  matrix[K,2] beta_free;                // difficulty of question k by group
  //matrix sigma_beta[2,2];       //cov martrix for beta
  corr_matrix[2] Omega; // correlation of difficulties
  vector<lower=0>[2] tau; // scale of difficulties
}
transformed parameters {
  matrix[K,2] beta;
  real tmp_sum;
  for (i in 1:2) {
    tmp_sum = 0;
    for(k in 1:K) {
      tmp_sum = tmp_sum+beta_free[k,i];
    }
    for(k in 1:K) {
      beta[k,i] = beta_free[k,i]-tmp_sum/K;
    }
  }
}      
model {
  //hyperpriors
  to_vector(mu_alpha) ~ normal(0,1);   
  to_vector(sigma_alpha) ~ cauchy(0,25);
  #see 66/67 of stan manual
  tau ~ cauchy(0,2.5);
  Omega ~ lkj_corr(2);
  for (j in 1:J) //prior for ability
    alpha[j] ~ normal(mu_alpha[g[j]],sigma_alpha[g[j]]);
  for (k in 1:K) 
    beta_free[k,] ~ multi_normal(zeros,quad_form_diag(Omega,tau));

  for (n in 1:N)
    y[n] ~ bernoulli_logit(alpha[jj[n]] - beta[ kk[n] , g[jj[n]] ] );

}


// transformed parameters {
//   vector[2] beta_means;
//   real tmp_sum;
//   for (i in 1:2) {
//     tmp_sum = 0;
//     for(k in 1:K) {
//       tmp_sum = tmp_sum+beta[k,i];
//     }
//     beta_means[i] = tmp_sum/K;
//     }      
//   for (i in 1:2) {
//     for (k in 1:K) {
//       beta[k,i]=beta[k,i]-beta_means[i];
//     }
//   }
//   //print(beta);
// }

