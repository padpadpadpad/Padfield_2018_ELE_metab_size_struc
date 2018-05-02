
data{
    int N;                              // number of samples
    int K;                              // length of masses 
    int mass_n[N];                      // position of last mass of each sample / number of organisms
    int first_mass[N];                  // position of first mass of each sample
    vector[N] ln_rate;                  // rate
    vector[N] temp;                     // temperature
    vector[K] mass;                     // masses (log)
    vector[N] df;                       // dilution factor (log)
}
parameters{  
    real<lower=0> sigma;      // standard deviation
    real ln_r0;               // community normalisation constant
    real bT;                  // effect of temperature
    real a;                   // size scaling exponent
    
}
transformed parameters{
    vector[N] ln_biomass;     // define ln_biomass
    vector[N] mu;             // define mu
    
    // calculate log biomass
    for(i in 1:N){
    ln_biomass[i] = log(df[i]) + log_sum_exp(a * segment(mass, first_mass[i], mass_n[i]));
    }  
    mu = ln_r0 + ln_biomass + bT*temp;

}
model{

    ln_r0 ~ normal(0,10);
    bT ~ normal(0,10);
    a ~ normal(0,10);
    sigma ~ cauchy(0,2);

    // likelihood regression for size dependence of rate
    
    ln_rate ~  normal(mu, sigma);

}
generated quantities{

    // log likelihood
    vector[N] log_lik;

    // create temperature corrected rate from data
    vector[N] temp_cor_mu;
    vector[N] temp_cor_sim;
    // create size corrected biomass
    vector[N] size_cor_biom_mu;
    vector[N] size_cor_biom_sim;
    // create individual metabolic rate
    vector[N] ind_mu;

    for(i in 1:N){
        temp_cor_mu[i] = ln_rate[i] - log(exp(bT*temp[i]));
        temp_cor_sim[i] = normal_rng(ln_rate[i] - log(exp(bT*temp[i])), sigma);
        size_cor_biom_mu[i] = ln_biomass[i] + ln_r0;
        size_cor_biom_sim[i] = normal_rng(ln_biomass[i] + ln_r0, sigma);
        log_lik[i] = normal_lpdf(ln_rate[i] | mu[i], sigma);
        ind_mu[i] = log(exp(mu[i])/(mass_n[i]*df[i]));
    }
}
