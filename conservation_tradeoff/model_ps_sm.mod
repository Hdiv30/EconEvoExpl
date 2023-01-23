param eps := 10^(-6);
param gamma_param := 0.00000722;  # length to weight factor SMALK 1971-2019
param alpha := 3.08695688;  # length to weight exponent SMALK 1971-2019
# 5. natural mortality
param mu = 9.767037645676052e-16;
param kappa = 28.632718788724997;
param nu = 1.0359528166246261;
# 6. fishing
param omega = 0.15;  # harvesting accuracy parameter
# 7. mean m see 'data'
param l_sigma = 21.5;  # standard deviation for maturity
# 8. Beverton-Holt:
param Nmax = 2848.715979;
param ssbhalf = 419.69765808;

# 10. cost function 
param chi := 0.239;  # cost dependence on accessible biomass
param c = 1.2097254364874803;  # stock independent costs
# param tau = -0.010701319368378155;
param tau;
# SET
# 11. expenditures
param psi = 16.08;  # not expenditures but scaling parameter of fish consumption
param eta = 3.896;  # elasticity of subsitution between fish consumption and other consumption

param M := 10;  # new number of genotypes
# -------- change only these --------
param noEvo;
param T; #:=50; # time horizon
param hyper;
param delta; # := 1; # discount factor
param smoother;  # for optimisation
param xi = 8;  # elasticity of substitution for different size classes

param theta := 0.4993370985452277; 
param beta_preferences {s in 1..6};  # size class preferences

# -------- pre-defined rates and parameters --------
param g_sm {s in 1..6, m in 1..M}; # size and maturation specific transition proportion (growth)
param n_0 {s in 1..6, m in 1..M};  # initial population in millions
param hs_guess {t in 0..100, s in 1..6}; # steady state ps harvest in million numbers
param n_final_ps {s in 1..6, m in 1..M};  # steady state population in millions for ps
param n_final_default {m in 1..M, s in 1..6};  # steady state final value
param n_tm_guess {t in 1..51, m in 1..M};  # guessed population with genotype m for 5 years (Evo, delta=1, tau=0)
param n_1963 {s in 1..6, m in 1..M};  # population in 1963 (start of data)
param H_B_guess {t in 0..1000};
param sigmas_guess {t in 0..1000};
param h_s_eq {s in 1..6};  # for ps
param n_eq {s in 1..6, m in 1..M};  # population in eq. roughly
param w_class {s in 1..6};  # minimium weight (kg) in that s-class
param w_mean {s in 1..6}; #.. corresponding weight
param l_mean {s in 1..6} = (w_mean[s] / gamma_param ) ^ (1 / alpha);  # mean length in that class and..
param l_m {m in 1..M};  # mean maturation length 40, 45, ..., 80
param reproduction_factor {s in 1..6, m in 1..M};  # cdf. proportion that is mature at size s for each m-class
param C {t in 0..T} = c * exp(t * tau);  # t=0 is 2019

# for evoDistribution:
param n_size{s in 1..6} = sum{m in 1..M} n_0[s, m];  # population by size class only
param n_evoDistribution {s in 1..6, m in 1..M};  # initial population for a given mean m

# for noEvo
param freq_m {m in 1..M};


# stock numbers
var n {t in 0..(T+1), s in 1..6, m in 1..M}>=0;  # > for eco-evo eq that is not 0
var B_t {t in 0..T} = sum{s in 1..6, m in 1..M} (w_mean[s]*n[t,s,m]);

# for target:
var mean_lm {t in 0..(T+1)} = (sum{m in 1..M} (sum{s in 1..6} n[t, s, m]) * l_m[m]) / (sum{s in 1..6, m in 1..M} n[t, s, m]);
param t_target;
param target_lm;

# --------- 1. reproduction  ---------
var ssb_s {t in 0..T, s in 1..6, m in 1..M} = (w_mean[s] * n[t,s,m] * reproduction_factor[s,m]); 
var ssb {t in 0..T} = sum{s in 1..6, m in 1..M} ssb_s[t, s, m]; 
var ssb_scaled_s {t in 0..T, s in 1..6, m in 1..M} = (w_mean[s]^hyper * n[t,s,m] * reproduction_factor[s,m]);
var ssb_scaled_total {t in 0..T} = sum{s in 1..6, m in 1..M} ssb_scaled_s[t, s, m];  
var bev_holt {t in 0..T} = Nmax * ssb[t] / (ssbhalf + ssb[t]);  

var r_sm_i {t in 0..T, s in 1..6, m in 1..M} = bev_holt[t] * ssb_scaled_s[t, s, m] / ssb_scaled_total[t];  # with hyperallometric scaling
var r_sm_Evo {t in 0..T, s in 1..6, m in 1..M} = (1 - theta) * r_sm_i[t,s,m] + (if m>1 then theta/2*r_sm_i[t,s,m-1]) + (if m<M then theta/2*r_sm_i[t,s,m+1]) + (if m=1 then theta/2*r_sm_i[t,s,1]) + (if m=M then theta/2*r_sm_i[t,s,M]);
# if noEvo:
var r_sm {t in 0..T, s in 1..6, m in 1..M} = if (noEvo = 1) then (freq_m[m] * (sum{j in 1..M} r_sm_i[t, s, j])) else r_sm_Evo[t, s, m];
var recruits {t in 0..T, m in 1..M} = sum{s in 1..6} r_sm[t,s,m];
var population {t in 0..T, s in 1..6, m in 1..M} = if s=1 then n[t, 1, m] + recruits[t, m] else n[t,s,m];

# --------- 2. growth ---------
var adults_Evo {t in 0..T, s in 2..6, m in 1..M} = population[t,s,m]*(1-g_sm[s,m]) + population[t,s-1,m]*(g_sm[s-1,m]);
var grownRecruits_Evo {t in 0..T, m in 1..M} = population[t, 1, m] * (1 - g_sm[1, m]);
# if noEvo:
var adults {t in 0..T, s in 2..6, m in 1..M} = if (noEvo = 1) then  (freq_m[m] * (sum{j in 1..M} adults_Evo[t, s, j])) else adults_Evo[t, s, m];
var grownRecruits {t in 0..T, m in 1..M} = if (noEvo = 1) then  (freq_m[m] * (sum{j in 1..M} grownRecruits_Evo[t, j])) else grownRecruits_Evo[t, m];

# --------- 3. natural mortality ---------
param d_s {s in 1..6} := exp(-(mu+kappa*(l_mean[s]^(-nu))));
var survivors_natural {t in 0..T, s in 1..6, m in 1..M} = if s=1 then grownRecruits[t,m]*d_s[1] else adults[t,s,m]*d_s[s];
var mean_weight_before_harvest {t in 0..T} = sum{s in 1..6} (w_mean[s] * (sum{m in 1..M}survivors_natural[t, s, m]) / (sum{j in 1..6, m in 1..M}survivors_natural[t, j, m]));

# --------- 4. fishing ---------
var h_sm {t in 0..T, s in 1..6, m in 1..M}>= 0;
var h_s {t in 0..T, s in 1..6}= sum{m in 1..M} h_sm[t, s, m];
var a_N {t in 0..T, s in 1..6, m in 1..M} = survivors_natural[t,s,m];
var a_B {t in 0..T, s in 1..6} = sum{m in 1..M} (a_N[t,s,m] * w_mean[s]);
var a_Bm {t in 0..T, m in 1..M} = sum{s in 1..6} (a_N[t,s,m] * w_mean[s]);
var A_B {t in 0..T} = sum{s in 1..6} a_B[t, s];  # in 1000 tons
# distribute size class harvest according to m-distribution in that class (all genotypes are harvested equally)
var h {t in 0..T, s in 1..6} = w_mean[s] * h_s[t,s];  # in 1000 tons
var h_m {t in 0..T, m in 1..M} = sum{s in 1..6} h_sm[t, s, m];
var h_mw {t in 0..T, m in 1..M} = sum{s in 1..6} (w_mean[s] * h_sm[t, s, m]);
var H_B {t in 0..T} = sum{s in 1..6} h[t,s];  # total harvest in 1000 tons
var survivors {t in 0..T, s in 1..6, m in 1..M} = survivors_natural[t,s,m] - h_sm[t,s,m];

# --------- 5. economics ---------
var costs {t in 0..T} = c * H_B[t] ;  # producer supply costs
var utility_from_fish {t in 0..T}= (sum{s in 1..6} (beta_preferences[s] * (h[t, s] + eps)^((xi - 1) / xi)))^(xi/(xi-1));
var utility {t in 0..T} = psi * (eta / (eta-1))  * utility_from_fish[t]^((eta-1)/eta); # consumer utility
var surplus {t in 0..T} = utility[t] - costs[t];

# --------- optimise ---------
maximize discounted_smoothed_surplus: sum{t in 0..T} (delta^t * surplus[t])^(smoother);  # - sum{t in 2..T} penalty[t];  
maximize nothing: 0;
maximize maximise_lm: mean_lm[t_target];

subject to population_dynamics {t in 1..(T+1), s in 1..6, m in 1..M}: n[t,s,m] = survivors[t-1,s,m];
# subject to steady_state {t in 1..(T+1), s in 1..6, m in 1..M}: n[t,s,m] = n[t-1,s,m];
# subject to steady_state_H {t in 1..T}: H_B[t] = H_B[t-1];
# subject to initial_population {s in 1..6, m in 1..M}: n[0,s,m] = n_0[s,m];
# subject to initial_population_evoDistribution {s in 1..6, m in 1..M}: n[0,s,m] = n_evoDistribution[s,m];
subject to evo_target {t in t_target..(T+1)}: mean_lm[t]+eps>=target_lm;  # all future lm after t_target shall be larger than a set target
subject to evo_target_nodecrease {t in t_target..(T+1)}: mean_lm[t]+eps>=mean_lm[t-1];  # all future lm after t_target shall be larger than a set target
