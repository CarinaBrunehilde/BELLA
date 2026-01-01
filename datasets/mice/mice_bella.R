library(geepack)
library(tidyverse)
library(abind)
library(posterior)
library(bayesplot)
library(kableExtra)

load("datasets/mice/miceweight.RData")

mice_w <- pivot_wider(miceweight, names_from = week, values_from = weight)
miceweight <- miceweight %>% mutate(week2 = week-4) %>%
  mutate(interTreat=if_else(week2 <0, "C", treat))

X_mice <- model.matrix(~interTreat+week+interTreat*week, data=miceweight)
id <- as.numeric(rep(1:52, each=11))
X_mice <- data.frame(id, X_mice)
y_mice <- (miceweight$weight)/10
y_mat <- matrix(y_mice,52,11, byrow = T)
cov(y_mat)
var_mice <- diag(cov(y_mat))
mean(var_mice)
var(var_mice)
class(y_mice)
n  <- length(unique(X_mice$id))   
ni <- as.vector(table(X_mice$id))   
N  <- sum(ni)  

                            
mice_BELLA <- BELLA_chain_thin(X = X_mice, y_mice, nsim = 12000, nwarm=2000, thin=10, 
                                nchain = 3, L = 20, a = 4.2, b = 0.73, epsilon = .002)


fit <- mice_BELLA  
fit$acceptance.rate

## For betas, change list for array (iter, par, chain)
betas_array <- abind(fit$beta, along = 3)
## Change for (iter, chain, par) - because this is how posterior package works
betas_array <- aperm(betas_array,c(1,3,2))
## Add sigma and rho
parametros <- abind(betas_array,fit$sigma,fit$rho, along=3)
par_matrix <- as_draws_matrix(parametros)
colnames(par_matrix) <- c("beta[0]", "beta[1]", "beta[2]", "beta[3]", "beta[4]",
                          "beta[5]","sigma^2", "rho")


## get results
draws <- as_draws_array(parametros)
summarise_draws(par_matrix)
summarize_draws(par_matrix, "mean", "mcse_mean", function(x) quantile(x, probs = c(0.025, 0.975)))
summarize_draws(par_matrix, "mean", "mcse_mean", "sd")
mice <- summarize_draws(par_matrix, "mean", "sd","mcse_mean", "rhat", "ess_bulk", function(x) quantile(x, probs = c(0.025, 0.975)))
xtable(mice)

##Plots
color_scheme_set(beyonce_palette(101))
labels_params <- c(
  "beta[0]"  = expression(beta[0]),
  "beta[1]"  = expression(beta[1]),
  "beta[2]"  = expression(beta[2]),
  "beta[3]"  = expression(beta[3]),
  "beta[4]"  = expression(beta[4]),
  "beta[5]"  = expression(beta[5]),
  "sigma^2"  = expression(sigma^2),
  "rho"      = expression(rho)
)

mcmc_dens_overlay(par_matrix)

p <- mcmc_acf_bar(par_matrix, lags=50)
p + facet_grid(Chain ~ Parameter, labeller = labeller(Parameter = label_parsed)) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "gray90", color = "gray40"),
        strip.text = element_text(face = "bold"))

trace <- mcmc_trace(par_matrix)
trace +  facet_wrap(~parameter, labeller = labeller(parameter = label_parsed), 
                    scales = "free_y")+
  theme_minimal() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "gray90", color = "gray40"),
        strip.text = element_text(face = "bold"))





