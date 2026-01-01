library(tidyverse)
library(posterior)
library(dplyr)
library(abind)
library(bayesplot)

dados <- read.table("datasets/cd4/actg398.dat", header = TRUE)

hiv <- dados %>%
  mutate(lcd4=log(cd4+1), 
         id = as.numeric(factor(patid)), 
         logrna = str_replace_all(logrna, "\\[|\\]", ""),
         logrna = as.numeric(logrna), 
         studyday = as.numeric(studyday),
         stday_scaled = scale(studyday),
         week = ceiling(studyday/7),
         weekc = week - mean(week),
         trtarm = factor(trtarm)) %>% 
  dplyr::select(id, everything()) %>%
  drop_na()

## Number of subjects and time-points of each subject
n  <- length(unique(hiv$id))  
ni <-  table(hiv$id)

## Taking subjects with only 1 observation
obs1 = names(ni[ni == 1])  
hiv = hiv %>% filter(!id %in% obs1)

hiv <- hiv %>% 
  mutate(id = as.integer(factor(id)))

## The new numbers of subjects and time points of each subjects 
n  <- length(unique(hiv$id))
ni <- as.vector(table(hiv$id))


var_by_time <- hiv %>% group_by(week) %>%  summarise(n(),
                                                     vartime=var(lcd4))


##---------------------------------------------------------------------------
##BELLA estimation for cd4 data from Biometrika paper - just for treatment D
##--------------------------------------------------------------------------
hivD <- hiv %>% filter(trtarm=="D") %>% 
  mutate(id = as.numeric(factor(patid)))

n  <- length(unique(hivD$id))
ni <- as.vector(table(hivD$id))



var_by_timeD <- hivD  %>% group_by(week) %>%  
                          summarise(n(),
                                    vartime=var(lcd4))

##Time as day
X_hiv <- model.matrix(~stday_scaled+logrna, data=hivD)
id    <- hivD$id
X_hiv <- cbind(id, X_hiv)

Y_hiv <- hivD %>% dplyr::select(lcd4) 
Y_hiv <- unlist(Y_hiv)

cd4_D_day<- BELLA_chain_thin(X=X_hiv,y=Y_hiv,nsim=6000,nwarm=1000, thin=5 , 
                             nchain=3, epsilon=0.005, L=10, a=a,b=b)
saveRDS(cd4_D_day, "data sets/cd4/cd4_qiu/cd4_D_day.rds")

##More iteractions, bigger epsilon and new hyperparameters
cd4_D_BELLA<- BELLA_chain_thin(X=X_hiv,y=Y_hiv, nsim=12000, nwarm=3000, thin=9 , 
                             nchain=3, epsilon=0.0067, L=10, a=2.1,b=1.1)


### --------------------------------------------------------
### Results for group D - cd4_D_BELLA
### --------------------------------------------------------
fit <- cd4_D_BELLA
fit$acceptance.rate
## For betas, change list for array (iter, par, chain)
betas_array <- abind(fit$beta, along = 3)
## Change for (iter, chain, par) - because this is how posterior package works
betas_array <- aperm(betas_array,c(1,3,2))
## Add sigma and rho
parametros <- abind(betas_array,fit$sigma,fit$rho, along=3)
par_matrix <- as_draws_matrix(parametros)
colnames(par_matrix) <- c("beta[0]", "beta[1]", "beta[2]","sigma^2", "rho")

## get results
draws <- as_draws_array(parametros)
summarise_draws(par_matrix)
summarize_draws(par_matrix, "mean", "sd","mcse_mean", "rhat", "ess_bulk", function(x) quantile(x, probs = c(0.025, 0.975)))


##Plots
color_scheme_set(beyonce_palette(101))
labels_params <- c(
  "beta[0]"  = expression(beta[0]),
  "beta[1]"  = expression(beta[1]),
  "beta[2]"  = expression(beta[2]),
  "sigma^2"  = expression(sigma^2),
  "rho"      = expression(rho)
)

p <- mcmc_acf_bar(par_matrix, lags = 50)
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

mcmc_dens(par_matrix)
mcmc_violin(par_matrix)




