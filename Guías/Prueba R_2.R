## ## VERIFICACION DE FUNCIONAMIENTO DE R PARA CURSO ## ##

# Notas ####
# - Script para verificar que esté instalado lo necesario para el curso

rm(list=ls())

suppressPackageStartupMessages(library(cmdstanr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggridges))

suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(bayesplot))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(DHARMa))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(loo))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))

file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  parallel_chains = 4,
  refresh = 500,
  show_messages = F
)

rsm_fit <- fit$summary()

draws_df <- as.data.frame(fit$draws(format = "df"))
draws_df_lng <- draws_df[, c("lp__", "theta")] %>%
  pivot_longer(c(lp__, theta))

p <- ggplot(draws_df_lng) +
  geom_density(aes(value, group = name)) +
  facet_wrap(~name, scales = "free_x") +
  labs(title = paste0("Verficación = ", round(sum(rsm_fit$mean),2)))
print(p)
