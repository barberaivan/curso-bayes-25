## ## TP 4B: MODELO JERÁRQUICO (PARTIAL POOLING) ## ##

# Notas ----
# - Modelo jerárquico con ordenada al origen aleatoria
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Librerías ----
library(tidyverse)

# Modelo jerarquico simple ----
## Definiciones modelo ----
x <- seq(-10, 20, 2)
N_jerarq <- 12
desvio_media <- 2
media_a <- 3
desvio_a <- 3
b <- .4
veces <- 50
N <- length(x)

## Generar datos ----
d_sim <- NULL
for(i_jerarq in 1:N_jerarq){
  a_jerarq_tmp <- rnorm(1, media_a, desvio_a)
  y_med_tmp <- a_jerarq_tmp + b * x
  
  for(i_vez in 1:veces){
    y_obs_tmp <- rnorm(N, y_med_tmp, desvio_media)
    
    d_sim <- rbind(
      d_sim,
      cbind.data.frame(
        i_jerarq = i_jerarq,
        i_vez = i_vez,
        x = x,
        y_media = y_med_tmp,
        y_obs = y_obs_tmp))
  }
}

## Resumir datos simulados ----
rsm_d_sim <- d_sim %>%
  group_by(i_jerarq, x) %>%
  summarise(
    y_media = median(y_media),
    y_obs_mediana = median(y_obs),
    y_obs_q_05 = quantile(y_obs, .05),
    y_obs_q_95 = quantile(y_obs, .95)
  )

## Visualizar datos simulados ----
ggplot(rsm_d_sim, aes(group = i_jerarq)) +
  geom_line(aes(x, y_media), color = "orange", alpha = 1)+
  geom_line(aes(x, y_obs_mediana), color = "blue", alpha = .5) +
  geom_ribbon(aes(x, ymin = y_obs_q_05, ymax = y_obs_q_95),
              fill = "blue", color = NA, alpha = .2)+
  facet_wrap(~i_jerarq)+
  labs(x = "x", y = "y", title = "TP4B: Modelo con ordenada jerárquica")

ggplot(d_sim) +
  geom_line(aes(x, y_media, group = i_jerarq),
            color = "orange", alpha = .5) +
  labs(x = "x", y = "y", title = "TP4B: Modelo con ordenada jerárquica")
