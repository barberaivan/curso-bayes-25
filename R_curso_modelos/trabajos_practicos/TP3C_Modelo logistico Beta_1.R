## ## TP 3C: MODELOS NO LINEALES - LOGISTICO CON FRACCIÓN ## ##

# Notas ----
# - Especifica modelo logístico para fracción con incertidumbre en parámetros
# - El modelo es: y ~ Beta(mu * kappa, (1 - mu) * kappa);
#     mu = 1 / (1 + e^(-(a + b * x))))
# - Simula datos muchas veces para representar incertidumbre
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Paquetes ----
library(tidyverse)

# Definiciones ----
x <- seq(-10, 20, 2)
media_a <- -2
desvio_a <-.2
media_b <- .3
desvio_b <- .01
kappa <- 10
N <- length(x)
veces <- 100

# Generar datos ----
d_sim <- NULL
for(i_vez in 1:veces){
  ## Calcular valores
  a_tmp <- rnorm(1, media_a, desvio_a)
  b_tmp <- rnorm(1, media_b, desvio_b)
  y_media_tmp <- 1 / (1 + exp(-(a_tmp + b_tmp * x)))
  
  alfa_tmp <- y_media_tmp * rgamma(N, kappa)
  beta_tmp <- (1 - y_media_tmp) * rgamma(N, kappa)
  y_obs_tmp <- rbeta(N, alfa_tmp, beta_tmp)
  
  ## Guardar valores simulados
  d_sim <- rbind(
    d_sim,
    cbind.data.frame(
      i_vez = i_vez,
      x = x,
      y_media = y_media_tmp,
      y_obs = y_obs_tmp))
}

# Resumir datos simulados ----
rsm_d_sim <- d_sim %>%
  group_by(x) %>%
  summarise(
    y_media_mediana = median(y_media),
    y_media_q_05 = quantile(y_media, .05),
    y_media_q_95 = quantile(y_media, .95),
    y_obs_mediana = median(y_obs),
    y_obs_q_05 = quantile(y_obs, .05),
    y_obs_q_95 = quantile(y_obs, .95)
  )

# Visualizar datos ----
p1 <- ggplot(rsm_d_sim) +
  geom_line(aes(x, y_media, group = i_vez), d_sim, color = "orange", alpha = .3) +
  geom_point(aes(x, y_obs), d_sim, color = "blue", alpha = .3) +
  labs(x = "x", y = "y", title = "Modelo lineal con error")
p1

p2 <- p1 +
  geom_line(aes(x, y_media_mediana), color = "orange")+
  geom_line(aes(x, y_obs_mediana), color = "blue")+
  geom_ribbon(aes(x, ymin = y_media_q_05, ymax = y_media_q_95),
              fill = "orange", color = NA, alpha = .4)+
  geom_ribbon(aes(x, ymin = y_obs_q_05, ymax = y_obs_q_95),
              fill = "blue", color = NA, alpha = .2)+
  labs(x = "x", y = "y", title = "Modelo lineal con error")
p2
