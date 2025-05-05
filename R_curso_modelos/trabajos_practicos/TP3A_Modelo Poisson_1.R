## ## TP 3A: MODELOS NO LINEALES - POISSON ## ##

# Notas ----
# - Especifica modelo exponencial Poisson con incertidumbre en parámetros
# - El modelo es: y ~ Poisson(e^(a + b * x))
# - Simula datos muchas veces para representar incertidumbre
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Paquetes ----
library(tidyverse)

# Definiciones ----
x <- seq(-10, 20, 2)
media_a <- .2
desvio_a <- .02
media_b <- .1
desvio_b <- .02
N <- length(x)
veces <- 100

# Simular datos ----
d_sim <- NULL
for (i_vez in 1:veces) {
  ## Calcular valores
  a_tmp <- rnorm(1, media_a, desvio_a)
  b_tmp <- rnorm(1, media_b, desvio_b)
  y_media_tmp <- exp(a_tmp + b_tmp * x)
  y_obs_tmp <- rpois(N, y_media_tmp)
  
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

# Visualizar datos simulados ----
## Puntos simulados
p1 <- ggplot(rsm_d_sim) +
  geom_point(data = d_sim, aes(x, y_obs), color = "blue", alpha = .1) +
  labs(x = "x", y = "y", title = "TP3A: Modelo de Poisson")
p1


## Líneas individuales ('spaghetti plot') (para subset de simulaciones)
p2 <- ggplot(d_sim[d_sim$i_vez %in% 1:20,], aes(group = i_vez)) +
  geom_line(aes(x, y_obs), color = "blue", alpha = .1) +
  geom_line(aes(x, y_media), color = "orange", alpha = .4) +
  labs(x = "x", y = "y", title = "TP3A: Modelo de Poisson")
p2

## Intervalos de incertidumbre
p3 <- p2 +
  geom_line(aes(x, y_media_mediana), rsm_d_sim, color = "orange")+
  geom_line(aes(x, y_obs_mediana), rsm_d_sim, color = "blue")+
  geom_ribbon(aes(x, ymin = y_media_q_05, ymax = y_media_q_95), rsm_d_sim,
              fill = "orange", color = NA, alpha = .4)+
  geom_ribbon(aes(x, ymin = y_obs_q_05, ymax = y_obs_q_95), rsm_d_sim,
              fill = "blue", color = NA, alpha = .2)
p3
