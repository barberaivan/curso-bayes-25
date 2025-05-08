## ## TP 4B: MODELO CON DETECCIÓN IMPERFECTA ## ##

# Notas ----
# - Especifica un modelo con detección imperfecta de la variable de respuesta
# - Simula datos muchas veces para representar incertidumbre
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Paquetes ----
library(tidyverse)

# Definiciones modelo ----
x <- seq(-10, 20, 2)
N <- length(x)
veces <- 100
a <- .2
b <- .1
p <- .6

## Generar datos ----
d_sim <- NULL
for(i_vez in 1:veces){
  ## Calcular valores
  lambda_tmp <- exp(a + b * x)
  z_tmp <- rpois(N, lambda_tmp)
  y_tmp <- rbinom(N, z_tmp, p)
  
  ## Guardar valores simulados
  d_sim <- rbind(
    d_sim,
    cbind.data.frame(
      i_vez = i_vez,
      x = x,
      lambda = lambda_tmp,
      z = z_tmp,
      y = y_tmp))
}

## Resumir datos simulados ----
rsm_d_sim <- d_sim %>%
  group_by(x) %>%
  summarise(
    lambda_mediana = median(lambda),
    z_mediana = median(z),
    z_q_05 = quantile(z, .05),
    z_q_95 = quantile(z, .95),
    y_mediana = median(y),
    y_q_05 = quantile(y, .05),
    y_q_95 = quantile(y, .95)
  )

## Visualizar datos simulados ----
ggplot(rsm_d_sim) +
  geom_line(aes(x, lambda_mediana), color = "orange") +
  geom_line(aes(x, z_mediana), color = "green") +
  geom_line(aes(x, y_mediana), color = "blue") +
  geom_ribbon(aes(x, ymin = z_q_05, ymax = z_q_95),
              fill = "green", color = NA, alpha = .1)+
  geom_ribbon(aes(x, ymin = y_q_05, ymax = y_q_95),
              fill = "blue", color = NA, alpha = .1)+
  labs(x = "x", y = "y", title = "TP4C: Modelo con detección imperfecta")
