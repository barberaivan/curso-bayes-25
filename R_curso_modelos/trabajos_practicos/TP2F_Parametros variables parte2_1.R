## ## TP 2F: MODELO LINEAL CON PARÁMETROS VARIABLES - PARTE 2 ## ##

# Notas ----
# - Especifica modelo lineal con incertidumbre en variable dependiente
#   y en parámetros
# - Simula datos muchas veces para representar incertidumbre
# - Repite en un loop 'for' la simulación usada en TP2E
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Paquetes ----
library(tidyverse)

# Definiciones ----
path.datos <- file.path("datos", "barbera_data_fire_total_climate.csv")

x <- seq(-5, 25, 2)
media_a <- 2    #media de 'a'
desvio_a <- .5  #desvío estándar de 'a'
media_b <- .7   #media de 'b'
desvio_b <- .2  #desvío estándar de 'b'
desvio_media <- 1.2 #desvío estándar de la media del modelo
N <- length(x)
veces <- 100 #cantidad de veces que se repite la simulación

# Abrir datos de fuego ----
d <- read.csv(path.datos)

# Simular datos ----
d_sim <- NULL
for (i_vez in 1:veces) { #para cada vez
  ## Calcular valores
  a_tmp <- rnorm(1, media_a, desvio_a) #genera muestra de valores de 'a' (es un vector)
  b_tmp <- rnorm(1, media_b, desvio_b) #genera muestra de valores de 'b' (es un vector)
  y_media_tmp <- a_tmp + b_tmp * x     #calcula media del modelo (es un vector)
  y_obs_tmp <- rnorm(N, y_media_tmp, desvio_media) #genera muestras de observaciones (es un vector)
  
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
rsm_d_sim <- d_sim %>%            #resume 'd_sim'
  group_by(x) %>%                 #para cada valor de 'x'
  summarise(                      #calcula...
    y_media_mediana = median(y_media),
    y_media_q_05 = quantile(y_media, .05),
    y_media_q_95 = quantile(y_media, .95),
    y_obs_mediana = median(y_obs),
    y_obs_q_05 = quantile(y_obs, .05),
    y_obs_q_95 = quantile(y_obs, .95)
  )

# Visualizar datos simulados ----
## Líneas individuales ('spaghetti plot')
ggplot(d_sim, aes(group = i_vez)) +
  geom_line(aes(x, y_obs), color = "blue", alpha = .15) +
  geom_line(aes(x, y_media), color = "orange", alpha = .25) +
  labs(x = "x", y = "y", title = "TP2F: Modelo lineal con parámetros variables")

## Intervalos intercuantiles
ggplot(rsm_d_sim) +
  geom_ribbon(aes(x, ymin = y_media_q_05, ymax = y_media_q_95),
              fill = "orange", color = NA, alpha = .4) +
  geom_ribbon(aes(x, ymin = y_obs_q_05, ymax = y_obs_q_95),
              fill = "blue", color = NA, alpha = .2) +
  # geom_point(data = d, aes(temp, log10(area_ha))) +
  labs(x = "x", y = "y", title = "TP2F: Modelo lineal con parámetros variables")
