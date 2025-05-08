## ## TP 3B: MODELOS NO LINEALES - LOGÍSTICO BINARIO ## ##

# Notas ----
# - Especifica modelo logístico con incertidumbre en parámetros
# - El modelo es: y ~ Bernoulli(1 / (1 + e^(-(a + b * x))))
# - Simula datos muchas veces para representar incertidumbre
# - Comparar simulación con datos de tabaquillo (Renison)
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Paquetes ----
library(tidyverse)

# Definiciones ----
path.datos <- file.path("datos", "renison_data_trees_identified.csv")

x <- seq(-100, 1200, 100)
b <- .01
a <- -6
N <- length(x)
veces <- 10

# Simular datos ----
y_media <- 1 / (1 + exp(-(a + b * x)))
d_sim <- NULL
for (i_vez in 1:veces) {
  ## Calcular valores
  y_obs_tmp <- rbinom(N * veces, 1, y_media)
  
  ## Guardar valores simulados
  d_sim <- rbind(
    d_sim,
    cbind.data.frame(
      i_vez = i_vez,
      x = x,
      y_media = y_media,
      y_obs = y_obs_tmp))
}

# Resumir datos simulados ----
rsm_d_sim <- d_sim %>%
  group_by(x) %>%
  summarise(
    y_media_mediana = median(y_media),
    y_obs_media = mean(y_obs),
    y_obs_q_05 = quantile(y_obs, .05),
    y_obs_q_95 = quantile(y_obs, .95)
  )

# Visualizar datos simulados ----
## Valores simulados
p1 <- ggplot(rsm_d_sim) +
  geom_line(aes(x, y_media, group = i_vez), d_sim, color = "orange", alpha = .4) +
  geom_jitter(aes(x, y_obs), d_sim,
              width = 5, height = 0, color = "blue", alpha = .3) +
  labs(x = "x", y = "y", title = "TP3B: Modelo logístico respuesta binaria")
p1

# Abrir datos de tabaquillos ----
d <- read.csv(path.datos)

## Graficar supervivencia en función de la altura de cada árbol
p3 <- ggplot(d) +
  geom_point(aes(height, surv), alpha = .1) +
  labs(x = "Altura inicial (cm)", y = "Supervivencia", size = "N obs",
       title = "TP3B: Supervivencia en función de altura")
p3

## Calcular fracción que sobrevive según categoría de altura de cada árbol
seq_breaks <- seq(25, max(d$height) + 50, 50)
rsm_d_height <- d %>%
  mutate(height_ctg = cut(height, breaks = seq_breaks, right = F)) %>%
  group_by(height_ctg) %>%
  summarise(cant = n(),
            n_surv = sum(surv)) %>%
  mutate(fr_surv = n_surv/cant)
levels(rsm_d_height$height_ctg) <- seq_breaks + 25
rsm_d_height$height_ctg2 <- as.numeric(as.character(rsm_d_height$height_ctg))

## Graficar fracción que sobrevive según categoría de altura
p4 <- p3 +
  geom_point(data = rsm_d_height, aes(height_ctg2, fr_surv, size = cant),
             color = "slategray4", alpha = .8)
p4

## Graficar con predicción de modelo logístico
p4 +
  geom_line(data = rsm_d_sim, aes(x, y_media_mediana), color = "orange")
