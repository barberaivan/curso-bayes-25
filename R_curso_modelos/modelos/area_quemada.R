rm(list = ls()) # Limpiamos el entorno

# Paquetes ----------------------------------------------------------------

library(tidyverse)  # gráficos et al
library(cmdstanr)   # ajustar modelos con Stan
library(posterior)  # resumir posteriores y diagnósticos de MCMC
library(bayesplot)  # visualizar posteriores
library(DHARMa)     # probabilidad acumulada
library(patchwork)  # unir gráficos

# Funciones ---------------------------------------------------------------

# ggplot custom theme
nice_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major =  element_line(),

    axis.line = element_line(linewidth = 0.35),

    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),

    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

theme_set(theme_classic())

# Funciones para resumir posteriores
ci <- function(x, width = 0.95, name = "mu") {
  # definimos probabilidad acumulada según la amplitud del intervalo (ci)
  low <- (1 - width) / 2
  high <- width + low

  # calculamos cuantiles asociados
  qq <- quantile(x, probs = c(low, high), method = 8)

  # damos nombres y ordenamos
  nn <- paste(name, c("_lwr", "_upr"), sep = "")
  names(qq) <- nn

  return(qq)
}

mean_ci <- function(x, width = 0.95, name = "mu") {
  # media
  mm <- mean(x)
  names(mm) <- paste(name, "_mean", sep = "")

  # cuantiles
  qq <- ci(x, width, name)

  return(c(mm, qq))
}

# Datos -------------------------------------------------------------------

datos <- read.csv(file.path("datos", "barbera_data_fire_vegetation.csv"))

# fracción quemada de lo quemable
datos$prop <- datos$area_ha / datos$area_ha_available

# Como hay muchos ceros, modelaremos con una beta cero inflada:
# probabilidad de que haya fuego y, dado que hay, cuánto se quema.
datos$fire <- as.numeric(datos$area_ha > 0)

# casteamos como factor la vegetación y como numérica
datos$vegetation_class <- factor(datos$vegetation_class)
datos$veg_num <- as.numeric(datos$vegetation_class)

# Importamos clima para toda el área de estudio
dclim <- read.csv(file.path("datos", "barbera_data_fire_total_climate.csv"))

# agregamos las variables climáticas
datos_full <- left_join(datos, dclim, by = "year")

# Exploración de efectos del clima climáticos -----------------------------

# Para elegir una predictora
clim_names <- c("pp", "temp", "vpd", "fwi")
dlong <- pivot_longer(
  datos_full, all_of(which(names(datos_full) %in% clim_names)),
  names_to = "variable", values_to = "clim_value"
)

ggplot(dlong, aes(clim_value, prop, color = vegetation_class)) +
  geom_smooth(method = "gam", se = F) +
  geom_point() +
  facet_grid(
    rows = vars(vegetation_class), cols = vars(variable),
    scales = "free", margins = "all"
  )

# Utilizamos el fwi. Lo agregamos a datos
datos$fwi <- datos_full$fwi

# Lo estandarizamos. Antes calculamos media y desvío para poder desestandarizar
# luego.
fwi_mean <- mean(datos$fwi)
fwi_sd <- sd(datos$fwi)
datos$fwiz <- (datos$fwi - fwi_mean) / fwi_sd



# Definición de modelos ---------------------------------------------------

# Dos modelos: GLM beta cero inflado con efectos del fwi por tipo de vegetación,
# considerando a la vegetación como efecto fijo (no pooling) o aleatorio
# (partial pooling)


# Datos para Stan ---------------------------------------------------------

# Creamos una lista nombrada para pasarle los datos. Los nombres de cada
# elemento de la lista tienen que ser exactamente los que definimos en
# la sección data {}

N <- nrow(datos)
K <- max(datos$veg_num)
B <- sum(datos$fire)
rows_fire <- which(datos$fire == 1)

stan_data <- list(
  N = N,
  K = K,
  B = B,

  prop = datos$prop,
  fire = datos$fire,
  fwiz = datos$fwiz,
  veg = datos$veg_num,

  rows_fire = rows_fire
)

# Modelo 1: no pooling ----------------------------------------------------

# MCMC --------------------------------------------------------------------

# Compilamos
model1 <- cmdstan_model(file.path("modelos", "area_quemada_01.stan"))

# Y muestreamos la posterior
fit1 <- model1$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Verificamos el HMC
fit1$cmdstan_diagnose()

# Verificación predictiva posterior  --------------------------------------

# Extraemos las muestras de p1, p2 y phi_vec, aprovechando que las calculamos en
# Stan
p1 <- fit1$draws("p1") |> as_draws_matrix() |> t()
p2 <- fit1$draws("p2") |> as_draws_matrix() |> t()
phi <- fit1$draws("phi_vec") |> as_draws_matrix() |> t()

# Subseteamos p2 y phi, tomando sólo las filas con presencia
p2 <- p2[rows_fire, ]
phi <- phi[rows_fire, ]

# Parámetros para la distribución beta
aa <- p2 * phi
bb <- (1 - p2) * phi

S <- ncol(p1)

# Simulamos de la predictiva posterior para ambos procesos por separado
# (presencia y proporción)
ysim_pres <- matrix(NA, N, S)
ysim_prop <- matrix(NA, B, S) # menos obs

for (s in 1:S) {
  ysim_pres[, s] <- rbinom(N, size = 1, prob = p1[, s])
  ysim_prop[, s] <- rbeta(B, aa[, s], bb[, s])
}

# DHARMa
res_pres <- createDHARMa(
  simulatedResponse = ysim_pres,
  observedResponse = datos$fire,
  integerResponse = TRUE # se da cuenta solo, pero por las dudas
)
plot(res_pres)
plotResiduals(res_pres, form = datos$fwi)
plotResiduals(res_pres, form = datos$vegetation_class)

res_prop <- createDHARMa(
  simulatedResponse = ysim_prop,
  observedResponse = datos$prop[rows_fire],
  integerResponse = FALSE # se da cuenta solo, pero por las dudas
)
plot(res_prop, rank = F)
plotResiduals(res_prop, form = datos$fwi[rows_fire])
plotResiduals(res_prop, form = datos$vegetation_class[rows_fire])

# Trasponemos la matriz de datos simulados porque a bayesplot le gustan las
# iteraciones de la posterior como filas (las teníamos en las columnas)
ysim_prop_t <- t(ysim_prop)
datos2 <- datos[rows_fire, ]

# posterior predictive intervals en función de algo
for (v in 1:K) {
  rr <- datos2$veg_num == v

  pp <-
  ppc_intervals(
    y = datos2$prop[rr],
    yrep = ysim_prop_t[, rr],
    x = datos2$fwi[rr],
    prob = 0.9
  ) +
  labs(x = "FWI", y = "Proporción quemada",
       title = levels(datos$vegetation_class)[v])

  print(pp)
}

# Predicciones ------------------------------------------------------------

# Graficaremos la probabilidad de que se queme, la proporción esperada dado que
# se quema, y la proporción esperada marginal a que se queme o no. Esta última
# es el producto de las anteriores.

nrep <- 200 # largo de la secuencia de FWI
fwi_seq <- seq(min(datos$fwi), max(datos$fwi), length.out = nrep)

# Datos de predicción incluyendo tipo de vegetación
pred <- expand.grid( # crea todas las combinaciones de lo que le demos
  fwi = fwi_seq,
  veg_num = 1:K
)
pred$fwiz <- (pred$fwi - fwi_mean) / fwi_sd
npred <- nrow(pred)

# Matrices para guardar los predichos
pres_mat <- prop_mat <- matrix(NA, npred, S)

# Extraemos muestras de Stan
a1 <- fit1$draws("a1") |> as_draws_matrix() |> t()
b1 <- fit1$draws("b1") |> as_draws_matrix() |> t()
a2 <- fit1$draws("a2") |> as_draws_matrix() |> t()
b2 <- fit1$draws("b2") |> as_draws_matrix() |> t()

for (s in 1:S) {   # loop sobre muestras
  for (v in 1:K) { # loop sobre tipos de vegetación
    # filas correspondientes al tipo de vegetación v
    rr <- pred$veg_num == v
    pres_mat[rr, s] <- plogis(a1[v, s] + b1[v, s] * pred$fwiz[rr])
    prop_mat[rr, s] <- plogis(a2[v, s] + b2[v, s] * pred$fwiz[rr])
  }
}

# El producto
mu_mat <- pres_mat * prop_mat

# Resumimos todas y las unimos (formato long)
summs0 <- rbind(
  apply(mu_mat, 1, mean_ci) |> t() |> as.data.frame(),
  apply(prop_mat, 1, mean_ci) |> t() |> as.data.frame(),
  apply(pres_mat, 1, mean_ci) |> t() |> as.data.frame()
)

# Agregamos nombre de variable
vnames <- c("Proporción\nincondicional", "Proporción | quema", "Prob. quema")
summs0$variable <- factor(rep(vnames, each = npred), levels = vnames)

# Agregamos predictoras
summs <- cbind(summs0, rbind(pred, pred, pred))
summs$vegetation_class <- levels(datos$vegetation_class)[summs$veg_num]

# Graficamos
pp1 <-
ggplot(summs[summs$variable == "Prob. quema", ],
       aes(fwi, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  geom_ribbon(alpha = 0.6, color = NA, fill = "orange") +
  geom_line() +
  labs(x = NULL, y = NULL, title = "Prob. quema") +
  facet_wrap(vars(vegetation_class), ncol = 1, strip.position = "right") +
  nice_theme()

pp2 <-
  ggplot(summs[summs$variable == "Proporción | quema", ],
         aes(fwi, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  geom_ribbon(alpha = 0.6, color = NA, fill = "orange") +
  geom_line() +
  labs(x = "FWI", y = NULL, title = "Proporción | quema") +
  facet_wrap(vars(vegetation_class), ncol = 1, strip.position = "right") +
  nice_theme() +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank()
  )

pp3 <-
  ggplot(summs[summs$variable == "Proporción\nincondicional", ],
         aes(fwi, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  geom_ribbon(alpha = 0.6, color = NA, fill = "orange") +
  geom_line() +
  labs(x = NULL, y = "Predicción", title = "Proporción incondicional") +
  facet_wrap(vars(vegetation_class), ncol = 1, strip.position = "right") +
  nice_theme() +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank()
  )

# merge
(pp3 | pp2 | pp1)

# R2 ----------------------------------------------------------------------

# En este caso, el R2 del modelo es un poco entreverado porque la distribución
# es una mezcla entre una bernoulli y una beta. Definir la media es sencillo,
# pero la varianza es más compleja. Preguntándole a ChatGPT, obtuve lo siguiente:

# Vector para almacenar el R2
R2 <- numeric(S)

# Loop sobre muestras de la posterior
for (s in 1:S) {
  # La media es el producto
  mu <- as.vector(p1[, s] * p2[, s])
  # varianza de los predichos = varianza explicada
  var_fit <- var(mu)

  # Varianza residual para cada observación.

  # Varianza de la beta
  beta_var <- p2[, s] * (1 - p2[, s]) / (1 + phi[, s])

  # Compute hurdle beta variance
  var_res <- p1[, s] * beta_var + p1[, s] * (1 - p1[, s]) * p2[, s] ^ 2

  # promediamos
  var_res_mean <- mean(var_res)

  # R2
  R2[s] <- var_fit / (var_fit + var_res_mean)
}

plot(density(R2, from = 0, to = 1), main = NA, xlab = "R2")


# Modelo 2: partial pooling -----------------------------------------------

# Un modelo jerárquico permite que los grupos con pocas observaciones o más
# ruido tengan estimaciones menos extremas. En este caso, podría ayudar a
# obtener estimaciones más conservativas para los bosques secos.

# MCMC --------------------------------------------------------------------

# Compilamos
model2 <- cmdstan_model(file.path("modelos", "area_quemada_02.stan"))

# Y muestreamos la posterior
fit2 <- model1$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Verificamos el HMC
fit2$cmdstan_diagnose()

# El resto del código sigue igual, porque la vegetación se trata, en la práctica,
# como un efecto fijo. Sólo usamos un modelo jerárquico para estimar mejor, no
# para obtener predicciones marginales al tipo de vegetación.