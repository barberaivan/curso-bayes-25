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

datos <- read.csv(file.path("datos", "renison_data_trees_identified.csv"))

# Calculamos crecimiento en 15 años
datos$g <- datos$height2 - datos$height
datos <- subset(datos, !is.na(g)) # nos quedamos con los árboles que sobrevivieron

# Predictoras a considerar:
# height, elev, pforest, manag, fire, plot (random effect)

# Exploramos datos

# Altura inicial
ggplot(datos, aes(height, g)) +
  geom_point(alpha = 0.5) +
  facet_grid(manag ~ fire, axes = "all")

# Altitud
ggplot(datos, aes(elev, g)) +
  geom_point(alpha = 0.5) +
  facet_grid(manag ~ fire, axes = "all")

# Cobertura de bosque
ggplot(datos, aes(pforest, g)) +
  geom_point(alpha = 0.5) +
  facet_grid(manag ~ fire, axes = "all")

# Estandarizamos predictoras continuas
datos$h <- scale(datos$height) |> as.numeric()
datos$a <- scale(datos$elev) |> as.numeric()
datos$f <- scale(datos$pforest) |> as.numeric()
# a de altitud, porque la e de elevation en Stan es el nro e :/

# Guardamos medias y desvíos para desestandarizar o estandarizar luego
mmm <- datos[, c("height", "elev", "pforest")] |> as.matrix()
means <- colMeans(mmm)
sds <- apply(mmm, 2, sd)

# Consideraremos interacción manejo * fuego. Creamos dummies, utilizando estos
# niveles de referencia:
# manejo = Conservation
# fuego = Unburned plots
# (creamos dummies para lo opuesto)
datos$r <- as.numeric(datos$manag == "Ranching")
datos$b <- as.numeric(datos$fire == "Burned plots")
datos$rb <- datos$r * datos$b
# combinación, para modelar interacción.

# parcela (plot), numérico
datos$p <- as.numeric(as.factor(as.character(datos$plot)))

N <- nrow(datos)
np <- max(datos$p)

# Datos para Stan
stan_data <- list(
  N = N,
  K = 7, # cantidad de efectos fijos
  np = np,

  g = datos$g,

  h = datos$h,
  a = datos$a,
  f = datos$f,

  r = datos$r,
  b = datos$b,
  rb = datos$rb,

  p = datos$p
)

# Modelo 01: distribución normal ------------------------------------------

# compilamos
m1 <- cmdstan_model(file.path("modelos", "tabaquillos_01.stan"))

# muestreamos
fit1 <- m1$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Verificamos el HMC
fit1$cmdstan_diagnose()

## Revisamos residuos DHARMa

mu <- fit1$draws("mu") |> as_draws_matrix() |> t()
sigma <- fit1$draws("sigma_g") |> as.vector()
S <- ncol(mu)

ysim <- matrix(NA, N, S)

for (s in 1:S) {
  ysim[, s] <- rnorm(N, mu[, s], sigma[s])
}

res <- createDHARMa(
  simulatedResponse = ysim,
  observedResponse = datos$g
)
plot(res)
# Se ve algo como subdispersión: gran parte de los datos se acumulan más al
# centro de la distribución en relación a lo que el modelo espera.
# Hay algunos outliers.

# Probemos t de student

# Modelo 02: distribución t de student -----------------------------------

# compilamos
m2 <- cmdstan_model(file.path("modelos", "tabaquillos_02.stan"))

# muestreamos
fit2 <- m2$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Verificamos el HMC
fit2$cmdstan_diagnose()

## Revisamos residuos DHARMa

mu <- fit2$draws("mu") |> as_draws_matrix() |> t()
sigma <- fit2$draws("sigma_g") |> as.vector()
nu <- fit2$draws("nu") |> as.vector()

S <- ncol(mu)

ysim <- matrix(NA, N, S)

for (s in 1:S) {
  ysim[, s] <- mu[, s] + sigma[s] * rt(N, df = nu[s]) # rt() está estandarizada
}

res <- createDHARMa(
  simulatedResponse = ysim,
  observedResponse = datos$g
)
plot(res)
# Mejoró un poco, pero sigue estando un poco errado.

# Resultados modelo 2 -----------------------------------------------------

# Resumen de las marginales (lo que le pidamos en "variables")
summ <- fit2$summary(variables = c("coef", "sigma_ranef", "sigma_g", "nu"))
print(summ)

# Extraemos las muestras para calculas más cosas
d <- fit2$draws(
  variables = c("coef", "sigma_ranef", "sigma_g", "nu"),
  format = "draws_df"
)

# Algunos enunciados probabilísticos que pueden interesar (dependiendo del
# contexto)

S <- nrow(d) # cantidad de muestras de la posterior

# Pr(media_ranching > media_conservación)
sum(d$'coef[2]' > 0) / S
# coef[2] es la diferencia ranching - conservación, porque conservación es el
# nivel de referencia

# Pr(media_quemado > media_no_quemad)
sum(d$'coef[3]' > 0) / S
# coef[3] es la diferencia quemado - no quemado, porque no quemado es el
# nivel de referencia.
# Esta prob da muy baja, casi cero. Eso indica que el enunciado opuesto
# tiene probabilidad casi 1: el crecimiento es mayor en parcelas no quemadas.

# ¿Cuán probable es que el efecto del manejo varíe según quemado/no quemado, o,
# equivalentemente, que el efecto del fuego varíe esgún el manejo?
# Esto se representa con el término de interacción.
sum(d$`coef[4]` > 0) / S
# Si la prob es cercana a 0.5 (lo es), significa que tanto un efecto positivo
# un efecto negativo son probables, lo cual indica que no hay mucha evidencia
# a favor de una interacción manejo * fuego.

# Probabilidad de que el crecimiento aumente en función de la altura inicial, la
# altitud y la cobertura de bosque:

sum(d$`coef[5]` > 0) / S
# si esto es cero, significa que la altura inicial disminuye el crecimiento, con
# probabilidad 1.

sum(d$`coef[6]` > 0) / S
# El efecto de la altitud es mayormente negativo

sum(d$`coef[7]` > 0) / S
# A mayor cobertura de bosque, el crecimiento es mayor, con probabilidad cercana
# a 1 (mucha certeza sobre el signo.)

# Comparación de previa contra posterior.
# Usamos density() para obtener la densidad empírica a partir de muestras,
# y curve() para dibujar la previa.

# intercept (coef[1])
plot(density(d$`coef[1]`), main = NA, xlab = "coef[1]: intercept",
     xlim = c(-100, 100))
curve(dnorm(x, 0, 1000), col = 4, add = T)

# sigma_g (residual) y sigma_ranef (variación entre plots)
plot(density(d$sigma_ranef), main = NA, xlab = "sigma_ranef", xlim = c(0, 200),
     ylim = c(0, 0.13))
lines(density(d$sigma_g), col = "red", ylim = c(0, 0.1))
curve(dnorm(x, 0, 1000) * 2, col = 4, add = T)
text(100, 0.06, "sigma_ranef")
text(25, 0.1, "sigma_g", col = "red")
text(150, 0.015, "previa", col = 4)

# nu, los grados de libertad
plot(density(d$nu), main = NA, xlab = expression(nu), xlim = c(0, 20))
curve(dgamma(x, 2, 0.1), col = 4, add = T)


# Predicciones (con modelo 2) ---------------------------------------------

# Graficaremos predicciones parciales para el crecimiento en función de la
# altura inicial del árbol, para todas las combinaciones de manejo y fuego.
# Las demás predictoras se fijarán en sus medias (0).

nrep <- 200
pred <- expand.grid( # crea todas las combinaciones
  h = seq(min(datos$h), max(datos$h), length.out = nrep),
  r = c(0, 1), # conservación y ganado
  b = c(0, 1)  # quemado y no quemado
  # ignoramos las demás predictoras porque asumimos que están en sus medias
  # (= 0)
)

# término de interacción
pred$rb <- pred$r * pred$b

# desestandarizamos para graficar luego
pred$height <- pred$h * sds["height"] + means["height"]

# recuperamos factores a partir de las dummy
pred$manag <- "Conservation"
pred$manag[pred$r == 1] <- "Ranching"
pred$fire <- "Unburned plots"
pred$fire[pred$b == 1] <- "Burned plots"

# Calculamos mu y simulamos datos para luego obtener los intervalos de predicción
npred <- nrow(pred)
mumat <- ymat <- matrix(NA, npred, S)

# Extraemos sigma_ranef para simular efectos aleatorios de los plots
sigma_ranef <- fit2$draws("sigma_ranef") |> as.vector()

# El sigma total es la raíz de la suma de cuadrados de los dos desvíos
sigma_total <- sqrt(sigma ^ 2 + sigma_ranef ^ 2)
# Así el intervalo de predicción incorporará variabilidad entre plots

# Extraemos los coeficientes de regresión
coefs <- fit2$draws("coef") |> as_draws_matrix() |> t()

for (s in 1:S) {
  mumat[, s] <-
    # Intercept:
    coefs[1, s] +
    # Efecto de pred. categóricas:
    coefs[2, s] * pred$r + coefs[3, s] * pred$b + coefs[4, s] * pred$rb +
    # Efecto de altura inicial
    coefs[5, s] * pred$h

  ymat[, s] <- rt(npred, nu[s]) * sigma_total[s] + mumat[, s]
}

# Resumimos
pred_full <- cbind(
  pred,
  apply(mumat, 1, mean_ci) |> t() |> as.data.frame(),
  apply(ymat, 1, ci, name = "y") |> t() |> as.data.frame()
)

# Graficamos
ggplot(pred_full, aes(height, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  # Intervalo de predicción
  geom_ribbon(
    aes(height, ymin = y_lwr, ymax = y_upr), inherit.aes = F,
    fill = "blue", alpha = 0.3, color = NA
  ) +
  # Intervalo de credibilidad para la media
  geom_ribbon(
    fill = "orange", alpha = 0.6, color = NA
  ) +
  # Media de la media
  geom_line() +
  # Datos
  geom_point(
    aes(height, g), datos, size = 1.5, alpha = 0.5, inherit.aes = F
  ) +
  # cosmética
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
  facet_grid(manag ~ fire, axes = "all") +
  labs(x = "Altura inicial (cm)", y = "Crecimiento (cm)") +
  nice_theme()

# Lo mismo se puede hacer con las otras dos predictoras continuas

# R2 ----------------------------------------------------------------------

# El problema con la t es que si nu <= 2,  la varianza no está definida,
# con lo cual, no podemos calcular el R2. Pero por suerte, en una proporción
# muy baja de las muetras pasa eso, así que ignoraremos el problema.
sum(nu <= 2)

# Calcularemos el R2 incondicional al efecto aleatorio del plot, es decir, sin
# contar su efecto en la variabilidad explicada (porque el plot no explica nada).
# Esto implica que la media no es mu, sino sólo la parte de efectos fijos, que
# en Stan llamamos fixef
fixef <- fit2$draws("fixef") |> as_draws_matrix() |> t()

# Vector para almacenar el R2
R2 <- numeric(S)

# Loop sobre muestras de la posterior
for (s in 1:S) {
  # varianza de los predichos = varianza explicada
  var_fit <- var(fixef[, s]) # calculada en Stan

  # Varianza residual para cada observación, usando sigma_total y como se define
  # para la t de student (pero es constante entre observaciones!)
  var_res <- ifelse(
    nu[s] > 2, # si nu > 2
    sigma_total[s] ^ 2 * nu[s] / (nu[s] - 2), # calculá la varianza
    NA # si no, ponele NA
  )

  if (is.na(var_res)) {
    R2[s] <- NA
  } else {
    R2[s] <- var_fit / (var_fit + var_res)
  }
}

# Ahora, considerando el efecto del plot como varianza explicada. Usamos
# sigma, no sigma_total, y la media es mu, que incluye tanto fixef como ranef.

R2_cond <- numeric(S) # cond por condicional a los ranef

# Loop sobre muestras de la posterior
for (s in 1:S) {
  # varianza de los predichos = varianza explicada
  var_fit <- var(mu[, s]) # calculada en Stan

  # Varianza residual para cada observación, usando sigma y como se define
  # para la t de student (pero es constante entre observaciones!)
  var_res <- ifelse(
    nu[s] > 2, # si nu > 2
    sigma[s] ^ 2 * nu[s] / (nu[s] - 2), # calculá la varianza
    NA # si no, ponele NA
  )

  if (is.na(var_res)) {
    R2_cond[s] <- NA
  } else {
    R2_cond[s] <- var_fit / (var_fit + var_res)
  }
}

par(mfrow = c(1, 2))
plot(density(na.omit(R2), from = 0, to = 1), main = NA,
     xlab = "R2 incondicional")
plot(density(na.omit(R2_cond), from = 0, to = 1), main = NA,
     xlab = "R2 condicional")
par(mfrow = c(1, 1))