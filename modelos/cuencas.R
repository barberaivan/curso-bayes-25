# Paquetes ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(viridis)
library(truncnorm)

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

# Cargamos y preparamos datos ---------------------------------------------

# Registros de caudal por cuenca y fecha
dcaudal <- read.csv(file.path("datos", "cingolani_data_caudales_fechas.csv"))

# Variables predictoras por cuenca
dvar <- read.csv(file.path("datos", "cingolani_data_caudales_variables.csv"))

# Extraemos las columnas con los registros por fecha (fechas en columnas)
date_names_0 <- names(dcaudal)[grep("F", names(dcaudal))]
date_names <- date_names_0[grep("mm", date_names_0, invert = T)]
dsub <- dcaudal[, c("cca", date_names)]

# número de semanas (t), cuencas (c) y observaciones (o)
nt <- length(date_names)
nc <- nrow(dsub)
nobs <- nt * nc

# Alargamos
dlong <- pivot_longer(
  dsub, all_of(date_names), names_to = "fecha_char", values_to = "caudal"
)

# Extraemos el día como numérico, arrancando en cero
dlong$dia <- sapply(1:nrow(dlong), function(i) {
  strsplit(dlong$fecha_char[i], split = "F")[[1]][2] |> as.numeric()
})

# Calculamos semana, para facilitar interpretación
dlong$semana <- dlong$dia / 7

# Cuenca como factor, servirá para graficar
dlong$cca_fac <- factor(as.character(dlong$cca),
                        levels = as.character(unique(dsub$cca)))

# Extraemos rocosidad total
dlong$roca <- rep(dvar$RT, each = nt)

# Gráfico exploratorio ----------------------------------------------------

ggplot(dlong, aes(semana, caudal, color = roca, group = cca_fac)) +
  # geom_line() +
  geom_smooth(se = F) +
  geom_point() +
  scale_color_viridis(option = "A") +
  nice_theme()



# Exploramos función de decaimiento de Weibul -----------------------------

a <- 1   # máximo
l <- 3   # cuánto tarda en bajar
k <- 0.5   # forma (de exponencial a sigmoidea)
curve(a * exp(-(x / l) ^ k), to = 10, ylim = c(0, 1))

# Simulamos k
cc <- rgb(0, 0, 0, 0.1)

a <- 1      # máximo
l <- 2      # cuánto tarda en bajar
k <- 0.05   # forma (de exponencial a sigmoidea)
curve(a * exp(-(x / l) ^ k), to = 10, col = cc, ylim = c(0, a))
for (i in 1:100) {
  k <- rtruncnorm(1, a = 0, mean = 1.5, sd = 3)
  curve(a * exp(-(x / l) ^ k), add = T, col = cc)
}

# Simulamos k y l
a <- 1      # máximo
l <- 2      # cuánto tarda en bajar
k <- 0.05   # forma (de exponencial a sigmoidea)
curve(a * exp(-(x / l) ^ k), to = 10, col = cc, ylim = c(0, a))
for (i in 1:100) {
  k <- rtruncnorm(1, a = 0.5, mean = 1.5, sd = 2)
  l <- rtruncnorm(1, a = 0, mean = 2.5, sd = 2)
  curve(a * exp(-(x / l) ^ k), add = T, col = cc)
}

# a es el valor en la primera fecha.
mean(dsub$F0)
range(dsub$F0)

# a puede estar centrada en
mu_a <- log(mean(dsub$F0))
sd_a <- 2
curve(dlnorm(x, mu_a, sd_a), to = 80)

# k debería estar medio lejos de cero, con 10 siendo un valor muy extremo
mu_k <- log(2)
sd_k <- 0.7
curve(dlnorm(x, mu_k, sd_k), to = 10)

# l debería estar medio lejos de cero, con 6 siendo un valor alto
mu_l <- log(2.5)
sd_l <- 1.5
curve(dlnorm(x, mu_l, sd_l), to = 10)

# A las medias les pongo estas previas,
# y los log-sds, les pongo normal truncada con sd_param * 2


# Previa para phi
mu <- 15
phi <- 5.5
curve(dgamma(x, shape = 1 / phi, rate = 1 / (mu * phi)), to = 80, col = cc)
for (i in 1:100) {
  log_phi <- rnorm(1, log(0.3), 1.5)
  phi <- exp(log_phi)

  phi <- rtruncnorm(1, 0, mean = 0.5, sd = 2)
  curve(dgamma(x, shape = 1 / phi, rate = 1 / (mu * phi)), col = cc, add = T)
}


# Datos para Stan ---------------------------------------------------------

# Reemplazamos ceros por algo pequeño
dlong$caudal2 <- dlong$caudal
dlong$caudal2[dlong$caudal == 0] <- 0.0001

stan_data <- list(
  nc = nc, nobs = nobs,
  semana = dlong$semana,
  caudal = dlong$caudal2,
  cuenca = dlong$cca
)

# Modelo 01 (efectos fijos) ----------------------------------------------

model1 <- cmdstan_model("cuencas_01.stan")

opt1 <- model1$optimize(data = stan_data)

fit1 <- model1$sample(
  data = stan_data,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  parallel_chains = 4
)



plot(density(fit$draws("phi") |> as.numeric()))
# Modelo 02 (jerárquico) -------------------------------------------------

