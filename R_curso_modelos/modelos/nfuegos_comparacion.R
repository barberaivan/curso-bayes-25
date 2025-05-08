rm(list = ls()) # Limpiamos el entorno

# Eliminamos modelos compilados de stan, para evitar incompatibilidades
file.remove(list.files("modelos", pattern = "^[^.]+$", full.names = TRUE))

# Paquetes ----------------------------------------------------------------

library(cmdstanr)  # ajustar modelos con Stan
library(loo)       # comparación de modelos

# Funciones ---------------------------------------------------------------

# Función para calcular el log del promedio de la verosimilitud de forma
# numéricamente estable. Toma input el vector de log_lik correspondiente a una
# observación para todas las muestras de la posterior.
log_mean_exp <- function(x) {
  M <- max(x)
  return(M + log(mean(exp(x - M))))
}

set.seed(23434) # Para reproducibilidad


# Datos -------------------------------------------------------------------

# Cargamos datos
datos <- read.csv(file.path("datos", "barbera_data_fire_total_climate.csv"))

# Estandarizamos las predictoras
datoz <- datos[, c("fires", "fwi", "pp", "temp")]
vars <- c("fwi", "pp", "temp")
for (v in vars) {
  datoz[, v] <- as.numeric(scale(datos[, v]))
}

# Datos para stan. Usaremos la misma lista para ambos modelos, ya que uno es una
# simplificación del otro. A Stan no le molesta que le pasemos datos de más.
stan_data <- list(
  N = nrow(datoz),
  y = datoz$fires,
  fwi = datoz$fwi,
  pp = datoz$pp,
  temp = datoz$temp
)

# Modelo 1 ----------------------------------------------------------------

model_path1 <- file.path("modelos", "modelo_nfuegos_comp1.stan")
model1 <- cmdstan_model(model_path1)

fit1 <- model1$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 10000
)

# Verificamos el HMC
fit1$cmdstan_diagnose()

# Calculamos el PSIS-LOO
loo1 <- fit1$loo()
print(loo1) # Todo OK.

# Modelo 2 ----------------------------------------------------------------

model_path2 <- file.path("modelos", "modelo_nfuegos_comp2.stan")
model2 <- cmdstan_model(model_path2)

fit2 <- model2$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 10000
)

# Verificamos el HMC
fit2$cmdstan_diagnose()

# Calculamos el PSIS-LOO
loo2 <- fit2$loo()
print(loo2) # Algunos no tan bien, pero zafa.

# Comparación con PSIS-LOO ------------------------------------------------

print(loo1); print(loo2)
loo_compare(loo1, loo2)

# Buscamos elpd_loo más alto, o looic más bajo.
# Es mejor el modelo 1, aunque por poco.
# Por las dudas, hacemos K-fold CV.
# looic = - 2 * elpd
# (para ser expresado en unidades de devianza, como el AIC o WAIC)

# K-fold CV ---------------------------------------------------------------

# Como el modelo corre muy rápido, podríamos definir K = N, pero el código
# con K < N es un poco más rebuscado, así que usaremos K = N / 2 a modo de
# ejemplo.

N <- nrow(datos)
O <- 2 # observaciones dejadas afuera por partición
K <- N / O # particiones
fold_id <- rep(1:K, each = O)
# aleatorizamos el id de cada observación
fold_id <- fold_id[sample(1:N, N)]
row_ids <- 1:N

# Lista para llenar con matrices de loglik de las observaciones
log_lik_list1 <- vector("list", length = K)

# Compilamos modelo (se usará un código genérico)
model1CV <- cmdstan_model(file.path("modelos", "modelo_nfuegos_comp1_KfoldCV.stan"))

# Loop para obtener la matriz log_lik de las observaciones dejadas fuera del ajuste
for (k in 1:K) {
  print(k)

  # Preparamos datos
  stan_data_k <- list(
    N = nrow(datoz),
    y = datoz$fires,
    fwi = datoz$fwi,
    pp = datoz$pp,
    temp = datoz$temp,
    O = O,
    train_ids = row_ids[fold_id != k],
    test_ids = row_ids[fold_id == k]
  )

  # Ajustamos modelo
  fit_k <- model1CV$sample(
    data = stan_data_k,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 10000
  )

  # Extraemos matriz de loglik y la guardamos en la lista
  log_lik_list1[[k]] <- fit_k$draws("log_lik", format = "draws_matrix")
}

# Juntamos todas las loglik matrices en una sola
log_lik1_cv <- do.call("cbind", log_lik_list1)

# Calculamos el elpd de cada observación
elpd_pointwise1 <- apply(log_lik1_cv, 2, log_mean_exp)
elpd1_cv <- sum(elpd_pointwise1)
elpd1_cv_se <- sqrt(N) * sd(elpd_pointwise1)


## Repetimos el procedimiento para el modelo 2

# Lista para llenar con matrices de loglik de las observaciones
log_lik_list <- vector("list", length = K)

# Compilamos modelo (se usará un código genérico)
model2CV <- cmdstan_model(file.path("modelos", "modelo_nfuegos_comp2_KfoldCV.stan"))

# Loop para obtener la matriz log_lik de las observaciones dejadas fuera del ajuste
for (k in 1:K) {
  print(k)

  # Preparamos datos
  stan_data_k <- list(
    N = nrow(datoz),
    y = datoz$fires,
    fwi = datoz$fwi,
    pp = datoz$pp,
    temp = datoz$temp,
    O = O,
    train_ids = row_ids[fold_id != k],
    test_ids = row_ids[fold_id == k]
  )

  # Ajustamos modelo
  fit_k <- model2CV$sample(
    data = stan_data_k,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 10000
  )

  # Extraemos matriz de loglik y la guardamos en la lista
  log_lik_list[[k]] <- fit_k$draws("log_lik", format = "draws_matrix")
}

# Juntamos todas las loglik matrices en una sola
log_lik2_cv <- do.call("cbind", log_lik_list)

# Calculamos el elpd de cada observación
elpd_pointwise2 <- apply(log_lik2_cv, 2, log_mean_exp)
elpd2_cv <- sum(elpd_pointwise2)
elpd2_cv_se <- sqrt(N) * sd(elpd_pointwise2)


## Comparación basada en K-fold CV

# Comparamos de a pares, igual que hace loo_compare. De esta manera, el
# error estándar de la diferencia en elpd es mucho menor.
# (ver ?loo_compare)

# restamos peor - mejor
elpd_diff_pointwise <- elpd_pointwise2 - elpd_pointwise1
elpd_diff <- sum(elpd_diff_pointwise)
elpd_diff_se <- sqrt(N) * sd(elpd_diff_pointwise)

elpd_diff    # el modelo 2 tiene una loglik menor
elpd_diff_se # pero el error estándar suele ser relativamente grande.

# Los resultados no difieren mucho de la comparación aproximada basada en
# PSIS-LOO, pero a este resultado le confiamos más.
# De utilizar algún modelo, usaríamos las versiones ajustadas al set de datos
# completo.