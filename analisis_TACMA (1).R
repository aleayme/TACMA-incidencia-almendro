set.seed(123)
library(jagsUI)
library(ggplot2)
library(mgcv)
library(gridExtra)
library(caret)

# Carga 
load("data/caso1.RData")

datos <- na.omit(caso1)
datos$isol      <- as.factor(datos$isol)
datos$replicate <- as.factor(datos$replicate)
datos$temp      <- as.numeric(datos$temp)
datos$wet       <- as.numeric(datos$wet)

# Identificador único por planta 
datos$id_planta <- as.factor(1:nrow(datos))

# Variable respuesta: proporción de hojas infectadas por planta.
datos$incidencia <- datos$damaged_leaves / datos$total_leaves

# Transformaciones para la exploración gráfica
inc_acot <- pmin(pmax(datos$incidencia, 0.001), 0.999)
datos$logit_incidencia   <- qlogis(inc_acot)
datos$cloglog_incidencia <- log(-log(1 - inc_acot))

tasa   <- "incidencia"
y_name <- "incidencia"
covars <- c("temp", "wet") 
facts  <- c("isol", "replicate") 

# Summary
summary_df <- as.data.frame.matrix(summary(datos[, c(tasa, covars, facts)]))
print(summary_df)

# Freq. table
freq_isol <- as.data.frame(table(Aislado = datos$isol))
freq_rep <- as.data.frame(table(Réplica = datos$replicate))

print(freq_isol)
print(freq_rep)

#Scatter plot
pairs(datos[, c(tasa, covars, "logit_incidencia")], 
      col = "darkblue", main = "Matriz de Dispersión")

# Corr. plot
matriz_cor <- cor(datos[, c(tasa, covars)])
print(matriz_cor, digits = 3)

# Tabla descriptiva
media_y <- mean(datos[[tasa]], na.rm = TRUE)
var_y   <- var(datos[[tasa]], na.rm = TRUE)

momentos_df <- data.frame(Métrica = c("Media", "Varianza"), Valor = c(media_y, var_y))
print(momentos_df, digits = 3)

# Tabla dispersión Pearson
# GLM Binomial preliminar, permite evaluar sobredispersión (residuos de Pearson).
# 'replicate' anidada en 'isol' (R1-R2 en phal_4, R3-R4 en phal_45)
mod_binomial <- glm(incidencia ~ temp + wet + isol + replicate, 
                    family = binomial, 
                    weights = total_leaves,
                    data = datos)

dispersion_resid <- sum(residuals(mod_binomial, type = "pearson")^2) / df.residual(mod_binomial)

estado_disp <- if (dispersion_resid > 1.5) {
  "Indicio claro de sobredispersión (Parámetro > 1.5)"
} else {
  "Parámetro cercano a 1. No se detecta sobredispersión"
}

tabla_disp <- data.frame(Estadístico = "Ratio Pearson / df", 
                         Valor = dispersion_resid, 
                         Diagnóstico = estado_disp)

print(tabla_disp, digits = 3)

# Histograma plot
hist(datos[[tasa]], main = "Distribución de la Incidencia", 
     xlab = "Incidencia (Proporción)", breaks = 20, col = "gray")


# Boxplot 
par(mfrow = c(1, 2))
boxplot(datos[[tasa]] ~ datos$isol, 
        main = "Incidencia vs Aislado", 
        xlab = "Aislado del patógeno", ylab = "Incidencia", col = "lightblue")

boxplot(datos[[tasa]] ~ datos$replicate, 
        main = "Incidencia vs Réplica", 
        xlab = "Réplica experimental", ylab = "Incidencia", col = "lightgreen")
par(mfrow = c(1, 1))

# Logit boxplot
par(mfrow = c(1, 2))
boxplot(datos[["logit_incidencia"]] ~ datos$isol, 
        main = "Logit vs Aislado", 
        xlab = "Aislado del patógeno", ylab = "logit(Incidencia)", col = "lightblue")

boxplot(datos[["logit_incidencia"]] ~ datos$replicate, 
        main = "Logit vs Réplica", 
        xlab = "Réplica experimental", ylab = "logit(Incidencia)", col = "lightgreen")
par(mfrow = c(1, 1))

# Clima plot
p3 <- ggplot(datos, aes(x = temp, y = logit_incidencia)) +
  geom_point(alpha = 0.5, aes(color = isol)) +
  geom_smooth(method = "loess", formula = y ~ x, color = "darkred") + 
  labs(title = "Logit de la Incidencia vs Temperatura",
       x = "Temperatura (ºC)",
       y = "Logit(Incidencia)") +
  theme_minimal()

p4 <- ggplot(datos, aes(x = wet, y = logit_incidencia)) +
  geom_point(alpha = 0.5, aes(color = isol)) +
  geom_smooth(method = "loess", formula = y ~ x, color = "darkblue") +
  labs(title = "Logit de la Incidencia vs Humedad",
       x = "Periodo de humectación (h)",
       y = "Logit(Incidencia)") +
  theme_minimal()

grid.arrange(p3, p4, nrow = 2)

# Funcion link 
# Comparación entre logit y cloglog.

# Modelo GAMM preliminar con enlace Logit
mod0_logit <- gam(incidencia ~ isol + s(temp, k = 5) + s(wet, k = 4) + 
                     s(replicate, bs = "re") + s(id_planta, bs = "re"),
                   family = binomial(link = "logit"),
                   weights = total_leaves, data = datos, method = "REML")

# Modelo GAMM preliminar con enlace Cloglog
mod0_cloglog <- gam(incidencia ~ isol + s(temp, k = 5) + s(wet, k = 4) + 
                     s(replicate, bs = "re") + s(id_planta, bs = "re"),
                   family = binomial(link = "cloglog"),
                   weights = total_leaves, data = datos, method = "REML")

# AIC
tabla_devianza <- data.frame(
  Enlace = c("Logit", "Cloglog"),
  `Devianza Residual` = c(deviance(mod0_logit), deviance(mod0_cloglog)),
  `df Resid.` = c(df.residual(mod0_logit), df.residual(mod0_cloglog)),
  check.names = FALSE
)

print(tabla_devianza, digits = 3)

# Modelos
# M0: efectos aditivos no paramétricos.
mod0 <- gam(incidencia ~ isol + s(temp, k = 5) + s(wet, k = 4) + 
              s(replicate, bs = "re") + s(id_planta, bs = "re"),
            family = binomial(link = "logit"),
            weights = total_leaves, data = datos, method = "REML")

# M1: interacción no lineal mediante tensor suavizado.
mod1 <- gam(incidencia ~ isol + te(temp, wet, k = c(5, 4)) + 
              s(replicate, bs = "re") + s(id_planta, bs = "re"),
            family = binomial(link = "logit"),
            weights = total_leaves, data = datos, method = "REML")

# M2: interacción paramétrica (lineal).
mod2 <- gam(incidencia ~ isol + temp * wet + 
              s(replicate, bs = "re") + s(id_planta, bs = "re"),
            family = binomial(link = "logit"),
            weights = total_leaves, data = datos, method = "REML")

# M3: modelo de referencia sin efectos aleatorios.
mod3 <- gam(incidencia ~ isol + temp * wet, 
            family = binomial(link = "logit"), 
            weights = total_leaves, data = datos, method = "REML")

# Diagnosticos
k.check(mod0); summary(mod0)
k.check(mod1); summary(mod1)
k.check(mod2); summary(mod2)
summary(mod3)

# Tabla seleccion 
lista_mod <- list(mod0, mod1, mod2, mod3)

tabla_seleccion <- data.frame(
  `Especificación Estructural` = c(
    "M0: Preliminar aditivo (s(temp) + s(wet))",
    "M1: Interacción tensor (te(temp, wet))",
    "M2: Interacción paramétrica (temp * wet)",
    "M3: Basal (GLM)"
  ),
  `edf Totales`          = sapply(lista_mod, function(m) sum(m$edf)),
  `AIC`                  = sapply(lista_mod, AIC),
  `% Devianza Explicada` = sapply(lista_mod, function(m) sprintf("%.1f%%", 100 * summary(m)$dev.expl)),
  check.names = FALSE
)

print(tabla_seleccion, digits = 5)

# ---- diag-plot ----
par(mfrow = c(2, 2))
gam.check(mod2)
par(mfrow = c(1, 1))

# ---- datosbay ----
# Estandarizacion covariables
datos$temp_std <- as.numeric(scale(datos$temp))
datos$wet_std  <- as.numeric(scale(datos$wet))

# El modelo JAGS asume un aislado con 2 niveles.
stopifnot(nlevels(datos$isol) == 2)

datos_jags <- list(
  n = nrow(datos),
  y_damaged = datos$damaged_leaves,
  n_total = datos$total_leaves,
  
  temp = datos$temp_std,
  wet  = datos$wet_std,
  
  isol = as.numeric(datos$isol) - 1,      # 0 = phal_4, 1 = phal_45 (orden alfabético de los niveles)
  rep  = as.numeric(datos$replicate),     # 1 a nlevels
  n_reps = nlevels(datos$replicate)
)

# JAGS
cat(file = "ModeloCaso1_Definitivo.jags", "model {
  for (i in 1:n) {
    y_damaged[i] ~ dbin(pi[i], n_total[i])
    logit(pi[i]) <- beta0 + 
                    beta_isol * isol[i] + 
                    beta_temp * temp[i] + 
                    beta_wet * wet[i] + 
                    beta_int * (temp[i] * wet[i]) + 
                    b_rep[rep[i]] + 
                    epsilon[i]
    
    epsilon[i] ~ dnorm(0, tau_obs)
  }
  
  for (j in 1:n_reps) { 
    b_rep[j] ~ dnorm(0, tau_rep) 
  }
  
  beta0     ~ dnorm(0, 0.1)      
  beta_isol ~ dnorm(0, 0.1)
  beta_temp ~ dnorm(0, 0.1)
  beta_wet  ~ dnorm(0, 0.1)
  beta_int  ~ dnorm(0, 0.1)
  
  tau_rep <- pow(sd_rep, -2)
  sd_rep  ~ dunif(0, 10)
  
  tau_obs <- pow(sd_obs, -2)
  sd_obs  ~ dunif(0, 10)
}")

# Parametros
parametros <- c("beta0", "beta_isol", "beta_temp", "beta_wet", "beta_int", "sd_rep", "sd_obs")

# Valores iniciales
set.seed(123)
cadena <- 0
iniciales <- function() {
  cadena <<- cadena + 1
  list(
    beta0     = rnorm(1, 0, 0.1),
    beta_isol = rnorm(1, 0, 0.1),
    beta_temp = rnorm(1, 0, 0.1),
    beta_wet  = rnorm(1, 0, 0.1),
    beta_int  = rnorm(1, 0, 0.1),
    sd_rep    = runif(1, 0.1, 1),
    sd_obs    = runif(1, 0.1, 1),
    .RNG.name = "base::Wichmann-Hill",
    .RNG.seed = 1000 + cadena
  )
}

Resul_Definitivo <- jags(
  data = datos_jags,
  inits = iniciales,
  parameters.to.save = parametros,
  model.file = "ModeloCaso1_Definitivo.jags",
  n.chains = 3,
  n.iter = 50000,     
  n.burnin = 15000,   
  n.thin = 10,
  parallel = TRUE
)

# Convergencia
# Diagnóstico de convergencia: R-hat < 1.05 y n.eff elevado. 
tabla_conv <- data.frame(
  Parámetro = parametros,
  `R-hat`   = Resul_Definitivo$summary[parametros, "Rhat"],
  `n.eff`   = Resul_Definitivo$summary[parametros, "n.eff"],
  check.names = FALSE, row.names = NULL
)
print(tabla_conv, digits = 3)

# Validación cruzada 5-Fold. Se evalúa la capacidad predictiva poblacional
# excluyendo los efectos aleatorios (plantas no observadas previamente).

set.seed(123) 
K <- 5
folds <- createFolds(datos$incidencia, k = K, list = TRUE)

mae_cv <- numeric(K)

for (i in 1:K) {
  idx_test   <- folds[[i]]
  train_data <- datos[-idx_test, ]
  test_data  <- datos[idx_test, ]
  
  mod_train <- gam(incidencia ~ isol + temp * wet + 
                     s(replicate, bs = "re") + s(id_planta, bs = "re"), 
                   family = binomial(link = "logit"), 
                   weights = total_leaves,
                   data = train_data, 
                   method = "REML")
  
  preds <- predict(mod_train, 
                   newdata = test_data, 
                   type = "response", 
                   exclude = c("s(id_planta)", "s(replicate)"),
                   newdata.guaranteed = TRUE)
  
  mae_cv[i] <- mean(abs(test_data$incidencia - preds), na.rm = TRUE)
}

mae_final <- mean(mae_cv)

# Tabla frecuentista
pt <- summary(mod2)$p.table
stopifnot(nrow(pt) == 5)   # Intercepto, isol, temp, wet, temp:wet
z  <- qnorm(0.975)

tabla_frecuentista <- data.frame(
  Predictor = c("Intercepto", "Aislado (phal_45)", "Temperatura", "Humedad", "Interacción (Temp:Wet)"),
  `Log-Odds (Beta)` = pt[, "Estimate"],
  `Odds Ratio (OR)` = exp(pt[, "Estimate"]),
  `IC 2.5%`  = exp(pt[, "Estimate"] - z * pt[, "Std. Error"]),
  `IC 97.5%` = exp(pt[, "Estimate"] + z * pt[, "Std. Error"]),
  check.names = FALSE, row.names = NULL
)

print(tabla_frecuentista, digits = 5)

# Tabla bayesiana
pb <- Resul_Definitivo$summary[c("beta0", "beta_isol", "beta_temp", "beta_wet", "beta_int"), ]

tabla_bayesiana <- data.frame(
  Predictor = c("Intercepto", "Aislado (phal_45)", "Temperatura (Std)", "Humedad (Std)", "Interacción (Std)"),
  `Media Log-Odds`   = pb[, "mean"],
  `Odds Ratio Medio` = exp(pb[, "mean"]),
  `IC 2.5%`  = exp(pb[, "2.5%"]),     # exp() es monótona: los cuantiles se transforman directamente
  `IC 97.5%` = exp(pb[, "97.5%"]),
  check.names = FALSE, row.names = NULL
)

print(tabla_bayesiana, digits = 3)


# Todas las combinaciones posibles de clima para predecir
newdat <- expand.grid(
  temp = seq(min(datos$temp), max(datos$temp), length = 100),
  wet  = seq(min(datos$wet),  max(datos$wet),  length = 100),
  isol = levels(datos$isol)[1],
  id_planta = datos$id_planta[1], 
  replicate = datos$replicate[1]  
)

# Predicción excluyendo el ruido de la planta y de la réplica
newdat$pred <- predict(mod2, newdata = newdat, type = "response", 
                       exclude = c("s(id_planta)", "s(replicate)"))

p_sup <- ggplot(newdat, aes(x = temp, y = wet, fill = pred)) +
  geom_raster() +
  scale_fill_viridis_c() +
  labs(title = "Superficie de probabilidad predicha",
       x = "Temperatura (ºC)",
       y = "Periodo de humectación (h)",
       fill = "Probabilidad") +
  theme_minimal()
print(p_sup)

# Figura 1 del informe (informe/Informe_TACMA.tex)
dir.create("informe/figuras", showWarnings = FALSE, recursive = TRUE)
ggsave("informe/figuras/figura1.png", p_sup, width = 6, height = 5, dpi = 300)

# Prueba de la herramienta INLA 
library(INLA)

datos$rep_id    <- as.integer(datos$replicate)
datos$planta_id <- seq_len(nrow(datos))

# Modelo M2 con priors PC: P(sigma > u) = alpha para sigma_u (réplica) y sigma_v (planta).
# Coeficientes: N(0, 10) (prec = 0.1), igual que en JAGS.
ajusta_inla <- function(u_rep = 1, u_planta = 1, alpha = 0.01) {
  pc <- function(u) list(prec = list(prior = "pc.prec", param = c(u, alpha)))
  inla(
    damaged_leaves ~ isol + temp_std * wet_std +
      f(rep_id,    model = "iid", hyper = pc(u_rep)) +
      f(planta_id, model = "iid", hyper = pc(u_planta)),
    family  = "binomial",
    Ntrials = datos$total_leaves,
    data    = datos,
    control.fixed   = list(mean = 0, prec = 0.1, mean.intercept = 0, prec.intercept = 0.1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
}

# sigma = 1/sqrt(precisión): mediana e intervalo del 95 %
resumen_sd <- function(fit, nombre) {
  m <- inla.tmarginal(function(x) 1 / sqrt(x), fit$marginals.hyperpar[[nombre]])
  q <- inla.qmarginal(c(0.025, 0.5, 0.975), m)
  c(mediana = q[2], q2.5 = q[1], q97.5 = q[3])
}

# --- Análisis principal (u = 1) ---
fit_inla <- ajusta_inla(u_rep = 1, u_planta = 1)

fit_inla$summary.fixed                                                  # log-odds
exp(fit_inla$summary.fixed[, c("mean", "0.025quant", "0.975quant")])    # odds ratios
rbind(sd_rep   = resumen_sd(fit_inla, "Precision for rep_id"),
      sd_planta = resumen_sd(fit_inla, "Precision for planta_id"))

# --- Sensibilidad al prior de sigma_u: u = 0.5, 1, 3, 10 ---
u_vals <- c(0.5, 1, 3, 10)

fmt <- function(fx, nombre)
  sprintf("%.3f [%.3f, %.3f]", fx[nombre, "mean"], fx[nombre, "0.025quant"], fx[nombre, "0.975quant"])

sens <- do.call(rbind, lapply(u_vals, function(u) {
  fit <- ajusta_inla(u_rep = u)
  fx  <- fit$summary.fixed
  s_u <- resumen_sd(fit, "Precision for rep_id")
  data.frame(
    u_prior      = u,
    Aislado      = fmt(fx, "isolphal_45"),
    Temperatura  = fmt(fx, "temp_std"),
    Humectacion  = fmt(fx, "wet_std"),
    Interaccion  = fmt(fx, "temp_std:wet_std"),
    sigma_u_med  = round(s_u["mediana"], 3),
    sigma_u_IC95 = sprintf("[%.2f, %.2f]", s_u["q2.5"], s_u["q97.5"])
  )
}))
print(sens, row.names = FALSE)
