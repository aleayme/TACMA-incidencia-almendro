# Modelización de la incidencia fúngica en almendro (Caso 1)

Análisis de la proporción de hojas dañadas por dos aislados de *Phaeoacremonium* (`phal_4`, `phal_45`) en función de la temperatura y del periodo de humectación, mediante modelos binomiales mixtos (GAMM, mgcv) y su aproximación bayesiana (JAGS; opcionalmente INLA y Stan/brms).

**Autores:** Alejandro Aymerich Pérez, Rafael Casajús Clemente. 

## Estructura

```
.
├── codigoTACMA.R 
├── data/
│   └── caso1.RData      # datos: 700 plantas (3 con valores perdidos)
├── .gitignore
└── README.md
```

## Diseño experimental

Factorial completo: 7 temperaturas × 5 periodos de humectación × 4 réplicas × 5 plantas. Las réplicas están **anidadas** en el aislado (R1-R2 en `phal_4`; R3-R4 en `phal_45`).

## Requisitos

- R (≥ 4.1) y RStudio (recomendado)
- [JAGS](https://mcmc-jags.sourceforge.io/) instalado en el sistema
- Una distribución de LaTeX con `xelatex` (p. ej. TinyTeX: `tinytex::install_tinytex()`)
- Paquetes de R:

```r
install.packages(c("rmarkdown", "knitr", "jagsUI", "ggplot2", "mgcv", "gridExtra", "caret"))
```

## Reproducibilidad

Se fijan las semillas de R (`set.seed(123)`) y la del generador de JAGS por cadena. Al final del informe se listan las versiones de R y de los paquetes.
