setwd("/Users/user/Desktop/Scripts Perdita maritima")
perdita<- read.csv(choose.files())
library(DHARMa)
library(lme4)
library(ggplot2)
library(lmtest)
library(pscl)
library(glmmTMB)
library(car)
library(performance)
library(sjPlot)
library(ggeffects)
library(lattice)
library(AER)
library(MASS)
library(see)
library(broom.mixed)
library(purrr)
library(htmltools)
library(webshot2)
library(MuMIn)

perdita$weather<-factor(perdita$weather,
                        levels = c("Cloudy","Rain","Sunny"))
perdita$weather<-relevel(perdita$weather, ref = "Sunny")

perdita$time_num<-as.numeric(gsub(":",".",perdita$time))
perdita$time_centered<-perdita$time_num-10
perdita$observer<-factor(perdita$observer)
perdita$day<-factor(perdita$day)

variables_respuesta<-c("males_activity",
                       "females_foraging","males_following",
                       "males_bumping","m_vs_m")
print(variables_respuesta)
length(variables_respuesta)
dput(variables_respuesta)

formula_fija<- ~weather+poly(time_num,2)+(1|observer)+(1|day)
class(formula_fija)

lista_modelos<-vector("list", length(variables_respuesta))
names(lista_modelos)<- variables_respuesta
str(lista_modelos)

for (i in 1:length(variables_respuesta)) {
  var_actual<- variables_respuesta [i]
  formula_completa<- as.formula(paste(var_actual, "~weather+poly(time_num,2)+(1|observer)+(1|day)"))
  modelo_actual<- glmmTMB(formula_completa, data = perdita, family = nbinom2)
  lista_modelos[[var_actual]]<- modelo_actual
}
str(lista_modelos, max.level = 1)
class(lista_modelos[[1]])
summary(modelo_actual)

clases<- sapply(lista_modelos, class)
print(clases)

modelo_ejemplo<-lista_modelos[["males_activity"]]
Sim_res_ejemplo<-simulateResiduals(modelo_ejemplo, n=250)
plot(Sim_res_ejemplo)
test_residual_ejemplo<- testResiduals(Sim_res_ejemplo)
print(test_residual_ejemplo)


resultados_diagnostico<-data.frame(
  variable=variables_respuesta,
  ks_p= NA,
  dispersion_p=NA,
  outliers_p=NA,
  adecuado=NA
)

for (i in 1:length(variables_respuesta)) {
  var_actual<-variables_respuesta[i]
  modelo_actual<-lista_modelos[[var_actual]]
  
  sim_res<-simulateResiduals(modelo_actual,n=250)
  
  resultados_diagnostico$ks_p[i]<-testUniformity(sim_res, plot = FALSE)$p.value
  resultados_diagnostico$dispersion_p[i]<-testDispersion(sim_res, plot = FALSE)$p.value
  resultados_diagnostico$outliers_p[i]<-testOutliers(sim_res,plot = FALSE)$p.value
  
  resultados_diagnostico$adecuado[i]<-all(
    c(resultados_diagnostico$ks_p[i],
      resultados_diagnostico$dispersion_p[i],
      resultados_diagnostico$outliers_p[i]) > 0.05
  )
  cat("Complete diagnosis for:",var_actual, "\n")
}

print(resultados_diagnostico)

######Tablas de efectos fijos######
tabla_fijos<-map_dfr(lista_modelos,
                     ~tidy(.x, effects= "fixed", conf.int = TRUE),
                     .id= "variable_respuesta")
print(tabla_fijos, n=30)
write.csv(tabla_fijos, "coeficientes fijos 5 modelos.csv",row.names=FALSE)

####tabla de efectos aleatorios####

tabla_aleatorios<-map_dfr(lista_modelos,
                           ~tidy(.x, effects= "ran_pars", conf.int = TRUE),
                           .id = "variable_respuesta")
print(tabla_aleatorios)
write.csv(tabla_aleatorios, "vaarianzas aleatorias de 5 modelos.csv", row.names = FALSE)

####tabla elegante####

# --- TABLA COMPLETA CON KABLEEXTRA (INCLUYENDO R²) ---
library(kableExtra)
library(broom.mixed)
library(MuMIn)
library(tidyverse)

# 1. Función para extraer TODO de un modelo
extraer_todo_modelo <- function(modelo, nombre_comportamiento) {
  
  # A. Coeficientes (efectos fijos)
  coefs <- tidy(modelo, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(effect == "fixed" & term %in% c("(Intercept)", "weatherCloudy", "weatherRain")) %>%
    mutate(
      valor_formateado = sprintf("%.2f (%.2f–%.2f)%s", 
                                 estimate, conf.low, conf.high,
                                 case_when(
                                   p.value < 0.001 ~ "***",
                                   p.value < 0.01 ~ "**",
                                   p.value < 0.05 ~ "*",
                                   TRUE ~ ""
                                 ))
    ) %>%
    select(term, valor_formateado)
  
  # B. R² (delta method)
  r2_vals <- r.squaredGLMM(modelo)["delta", ]
  r2_text <- sprintf("%.3f / %.3f", r2_vals["R2m"], r2_vals["R2c"])
  
  # C. Estadísticos de ajuste
  data.frame(
    Comportamiento = nombre_comportamiento,
    Intercepto = coefs$valor_formateado[coefs$term == "(Intercept)"],
    `Clima: Nublado` = coefs$valor_formateado[coefs$term == "weatherCloudy"],
    `Clima: Lluvia` = coefs$valor_formateado[coefs$term == "weatherRain"],
    `R² (marginal/condicional)` = r2_text,
    AIC = sprintf("%.1f", AIC(modelo)),
    Deviance = sprintf("%.1f", deviance(modelo)),
    check.names = FALSE
  )
}

# 2. Nombres de comportamientos
nombres_comportamientos <- c(
  "Males activity",
  "Females foraging", 
  "Males following females",
  "Males bumping", 
  "male vs male"
)

# 3. Extraer datos de los 5 modelos
datos_completos <- map2_dfr(
  lista_modelos,
  nombres_comportamientos,
  extraer_todo_modelo
)

# 4. Crear tabla profesional
tabla_final <- datos_completos %>%
  rename(
    Comportamiento = Comportamiento,
    `(Intercept)` = Intercepto,
    `Clima: Nublado vs. Soleado` = `Clima: Nublado`,
    `Clima: Lluvia vs. Soleado` = `Clima: Lluvia`
  ) %>%
  kbl(
    caption = "Results of GLMM models for Perdita maritima behavior",
    align = c("l", "c", "c", "c", "c", "c", "c"),
    col.names = c("Behavior", "Intercept", "Cloudy", "Rainy", 
                  "R² (marginal/conditional)", "AIC", "Deviance"),
    booktabs = TRUE
  ) %>%
  kable_classic(
    full_width = FALSE,
    font_size = 11,
    html_font = "Arial"
  ) %>%
  add_header_above(c(" " = 1, "Coefficients (IRR with IC 95%)" = 3, 
                     "Fit statistics" = 3)) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, bold = TRUE) %>%  # Columna de comportamientos en negrita
  column_spec(5, bold = TRUE) %>%  # Columna de R² en negrita
  footnote(
    general = "IRR = Incidence Rate Ratios. 95% confidence intervals. *** p<0.001, ** p<0.01, * p<0.05. R² calculados con método delta (MuMIn::r.squaredGLMM). Modelos: nbinom2 ~ weather + poly(time_num,2) + (1|observer) + (1|day).",
    general_title = "Note: ",
    footnote_as_chunk = TRUE
  )

# 5. Mostrar en RStudio
tabla_final

# 6. Guardar en múltiples formatos
# HTML (se ve mejor)
save_kable(tabla_final, "tabla_completa_perdita.html")

# PNG (imagen)
save_kable(tabla_final, "table_perdita.png", 
           zoom = 2, vwidth = 1000)

# PDF (para publicaciones)
# save_kable(tabla_final, "tabla_completa_perdita.pdf")

# 7. Mensaje de confirmación
cat("\n✅ TABLAS CREADAS CON ÉXITO:\n")
cat("   📄 HTML: tabla_completa_perdita.html\n")
cat("   🖼️  PNG:  tabla_completa_perdita.png\n")
cat("\n📊 RESUMEN ESTADÍSTICO:\n")
print(datos_completos)


#######grafico de efectos marginales#####

library(ggeffects)
library(ggplot2)
library(patchwork)

cat("Valores únicos de time_num en tus datos:\n")
print(sort(unique(perdita$time_num)))
# 2. Definir la secuencia CORRECTA para predicciones
# Desde el mínimo hasta el máximo, con incremento de 0.1 para suavizar
time_seq <- seq(min(perdita$time_num, na.rm = TRUE), 
                max(perdita$time_num, na.rm = TRUE), 
                by = 0.1)
print(round(time_seq, 1))

predicciones_lista <- list()
nombres_bonitos <- c("Males activity",
                     "Females foraging", 
                     "Males following females",
                     "Males bumping", 
                     "male vs male")
for (i in 1:5) {
  modelo <- lista_modelos[[i]]
  nombre <- nombres_bonitos[i]
  
  cat("Calculando predicciones para:", nombre, "\n")
  
  # OPCIÓN A: Predicciones con efectos fijos solamente (recomendado para patrones generales)
  # Usamos type = "fixed" o type = "fe" (fixed effects)
  tryCatch({
    # Intentar con type = "fixed" primero
    pred <- ggpredict(modelo, 
                      terms = "time_num [all]",  # "all" usa todo el rango
                      type = "fixed")            # <-- CORRECTO
    
  }, error = function(e) {
    cat("Intento 1 falló:", e$message, "\n")
    # Intentar alternativa
    pred <<- ggpredict(modelo, 
                       terms = paste0("time_num [", time_min, ":", time_max, " by=0.1]"))
  })
  
  # Añadir nombre
  pred$comportamiento <- nombre
  predicciones_lista[[i]] <- pred
  
  cat("   ✓ Completado\n")
}
predicciones_todas <- bind_rows(predicciones_lista)
cat("\n=== RESUMEN DE PREDICCIONES ===\n")
print(head(predicciones_todas))
cat("Número de filas:", nrow(predicciones_todas), "\n")

crear_etiqueta_hora <- function(x) {
  sapply(x, function(t) {
    horas <- floor(t)
    minutos <- round((t - horas) * 100)
    if (minutos == 30) {
      paste0(horas, ":30")
    } else {
      paste0(horas, ":00")
    }
  })
}

# 7. Gráfico MULTIPANEL

etiquetas_seleccionadas <- sort(unique(perdita$time_num))
indices_mostrar<-seq(1,length(etiquetas_seleccionadas), by= 2)

if (nrow(predicciones_todas) > 0) {
  grafico_completo <- ggplot(predicciones_todas, 
                             aes(x = x, y = predicted, 
                                 ymin = conf.low, ymax = conf.high,
                                 fill = comportamiento, color = comportamiento)) +
    geom_ribbon(alpha = 0.2, linetype = 0) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ comportamiento, scales = "free_y", ncol = 2) +
    labs(
      title = "Daily activity patterns in Perdita maritima",
      subtitle = "GLMM model predictions (fixed effects)",
      x = "Time",
      y = "Predicted count for 10 minutes",
      caption = "Shaded bands = 95% confidence intervals"
    ) +
    scale_x_continuous(
      breaks = etiquetas_seleccionadas[indices_mostrar],
      labels = crear_etiqueta_hora(etiquetas_seleccionadas[indices_mostrar])
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 10),
      panel.spacing = unit(1, "lines"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
    )}
print(grafico_completo)  


####grafico ultra simple#####
# --- VERSIÓN SIMPLE que SIEMPRE funciona ---
library(ggplot2)

# 1. Crear datos para predicción manual
nuevos_datos <- expand.grid(
  time_num = seq(min(perdita$time_num), max(perdita$time_num), by = 0.1),
  weather = "Sunny",
  observer = levels(perdita$observer)[1],  # Primer observador como referencia
  day = levels(perdita$day)[1]             # Primer día como referencia
)

# 2. Predecir para cada modelo
nombres_bonitos <- c("Males activity",
                     "Females foraging", 
                     "Males following females",
                     "Males bumping", 
                     "male vs male")

graficos_individuales <- list()

for (i in 1:5) {
  modelo <- lista_modelos[[i]]
  
  # Predicción manual
  pred <- predict(modelo, newdata = nuevos_datos, type = "response", se.fit = TRUE)
  
  nuevos_datos$predicho <- pred$fit
  nuevos_datos$conf.low <- pred$fit - 1.96 * pred$se.fit
  nuevos_datos$conf.high <- pred$fit + 1.96 * pred$se.fit
  nuevos_datos$comportamiento <- nombres_bonitos[i]
  
  # Crear gráfico individual
  g <- ggplot(nuevos_datos, aes(x = time_num, y = predicho)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
    geom_line(color = "darkblue", size = 1) +
    labs(
      title = nombres_bonitos[i],
      x = "Time",
      y = "Predicted count"
    ) +
    scale_x_continuous(
      breaks = sort(unique(perdita$time_num)),
      labels = function(x) {
        horas <- floor(x)
        minutos <- ifelse(round((x - horas) * 100) == 30, "30", "00")
        paste0(horas, ":", minutos)
      }
    ) +
    theme_minimal()
  
  graficos_individuales[[i]] <- g
}

# 3. Combinar gráficos
library(patchwork)
grafico_final <- wrap_plots(graficos_individuales, ncol = 2) +
  plot_annotation(
    title = "Daily activity patterns in Perdita maritima",
    subtitle = "GLMM model predictions (sunny days)"
  )

print(grafico_final)

ggsave("patrones_manual.png", grafico_final, width = 10, height = 8)

####### weather effects######
library(dotwhisker)

# Extraer coeficientes de clima
coef_clima <- tabla_fijos %>%
  filter(term %in% c("weatherCloudy", "weatherRain")) %>%
  mutate(
    term = recode(term,
                  "weatherCloudy" = "Cloudy vs Sunny",
                  "weatherRain" = "Rainy vs Sunny"),
    variable_respuesta = recode(variable_respuesta,
                                "males_activity" = "Males activity",
                                "females_foraging" = "Females foraging",
                                "males_following" = "Males following",
                                "males_bumping" = "Males bumping",
                                "m_vs_m" = "Male-male fights")
  )

ggplot(coef_clima, aes(x = estimate, y = variable_respuesta, color = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                 height = 0.2, position = position_dodge(width = 0.5)) +
  labs(title = "Weather effects on Perdita maritima behaviors",
       subtitle = "Incidence Rate Ratios (IRR) with 95% CI",
       x = "IRR (log scale)", y = "Behavior") +
  scale_x_log10() +
  theme_minimal() +
  theme(legend.title = element_blank())

#########grafico 2x2##########
# ============================================
# PANEL 2×2 COMPLETO - CÓDIGO EJECUTABLE
# ============================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(ggeffects)

# -----------------------------------------------------
# 1. PREPARAR DATOS NECESARIOS
# -----------------------------------------------------

# A) Datos para gráfico temporal (P1)
# (usamos los que ya tienes de predicciones_todas)

# B) Datos para forest plot (P2)
# Extraer coeficientes de clima de tu tabla_fijos
coef_clima <- tabla_fijos %>%
  filter(term %in% c("weatherCloudy", "weatherRain")) %>%
  mutate(
    term = recode(term,
                  "weatherCloudy" = "Cloudy vs Sunny",
                  "weatherRain" = "Rainy vs Sunny"),
    variable_respuesta = recode(variable_respuesta,
                                "males_activity" = "Males activity",
                                "females_foraging" = "Females foraging",
                                "males_following" = "Males following",
                                "males_bumping" = "Males bumping",
                                "m_vs_m" = "Male-male fights"),
    # Calcular IRR (exponenciar)
    IRR = exp(estimate),
    CI_low = exp(conf.low),
    CI_high = exp(conf.high)
  )

# C) Datos para efectos de observador (P3)
# Para el primer modelo (males_activity)
modelo_actual <- lista_modelos[[1]]
ranef_obs <- ranef(modelo_actual)$cond$observer
ranef_obs_df <- data.frame(
  observer = rownames(ranef_obs),
  efecto = ranef_obs[,"(Intercept)"],
  stringsAsFactors = FALSE
)

# Calcular error estándar aproximado (VERSIÓN CORREGIDA para glmmTMB)
ranef_obj <- ranef(modelo_actual, condVar = TRUE)
if (!is.null(attr(ranef_obj$cond$observer, "postVar"))) {
  var_matrix <- attr(ranef_obj$cond$observer, "postVar")
  se_obs <- sqrt(diag(var_matrix)[1])  # Solo el primer elemento (intercept)
} else {
  # Alternativa: usar varianza del modelo
  summ <- summary(modelo_actual)
  var_obs <- summ$varcor$cond$observer[1,1]
  se_obs <- sqrt(var_obs)
}
ranef_obs_df$se <- se_obs

# Calcular IRR en ranef_obs_df ANTES de usarlo en p3
ranef_obs_df$IRR <- exp(ranef_obs_df$efecto)  # <<< MOVER ESTA LÍNEA AQUÍ

# D) Datos para heatmap (P4)
# Crear grid de predicciones para heatmap
grid_data <- expand.grid(
  time_num = seq(7, 13, by = 0.25),  # Cada 15 minutos
  weather = c("Sunny", "Cloudy", "Rain"),
  observer = levels(perdita$observer)[1],
  day = levels(perdita$day)[1]
)

# Predecir para males_activity
grid_data$predicted <- predict(modelo_actual, 
                               newdata = grid_data, 
                               type = "response")

# -----------------------------------------------------
# 2. CREAR LOS 4 SUBGRÁFICOS
# -----------------------------------------------------

# ---- P1: Gráfico temporal (simplificado) ----
p1 <- ggplot(predicciones_todas, 
             aes(x = x, y = predicted, color = comportamiento)) +
  geom_line(size = 0.8) +
  labs(title = "A) Daily activity patterns",
       subtitle = "All behaviors over time",
       x = "Time of day", 
       y = "Predicted count") +
  scale_x_continuous(
    breaks = c(7, 9, 11, 13),
    labels = c("7:00", "9:00", "11:00", "13:00")
  ) +
  scale_color_brewer(palette = "Set1", name = "Behavior") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

# ---- P2: Forest plot de efectos climáticos ----
p2 <- ggplot(coef_clima, 
             aes(x = IRR, y = variable_respuesta, color = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), 
                 height = 0.2, 
                 position = position_dodge(width = 0.5),
                 size = 0.7) +
  labs(title = "B) Weather effects (IRR)",
       subtitle = "Cloudy/Rainy vs Sunny",
       x = "Incidence Rate Ratio (log scale)", 
       y = "Behavior") +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 2)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"), 
                     name = "Comparison") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")

# ---- P3: Efectos de observador ----
p3 <- ggplot(ranef_obs_df, 
             aes(x = reorder(observer, IRR), y = IRR)) +
  geom_bar(stat = "identity", fill = "#4DAF4A", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = exp(efecto - 1.96*se), 
                    ymax = exp(efecto + 1.96*se)), 
                width = 0.2, color = "darkgreen") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_text(aes(label = sprintf("%.2f", IRR)), 
            vjust = -0.5, size = 3) +
  labs(title = "C) Observer random effects",
       subtitle = "Male activity (IRR relative to mean)",
       x = "Observer", 
       y = "IRR (1 = population mean)") +
  theme_minimal(base_size = 10) +
  coord_flip()

# MODIFICA la sección del heatmap (p4) así:

# Crear predicciones para TODOS los comportamientos
grid_data_all <- expand.grid(
  time_num = seq(7, 13, by = 0.25),
  weather = c("Sunny", "Cloudy", "Rain"),
  observer = levels(perdita$observer)[1],
  day = levels(perdita$day)[1],
  comportamiento = nombres_bonitos,  # <<< AÑADIR
  stringsAsFactors = FALSE
)

# Predecir para cada comportamiento
predictions_list <- list()
for (i in 1:5) {
  modelo <- lista_modelos[[i]]
  nombre <- nombres_bonitos[i]
  
  sub_grid <- grid_data_all[grid_data_all$comportamiento == nombre, ]
  sub_grid$predicted <- predict(modelo, newdata = sub_grid, type = "response")
  predictions_list[[i]] <- sub_grid
}

grid_data_all <- do.call(rbind, predictions_list)

# Crear heatmap con facet
p4 <- ggplot(grid_data_all, aes(x = time_num, y = weather, fill = predicted)) +
  geom_tile() +
  facet_wrap(~comportamiento, ncol = 2, scales = "free") +  # <<< FACET aquí
  scale_fill_viridis_c(name = "Predicted\ncount", 
                       option = "C", 
                       direction = -1) +
  labs(title = "D) Heatmap: Time × Weather × Behavior",
       subtitle = "All behavioral predictions",
       x = "Time of day", 
       y = "Weather condition") +
  scale_x_continuous(
    breaks = c(7, 9, 11, 13),
    labels = c("7:00", "9:00", "11:00", "13:00")
  ) +
  theme_minimal(base_size = 9) +  # Reducir tamaño base por los facets
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text = element_text(size = 8, face = "bold"),  # Texto de facets
    panel.spacing = unit(0.5, "lines")  # Espacio entre facets
  )
# -----------------------------------------------------
# 3. COMBINAR EN PANEL 2×2
# -----------------------------------------------------

# Versión 1: Todos del mismo tamaño
panel_final <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Comprehensive analysis of Perdita maritima behavior",
    subtitle = "GLMM results: Temporal patterns, weather effects, and observer variability",
    caption = "Data: Beach observations (9-12 Oct) | Models: glmmTMB nbinom2",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(size = 9, color = "gray40")
    )
  )

# Mostrar panel
print(panel_final)
