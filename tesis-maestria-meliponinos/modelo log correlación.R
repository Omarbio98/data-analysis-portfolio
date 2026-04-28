setwd("/Users/user/Desktop/Gráficas y scripts maestria")
corre<- read.csv(choose.files(),na.strings = c("NA",""))
library(lmtest)
library(ggplot2)
library(tidyr)
library(dplyr)
####Comprobamos normalidad en los datos de celdas de cría y volumen de las cavidades###
####We verified normality in the data for brood cells and cavity volume###
tapply(corre$Volumen, corre$Especie, shapiro.test)
tapply(corre$Celdas, corre$Especie, shapiro.test)
###En base a si los datos tuvieron una distribución normal o no, aplicamos diferentes tests para evaluar diferencias entre especies###
###Based on whether the data had a normal distribution or not, we applied different tests to evaluate differences between species###
t.test(Celdas~Especie, data = corre)
wilcox.test(Volumen~Especie, data = corre)


######modelo logaritmico#######
###logarithmic model###
modelo_log <- lm(Celdas ~ log(Volumen), data = corre)
summary(modelo_log)
shapiro.test(residuals(modelo_log))
bptest(modelo_log)

# Extraer coeficientes del modelo log
# Extract coefficients from the log model
coef_log <- coef(modelo_log)
intercept_log <- round(coef_log[1], 2)
slope_log <- round(coef_log[2], 2)  # pendiente para log(Volumen)

# Calcular R² y p-value
r2_log <- round(summary(modelo_log)$r.squared, 3)
p_val_log <- formatC(summary(modelo_log)$coefficients[2,4], format = "e", digits = 2)

# Crear secuencia suave de volúmenes para la curva
# Create a smooth sequence of volumes for the curve
volumen_seq <- seq(min(corre$Volumen), max(corre$Volumen), length.out = 100)

# Predecir valores con intervalo de confianza al 95%
# Predict values with a 95% confidence interval
# interval = "confidence" da el IC para la media (la curva)
predicciones <- predict(modelo_log, 
                        newdata = data.frame(Volumen = volumen_seq), 
                        interval = "confidence", 
                        level = 0.95)

# Convertir a data frame para ggplot
curva_log <- data.frame(
  Volumen = volumen_seq,
  Celdas = predicciones[, 1],  # valor ajustado (fit)
  LI = predicciones[, 2],       # límite inferior (lwr)
  LS = predicciones[, 3]        # límite superior (upr)
)

# GRÁFICO CON BANDA DE CONFIANZA
# GRAPH WITH CONFIDENCE BAND
ggplot(corre, aes(x = Volumen, y = Celdas, color = Especie)) +
  
  # Banda de confianza (primero para que quede detrás de los puntos)
  # Trust band (first so that it is behind the dots)
  geom_ribbon(data = curva_log, 
              aes(x = Volumen, ymin = LI, ymax = LS), 
              fill = "#74a9cf", alpha = 0.3, inherit.aes = FALSE) +
  
  # Puntos originales
  # Original dots
  geom_point(size = 3, alpha = 0.8) +
  
  # Curva del modelo logarítmico
  # Logarithmic model curve
  geom_line(data = curva_log, aes(x = Volumen, y = Celdas), 
            color = "#2c3e50", size = 1.5, inherit.aes = FALSE) +
  
  # Línea horizontal en cero (opcional, para referencia)
  # Horizontal line at zero (optional, for reference)
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.5) +
  
  # Etiquetas y títulos
  # Labels ant titles
  labs(x = "Volumen de cavidad (L)", 
       y = "Número de celdas de cría",
       title = "Relación entre volumen de cavidad y celdas de cría",
       subtitle = paste0("Modelo logarítmico con IC 95%",
                         "\nCeldas = ", intercept_log, " + ", slope_log, " × log(Volumen)",
                         "\nR² = ", r2_log, " | p = ", p_val_log),
       color = "Especie",
       caption = "El área sombreada representa el intervalo de confianza al 95% para la media") +
  
  
  # Escala de colores para especies
 # rojo y azul daltónicos-amigables
  scale_color_manual(values = c("#e41a1c", "#377eb8"),
                     labels= c(expression(italic("Scaptotrigona pectoralis")),
                               expression(italic("Plebeia fulvopilosa"))))+
  # Tema mejorado
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30", size = 10),
        plot.caption = element_text(hjust = 0.5, color = "gray60", size = 8, face = "italic"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "gray30"),
        axis.title = element_text(face = "bold"))



