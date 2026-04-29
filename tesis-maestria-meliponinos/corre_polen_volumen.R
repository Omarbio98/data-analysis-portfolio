setwd("/Users/user/Desktop/Gráficas y scripts maestria")
corre<- read.csv(choose.files(),na.strings = c("NA",""))
library(lmtest)
library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
corre$Especie <- as.factor(corre$Especie)
### En este caso, varios modelos de regresión fueron aplicados para mostrar como el volumen de la cavidad afecta el total de polen y y que curva describe mejor la relación entre las 2 variables### 
###In this case, several regression models were applied to show how cavity volume affects total pollen and which curve best describes the relationship between the two variables###
modelo_lm <- lm(Polen~Cavidad, data = corre)
### En los modelos log(Polen), se usó +1 debido a que en algunos nidos se encontró 0 polen###
### In the log(Pollen) models, +1 was used because 0 pollen was found in some nests###
modelo_pot_polen <- lm(log(Polen+1)~log(Cavidad),data = corre)

modelo_exp_polen <- lm(log(Polen+1)~ Cavidad, data=corre)

modelo_log <- lm(Polen ~ log(Cavidad), data = corre)

modelo_mixto <- lmer(log(Polen+1) ~ log(Cavidad) + (1 | Especie), 
                     data = corre)
###Comparé el AIC de todos los modelos para elegir el que muestra un mejor ajuste###
###I compared the AIC of all models to choose the one that shows the best fit###
AIC(modelo_lm, modelo_pot_polen, modelo_exp_polen, modelo_log, modelo_mixto)

###Finalmente, a comprobar los supuestos del modelo con el mejor ajuste###
###Finally, to verify the assumptions of the model with the best fit###
summary(modelo_pot_polen)
shapiro.test(residuals(modelo_pot_polen))
bptest(modelo_pot_polen)
### Sin embargo, los supuestos mostraron heterocedasticidad y falta de normalidad en los residuales del modelo con mejor ajuste###
### However, the assumptions showed heteroscedasticity and a lack of normality in the residuals of the best-fitting model###

# Esto se debe a los 0 en algunos nidos###
# This is due the 0 in some nests###
sum(corre$Polen == 0, na.rm = TRUE)

###Por ende, los nidos con polen=0 fueron excluídos###
### Therefore, the nests with pollen=0 were excluded###

corre_sin0 <- subset(corre, Polen>0)
modelo_pot_polen_positivo <- lm(log(Polen)~log(Cavidad),data = corre_sin0)
summary(modelo_pot_polen_positivo)

shapiro.test(residuals(modelo_pot_polen_positivo))
bptest(modelo_pot_polen_positivo)
confint(modelo_pot_polen_positivo, parm = "log(Cavidad)")



# 2. EXTRAER COEFICIENTES DEL MODELO
coef_pot <- coef(modelo_pot_polen_positivo)
intercept_pot <- round(coef_pot[1], 2)
slope_pot <- round(coef_pot[2], 2)  # pendiente para log(Volumen) = exponente b

# 3. CALCULAR R² Y P-VALUE
r2_pot <- round(summary(modelo_pot_polen_positivo)$r.squared, 3)
p_val_pot <- formatC(summary(modelo_pot_polen_positivo)$coefficients[2,4], format = "e", digits = 2)

# 4. CREAR SECUENCIA DE VOLÚMENES PARA LA CURVA
volumen_seq <- seq(min(corre$Cavidad), max(corre$Cavidad), length.out = 100)

# 5. PREDECIR VALORES CON INTERVALO DE CONFIANZA AL 95%
predicciones <- predict(modelo_pot_polen_positivo, 
                        newdata = data.frame(Cavidad = volumen_seq), 
                        interval = "confidence", 
                        level = 0.95)

# 6. CONVERTIR PREDICCIONES A ESCALA ORIGINAL (deshacer log)
# IMPORTANTE: exp() invierte el logaritmo
curva_pot <- data.frame(
  Cavidad = volumen_seq,
  Polen = exp(predicciones[, 1]),        # valor ajustado (fit) en escala original
  LI = exp(predicciones[, 2]),           # límite inferior en escala original
  LS = exp(predicciones[, 3])            # límite superior en escala original
)

# 7. GRÁFICO
ggplot(corre, aes(x = Cavidad, y = Polen, color = Especie)) +
  
  # Banda de confianza (primero para que quede detrás)
  geom_ribbon(data = curva_pot, 
              aes(x = Cavidad, ymin = LI, ymax = LS), 
              fill = "#2ecc71", alpha = 0.2, inherit.aes = FALSE) +
  
  # Puntos originales
  geom_point(size = 3, alpha = 0.8) +
  
  # Curva del modelo potencial
  geom_line(data = curva_pot, aes(x = Cavidad, y = Polen), 
            color = "#27ae60", size = 1.5, inherit.aes = FALSE) +
  
  # Etiquetas y títulos
  labs(x = "Volumen de cavidad (L)", 
       y = "Polen almacenado (gr)",
       title = "Relación entre volumen de cavidad y almacenamiento de polen",
       subtitle = paste0("Modelo potencial: Polen = exp(", intercept_pot, " + ", slope_pot, " × log(Volumen))",
                         "\nR² = ", r2_pot, " | p = ", p_val_pot,
                         "\nExponente b = ", slope_pot),
       color = "Especie",
       caption = "El área sombreada representa el IC 95% para la media") +
  
  # Escala de colores para especies (ajusta según tengas)
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a"),
                     labels= c(expression(italic("Scaptotrigona pectoralis")),
                               expression(italic("Plebeia fulvopilosa")),
                               expression(italic("Cephalotrigona zexmeniae"))))+ # rojo, azul, verde
  
  # Tema
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30", size = 10),
        plot.caption = element_text(hjust = 0.5, color = "gray60", size = 8))

