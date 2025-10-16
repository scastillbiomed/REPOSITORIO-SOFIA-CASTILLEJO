# SofíaCastillejoSánchez_Trabajo2.R
# Trabajo final Bioinformática - Curso 25/26
# Análisis de parámetros biomédicos por tratamiento

# 1. Cargar librerías (si necesarias) y datos del archivo "datos_biomed.csv". (0.5 pts)

datos = read.csv("datos_biomed.csv")
# Instalar librerías
if(!require("ggplot2")) install.packages("ggplot2",dependencies=TRUE)
if(!require("dplyr")) install.packages("dplyr",dependencies=TRUE)

# Cargar librerías
library(ggplot2)
library(dplyr)

# 2. Exploración inicial con las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos? (0.5 pts)

#La función head() enseña las primeras filas
head(datos)
#La función summary() devuelve un resumen estadístico 
summary(datos)
#La función dim() muestra las dimensiones de la tabla (filas y columnas)
dim(datos)
#La función str() indica el tipo de datos que tienen las columnas 
str(datos)

#Hay 5 variables ("ID", "Tratamiento", "Glucosa", "Presion", "Colesterol") y 3 tratamientos (FarmacoA, FarmacoB, Placebo)

# 3. Una gráfica que incluya todos los boxplots por tratamiento. (1 pt)

#cada boxplot representa los niveles de glucosa (eje y) en base al tratamiento (eje x).

ggplot(datos, aes(x=Tratamiento, y=Glucosa, fill=Tratamiento)) +
geom_boxplot(alpha=0.7) + 
labs(title="NIVELES DE GLUCOSA POR TRATAMIENTO", 
x= "Tratamiento", 
y="Glucosa (mg/dL)")+
theme_minimal(base_size=14)

# 4. Realiza un violin plot (investiga qué es). (1 pt)
#El violin plot es un tipo de gráfica que permite ver la distribución de los datos y su densidad de población. 
#Para realizar un violin plot, se siguen los mismos pasos que para el BoxPlot pero usaando la función geom_violin. 

ggplot(datos, aes(x=Tratamiento, y=Glucosa, fill=Tratamiento)) +
geom_violin(alpha=0.7) + 
labs(title="NIVELES DE GLUCOSA POR TRATAMIENTO", 
x= "Tratamiento", 
y="Glucosa (mg/dL)")+
theme_minimal(base_size=14)
names(datos)

# 5. Realiza un gráfico de dispersión "Glucosa vs Presión". Emplea legend() para incluir una leyenda en la parte inferior derecha. (1 pt)
datos$Tratamiento = as.factor(datos$Tratamiento)
plot(datos$Glucosa, datos$Presion, col=datos$Tratamiento, pch=19, 
main="Relación entre glucosa y presión arterial", xlab="Glucosa(mg/dl)", ylab="Presión (mmHg)")
legend("bottomright", legend =levels(datos$Tratamiento), 
col=1:length(levels(datos$Tratamiento)), pch=19, title = "Tratamiento")


# 6. Realiza un facet Grid (investiga qué es): Colesterol vs Presión por tratamiento. (1 pt)
#Un face Grid permite hacer y visualizar varias gráficas a la vez según la variable. Para ello es necesario introducir la función facet_grid.
ggplot(datos, aes(x=Tratamiento, y=Glucosa, fill=Tratamiento)) +
geom_point(size=3, alpha=0.7) + 
facet_grid(.~Tratamiento)+
labs(title="COLESTEROL VS PRESIÓN POR TRATAMIENTO", 
x= "Colesterol(mg/dL)", 
y="Presión Arterial (mmHg)")+
theme_minimal(base_size=14)


# 7. Realiza un histogramas para cada variable. (0.5 pts)
#Hago un histograma para cada variable estudiada. Con hist se crea la gráfica, con $ se indica la columna de la variable. 

#Histograma de glucosa
hist(datos$Glucosa, main = "Histograma de Glucosa",
     col = "orange", xlab = "Glucosa (mg/dL)")

#Histograma de presión arterial
hist(datos$Presion, main = "Histograma de Presion Arterial",
     col = "lightgreen", xlab = "Presión (mmHg)")

#Histograma de colesterol
hist(datos$Colesterol, main = "Histograma de Colesterol",
     col = "blue", xlab = "Colesterol (mg/dL)")


# 8. Crea un factor a partir del tratamiento. Investiga factor(). (1 pt)
#El código factor sirve para ordenar los datos según los grupos, en este caso, de tratamiento. 
datos$Tratamiento <- factor(datos$Tratamiento)
datos$Tratamiento

#9. Obtén la media y desviación estándar de los niveles de glucosa por tratamiento. Emplea aggregate() o apply(). (0.5pts)

aggregate (Glucosa~Tratamiento, data=datos, FUN=mean)
aggregate (Glucosa~Tratamiento, data=datos, FUN=sd)

# 10. Extrae los datos para cada tratamiento y almacenalos en una variable. Ejemplo todos los datos de Placebo en una variable llamada placebo. (1 pt)

#Guardamos cada grupo en variables diferentes, cogiendo las filas que contengan ese tratamiento. 
farmacoA = datos[datos$Tratamiento == "FarmacoA",]
farmacoB = datos[datos$Tratamiento == "FarmacoB",]
placebo = datos [datos$Tratamiento == "Placebo",]

# 11. Evalúa si los datos siguen una distribución normal y realiza una comparativa de medias acorde. (1 pt)

#Los datos siguen una distribución normal si p>0.05. Para comparar las medias, si siguen distribución normal se hace ANOVA. Si no, Kruskal-Wallis. 
#Según los datos de ANOVA, si pr>0.05, no hay diferencias significativas entre los diferentes tratamientos. 

#GLUCOSA
shapiro_glu = shapiro.test(datos$Glucosa)
shapiro_glu
if (shapiro_glu$p.value>0.05) {
anova_glu = aov(Glucosa~Tratamiento, datos)
summary(anova_glu)} else {kruscal.test(Glucosa~Tratamiento, datos)}
#Los datos de glucosa siguen una distribución normal (p-value = 0.91). 
#No hay diferencias significativas entre los tratamientos (pr = 0.1).

#PRESIÓN
shapiro_pre = shapiro.test(datos$Presion)
shapiro_pre
if (shapiro_pre$p.value>0.05) {
anova_pre = aov(Presion~Tratamiento, datos)
summary(anova_pre)} else {kruscal.test(Presion~Tratamiento, datos)}
#Los datos de presión sí siguen una distribución normal (p-value = 0.6399). 
#Sí hay diferencias significativas entre los tratamientos (pr = 1.2e-05).

#COLESTEROL
shapiro_col = shapiro.test(datos$Colesterol)
shapiro_col
if (shapiro_col$p.value>0.05) {
anova_col = aov(Colesterol~Tratamiento, datos)
summary(anova_col)} else {kruscal.test(Colesterol~Tratamiento, datos)}
#Los datos de presión sí siguen una distribución normal (p-value = 0.9462). 
#Sí hay diferencias significativas entre los tratamientos (pr = 0.00429).

# 12. Realiza un ANOVA sobre la glucosa para cada tratamiento. (1 pt)
#He hecho manualmente un archivo de texto con los datos para poder aplicar el anova. 
datos_anova = read.table("datos_r.txt", sep="", header=TRUE)
anova = aov(Glucosa~Tratamiento, data = datos_anova)
summary(anova)



