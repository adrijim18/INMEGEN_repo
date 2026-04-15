# ============================================================
# PCA DE METILACIÓN - DATOS TCGA (LinkedOmics)
# ============================================================
# Proyecto: Detección de cáncer por metilación de cfDNA
# Institución: INMEGEN
# Datos: M-values de arrays Illumina 450K
# Tipos de cáncer: BRCA, LUAD, LUSC, PRAD, THCA
# Todos los datos son muestras tumorales (LinkedOmics)
# Nota: Datos de tejido tumoral primario, no cfDNA.
# El PCA es un paso exploratorio para validar que los perfiles
# de metilación son tejido-específicos, supuesto central del
# proyecto. Respaldado por Moss et al. 2018, Katsman et al. 2022.
# ============================================================

library(ggplot2)
library(data.table)

# --- 1. PARÁMETROS ---
N_CPGS <- 20000
set.seed(42)

# --- 2. LEER SOLO LOS IDs DE CpG DE CADA ARCHIVO ---
# Antes de cargar los datos completos, identificamos qué CpGs
# están presentes en los 5 datasets para poder compararlos
cat("Leyendo IDs de CpGs de cada archivo...\n")

ids_brca <- fread("BRCA.cct", select = 1)[[1]]
ids_luad <- fread("LUAD.cct", select = 1)[[1]]
ids_lusc <- fread("LUSC.cct", select = 1)[[1]]
ids_prad <- fread("PRAD.cct", select = 1)[[1]]
ids_thca <- fread("THCA.cct", select = 1)[[1]]

# --- 3. ENCONTRAR CpGs EN COMÚN ---
# Solo podemos comparar muestras en sitios CpG medidos en todos
cpgs_comunes <- Reduce(intersect, list(
  ids_brca, ids_luad, ids_lusc, ids_prad, ids_thca
))
cat("CpGs en común entre los 5 datasets:", length(cpgs_comunes), "\n")

# Muestreo aleatorio DE LOS CpGs COMUNES
cpgs_seleccionados <- sample(cpgs_comunes, min(N_CPGS, length(cpgs_comunes)))
cat("CpGs seleccionados para el PCA:", length(cpgs_seleccionados), "\n\n")

# Limpiar memoria
rm(ids_brca, ids_luad, ids_lusc, ids_prad, ids_thca)
gc()

# --- 4. FUNCIÓN PARA CARGAR SOLO LOS CpGs SELECCIONADOS ---
cargar_cct <- function(archivo, cpgs) {
  cat("Cargando", archivo, "...\n")
  datos <- fread(archivo, sep = "\t", header = TRUE)
  ids   <- datos[[1]]
  datos <- as.matrix(datos[, -1, with = FALSE])
  rownames(datos) <- ids
  datos <- datos[cpgs, ]
  cat("  Dimensiones:", nrow(datos), "CpGs x", ncol(datos), "muestras\n\n")
  return(datos)
}

# --- 5. CARGAR LOS 5 ARCHIVOS ---
brca <- cargar_cct("BRCA.cct", cpgs_seleccionados)
luad <- cargar_cct("LUAD.cct", cpgs_seleccionados)
lusc <- cargar_cct("LUSC.cct", cpgs_seleccionados)
prad <- cargar_cct("PRAD.cct", cpgs_seleccionados)
thca <- cargar_cct("THCA.cct", cpgs_seleccionados)

# --- 6. UNIR EN UNA SOLA MATRIZ ---
datos <- cbind(brca, luad, lusc, prad, thca)
cat("Matriz final:", nrow(datos), "CpGs x", ncol(datos), "muestras\n\n")

# Guardar etiquetas ANTES de borrar los objetos individuales
tipo_cancer <- c(
  rep("BRCA", ncol(brca)),
  rep("LUAD", ncol(luad)),
  rep("LUSC", ncol(lusc)),
  rep("PRAD", ncol(prad)),
  rep("THCA", ncol(thca))
)

rm(brca, luad, lusc, prad, thca)
gc()

# --- 7. ELIMINAR NAs ---
filas_completas <- apply(datos, 1, function(x) !any(is.na(x)))
cat("CpGs sin ningún NA:", sum(filas_completas), "de", nrow(datos), "\n\n")
datos <- datos[filas_completas, ]

# --- 8. CALCULAR PCA ---
cat("Calculando PCA... (puede tardar unos minutos)\n")
pca_resultado <- prcomp(t(datos), center = TRUE, scale. = FALSE)

varianza <- summary(pca_resultado)$importance[2, 1:5] * 100
cat("Varianza explicada por PC1-PC5:\n")
print(round(varianza, 1))

# --- 9. PREPARAR DATOS PARA GRAFICAR ---
pca_df <- data.frame(
  PC1    = pca_resultado$x[, 1],
  PC2    = pca_resultado$x[, 2],
  cancer = tipo_cancer
)

# --- 10. GRÁFICA: PCA coloreado por tipo de cáncer ---
colores <- c(
  "BRCA" = "#E74C3C",
  "LUAD" = "#3498DB",
  "LUSC" = "#2ECC71",
  "PRAD" = "#9B59B6",
  "THCA" = "#F39C12"
)

g1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cancer)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = colores) +
  labs(
    title    = "PCA de perfiles de metilación por tipo de cáncer",
    subtitle = "Datos TCGA vía LinkedOmics | M-values | Arrays Illumina 450K",
    x        = paste0("PC1 (", round(varianza[1], 1), "% varianza)"),
    y        = paste0("PC2 (", round(varianza[2], 1), "% varianza)"),
    color    = "Tipo de cáncer"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(g1)
ggsave("PCA_por_cancer.png", g1, width = 8, height = 6, dpi = 300)
cat("Guardada: PCA_por_cancer.png\n")
cat("\n¡Análisis completado!\n")
