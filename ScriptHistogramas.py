import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Ruta al archivo .cct directo, me cuesta trabajo encontrar la ruta más adelante
ruta = r"C:\Users\adyji\Documents\INMEGEN\Datos LinkedOmics\Human__TCGA_THCA__JHU_USC__Methylation__Meth450__01_28_2016__BI__CpG__Firehose_Methylation_Prepocessor.cct"

# Cargar archivo usando el tabulador como separación
datos = pd.read_csv(ruta, sep="\t")

# La primera columna es el nombre del gen, las demás son muestras
muestras = datos.columns[1:]

# Crear los histogramas para cada muestra
for muestra in muestras:
    plt.figure(figsize=(6,4))
    
    # Ajustar bins: aquí de -3 a 3 en incrementos de 0.1
    plt.hist(datos[muestra].dropna(), bins=np.arange(-3, 3.1, 0.1),
             color="palegreen", edgecolor="black")
    
    # Limitar el rango del eje X, luego investigar como alterar lo de los eje Y, pero por lo pronto funciona

    plt.xlim(-3, 3)
    
    plt.xlabel("Valor de metilación")
    plt.ylabel("Frecuencia")
    plt.title(f"Histograma de {muestra}")
    plt.tight_layout()
    plt.show()
    
# Número de filas que en realidad son numero de sitios CpG

num_valores = datos.shape[0]
print(f"Cada muestra tiene {num_valores} valores de metilación")