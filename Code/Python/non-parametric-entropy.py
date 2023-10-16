# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 17:49:46 2023

@author: luiso
"""

import numpy as np
import matplotlib.pyplot as plt
#from scipy import stats
import math
#import numpy as np
import scipy.special

# Parámetros de la distribución gamma SAR
mu = 2.0
L = 3.0

# Tamaño de la muestra
n = 5

# Número de muestras a generar
r = 2 # Puedes cambiar este valor según la cantidad de muestras que desees generar

def calcular_entropia(L, mu):
    # Calcular los términos necesarios
    termino1 = L 
    termino2 = np.log(scipy.special.gamma(L))
    termino3 = (1 - L) * scipy.special.psi(L)
    termino4 = np.log(L / mu)
    
    # Calcular la entropía H
    entropiaH = termino1 + termino2 + termino3 - termino4 
    
    return entropiaH

def calcular_entropia2(L, mu):
    # Calcular los términos necesarios
    termino1 = L 
    termino2 = np.log(scipy.special.gamma(L))
    termino3 = (1 - L) * scipy.special.psi(L)
    termino4 = np.log(mu / L)
    
    # Calcular la entropía H
    entropiaH = termino1 + termino2 + termino3 + termino4 
    
    return entropiaH


# Función para generar una muestra de la distribución gamma SAR
def gamma_sar_sample(L, mu, n):
    samples = np.random.gamma(shape=L, scale=mu / L, size=n)
    return samples

# Función para calcular la entropía no paramétrica utilizando el estimador Van Es
def van_es_estimator(data):
    n = len(data)
    m = int(np.sqrt(n) + 0.5)  # m-spacing
    data_sorted = sorted(data)
    sum_term1 = 0
    sum_term2 = 0

    for i in range(1, n - m):
        sum_term1 += math.log(((n + 1) / m) * (data_sorted[i + m] - data_sorted[i]))

    for k in range(m, n + 1):
        sum_term2 += 1 / k

    sum_term3 = math.log(m / (n + 1))

    return (sum_term1 / (n - m)) + sum_term2+ sum_term3

# Generar r muestras y calcular la entropía para cada muestra
muestras = []
entropias = []
for _ in range(r):
    muestra = gamma_sar_sample(L, mu, n)
    entropia = van_es_estimator(muestra)
    entropias.append(entropia)
    muestras.append(muestra)
    

# Visualizar los valores de entropía en un histograma
plt.hist(entropias, bins=30, density=True, alpha=0.6, color='b')
plt.title('Histograma de Entropías No Paramétricas')
plt.xlabel('Entropía')
plt.ylabel('Densidad de probabilidad')
plt.grid(True)
plt.show()


entropia_resultante = calcular_entropia(L, mu)
print("Entropía analitica:", entropia_resultante)

entropia_resultante = calcular_entropia2(L, mu)
print("Entropía analitica q:", entropia_resultante)
print("muestras:", muestra)

print("muestras2:",  muestras)