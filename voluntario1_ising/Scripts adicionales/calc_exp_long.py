"//////////////////////////////////////////////////////////////"
# Programa que obtiene el exponente crítico de la longitud de
# correlación para el modelo de Ising a partir de los datos de
# correlación simulados.
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit

# Función de ajuste
def func(x,a,b):
    return a*x + b

# Función que encuentra el índice del valor más cercano a uno dado en un array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# Lee el fichero de función de correlación SIN barras de error
# y la grafica para cada N

fileName = "correl_def.dat"
tmin = 2.1666666666666665
tmax = 2.611111111111111
tnum = 10

t_crit = 2.2353730544565877 # Valor obtenido en nuestra simulación 

with open(fileName, "r") as file:

    valuesLen = np.zeros(10-5)

    temp=np.linspace(tmin, tmax, tnum)
    temp_normalized = (temp-tmin)/(tmax-tmin)
    colors = [ cm.jet(t) for t in temp_normalized ]
    data=file.readlines()
    stop_old = 0

    for N in [16, 32, 64, 128]:

        # Leo el archivo
        stop_new = stop_old + (N+2)*10 + 1 # (N valores + 1 temp) x 10 temps + 1 salto de línea
        lines = data[stop_old:stop_new]

        # Elimino los valores de temperatura de los datos leídos
        del lines[::N+2]

        # Vacío variables que almacenan x e y
        temp_count = 0
        dist = 0
        x = []
        y = []

        # Leo por líneas
        for line_count in range(len(lines)):

            line = lines[line_count]
            # A cada línea (string) le quito el salto de línea final
            line = line.strip("\n")

            # Si no he llegado hasta el final del bloque de esta temp
            if line != "":
                # ... guardo datos
                x.append(dist)
                y.append(float(line))

                dist += 1
            # Si he llegado al final del bloque de esta temp...
            else:

                if N == 128: # Solo para 128 porque el resto no suelen llegar a 1/e
                    x,y = np.array(x), np.array(y)
                    x,y = x[x <= N//2], y[x <= N//2]
                    
                    if temp_count >= 3:
                        idx = find_nearest(y, 1/math.e)
                        valuesLen[temp_count-3] = x[idx]

                        if temp_count == 7: # Máxima temp en este rango
                            valuesTemp = temp[3:8] - t_crit

                            valueslogTemp = np.log(np.abs(valuesTemp))
                            valueslogLen = np.log(valuesLen)

                            param, cov = curve_fit(func, valueslogTemp, valueslogLen)
                            a, b = param
                            ua, ub = np.sqrt(np.diag(cov))

                            nu = a
                            u_nu = ua

                            print(f"Exponente critico: {nu} +/- {u_nu}")

                            x = valueslogTemp
                            y = valueslogLen

                            fig, ax = plt.subplots()
                            ax.set_xlabel(r"$\log{(|T-T_c|)}$")
                            ax.set_ylabel(r"$\log{\xi}$")

                            y_new_value = func(x, a, b)

                            ax.scatter(x, y, marker="^", color = "#0A369D")
                            ax.plot(x,y_new_value,color="#4472CA")

                            ax.set_xlim(-2.6, -1.2)
                            ax.set_ylim(0.6,)

                            plt.title(r"Exp. crit. $\nu$ para $N=128$", fontname = "Arial Rounded MT Bold")

                            plt.grid()

                            plt.show()

                            break

                # Paso al siguiente bloque de temp
                temp_count += 1
                dist = 0
                x = []
                y = []

        stop_old = stop_new