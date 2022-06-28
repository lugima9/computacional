"//////////////////////////////////////////////////////////////"
# Programa que obtiene el exponente crítico de la magnetización
# para el modelo de Ising a partir de los datos de magnetización
# simulados.
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Función de ajuste
def func(x,a,b):
    return a*x + b

# Función que encuentra el índice del valor más cercano a uno dado en un array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# Lee los ficheros de magnetización, energía y calor específico CON barras de error
# y los grafica para cada N
file = "magnet_def.dat"

beta, u_beta = np.array([]), np.array([])

with open(file, "r") as file:

    count = 0

    # Leo el archivo
    lines=file.readlines()

    # Empiezo con N = 16
    N = 16
    # Vacío variables que almacenan x, y, error de y y dimensiones
    x=[]
    y=[]
    yerror=[]
    dims = []

    # Leo por líneas
    for line_count in range(len(lines)):

        line = lines[line_count]
        # A cada línea (string) le quito el salto de línea final
        line = line.strip("\n")

        # Si no he llegado hasta el final del bloque de este N...
        if line != "":
            # ... guardo datos de temperatura, magnitud y error
            data = line.split(", ")
            x.append(float(data[0]))
            y.append(float(data[1]))
            yerror.append(float(data[2]))
        # Si he llegado al final del bloque de este N...
        else:
            # Calculamos exponente crítico para estas dimensiones

            x = np.array(x)
            y = np.array(y)
            idx = find_nearest(y, 0.5)

            t_crit = x[idx]

            # Tomamos valores a la izquierda del corte con m = 0.5
            valuesTemp = np.abs(t_crit - x[0:idx])
            valuesMagn = y[0:idx]

            valueslogTemp = np.log(valuesTemp)
            valueslogMagn = np.log(valuesMagn)

            # Ajuste
            param, cov = curve_fit(func, valueslogTemp, valueslogMagn)
            a, b = param
            ua, ub = np.sqrt(np.diag(cov))

            beta = np.append(beta, a)
            u_beta = np.append(u_beta, ua)

            print(f"Exponente critico para N = {N}: {a} +/- {ua}")
            
            count += 1

            # Paso al siguiente bloque de N
            N *= 2
            x=[]
            y=[]
            yerror=[]

dims = np.asarray([16,32,64,128]).astype(float)

fig, ax = plt.subplots()

# Graficamos alpha en función de 1/N

x = np.reciprocal(dims)
y = beta
yerror = u_beta

ax.set_xlabel(r"$1/N$")
ax.set_ylabel(r"$\beta$")

# Eliminamos el beta de N = 16

# x = np.delete(x, 0)
# y = np.delete(y, 0)
# yerror = np.delete(yerror, 0)

# Ajuste lineal para extrapolar

param, cov = curve_fit(func, x, y)
a, b = param
ua, ub = np.sqrt(np.diag(cov))

print(f"Exponente critico extrapolado: {b} +/- {ub}")

x_new_value = np.append(x, 0.)
y_new_value = func(x_new_value, a, b)

# Graficamos betas obtenidos y ajuste para extrapolar

ax.errorbar(x, y, yerror, linestyle="", marker="^", color = "#0A369D")
ax.plot(x_new_value,y_new_value,color="#4472CA")

ax.set_xlim(0,)

plt.title(r"Exp. crit. $\beta$ vs $1/N$", fontname = "Arial Rounded MT Bold")

plt.grid()

plt.show()