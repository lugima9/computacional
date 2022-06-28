"//////////////////////////////////////////////////////////////"
# Programa que obtiene la temperatura crítica del modelo de Ising
# a partir de calor específico simulados.
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Lee el fichero de calor específico con barras de error
# y los grafica para cada N
file = "specheat_def.dat"
ylabels = "Calor esp. prom."

with open(file, "r") as file:

    # Leo el archivo
    lines=file.readlines()

    # Empiezo con N = 16
    N = 16

    # Vacío variables que almacenan x, y, error de y
    x,y = [],[]
    dims, max_specHeat, max_temp = [],[],[]

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
        # Si he llegado al final del bloque de este N...
        else:
            # Determinamos el máximo
            dims.append(N)
            max_specHeat.append(max(y))
            posmax = np.argmax(y)
            max_temp.append(x[posmax])

            # Paso al siguiente bloque de N
            N *= 2
            x=[]
            y=[]


dims = np.asarray(dims).astype(float)
max_temp = np.asarray(max_temp)

fig, ax = plt.subplots()

x = np.reciprocal(dims)
y = max_temp

ax.set_xlabel(r"$1/N$")
ax.set_ylabel(r"$T_c ~[J/k_B]$")

# Función de ajuste
def func(x,a,b):
    return a*x + b

# Ajuste
param, cov = curve_fit(func, x, y)
a, b= param
ua, ub = np.sqrt(np.diag(cov))

for i in range(len(dims)):
    print(f"Temperatura crítica para N = {dims[i]}: {max_temp[i]}")

print(f"Temperatura crítica extrapolada: {b} +/- {ub}")

x_new_value = np.append(x, 0.)
y_new_value = func(x_new_value, a, b)

# Graficamos la T_c en función de 1/N y la regresión lineal para extrapolar
ax.scatter(x, y, marker="^", color = "#0A369D")
ax.plot(x_new_value,y_new_value,color="#4472CA")

ax.set_xlim(0,)
ax.set_ylim(2.2,2.6)

plt.suptitle(r"Temperatura crítica $T_c$ vs $1/N$", fontname = "Arial Rounded MT Bold")

plt.grid()

plt.show()