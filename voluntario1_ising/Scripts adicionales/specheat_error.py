"//////////////////////////////////////////////////////////////"
# Programa que calcula el error en el calor específico obtenido en
# en el modelo de Ising con un método Monte Carlo.      
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

import numpy as np
from numba import njit # Compilador de Python que mejora que optimiza la ejecución del código

"FUNCIONES"

# Función que calcula el calor específico a partir de una distribución gaussiana de <E> y <E^2>
# con medias y desviaciones estándar dadas
# Recibe: <E>, std(E), <E^2>, std(E^2), dimensiones red, temperatura
# Devuelve: <C>, 2sigma(C)  
@njit
def gaussianSpecHeat(e, e_err, e2, e2_err, N, temp):

    iters = 10000

    e_sample = np.random.normal(e, 2*e_err/(iters**0.5), int(1e6))
    e2_sample = np.random.normal(e2, 2*e2_err/(iters**0.5), int(1e6))

    specHeat_sample = (e2_sample - np.square(e_sample))/(N*N*temp)

    specHeat = np.mean(specHeat_sample)
    specHeat_err = 2*np.std(specHeat_sample)/(iters**0.5)

    return specHeat, specHeat_err


# Función que obtiene los datos del fichero de energías del script de simulación del modelo de Ising.
# Recibe: nombre del archivo
# Devuelve: nº de datos, temp, dimensiones red, <E>, std(E), <E^2>, std(E^2)
def gatherData(fileName):

    # Extraemos los datos
    file = open(fileName, "r")
    datos_raw = file.read().replace(" ", "").split("\n")
    datos = [datos_raw[i] for i in range(len(datos_raw)) if datos_raw[i] != ""]
    file.close()

    # Creamos arrays para devolver
    temp, N, e, e_err, e2, e2_err = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

    N = np.append(N, np.full(10, 16))
    N = np.append(N, np.full(10, 32))
    N = np.append(N, np.full(10, 64))
    N = np.append(N, np.full(10, 128))

    numeroDatos = len(datos)

    # Asignamos valores a los arrays
    for i in range(numeroDatos):
        d = datos[i].split(",")
        temp = np.append(temp, float(d[0]))
        e = np.append(e, float(d[1]))
        e_err = np.append(e_err, float(d[2]))
        e2 = np.append(e2, float(d[3]))
        e2_err = np.append(e2_err, float(d[4]))

    return numeroDatos, temp, N, e, e_err, e2, e2_err


"FUNCIÓN MAIN"

# Obtenemos los datos
fileName = "specheat_wip.dat"
numeroDatos, temp, N, e, e_err, e2, e2_err = gatherData(fileName)

# Calculamos el calor específico
specHeat, specHeat_err = np.zeros(numeroDatos), np.zeros(numeroDatos)

for i in range(numeroDatos):
    specHeat[i], specHeat_err[i] = gaussianSpecHeat(e[i], e_err[i], e2[i], e2_err[i], N[i], temp[i])

# Escribimos los resultados en un fichero
fSpecHeat = open("specheat_def.dat", "w")

for j in range(numeroDatos):
    fSpecHeat.write(f"{temp[j]}, {specHeat[j]}, {specHeat_err[j]}\n")
    if (j+1)%10==0:
        fSpecHeat.write("\n")

fSpecHeat.close()