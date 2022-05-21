"//////////////////////////////////////////////////////////////"
# Programa que simular el Modelo de Ising para medio magnéticos
# con el método Monte Carlo para distintas temperaturas.     
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

import numpy as np
import random as rd
import math
from numba import njit


"-------------------------------------FUNCIONES------------------------------------------"

# Función que genera una configuración inicial de espines en una matriz cuadrada
# Recibe: configuración elegida [int], tamaño fila [int]
# Devuelve: matriz de espines
def iniSim(choice, N):

    if choice == 1:
        # Genero configuración inicial de espines en una matriz de N x N
        spins = np.random.choice([-1,1],size=N*N).reshape(N, N)
    elif choice == 2:
        spins = -np.ones((N, N))    # Todos abajo
    else:
        spins = np.ones((N, N))     # Todos arriba

    print(f"Configuracion inicial: \n {spins}")

    return spins

# Función que simula el modelo de Ising (ejercicio obligatorio) y escribe los
# los resultados en un fichero en formato adecuada para el script de animación
# Recibe: temperatura [float], tamaño fila [int], matriz de espines [numpy array]
#         numero de pasos Monte Carlo [int]
# Devuelve: array de matriz de spins tras cada paso MC
@njit
def simulation(temp, N, spins, stepsMC):

    bigSpins = np.empty((stepsMC+1, N, N))

    # Guardo la config. inicial
    bigSpins[0] = spins

    # Damos varios paso de Monte Carlos (intentamos cambiar NxN espines)
    N2 = N*N

    for step in range(stepsMC):

        for i in range(N2):

            # Escogemos un elem aleat de la matriz de espines
            n = rd.randint(0, N-1)
            m = rd.randint(0, N-1)

            # Calculamos la energía teniendo en cuenta condiciones de contorno
            if n != 0 and n != N-1 and m != 0 and m != N-1:
                deltaE = 2*spins[n,m]*(spins[n+1, m] + spins[n-1, m] + spins[n, m+1] + spins[n, m-1])
            else:
                deltaE = 2*spins[n,m]*(spins[(n+1)%N, m] + spins[(n-1)%N, m] + spins[n, (m+1)%N] + spins[n, (m-1)%N])

            # Generamos número aleat unif entre 0 y 1
            # Como uniform genera [0,1), forzamos que se coja también el 1
            aleatUnif = np.random.uniform(0., 1.)

            # Comparamos con el mínimo entre 1 y exp(-deltaE/temp)
            # Si es menor, cambiamos el signo del espín
            exponential = math.exp(-float(deltaE)/temp)
            p = min(1., exponential)

            if aleatUnif < p:
                spins[n, m] = -spins[n, m]

        bigSpins[step+1] = spins

    return bigSpins


# Función que escribe en un fichero .dat el array de matrices de espín para cada paso MC
# Recibe: array de matrices para cada paso MC, pasos MC
# Devuelve: -
def writeToFile(bigSpins, stepsMC):
    file = open("ising_data.dat", "w")

    for i in range(stepsMC+1):
        np.savetxt(file, bigSpins[i].astype(int), fmt = "%s", delimiter=", ")        
        file.write("\n")
    
    file.close()

    print("El archivo ising_data.dat se ha creado con exito en la localizacion del script.")


"-------------------------------------MAIN------------------------------------------"

# Generamos la configuración inicial

# Fijamos la temperatura y el número de espines del imán 2D
temp = float(input("Temperatura a la que simular el sistema: \n"))
N = int(input("Numero de elementos en una fila de la matriz 2D de espines: \n"))

# Fijamos configuración inicial
choice = 0
while choice != 1 and choice != 2 and choice != 3:
    print("Opciones de config. inicial: Aleatorio (1), todos abajo (2), todos arriba (3)")
    choice = int(input("Elige una configuracion: \n"))

spins = iniSim(choice, N)

# Realizamos la simulación

stepsMC = 0
while stepsMC < 1:
    stepsMC = int(input("Introduzca el numero de pasos Monte Carlo que desea realizar: \n"))

bigSpins = simulation(temp, N, spins, stepsMC)

# Escribimos los resultados en un fichero
writeToFile(bigSpins, stepsMC)







    



