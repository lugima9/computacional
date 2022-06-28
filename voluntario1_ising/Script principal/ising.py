"#//////////////////////////////////////////////////////////////"
# Programa que simula el Modelo de Ising para medios magnéticos
# mediante Monte Carlo para distintas temperaturas y tamaños de red.
# Autor: Luis Gil Martín
# Universidad de Granada
"#//////////////////////////////////////////////////////////////"

import numpy as np
import math
from numba import njit # Compilador de Python que mejora el rendimiento en la ejecución del código
import time

"FUNCIONES"

# Función que genera una configuración inicial de espines en una matriz cuadrada
# Recibe: configuración elegida [int], tamaño fila[int]
# Devuelve: matriz de espines
def iniSim(choice, N):

    N = int(N)

    if choice == 1:
        # Genero configuración inicial de espines en una matriz de N x N
        spins = np.random.choice(np.array([-1,1]),size=N*N).reshape(N, N)
    elif choice == 2:
        spins = -np.ones((N, N))    # Todos abajo
    else:
        spins = np.ones((N, N))     # Todos arriba

    return spins.astype(np.int64)


# Función que cambia N^2 espines de la configuración de espines
# Recibe: configuración de espines, N^2
# Devuelve: configuración de espines actualizada
@njit
def changeSpins(spins, N, temp):

    N2 = N*N
        
    for i in range(N2):
        # Escogemos un elem aleat de la matriz de espines efectiva
        n = np.random.randint(0, N-1)
        m = np.random.randint(0, N-1)

        # Calculamos la energía teniendo en cuenta condiciones de contorno
        if n != 0 and n != N-1 and m != 0 and m != N-1:
            deltaE = 2.*spins[n,m]*(spins[n+1, m] + spins[n-1, m] + spins[n, m+1] + spins[n, m-1])
        else:
            deltaE = 2.*spins[n,m]*(spins[(n+1)%N, m] + spins[(n-1)%N, m] + spins[n, (m+1)%N] + spins[n, (m-1)%N])

        # Generamos número aleat unif entre 0 y 1
        aleatUnif = np.random.uniform(0., 1.)

        # Comparamos con el mínimo entre 1 y exp(-deltaE/temp)
        # Si es menor, cambiamos el signo del espín
        exponential = math.exp(-deltaE/temp)
        p = min(1., exponential)

        if aleatUnif < p:
            spins[n, m] = -spins[n, m]

    return spins


# Función que simula la red de espines con N y temperatura dados.
# Recibe: temperatura (unidades de k_B), N, config inicial, pasos MC
# Devuelve: (<m>, 2sigma(m)), (<e>, 2sigma(e)), (<E>, std(E), <E^2>, std(E^2)),
#            array de <correl>
@njit
def simulationExtra(temp, N, spins, stepsMC):

    N = int(N)
    stepsMC = int(stepsMC)

    # Calculamos los siguientes parámetros, que aparecen varias veces
    N2 = N*N
    iters = stepsMC//100
    iters_root = iters**0.5

    # Creamos los arrays para los valores cada paso MC
    magnet, energy, energySq = np.zeros(iters), np.zeros(iters), np.zeros(iters)
    correlArr = np.zeros((iters, N))

    # Creamos el array de promedios de función de correlación
    correlArrAvg = np.zeros(N)

    k = 0
    for step in range(stepsMC):

        # Cada 100 pMC se calculan las magnitudes
        if (step+1)%100 == 0:

            # Calcula la suma en valor absoluto de los espines
            magnet[k] = abs(np.sum(spins))

            # Calcula la energía de los espines y la función de correlación
            for i in range(N):
                for j in range(N):
                    # Se tienen en cuenta condiciones de contorno
                    if i != 0 and i != N-1 and j != 0 and j != N-1:
                        energy[k] += -0.5*spins[i,j]*(spins[i+1, j] + spins[i-1, j] + spins[i, j+1] + spins[i, j-1])
                    else:
                        energy[k] += -0.5*spins[i,j]*(spins[(i+1)%N, j] + spins[(i-1)%N, j] + spins[i, (j+1)%N] + spins[i, (j-1)%N])
                    
                    # *Es lo mismo hacer la suma del promedio que el promedio de la suma
                    for dist in range(N):       
                        correlArr[k, dist] += spins[i,j]*spins[(i+dist)%N, j]

            # Calcula el cuadrado de la energía de los espines 
            energySq[k] = energy[k]*energy[k]

            k += 1

        # Se cambian N^2 espines aleatorios
        spins = changeSpins(spins, N, temp)

    #----RESULTADOS PROMEDIO----

    # Divide todo los elementos por el número de medidas (para que sean promedios)
    # Multiplica por los factores necesarios para calcular las magnitudes de interés

    magnet /= N2 

    magnetAvg = np.sum(magnet)/iters                                   # Magnetización promedio ->           Resultado

    energyAvg = np.sum(energy)/iters                                   # Energía media
    energySqAvg = np.sum(energySq)/iters                               # Cuadrado de energía medio
    eAvg = energyAvg/(2*N2)                                            # Energía media e ->                  Resultado

    specHeatAvg = (energySqAvg - energyAvg*energyAvg)/(N2*temp)        # Calor específico ->                 Resultado

    correlArrAvg = np.sum(correlArr, axis=0)/(N2*iters)                # Función de correlación (array) ->   Resultado


    #----DESVIACIÓN TÍPICA----

    magnetStd, energyStd, energySqStd = 0., 0., 0.

    # Calculo de forma cumulativa la varianza*N de cada variable
    for k in range(len(magnet)):
        magnetStd += (magnet[k] - magnetAvg)*(magnet[k] - magnetAvg)
        energyStd += (energy[k] - energyAvg)*(energy[k] - energyAvg)
        energySqStd += (energySq[k] - energySqAvg)*(energySq[k] - energySqAvg)
    
    # Calculo la desviación típica
    magnetStd = (magnetStd)**(0.5) / iters_root
    energyStd = (energyStd)**(0.5) / iters_root
    energySqStd = (energySqStd)**(0.5) / iters_root

    #----BARRAS DE ERROR 2SIGMA----
    magnetError = 2*magnetStd/iters_root         

    energyError = energyStd/(2*N2)               # Solo hace falta multiplicar la cte
    energyError *= 2/iters_root

    # ***Calculamos la incertidumbre del calor específico en un script aparte mediante un método Monte Carlo
    # llamado "specheat_error.py"

    resultMagnet = np.array([magnetAvg, magnetError])
    resultEnergy = np.array([eAvg, energyError])
    resultSpecHeat = np.array([energyAvg, energyStd, energySqAvg, energySqStd, specHeatAvg])

    return resultMagnet, resultEnergy, resultSpecHeat, correlArrAvg


"FUNCIÓN PRINCIPAL"

#---PARÁMETROS---
stepsMC = int(1e6)                  # Se van a realizar 10^6 pasos Monte Carlo

# Rango de temperaturas (elegir rango original o cercano al punto crítico)
tmin = 2.1666666666666665
tmax = 2.611111111111111            
tnum = 10

temp = np.linspace(tmin, tmax, tnum)

N = np.array((16, 32, 64, 128))     # 4 tamaños distintos para la configuración de espines

choice = 2                          # Fijamos configuracion inicial ordenada

#---SIMULACIÓN---

sizes = len(N)
temps = len(temp)

fMagnet = open("magnet.dat", "w")
fEnergy = open("energy.dat", "w")
fSpecHeat1 = open("specheat1.dat", "w") # Aquí energía, energía^2 y desv. est.
fSpecHeat2 = open("specheat2.dat", "w") # Aquí solo el promedio
fCorrel = open("correl.dat", "w")

count = 0

# Para cada N
for i in range(sizes):

    # Se genera la config. inicial
    spins = iniSim(choice, N[i])

    # Para cada T
    for j in range(temps):

        start = time.time()

        resultMagnet, resultEnergy, resultSpecHeat, correlArrAvg = simulationExtra(temp[j], N[i], spins, stepsMC)

        # Se guardan en formato:
        # T_0,  m_0_0,  merror_0_0
        # T_1,  m_0_1,  merror_0_1
        # ...
        # T_9,  m_0_9,  merror_0_9
        # 
        # T_0,  m_1_0,  merror_1_0
        # T_1,  m_1_1,  merror_1_1
        # ...

        fMagnet.write(f"{temp[j]}, {resultMagnet[0]}, {resultMagnet[1]}\n")
        fEnergy.write(f"{temp[j]}, {resultEnergy[0]}, {resultEnergy[1]}\n")

        fSpecHeat1.write(f"{temp[j]}, {resultSpecHeat[0]}, {resultSpecHeat[1]}, {resultSpecHeat[2]}, {resultSpecHeat[3]}\n")
        fSpecHeat2.write(f"{temp[j]}, {resultSpecHeat[4]}\n")

        # Se guarda en formato:
        # T_0
        # f_0(0)
        # f_0(1)
        # ...
        # f_0(N_i-1)
        #
        # T_1
        # f_1(0)
        # f_1(1)
        # ...
        # f_1(N_i-1)
        #
        # ...

        fCorrel.write(f"{temp[j]}\n")
        for elem in range(len(correlArrAvg)):
            fCorrel.write(f"{correlArrAvg[elem]}\n")
        fCorrel.write("\n")

        count += 1

        end = time.time()

        print(f"Se han terminado ({count}) iteraciones completas.")
        print(f"Tiempo empleado en esta iteracion: {end - start} segundos")

    fMagnet.write("\n")
    fEnergy.write("\n")
    fSpecHeat1.write("\n")
    fSpecHeat2.write("\n")
    fCorrel.write("\n")

fMagnet.close()
fEnergy.close()
fSpecHeat1.close()
fSpecHeat2.close()
fCorrel.close()  