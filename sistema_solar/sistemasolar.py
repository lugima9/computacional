"#//////////////////////////////////////////////////////////////"
# Programa que simula la dinámica del Sistema Solar.
# OUTPUT: Periodo de traslación de cada planeta y comparación el tabulado,
#         energía cinética y momento angular total del sistema (conservado?), 
#         animación de las órbitas en modelo heliocéntrico y geocéntrico
#         (restando la posición de la Tierra a la de cada planeta)         
# Autor: Luis Gil Martín
# Universidad de Granada
"#//////////////////////////////////////////////////////////////"

import math
import matplotlib.pyplot as plt
import numpy as np
import random
import itertools

"------------------------------------CONSTANTES------------------------------------------"

c = 1.496e11            # [m] (distancia Tierra-Sol)
masaSol = 1.98847e30    # [kg]
G = 6.67430e-11         # [m3 kg-1 s-2]

"------------------------------------FUNCIONES------------------------------------------"

# Función que calcula la aceleración (x,y) de un planeta en un instante por acción del Sol y el resto de planetas
# Recibe: índice del planeta, array de posiciones de cada planeta en (x,y), vector de masas de cada planeta
# Devuelve: array 2D con aceleración (x,y) del planeta
def calcAcelPlaneta(estePlaneta, pos, masas):

    numeroPlanetas = len(pos)
    # Posición (x,y) con respecto al Sol del planeta y distancia cubo
    difPos = pos[estePlaneta]

    distCubo = (difPos[0]*difPos[0] + difPos[1]*difPos[1])**1.5

    # Aceleración (x,y) del planeta debido al Sol
    acel = -difPos/distCubo

    # Aceleración (x,y) del planeta debido al resto de planetas
    for i in range(numeroPlanetas):
        if i != estePlaneta:
            difPos = pos[estePlaneta] - pos[i]     #Posición relativa entre planetas
            distCubo = (difPos[0]*difPos[0] + difPos[1]*difPos[1])**1.5
            acel -= masas[i]*difPos/distCubo                     #Sumamos la aceleración por el resto de planetas a la del Sol

    return acel


# Función que calcula un paso del algoritmo de Verlet para obtener la posición y velocidad de cada planeta en el tiempo
# Recibe: paso, masas, posiciones, velocidades
# Devuelve: array de posiciones, array de velocidades
def pasoVerlet(h, m, r, v):

    numeroPlanetas = len(r)
    #Creo array de acel. 3D de 8 (nº planetas) x 2 (x,y)
    acel = np.zeros((numeroPlanetas, 2))

    # Para cada uno de los 8 planetas...
    for i in range(numeroPlanetas):
        #Calcula la aceleración (x,y) del planeta       
        acel[i] = calcAcelPlaneta(i, r, m) #OJO: se mete el array r entero
    
    #Cálculo intermedio (x,y) del planeta
    w = v + 0.5*h*acel

    #Calcula la nueva posición (x,y) del planeta
    r += h*w
    
    for i in range(numeroPlanetas):
        #Calcula la nueva aceleración (x,y) del planeta
        acel[i] = calcAcelPlaneta(i, r, m)

    #Calcula la nueva velocidad (x,y) del planeta
    v = w + 0.5*h*acel

    return r, v


# Función que obtiene los datos de los planetas de un fichero de texto
# Recibe:  paso, t_max, nombre del fichero con una fila por planeta, bool de alineación
# Devuelve: paso, t_max, nombres, masas, posiciones, velocidades (reescalados)
def condIni(nombreFichero, alinear):
    # Abrimos el fichero en modo lectura
    f = open(nombreFichero, "r")

    # Leemos el archivo saltando el header y separamos en lista por planeta
    next(f)
    datos = f.read().split('\n')
    f.close()

    # Guardamos cada dato de cada planeta en un elemento distinto de cada np.array
    nombres, masas, distancias, velocidades = np.array([]), np.array([]), np.array([]), np.array([])

    for i in range(len(datos)):
        d = datos[i].split()
        nombres = np.append(nombres, d[0])
        masas = np.append(masas, float(d[1]))
        distancias = np.append(distancias, float(d[2]))
        velocidades = np.append(velocidades, float(d[3]))
    

    # Reescalamos los datos
    masas, distancias, velocidades = reescala(masas, distancias, velocidades)

    # Iniciamos posiciones y velocidades
    if alinear == True:
        r, v = inicioAlineado(distancias, velocidades)
    else:
        r, v = inicioAleatorio(distancias, velocidades)

    return nombres, masas, r, v


# Función que realiza el cambio de unidades de reescalado
# Recibe: paso, t_max, array de masas, array de posiciones iniciales y array de velocidades iniciales
# Devuelve: paso, t_max, array de masas, array de posiciones iniciales y array de velocidades iniciales (reescalado)
def reescala(masas, distancias, velocidades):

    c = 1.496e11            # [m] (distancia Tierra-Sol)
    masaSol = 1.98847e30    # [kg]
    G = 6.67430e-11         # [m3 kg-1 s-2]

    masas = masas/masaSol
    distancias = distancias/c
    velocidades = velocidades*(c/(G*masaSol))**(0.5)

    return masas, distancias, velocidades


# Función que realiza el cambio de unidades de desescalado
# Recibe: array de masas, array de posiciones y array de velocidades
# Devuelve: array de masas, array de posiciones y array de velocidades iniciales (desescalado)
def desescala(masas, distancias, velocidades):

    c = 1.496e11            # [m] (distancia Tierra-Sol)
    masaSol = 1.98847e30    # [kg]
    G = 6.67430e-11         # [m3 kg-1 s-2]

    masas = masas*masaSol
    distancias = distancias*c
    velocidades = velocidades*(G*masaSol/c)**(0.5)

    return masas, distancias, velocidades


# Función que aleatoriza el inicio de la simulación
# Recibe: distancias y velocidades iniciales
# Devuelve: array de posiciones y array de velocidades con formato (x,y) para cada elemento-planeta
def inicioAleatorio(distancias, velocidades):

    numeroPlanetas = len(distancias)
    # Para cada planeta (tantos como elementos en distancias[])
    r, v = np.zeros((numeroPlanetas, 2)), np.zeros((numeroPlanetas, 2))
    for i in range(numeroPlanetas):
        # Genero una coordenada polar aleatoria
        angulo = random.random()*2*np.pi
        coseno = math.cos(angulo)
        seno = math.sin(angulo)

        # Calculo una posición inicial en cartesianas
        r[i,0] = distancias[i]*coseno
        r[i,1] = distancias[i]*seno

        # Calculo una velocidad inicial perpendicular al vector r en cartesianas
        v[i,0] = -velocidades[i]*seno
        v[i,1] = velocidades[i]*coseno        

    return r, v


# Función que alinea ápsides planetas al inicio de la simulación
# Recibe: distancias y velocidades iniciales
# Devuelve: array de posiciones y array de velocidades con formato (x,y) para cada elemento-planeta
def inicioAlineado(distancias, velocidades):

    numeroPlanetas = len(distancias)
    # Para cada planeta (tantos como elementos en distancias[])
    r, v = np.zeros((numeroPlanetas, 2)), np.zeros((numeroPlanetas, 2))
    for i in range(numeroPlanetas):

        # Fijo posición inicial en (x,0)
        r[i,0] = distancias[i]
        r[i,1] = 0.

        # Fijo velocidad inicial en (0, v)
        v[i,0] = 0.
        v[i,1] = velocidades[i]  

    return r, v


# Función que calcula el periodo en la primera vuelta si los planetas comienzan alineados (NO NECESITA DESESCALADO) 
# Recibe: instante, arrays con posición (x,y) y velocidad (x,y) instantánea de cada planeta, paso
# Devuelve: array con periodos simulados de cada planeta
def calcPeriodoAlin(t, r, v, paso, obtenido):

    numeroPlanetas = len(r)

    # A partir de cierto número de iteraciones de arranque, empezamos a obtener el periodo
    if t > 100.*paso:
        # Para cada planeta...
        for i in range(numeroPlanetas):
            # Si su coordenada y se encuentra cerca de 0, siendo x positiva, se obtiene el tiempo que ha pasado
            if obtenido[i] == False:
                if (r[i, 1] < 0.01 and r[i, 1] > -0.01) and (r[i, 0] > 0.):
                    obtenido[i] = True
                    break
                else:
                    periodos[i] += h

    return periodos, obtenido


# Función que calcula el error relativo de los periodos con respecto a los valores tabulados
# Recibe: array de periodos calculados, fichero con periodos tabulados
# Devuelve: array con periodos tabulados, array de error relativo por planeta
def errorRelPeriodo(periodos, ficheroTab):
        
    numeroPlanetas = len(periodos)

    # Inicializamos el array de periodos
    periodosTab = np.zeros(numeroPlanetas)

    #Leemos el fichero con periodo de cada planetas
    with open(ficheroTab, "r") as fTab:
        next(fTab)
        datosTab = fTab.read().split("\n")

    # Introducimos los periodos tabulados en un array
    for i in range(numeroPlanetas):
        d = datosTab[i].split()
        periodosTab[i] = float(d[1])

    # Calculamos el array de errores relativos
    errorRel = np.array([])
    for i in range(numeroPlanetas):
        errorRel = np.append(errorRel, 100*abs(periodos[i] - periodosTab[i])/periodosTab[i])

    return periodosTab, errorRel


# Función que calcula el momento angular de los planetas
# Recibe: array de masas, array posiciones (x,y) y array de velocidad (x,y) de cada planeta
# Devuelve: módulo del momento angular total del Sistema Solar
def calcMomAng(m, r, v):
    numeroPlanetas = len(r)
    momAng = np.zeros(3)

    # Desescalado
    m, r, v = desescala(m, r, v)

    for i in range(numeroPlanetas):
        momAng += m[i]*np.cross(r[i], v[i])

    momAngMod = momAng[0]*momAng[0] + momAng[1]*momAng[1] + momAng[2]*momAng[2]    
    momAngMod = math.sqrt(momAngMod)
    
    return momAngMod


# Función que calcula la energía de los planetas
# Recibe: array de posiciones (x,y) y velocidad (x,y) de cada planeta
# Devuelve: energía total del Sistema Solar
def calcEnergia(m, r, v):

    numeroPlanetas = len(r)
    energia = 0.
    energiaCin = 0.
    energiaPot = 0.
    G = 6.67430e-11         # [m3 kg-1 s-2]

    # Desescalado
    m, r, v = desescala(m, r, v)
    
    for i in range(numeroPlanetas):
        energiaCin += m[i]*(v[i,0]*v[i,0] + v[i,1]*v[i,1])
        for j in range(numeroPlanetas):
            if j != i:
                dist = r[i] - r[j]
                distCuad = dist[0]*dist[0] + dist[1]*dist[1]
                energiaPot += m[j]/distCuad
        energia += 0.5*energiaCin + G*m[i]*energiaPot

    return energia

"------------------------------------FUNCIÓN PRINCIPAL------------------------------------------"

# Preguntamos qué se desea hacer

eleccion = -1
while eleccion != 1 and eleccion != 2 and eleccion != 3:
    print("Si desea hacer solo la simulacion, introduzca 1")
    print("Si desea además calcular los periodos de los planetas y energía promedio del sistema en el tiempo, introduzca 2")
    print("Si desea visualizar el modelo geocentrico, introduzca 3")
    eleccion = int(input("Que desea hacer? \n"))

# Preguntamos condiciones de simulación
t=0.

t_max = -1.
while t_max < 0.:
    t_max = float(input("Introduzca el tiempo máximo de simulación: \n"))

h = -1.
while h < 0 or h > t_max:
    h = float(input("Introduzca el paso: \n"))

alinear = ""
while alinear != "Y" and alinear != "N":
    alinear = input("¿Desea que los planetas comiencen alineados? (Y/N): \n")

if alinear == "Y":
    alinear = True
else:
    alinear = False

# Leemos condiciones iniciales
nombreArchivo = "datos_planetas.txt"
nombres, masas, r, v = condIni(nombreArchivo, alinear)

# Abrimos un fichero donde se van a guardar los resultados de la simulación
fPos = open("planets_data.dat", "w")

# Iniciamos la simulación
numeroPlanetas = len(nombres)

# Si se ha elegido obtener también periodos, momento angular y energía (inicio alineado)
if eleccion == 2 and alinear == True:

    #Inicializamos el array que guarda cada instante
    tiempo = np.array([])

    #Inicializamos el array con periodos con el instante en el que se empieza a contar
    periodos = np.full(numeroPlanetas, 100.*h)
    obtenido = np.full(numeroPlanetas, False)
    vueltas = np.zeros(numeroPlanetas)

    #Inicializamos la energía y momento angular (módulo) totales
    momAngArr, energiaArr = np.array([]), np.array([])

    mostrar = 50
    while t < t_max:

        # Escribimos las posiciones en el fichero (1 de cada 50)
        if mostrar == 50:
            for i in range(numeroPlanetas): # Para cada planeta
                fPos.write(f"{r[i,0]},  {r[i,1]}\n")
            fPos.write("\n")

            mostrar = 0
        
        mostrar += 1

        r, v = pasoVerlet(h, masas, r, v)

        #HACER ESTO CON FUNCIONES APARTE, QUE SI NO ES UN PUTO LÍO

        #--CÁLCULO DEL PERIODO--
        periodos, obtenido= calcPeriodoAlin(t, r, v, h, obtenido)
        
        #--CÁLCULO DEL MOMENTO ANGULAR--
        momAngMod = calcMomAng(masas, r, v)
        momAngArr = np.append(momAngArr, momAngMod)

        #--CÁLCULO DE LA ENERGÍA--
        energia = calcEnergia(masas, r, v)
        energiaArr = np.append(energiaArr, energia)

        # Pasamos al siguiente instante
        t += h
        tiempo = np.append(tiempo, t)

    # Mostramos en un fichero los periodos (desescalados) y sus errores relativos

    pFich = open("resultados_periodos.dat", "w")

    pFich.write("El periodo de cada planeta es: \n")
    periodos = 58.1*periodos
    #*(c**3/(G*masaSol))**0.5
    for i in range(numeroPlanetas):
        pFich.write(f"{nombres[i]},  {periodos[i]}\n")

    periodosTab, errorRel = errorRelPeriodo(periodos, "periodos_tab.dat")
    pFich.write("Los periodos tabulados son: \n")
    for i in range(numeroPlanetas):
        pFich.write(f"{nombres[i]},  {periodosTab[i]}\n")

    pFich.write("Los errores relativos son: \n")
    for i in range(numeroPlanetas):
        pFich.write(f"{nombres[i]},  {errorRel[i]}\n")

    pFich.close()

    # Mostramos dos gráficas con momento angular y energía en el tiempo
    fig, axs = plt.subplots(2)
    axs[0].plot(tiempo, momAngArr)
    axs[1].plot(tiempo, energiaArr)

    axs[0].set_xlabel("Tiempo")
    axs[1].set_xlabel("Tiempo")
    axs[0].set_ylabel("Momento angular total(kg m s-1)")
    axs[1].set_ylabel("Energía total (J)")

    plt.show() 

# Si se ha seleccionado calcular periodos con inicio desalineado, para el programa.
elif eleccion == 2 and alinear == False:
    print("Para calcular periodos, los planetas deben iniciarse alineados.")

# Si solo se desea simular...
elif eleccion == 1:

    mostrar = 50
    while t < t_max:

        # Escribimos las posiciones en el fichero (1 de cada 50)
        if mostrar == 100:
            for i in range(numeroPlanetas): # Para cada planeta
                fPos.write(f"{r[i,0]},  {r[i,1]}\n")
            fPos.write("\n")

            mostrar = 0
        
        mostrar += 1
        
        r, v = pasoVerlet(h, masas, r, v)

        t += h

else:

    mostrar = 50
    while t < t_max:

        # Escribimos las posiciones en el fichero referidas a la Tierra(1 de cada 50)
        if mostrar == 50:
            for i in range(numeroPlanetas): # Para cada planeta
                fPos.write(f"{r[i,0] - r[2,0]},  {r[i,1] - r[2,1]}\n")
            fPos.write("\n")

            mostrar = 0
        
        mostrar += 1
        
        r, v = pasoVerlet(h, masas, r, v)

        t += h
    

fPos.close()
print("Los archivos se han creado con exito en la localizacion del script.")