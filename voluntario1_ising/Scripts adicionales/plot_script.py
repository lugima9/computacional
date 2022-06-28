"//////////////////////////////////////////////////////////////"
# Programa que grafica los resultados del modelo de Ising para 
# distintas temperaturas y tamaños
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools as it

# Lee los ficheros de magnetización, energía y calor específico CON barras de error
# y los grafica para cada N
files = ("magnet_def.dat", "energy_def.dat", "specheat_def.dat")
ylabels = ("Magn. prom", "Energía prom.", "Calor esp. prom.")

"MAGNETIZACIÓN"

markers = it.cycle(('x', 's', 'v', 'o'))
fig, ax = plt.subplots(1)
with open(files[0], "r") as file:

    # Leo el archivo
    lines=file.readlines()

    # Empiezo con N = 16
    N = 16
    # Vacío variables que almacenan x, y, error de y
    x=[]
    y=[]
    yerror=[]
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
            # ... plotteo magnitud vs temperatura para este N
            ax.errorbar(x, y, yerror, label = f"N = {N}", marker = next(markers), ms = 3, lw = 0.8)
            ax.axhline(y=0.5, color='gray', linestyle='dashed', linewidth = 0.5)
            ax.legend()
            ax.set_ylim(0,1)

            # Paso al siguiente bloque de N
            N *= 2
            x=[]
            y=[]
            yerror=[]

    ax.set_xlabel("Temperatura [$J/k_B$]")
    ax.set_ylabel(r"Magn. prom. $m_N$")
    plt.title("Magnetización", fontname = "Arial Rounded MT Bold")
    plt.show()

"ENERGÍA"

markers = it.cycle(('x', 's', 'v', 'o'))
fig, ax = plt.subplots(1)
with open(files[1], "r") as file:

    # Leo el archivo
    lines=file.readlines()

    # Empiezo con N = 16
    N = 16
    # Vacío variables que almacenan x, y, error de y
    x=[]
    y=[]
    yerror=[]
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
            # ... plotteo magnitud vs temperatura para este N
            ax.errorbar(x, y, yerror, label = f"N = {N}", marker = next(markers), ms = 3, lw = 0.8)
            ax.legend()
            ax.set_ylim(-1,-0.3)

            # Paso al siguiente bloque de N
            N *= 2
            x=[]
            y=[]
            yerror=[]

    ax.set_xlabel("Temperatura [$J/k_B$]")
    ax.set_ylabel(r"Energía prom. $e_N$")
    plt.title("Energía interna", fontname = "Arial Rounded MT Bold")
    plt.show()

"CALOR ESPECÍFICO"

markers = it.cycle(('x', 's', 'v', 'o'))
fig, ax = plt.subplots(1)
with open(files[2], "r") as file:

    # Leo el archivo
    lines=file.readlines()

    # Empiezo con N = 16
    N = 16
    # Vacío variables que almacenan x, y, error de y
    x=[]
    y=[]
    yerror=[]
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
            # ... plotteo magnitud vs temperatura para este N
            ax.errorbar(x, y, yerror, label = f"N = {N}", marker = next(markers), ms = 3, lw = 0.8)
            ax.legend()
            ax.set_ylim(0,5)

            # Paso al siguiente bloque de N
            N *= 2
            x=[]
            y=[]
            yerror=[]

    ax.set_xlabel("Temperatura [$J/k_B$]")
    ax.set_ylabel(r"Calor espec. $c_N$")
    plt.title("Calor específico", fontname = "Arial Rounded MT Bold")
    plt.show()


"FUNCIÓN DE CORRELACIÓN"

fileName = "correl_def.dat"
tmin = 2.1666666666666665 # Rango de temperaturas cercano al punto crítico
tmax = 2.611111111111111
tnum = 10

fig, axs = plt.subplots(2, 2)
with open(fileName, "r") as file:

    temp=np.linspace(tmin, tmax, tnum)
    temp_normalized = (temp-tmin)/(tmax-tmin)
    colors = [ cm.jet(t) for t in temp_normalized ]
    data=file.readlines()
    stop_old = 0

    for N in [16, 32, 64, 128]:

        # Indico en qué subplot mostrar cada conjunto de datos
        if N == 16:
            i, j = 0, 0
        elif N == 32:
            i, j = 0, 1
        elif N == 64:
            i, j = 1, 0
        else:
            i, j = 1, 1

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
                x,y = np.array(x),np.array(y)
                # ... plotteo correl vs dist para esta temp
                axs[i,j].plot(x[x <= N//2], y[x <= N//2], label = f"T = {temp[temp_count]}", color = colors[temp_count], lw = 1)
                axs[i,j].set_xlim(0, N//2)
                axs[i,j].set_ylim(0,1)

                # Paso al siguiente bloque de temp
                temp_count += 1
                dist = 0
                x = []
                y = []

        axs[0,0].set_ylabel(r"Correl. $f(i)$")
        axs[1,0].set_ylabel(r"Correl. $f(i)$")
        axs[i,j].set_xlabel(r"Distancia ($i$)")
        axs[i,j].set_title(fr"$N = {N}$")

        stop_old = stop_new

    plt.suptitle("Función de correlación", fontname = "Arial Rounded MT Bold")
    plt.tight_layout()
    plt.show()


"FUNCIÓN DE CORRELACIÓN N = 128"

fig, ax = plt.subplots()
with open(fileName, "r") as file:

    temp=np.linspace(tmin, tmax, tnum)
    temp_normalized = (temp-tmin)/(tmax-tmin)
    colors = [ cm.jet(t) for t in temp_normalized ]
    data=file.readlines()
    stop_old = 0

    for N in [16, 32, 64, 128]:

        # Indico en qué subplot mostrar cada conjunto de datos
        if N == 16:
            i, j = 0, 0
        elif N == 32:
            i, j = 0, 1
        elif N == 64:
            i, j = 1, 0
        else:
            i, j = 1, 1

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

                # ... plotteo correl vs dist para esta temp
                if N == 128:
                    x,y = np.array(x),np.array(y)
                    ax.plot(x[x <= N//2], y[x <= N//2], label = f"T = {temp[temp_count]}", color = colors[temp_count], lw = 1)

                # Paso al siguiente bloque de temp
                temp_count += 1
                dist = 0
                x = []
                y = []
                
        ax.set_xlabel(r"Distancia ($i$)")
        ax.set_ylabel(r"Correl. $f(i)$")

        stop_old = stop_new

    # # Barra de color con temperaturas
    # cmap = plt.get_cmap("jet", 10)
    # norm = mpl.colors.Normalize(vmin=tmin, vmax=tmax)
    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])

    # cbar = plt.colorbar(sm, ticks=np.linspace(tmin, tmax, 10))
    # cbar.ax.set_ylabel(r"Temperatura [$J/k_B$]", rotation=270)
    # cbar.ax.yaxis.set_label_coords(0.75, 0.5)

    # plt.suptitle(r"Función de correlación ($N=128$)", fontname = "Arial Rounded MT Bold")
    # plt.show()
