"//////////////////////////////////////////////////////////////"
# Programa que grafica la distancia Tierra-Marte en función del
# tiempo.
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

import matplotlib.pyplot as plt
import numpy as np

file = "in_orbit_dist.dat"

fig, ax = plt.subplots()
with open(file, "r") as file:

    # Leo el archivo
    lines=file.readlines()

    t = [float(line.split()[0]) for line in lines]
    r = [float(line.split()[1]) for line in lines]
    t=np.array(t)/10**6
    d=np.array(r)

    t, d = zip(*sorted(zip(t, d)))

    # Graficamos
    ax.plot(t, d, marker = "^", ms = 3, lw = 0.8)

    # Fijamos para la simulación actual (9) cuál fue el paso temporal y el momento de llegada a Marte
    h = 15  # [s]
    duration = 22603065 # [s]

    ax.axvline(x=duration/(h*10**6), linestyle="dashed", color = "blue", label="Llegada a Marte")

    ax.set_ylabel(r"Distancia nave-Marte normalizada")
    ax.set_xlabel(r"Tiempo normalizado $t/h ~(\times 10^6)$")

    ax.set_xlim(0,)

    ax.legend()

    plt.title("Evolución de la distancia nave-Marte", fontname = "Arial Rounded MT Bold")
    plt.show()