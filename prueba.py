import numpy as np
# Abrimos el fichero en modo lectura
f = open("datos_planetas.txt", "r")

#Leemos el archivo saltando el header y separamos en lista por planeta
next(f)
datos = f.read().split('\n')



f.close()

print(datos)
print(type(datos))