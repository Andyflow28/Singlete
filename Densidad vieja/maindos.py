import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

# Lista de archivos .dat de entrada y sus correspondientes valores de zeta
input_files = [
    ('DOSc00g000E04t02.dat', r'0.0001'),
    ('DOSc00g005E04t02.dat', r'0.0005'),
    ('DOSc00g010E04t02.dat', r'0.0010'),
    ('DOSc00g015E04t02.dat', r'0.0015'),
    ('DOSc00g020E04t02.dat', r'0.0020')
]

# Crear una figura para las graficaciones
plt.figure(figsize=(10, 6))

# Procesar cada archivo
for input_file, zeta_value in input_files:
    # Verificar si el archivo existe
    if os.path.isfile(input_file):
        # Leer el archivo .dat
        data = pd.read_csv(input_file, sep='\s+', header=None)

        # Convertir y guardar como .csv
        output_file = input_file.replace('.dat', '.csv')
        data.to_csv(output_file, index=False, header=False)

        # Graficar los datos
        energy = data[0]  # Suponiendo que la primera columna es la energía
        frequency = data[1]  # Suponiendo que la segunda columna es la frecuencia (gamma)

        # Comprobar si hay valores NaN y eliminarlos
        mask = ~np.isnan(energy) & ~np.isnan(frequency)
        energy = energy[mask]
        frequency = frequency[mask]

        # Interpolar para obtener valores entre el rango deseado
        f_interpolate = interp1d(energy, frequency, kind='linear', fill_value="extrapolate")
        energy_new = np.linspace(0, 3, num=100)  # Crear 100 puntos entre 0 y 1.90
        frequency_new = f_interpolate(energy_new)

        # Graficar la interpolación
        plt.plot(energy_new, frequency_new, marker='o', linestyle='-', label=r'$\zeta = ' + zeta_value + '$')
    else:
        print(f"Advertencia: {input_file} no se encuentra.")

# Configurar el gráfico
plt.title('Gráfico de Frecuencia vs. Energía')
plt.xlabel('Energía (eV)')
plt.ylabel(r'Frecuencia ($\omega$)')
plt.xlim(0, 3)  # Establecer límites del eje X entre 0 y 1.90
plt.grid()
plt.legend()  # Mostrar la leyenda
plt.savefig('grafico_frecuencia_vs_energia_dos.png')  # Guarda el gráfico como un archivo PNG
plt.show()  # Muestra el gráfico