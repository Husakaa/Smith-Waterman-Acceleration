import re
import matplotlib.pyplot as plt

# Inicializar listas de tiempos
tiempos_secuencial = []
tiempos_openmp = []
tiempos_cuda = []

# Leer archivo
with open(r"times.txt", "r", encoding="utf-8") as f:
    contenido = f.read()

# Usar expresiones regulares para capturar los tiempos
secuencial_matches = re.findall(r"Tiempo \d+ para Secuencial.*?Tiempo de CPU total:\s*([\d.]+)", contenido, re.DOTALL)
openmp_matches = re.findall(r"Tiempo \d+ para OpenMP.*?Tiempo de CPU total:\s*([\d.]+)", contenido, re.DOTALL)
cuda_matches = re.findall(r"Tiempo \d+ para CUDA.*?Tiempo GPU \+ transferencias:\s*([\d.]+)", contenido, re.DOTALL)

# Convertir a float
tiempos_secuencial = [float(t) for t in secuencial_matches]
tiempos_openmp = [float(t) for t in openmp_matches]
tiempos_cuda = [float(t)/1000 for t in cuda_matches]

# Calcular medias
def media(lista):
    return sum(lista) / len(lista) if lista else 0.0

promedios = {
    "Secuencial": media(tiempos_secuencial),
    "OpenMP": media(tiempos_openmp),
    "CUDA": media(tiempos_cuda)
}

# Imprimir promedios
print("Promedios de ejecución (segundos):")
for nombre, valor in promedios.items():
    print(f"{nombre}: {valor:.3f}")

# Graficar
plt.figure(figsize=(8, 5))
plt.bar(promedios.keys(), promedios.values(), color=["gray", "orange", "teal"])
plt.title("Tiempo medio de ejecución por aproximación")
plt.ylabel("Tiempo (segundos)")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("tiempos_promedio.png")
plt.show()

