"""
main.py
=======
Ejecuta en orden todos los ejercicios computacionales de la práctica.

Uso:
    python main.py              -> corre todos los ejercicios
    python main.py cap8_ej03    -> corre sólo los scripts indicados

Cada script también se puede ejecutar por separado, por ejemplo:
    python cap7_ej16.py
"""

import importlib
import sys
import time

# Orden de ejecución: primero las gráficas complementarias del Cap. 7,
# luego los ejercicios del Cap. 7 y al final los del Cap. 8.
SCRIPTS = [
    "cap7_ej01_grafica",
    "cap7_ej02",
    "cap7_ej03_grafica",
    "cap7_ej04",
    "cap7_ej05",
    "cap7_ej16",
    "cap7_ej17",
    "cap7_ej18_grafica",
    "cap7_ej19",
    "cap7_ej20",
    "cap7_ej21_grafica",
    "cap7_ej22",
    "cap8_ej03",
    "cap8_ej04",
]


def ejecutar(nombre):
    """Importa el módulo del ejercicio y llama a su función main()."""
    modulo = importlib.import_module(nombre)
    t0 = time.perf_counter()
    modulo.main()
    return time.perf_counter() - t0


def main():
    seleccion = sys.argv[1:] or SCRIPTS
    desconocidos = [s for s in seleccion if s not in SCRIPTS]
    if desconocidos:
        print(f"Scripts desconocidos: {desconocidos}\nDisponibles: {SCRIPTS}")
        sys.exit(1)

    tiempos, fallos = {}, []
    for nombre in seleccion:
        print("\n\n" + "#" * 78)
        print(f"#  {nombre}.py")
        print("#" * 78)
        try:
            tiempos[nombre] = ejecutar(nombre)
        except Exception as e:          # se sigue con los demás ejercicios
            fallos.append(nombre)
            print(f"\n*** ERROR en {nombre}: {type(e).__name__}: {e}")

    print("\n\n" + "=" * 50)
    print("RESUMEN DE EJECUCIÓN")
    print("=" * 50)
    for nombre, t in tiempos.items():
        print(f"  OK     {nombre:<22} {t:8.2f} s")
    for nombre in fallos:
        print(f"  FALLÓ  {nombre}")
    print("Las gráficas quedaron en la carpeta graficas/")
    sys.exit(1 if fallos else 0)


if __name__ == "__main__":
    main()
