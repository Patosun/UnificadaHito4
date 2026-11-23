"""
Capítulo 7 - Ejercicio 19
=========================
Enunciado: Resolver la EDO del Problema 18 con los programas de los
métodos, usando h = 0.05. Tabla con soluciones, exacta y errores
verdaderos; comentar precisión y TIEMPO DE CÓMPUTO de cada método.

EDO del Problema 18:
        dy/dx + y^2 + y/x = 1/x^2  ->  dy/dx = 1/x^2 - y^2 - y/x
        1 <= x <= 2 ,   y(1) = -1 ,   exacta: y(x) = -1/x

RESTRICCIÓN DE LA PRÁCTICA: el libro pide Euler, Heun, Euler modificado y
RK4; aquí sólo se usan EULER y RK4.

Medición del tiempo: una sola ejecución tarda microsegundos, por debajo de la
resolución práctica del reloj, así que cada método se ejecuta N_REPETICIONES
veces y se reporta el tiempo promedio por ejecución (time.perf_counter).
"""

import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import euler, rk4, imprimir_tabla, ruta_grafica

X0, XF, Y0 = 1.0, 2.0, -1.0
H = 0.05
N_REPETICIONES = 2000
METODOS = {"Euler": euler, "RK4": rk4}
# Evaluaciones de f por paso: Euler usa 1 pendiente, RK4 usa 4 (k1..k4)
EVALS_POR_PASO = {"Euler": 1, "RK4": 4}


def f(x, y):
    """dy/dx = 1/x^2 - y^2 - y/x."""
    return 1.0 / x ** 2 - y ** 2 - y / x


def exacta(x):
    """y(x) = -1/x."""
    return -1.0 / x


def medir_tiempo(metodo):
    """Tiempo promedio (s) de UNA resolución completa con el método dado."""
    t0 = time.perf_counter()
    for _ in range(N_REPETICIONES):
        metodo(f, X0, Y0, H, XF)
    return (time.perf_counter() - t0) / N_REPETICIONES


def resolver():
    """Devuelve x, {nombre: y}, {nombre: tiempo}."""
    soluciones, tiempos = {}, {}
    x = None
    for nombre, metodo in METODOS.items():
        x, soluciones[nombre] = metodo(f, X0, Y0, H, XF)
        tiempos[nombre] = medir_tiempo(metodo)
    return x, soluciones, tiempos


def imprimir_resumen(x, soluciones, tiempos):
    n_pasos = len(x) - 1
    y_ex = exacta(x)
    print(f"\nResumen (h = {H}, {n_pasos} pasos, promedio de "
          f"{N_REPETICIONES} ejecuciones)")
    cab = (f"{'Método':>7} | {'evals f':>8} | {'tiempo [ms]':>12} | "
           f"{'err máx abs':>12} | {'err máx %':>10}")
    print("-" * len(cab))
    print(cab)
    print("-" * len(cab))
    for nombre, y in soluciones.items():
        err = np.abs(y_ex - y)
        print(f"{nombre:>7} | {EVALS_POR_PASO[nombre] * n_pasos:8d} | "
              f"{tiempos[nombre] * 1e3:12.4f} | {err.max():12.4e} | "
              f"{(err / np.abs(y_ex)).max() * 100:10.4e}")
    print("-" * len(cab))
    razon_t = tiempos["RK4"] / tiempos["Euler"]
    razon_e = (np.abs(y_ex - soluciones["Euler"]).max()
               / np.abs(y_ex - soluciones["RK4"]).max())
    print(f"RK4 tarda {razon_t:.1f} veces más que Euler, pero su error es "
          f"{razon_e:.1e} veces menor.")


def graficar(x, soluciones):
    xx = np.linspace(X0, XF, 300)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(xx, exacta(xx), "k-", lw=2, label="Exacta: $-1/x$")
    ax1.plot(x, soluciones["Euler"], "o", ms=4, color="tab:red", label="Euler")
    ax1.plot(x, soluciones["RK4"], "s", ms=4, mfc="none", color="tab:blue",
             label="RK4")
    for nombre, color in (("Euler", "tab:red"), ("RK4", "tab:blue")):
        ax2.semilogy(x[1:], np.abs(exacta(x[1:]) - soluciones[nombre][1:]),
                     "o-", ms=3, color=color, label=nombre)
    ax1.set(xlabel="x", ylabel="y", title=f"Ej. 7.19: h = {H}")
    ax2.set(xlabel="x", ylabel="|error|", title="Error verdadero (log)")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3, which="both")
        ax.legend()
    fig.tight_layout()
    ruta = ruta_grafica("cap7_ej19_euler_rk4.png")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 19: dy/dx = 1/x^2 - y^2 - y/x, Euler y RK4, "
          "h = 0.05")
    print("=" * 70)
    x, soluciones, tiempos = resolver()
    imprimir_tabla(x, soluciones, exacta(x), titulo=f"Euler y RK4 con h = {H}")
    imprimir_resumen(x, soluciones, tiempos)
    graficar(x, soluciones)


if __name__ == "__main__":
    main()
