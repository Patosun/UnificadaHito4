"""
Capítulo 7 - Ejercicio 16
=========================
Enunciado: "Repetir el Problema 12 con programas de computadora, usando
h = 0.25, 0.1 y 0.01. Graficar para comparar las soluciones con distintos h
con la exacta, y armar una tabla con los errores verdaderos de las tres
soluciones de cada método."

EDO del Problema 12:
        dy/dx - (1 - y)(4 - y) = 0   ->   dy/dx = (1 - y)(4 - y)
        0 <= x <= 1 ,   y(0) = 0
        Exacta: y(x) = 4 (e^{3x} - 1) / (4 e^{3x} - 1)

RESTRICCIÓN DE LA PRÁCTICA: el Problema 12 pide Heun y Euler modificado;
aquí se reemplazan por EULER y RUNGE-KUTTA 4.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import (euler, rk4, imprimir_tabla, imprimir_convergencia,
                     ruta_grafica)

X0, XF, Y0 = 0.0, 1.0, 0.0
PASOS = [0.25, 0.1, 0.01]
METODOS = {"Euler": euler, "RK4": rk4}


def f(x, y):
    """dy/dx = (1 - y)(4 - y)   (no depende de x: EDO autónoma)."""
    return (1.0 - y) * (4.0 - y)


def exacta(x):
    """y(x) = 4 (e^{3x} - 1) / (4 e^{3x} - 1)."""
    e3x = np.exp(3.0 * x)
    return 4.0 * (e3x - 1.0) / (4.0 * e3x - 1.0)


def resolver():
    """Devuelve {h: (x, {"Euler": y_eu, "RK4": y_rk})}."""
    resultados = {}
    for h in PASOS:
        sol = {}
        x = None
        for nombre, metodo in METODOS.items():
            x, sol[nombre] = metodo(f, X0, Y0, h, XF)
        resultados[h] = (x, sol)
    return resultados


def imprimir_resultados(resultados):
    """Tabla de errores verdaderos por h (Euler y RK4 lado a lado) y
    resumen de convergencia de ambos métodos."""
    filas = []
    for h, (x, sol) in resultados.items():
        paso_tabla = None if h >= 0.1 else 0.1
        imprimir_tabla(x, sol, exacta(x), titulo=f"h = {h}",
                       paso_tabla=paso_tabla)
    for nombre in METODOS:
        for h, (x, sol) in resultados.items():
            err = np.abs(exacta(x) - sol[nombre])
            filas.append({"metodo": nombre, "h": h, "err_max": err.max(),
                          "err_final": err[-1]})
    imprimir_convergencia(filas, "Resumen de errores verdaderos (Euler vs RK4)")


def graficar(resultados):
    """Fila 1: soluciones de cada método para los tres h.
    Fila 2: error absoluto en escala log."""
    xx = np.linspace(X0, XF, 300)
    fig, ejes = plt.subplots(2, 2, figsize=(12, 8))
    for col, nombre in enumerate(METODOS):
        ax_sol, ax_err = ejes[0, col], ejes[1, col]
        ax_sol.plot(xx, exacta(xx), "k-", lw=2, label="Exacta")
        for h, (x, sol) in resultados.items():
            estilo = "o--" if len(x) < 30 else "-"
            ax_sol.plot(x, sol[nombre], estilo, ms=4, label=f"h = {h}")
            err = np.abs(exacta(x[1:]) - sol[nombre][1:])
            ax_err.semilogy(x[1:], np.maximum(err, 1e-17), estilo, ms=4,
                            label=f"h = {h}")
        ax_sol.set(xlabel="x", ylabel="y", title=f"{nombre}: solución")
        ax_err.set(xlabel="x", ylabel="|error|",
                   title=f"{nombre}: error verdadero (log)")
        for ax in (ax_sol, ax_err):
            ax.grid(alpha=0.3, which="both")
            ax.legend()
    fig.suptitle("Ej. 7.16:  dy/dx = (1 - y)(4 - y),  y(0) = 0")
    fig.tight_layout()
    ruta = ruta_grafica("cap7_ej16_euler_rk4.png")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 16: dy/dx = (1-y)(4-y), Euler y RK4, "
          "h = 0.25, 0.1, 0.01")
    print("=" * 70)
    resultados = resolver()
    imprimir_resultados(resultados)
    graficar(resultados)


if __name__ == "__main__":
    main()
