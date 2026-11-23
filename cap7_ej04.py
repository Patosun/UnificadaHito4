"""
Capítulo 7 - Ejercicio 4
========================
Enunciado: Resolver la EDO del Problema 3 con un programa, usando
h = 0.2, 0.1 y 0.01, y comparar con la solución exacta.

        dy/dx = 1 + y/x + (y/x)^2 ,   1 <= x <= 3 ,   y(1) = 0
        Exacta: y(x) = x tan(ln x)

Método: Euler. Es una EDO NO LINEAL (término (y/x)^2), lo que no cambia el
método: Euler sólo necesita poder evaluar f(x, y).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import euler, imprimir_tabla, imprimir_convergencia, ruta_grafica

X0, XF, Y0 = 1.0, 3.0, 0.0
PASOS = [0.2, 0.1, 0.01]


def f(x, y):
    """dy/dx = 1 + y/x + (y/x)^2."""
    return 1.0 + y / x + (y / x) ** 2


def exacta(x):
    """y(x) = x tan(ln x)."""
    return x * np.tan(np.log(x))


def resolver_todos_los_pasos():
    """Ejecuta Euler para cada h y devuelve {h: (x, y)}."""
    return {h: euler(f, X0, Y0, h, XF) for h in PASOS}


def imprimir_resultados(resultados):
    filas = []
    for h, (x, y) in resultados.items():
        y_ex = exacta(x)
        paso_tabla = 0.2      # se imprimen los mismos x en las tres tablas
        imprimir_tabla(x, {"Euler": y}, y_ex, titulo=f"Euler con h = {h}",
                       paso_tabla=paso_tabla)
        err = np.abs(y_ex - y)
        filas.append({"metodo": "Euler", "h": h, "err_max": err.max(),
                      "err_final": err[-1]})
    imprimir_convergencia(filas)


def graficar(resultados):
    xx = np.linspace(X0, XF, 300)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(xx, exacta(xx), "k-", lw=2, label="Exacta")
    for h, (x, y) in resultados.items():
        estilo = "o--" if len(x) < 30 else "-"
        ax1.plot(x, y, estilo, ms=4, label=f"Euler h = {h}")
        ax2.semilogy(x[1:], np.abs(exacta(x[1:]) - y[1:]), estilo, ms=4,
                     label=f"h = {h}")
    ax1.set(xlabel="x", ylabel="y", title="Ej. 7.4: Euler vs exacta")
    ax2.set(xlabel="x", ylabel="|error|", title="Error absoluto (escala log)")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3, which="both")
        ax.legend()
    fig.tight_layout()
    ruta = ruta_grafica("cap7_ej04_euler.png")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 4: Euler, dy/dx = 1 + y/x + (y/x)^2, "
          "h = 0.2, 0.1, 0.01")
    print("=" * 70)
    resultados = resolver_todos_los_pasos()
    imprimir_resultados(resultados)
    graficar(resultados)


if __name__ == "__main__":
    main()
