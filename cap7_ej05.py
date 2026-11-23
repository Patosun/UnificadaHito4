"""
Capítulo 7 - Ejercicio 5
========================
Enunciado: Usar Euler (con un programa) para resolver

        dy/dx = -y + x + 1 ,   0 <= x <= 5 ,   y(0) = 1

con h = 0.5, 0.1 y 0.01, y comparar con la exacta  y(x) = e^{-x} + x.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import euler, imprimir_tabla, imprimir_convergencia, ruta_grafica

X0, XF, Y0 = 0.0, 5.0, 1.0
PASOS = [0.5, 0.1, 0.01]


def f(x, y):
    """dy/dx = -y + x + 1."""
    return -y + x + 1.0


def exacta(x):
    """y(x) = e^{-x} + x."""
    return np.exp(-x) + x


def resolver_todos_los_pasos():
    return {h: euler(f, X0, Y0, h, XF) for h in PASOS}


def imprimir_resultados(resultados):
    filas = []
    for h, (x, y) in resultados.items():
        y_ex = exacta(x)
        imprimir_tabla(x, {"Euler": y}, y_ex, titulo=f"Euler con h = {h}",
                       paso_tabla=0.5)   # mismos x (cada 0.5) en todas las tablas
        err = np.abs(y_ex - y)
        filas.append({"metodo": "Euler", "h": h, "err_max": err.max(),
                      "err_final": err[-1]})
    imprimir_convergencia(filas)


def graficar(resultados):
    xx = np.linspace(X0, XF, 400)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(xx, exacta(xx), "k-", lw=2, label="Exacta")
    for h, (x, y) in resultados.items():
        estilo = "o--" if len(x) < 30 else "-"
        ax1.plot(x, y, estilo, ms=4, label=f"Euler h = {h}")
        ax2.semilogy(x[1:], np.abs(exacta(x[1:]) - y[1:]), estilo, ms=4,
                     label=f"h = {h}")
    ax1.set(xlabel="x", ylabel="y", title="Ej. 7.5: Euler vs exacta")
    ax2.set(xlabel="x", ylabel="|error|", title="Error absoluto (escala log)")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3, which="both")
        ax.legend()
    fig.tight_layout()
    ruta = ruta_grafica("cap7_ej05_euler.png")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 5: Euler, dy/dx = -y + x + 1, h = 0.5, 0.1, 0.01")
    print("=" * 70)
    resultados = resolver_todos_los_pasos()
    imprimir_resultados(resultados)
    graficar(resultados)


if __name__ == "__main__":
    main()
