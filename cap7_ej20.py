"""
Capítulo 7 - Ejercicio 20
=========================
Enunciado: Resolver

        dy/dx - (2/x) y = x^2 e^x   ->   dy/dx = x^2 e^x + 2y/x
        1 <= x <= 2 ,   y(1) = 0

con h = 0.1 y comparar en una tabla con la exacta  y(x) = x^2 (e^x - e).

RESTRICCIÓN DE LA PRÁCTICA: el libro pide RK2, RK3 y RK4; aquí sólo RK4.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import rk4, imprimir_tabla, ruta_grafica

X0, XF, Y0 = 1.0, 2.0, 0.0
H = 0.1


def f(x, y):
    """dy/dx = x^2 e^x + 2y/x."""
    return x ** 2 * np.exp(x) + 2.0 * y / x


def exacta(x):
    """y(x) = x^2 (e^x - e)."""
    return x ** 2 * (np.exp(x) - np.e)


def graficar(x, y):
    xx = np.linspace(X0, XF, 300)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(xx, exacta(xx), "k-", lw=2, label="Exacta: $x^2(e^x - e)$")
    ax1.plot(x, y, "s", color="tab:blue", label=f"RK4, h = {H}")
    ax2.semilogy(x[1:], np.abs(exacta(x[1:]) - y[1:]), "o-", color="tab:blue")
    ax1.set(xlabel="x", ylabel="y", title="Ej. 7.20: RK4 vs exacta")
    ax2.set(xlabel="x", ylabel="|error|", title="Error verdadero (log)")
    ax1.legend()
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    ruta = ruta_grafica("cap7_ej20_rk4.png")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 20: RK4, dy/dx = x^2 e^x + 2y/x, h = 0.1")
    print("=" * 70)
    x, y = rk4(f, X0, Y0, H, XF)
    imprimir_tabla(x, {"RK4": y}, exacta(x), titulo=f"RK4 con h = {H}")
    graficar(x, y)


if __name__ == "__main__":
    main()
