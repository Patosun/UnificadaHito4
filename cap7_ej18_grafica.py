"""
Capítulo 7 - Ejercicio 18 (complemento: sólo la gráfica)
========================================================
Enunciado: Usar Runge-Kutta de 4º orden para resolver

        dy/dx + y^2 + y/x = 1/x^2 ,     1 <= x <= 2 ,     y(1) = -1

con h = 0.25. Graficar junto con la solución exacta  y(x) = -1/x.

Primero se despeja la derivada:  dy/dx = 1/x^2 - y^2 - y/x.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import rk4, imprimir_tabla, ruta_grafica

X0, XF, Y0 = 1.0, 2.0, -1.0
H = 0.25


def f(x, y):
    """dy/dx = 1/x^2 - y^2 - y/x."""
    return 1.0 / x ** 2 - y ** 2 - y / x


def exacta(x):
    """y(x) = -1/x."""
    return -1.0 / x


def graficar(x, y):
    xx = np.linspace(X0, XF, 300)
    plt.figure(figsize=(7, 4.5))
    plt.plot(xx, exacta(xx), "k-", label="Exacta: $-1/x$")
    plt.plot(x, y, "s--", color="tab:blue", label=f"RK4, h = {H}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Ej. 7.18:  dy/dx = 1/x$^2$ - y$^2$ - y/x,  y(1) = -1")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    ruta = ruta_grafica("cap7_ej18_rk4.png")
    plt.savefig(ruta, dpi=150)
    plt.close()
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 18 (gráfica): RK4 con h = 0.25")
    print("=" * 70)
    x, y = rk4(f, X0, Y0, H, XF)
    imprimir_tabla(x, {"RK4": y}, exacta(x), titulo="Tabla de verificación:")
    graficar(x, y)


if __name__ == "__main__":
    main()
