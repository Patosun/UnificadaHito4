"""
Capítulo 7 - Ejercicio 3 (complemento: sólo la gráfica)
=======================================================
Enunciado: Usar el método de Euler para resolver

        dy/dx = 1 + y/x + (y/x)^2 ,     1 <= x <= 3 ,     y(1) = 0

con h = 0.2. Graficar junto con la solución exacta  y(x) = x tan(ln x).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import euler, imprimir_tabla, ruta_grafica

X0, XF, Y0 = 1.0, 3.0, 0.0
H = 0.2


def f(x, y):
    """dy/dx = 1 + y/x + (y/x)^2  (EDO no lineal)."""
    return 1.0 + y / x + (y / x) ** 2


def exacta(x):
    """y(x) = x tan(ln x)."""
    return x * np.tan(np.log(x))


def graficar(x, y):
    xx = np.linspace(X0, XF, 300)
    plt.figure(figsize=(7, 4.5))
    plt.plot(xx, exacta(xx), "k-", label="Exacta: $x\\tan(\\ln x)$")
    plt.plot(x, y, "o--", color="tab:red", label=f"Euler, h = {H}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Ej. 7.3:  dy/dx = 1 + y/x + (y/x)$^2$,  y(1) = 0")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    ruta = ruta_grafica("cap7_ej03_euler.png")
    plt.savefig(ruta, dpi=150)
    plt.close()
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 3 (gráfica): Euler con h = 0.2")
    print("=" * 70)
    x, y = euler(f, X0, Y0, H, XF)
    imprimir_tabla(x, {"Euler": y}, exacta(x), titulo="Tabla de verificación:")
    graficar(x, y)


if __name__ == "__main__":
    main()
