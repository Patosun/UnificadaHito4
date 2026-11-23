"""
Capítulo 7 - Ejercicio 1 (complemento: sólo la gráfica)
=======================================================
Enunciado: Usar el método de Euler para resolver

        dy/dx = 1 + y/x ,     1 <= x <= 2 ,     y(1) = 2

con h = 0.25. Graficar la solución calculada junto con la exacta
        y(x) = x ln(x) + 2x.

La tabla se resolvió a mano; aquí se genera la gráfica comparativa y se
imprime la tabla para poder verificar los cálculos manuales.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")                 # backend sin ventanas: sólo guarda PNG
import matplotlib.pyplot as plt

from metodos import euler, imprimir_tabla, ruta_grafica

# Datos del problema
X0, XF, Y0 = 1.0, 2.0, 2.0
H = 0.25


def f(x, y):
    """Lado derecho de la EDO: dy/dx = f(x, y)."""
    return 1.0 + y / x


def exacta(x):
    """Solución exacta y(x) = x ln x + 2x."""
    return x * np.log(x) + 2.0 * x


def graficar(x, y):
    """Curva exacta (continua) vs puntos de Euler."""
    xx = np.linspace(X0, XF, 300)
    plt.figure(figsize=(7, 4.5))
    plt.plot(xx, exacta(xx), "k-", label="Exacta: $x\\ln x + 2x$")
    plt.plot(x, y, "o--", color="tab:red", label=f"Euler, h = {H}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Ej. 7.1:  dy/dx = 1 + y/x,  y(1) = 2")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    ruta = ruta_grafica("cap7_ej01_euler.png")
    plt.savefig(ruta, dpi=150)
    plt.close()
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 1 (gráfica): Euler con h = 0.25")
    print("=" * 70)
    x, y = euler(f, X0, Y0, H, XF)
    imprimir_tabla(x, {"Euler": y}, exacta(x), titulo="Tabla de verificación:")
    graficar(x, y)


if __name__ == "__main__":
    main()
