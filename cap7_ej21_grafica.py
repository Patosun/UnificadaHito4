"""
Capítulo 7 - Ejercicio 21 (complemento: sólo la gráfica)
========================================================
Enunciado: Usar Euler y Runge-Kutta de 4º orden para resolver

        dy/dx = e^x / y ,     0 <= x <= 2 ,     y(0) = 1

con h = 0.5. Graficar ambas soluciones junto con la exacta
        y(x) = sqrt(2 e^x - 1).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import euler, rk4, imprimir_tabla, ruta_grafica

X0, XF, Y0 = 0.0, 2.0, 1.0
H = 0.5


def f(x, y):
    """dy/dx = e^x / y."""
    return np.exp(x) / y


def exacta(x):
    """y(x) = sqrt(2 e^x - 1)."""
    return np.sqrt(2.0 * np.exp(x) - 1.0)


def graficar(x, y_eu, y_rk):
    xx = np.linspace(X0, XF, 300)
    plt.figure(figsize=(7, 4.5))
    plt.plot(xx, exacta(xx), "k-", label="Exacta: $\\sqrt{2e^x-1}$")
    plt.plot(x, y_eu, "o--", color="tab:red", label=f"Euler, h = {H}")
    plt.plot(x, y_rk, "s--", color="tab:blue", label=f"RK4, h = {H}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Ej. 7.21:  dy/dx = e$^x$ / y,  y(0) = 1")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    ruta = ruta_grafica("cap7_ej21_euler_rk4.png")
    plt.savefig(ruta, dpi=150)
    plt.close()
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 21 (gráfica): Euler y RK4 con h = 0.5")
    print("=" * 70)
    x, y_eu = euler(f, X0, Y0, H, XF)
    _, y_rk = rk4(f, X0, Y0, H, XF)
    imprimir_tabla(x, {"Euler": y_eu, "RK4": y_rk}, exacta(x),
                   titulo="Tabla de verificación:")
    graficar(x, y_eu, y_rk)


if __name__ == "__main__":
    main()
