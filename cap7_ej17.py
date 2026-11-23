"""
Capítulo 7 - Ejercicio 17
=========================
Enunciado: Usar Runge-Kutta de 4º orden (programa) para resolver

        dy/dx - y = x^2   ->   dy/dx = y + x^2
        1 <= x <= 2 ,   y(1) = 1

con h = 0.1, 0.01, 0.001 y 0.0001, y comparar en una tabla con la exacta
        y(x) = 6 e^{x-1} - x^2 - 2x - 2.

Observación importante para la defensa: RK4 es O(h^4), así que al dividir h
por 10 el error debería caer ~10^4 veces. Pero la aritmética de doble
precisión tiene ~16 cifras significativas: cuando el error de truncamiento
llega a ~1e-13 / 1e-14, el ERROR DE REDONDEO (que crece con el número de
pasos) domina y el error deja de bajar o incluso sube.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import rk4, imprimir_convergencia, indices_tabla, ruta_grafica

X0, XF, Y0 = 1.0, 2.0, 1.0
PASOS = [0.1, 0.01, 0.001, 0.0001]


def f(x, y):
    """dy/dx = y + x^2."""
    return y + x ** 2


def exacta(x):
    """y(x) = 6 e^{x-1} - x^2 - 2x - 2."""
    return 6.0 * np.exp(x - 1.0) - x ** 2 - 2.0 * x - 2.0


def resolver():
    """Devuelve {h: (x, y)} usando RK4."""
    return {h: rk4(f, X0, Y0, h, XF) for h in PASOS}


def imprimir_tabla_comparativa(resultados):
    """Una sola tabla: en cada x = 1.0, 1.1, ..., 2.0 se muestra la exacta
    y el error relativo (%) de RK4 con cada h (los cuatro h comparten
    esos nodos)."""
    cab = f"{'x':>5} | {'y exacta':>16} |"
    for h in PASOS:
        cab += f" {'y RK4 h=' + str(h):>18} |"
    print("\nSoluciones RK4 vs exacta")
    print("-" * len(cab))
    print(cab)
    print("-" * len(cab))
    # Se toman los nodos de cada solución que coinciden con x = 1.0, 1.1, ...
    idx = {h: indices_tabla(x, 0.1) for h, (x, _) in resultados.items()}
    x_tabla = resultados[PASOS[0]][0]
    for k in range(len(x_tabla)):
        xk = x_tabla[k]
        fila = f"{xk:5.2f} | {exacta(xk):16.12f} |"
        for h in PASOS:
            _, y = resultados[h]
            fila += f" {y[idx[h][k]]:18.12f} |"
        print(fila)
    print("-" * len(cab))

    cab = f"{'x':>5} |" + "".join(f" {'err% h=' + str(h):>14} |" for h in PASOS)
    print("\nError relativo verdadero (%) = |exacta - RK4| / |exacta| * 100")
    print("-" * len(cab))
    print(cab)
    print("-" * len(cab))
    for k in range(len(x_tabla)):
        xk = x_tabla[k]
        fila = f"{xk:5.2f} |"
        for h in PASOS:
            _, y = resultados[h]
            yk = y[idx[h][k]]
            e = abs(exacta(xk) - yk) / abs(exacta(xk)) * 100
            fila += f" {e:14.4e} |"
        print(fila)
    print("-" * len(cab))


def imprimir_resumen(resultados):
    filas = []
    for h, (x, y) in resultados.items():
        err = np.abs(exacta(x) - y)
        filas.append({"metodo": "RK4", "h": h, "err_max": err.max(),
                      "err_final": err[-1]})
    imprimir_convergencia(filas)


def graficar(resultados):
    xx = np.linspace(X0, XF, 300)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(xx, exacta(xx), "k-", lw=2, label="Exacta")
    x01, y01 = resultados[0.1]
    ax1.plot(x01, y01, "s", color="tab:blue", label="RK4 h = 0.1")
    ax1.set(xlabel="x", ylabel="y", title="Ej. 7.17: RK4 vs exacta")
    # Error máximo vs h en escala log-log: la pendiente es el orden
    hs = np.array(PASOS)
    errores = np.array([np.abs(exacta(x) - y).max() for x, y in resultados.values()])
    ax2.loglog(hs, errores, "o-", label="error máximo RK4")
    ax2.loglog(hs, errores[0] * (hs / hs[0]) ** 4, "k--", label="referencia $O(h^4)$")
    ax2.set(xlabel="h", ylabel="error máximo absoluto",
            title="Convergencia (log-log)")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3, which="both")
        ax.legend()
    fig.tight_layout()
    ruta = ruta_grafica("cap7_ej17_rk4.png")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 17: RK4, dy/dx = y + x^2, "
          "h = 0.1, 0.01, 0.001, 0.0001")
    print("=" * 70)
    resultados = resolver()
    imprimir_tabla_comparativa(resultados)
    imprimir_resumen(resultados)
    graficar(resultados)


if __name__ == "__main__":
    main()
