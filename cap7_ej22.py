"""
Capítulo 7 - Ejercicio 22
=========================
Enunciado: Usar Euler y RK4 para resolver

        x dy/dx - 4y = x^5 e^x   ->   dy/dx = x^4 e^x + 4y/x
        1 <= x <= 2 ,   y(1) = 0

ELIGIENDO un paso h apropiado, y comparar en una tabla con la exacta
        y(x) = x^4 (e^x - e).

¿Cómo se elige h? (justificación)
---------------------------------
Se fija un criterio de precisión: error relativo máximo < TOL_PORCENTAJE (1 %)
en todo el intervalo. Se prueban pasos cada vez más chicos (lista
CANDIDATOS) y para cada método se toma el MAYOR h que cumple el criterio
(el más barato que alcanza la precisión pedida). Como la tabla compara ambos
métodos en los mismos x, se usa el h que exige el método menos preciso
(Euler); con ese h RK4 queda prácticamente exacto.

Nota: el error relativo es grande cerca de x = 1 porque allí y(1) = 0 y la
solución es muy pequeña (dividir por un número chico amplifica el error).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import euler, rk4, imprimir_tabla, error_porcentual, ruta_grafica

X0, XF, Y0 = 1.0, 2.0, 0.0
TOL_PORCENTAJE = 1.0
CANDIDATOS = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]
METODOS = {"Euler": euler, "RK4": rk4}


def f(x, y):
    """dy/dx = x^4 e^x + 4y/x."""
    return x ** 4 * np.exp(x) + 4.0 * y / x


def exacta(x):
    """y(x) = x^4 (e^x - e)."""
    return x ** 4 * (np.exp(x) - np.e)


def error_max_porcentual(metodo, h):
    """Máximo error relativo (%) en todo [X0, XF] (se ignora x = 1, donde
    la exacta vale 0 y el error relativo no está definido)."""
    x, y = metodo(f, X0, Y0, h, XF)
    return np.nanmax(error_porcentual(y, exacta(x)))


def estudiar_pasos():
    """Tabla 'error máximo % vs h' para ambos métodos y elección de h."""
    print(f"\nEstudio de convergencia para elegir h "
          f"(criterio: error relativo máximo < {TOL_PORCENTAJE} %; * = cumple)")
    cab = f"{'h':>8} | {'err máx % Euler':>16} | {'err máx % RK4':>16}"
    print("-" * len(cab))
    print(cab)
    print("-" * len(cab))
    h_elegido = {}
    for h in CANDIDATOS:
        errores = {n: error_max_porcentual(m, h) for n, m in METODOS.items()}
        marca = {n: (" *" if e < TOL_PORCENTAJE else "  ") for n, e in errores.items()}
        print(f"{h:8.4g} | {errores['Euler']:14.4e}{marca['Euler']} | "
              f"{errores['RK4']:14.4e}{marca['RK4']}")
        for n, e in errores.items():
            if e < TOL_PORCENTAJE and n not in h_elegido:
                h_elegido[n] = h        # el primero (mayor h) que cumple
    print("-" * len(cab))
    for n in METODOS:
        print(f"  Mayor h que cumple el criterio con {n:5}: "
              f"{h_elegido.get(n, 'ninguno de los candidatos')}")
    h_tabla = min(h_elegido.values())
    print(f"  => h elegido para la tabla comparativa: h = {h_tabla}")
    return h_tabla, h_elegido


def graficar(h, x, soluciones, h_elegido):
    xx = np.linspace(X0, XF, 300)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(xx, exacta(xx), "k-", lw=2, label="Exacta")
    # RK4 con SU h óptimo (mucho más grande) para mostrar su ventaja
    h_rk = h_elegido["RK4"]
    x_rk, y_rk = rk4(f, X0, Y0, h_rk, XF)
    ax1.plot(x, soluciones["Euler"], "--", color="tab:red", label=f"Euler h = {h}")
    ax1.plot(x_rk, y_rk, "s", color="tab:blue", label=f"RK4 h = {h_rk}")
    ax1.set(xlabel="x", ylabel="y", title="Ej. 7.22: soluciones")
    # Error máximo % vs h (log-log)
    hs = np.array(CANDIDATOS)
    for n, m in METODOS.items():
        ax2.loglog(hs, [error_max_porcentual(m, hh) for hh in hs], "o-", label=n)
    ax2.axhline(TOL_PORCENTAJE, color="gray", ls=":", label=f"criterio {TOL_PORCENTAJE} %")
    ax2.set(xlabel="h", ylabel="error relativo máximo [%]",
            title="Elección de h")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3, which="both")
        ax.legend()
    fig.tight_layout()
    ruta = ruta_grafica("cap7_ej22_euler_rk4.png")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 7 - EJERCICIO 22: dy/dx = x^4 e^x + 4y/x, Euler y RK4, "
          "h elegido")
    print("=" * 70)
    h, h_elegido = estudiar_pasos()
    soluciones = {}
    x = None
    for nombre, metodo in METODOS.items():
        x, soluciones[nombre] = metodo(f, X0, Y0, h, XF)
    imprimir_tabla(x, soluciones, exacta(x),
                   titulo=f"Tabla comparativa con h = {h} (se muestran x cada 0.1)",
                   paso_tabla=0.1)
    graficar(h, x, soluciones, h_elegido)


if __name__ == "__main__":
    main()
