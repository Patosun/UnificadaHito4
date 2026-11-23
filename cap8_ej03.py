"""
Capítulo 8 - Ejercicio 3
========================
Enunciado: Resolver el Problema 2 de nuevo, pero dividiendo el cuadrado
unitario en 100 x 100 intervalos y con un programa. Comparar las soluciones
en los nodos con la exacta.

Problema 2 (ecuación de Laplace, elíptica):
        u_xx + u_yy = 0 ,    0 <= x <= 1 ,  0 <= y <= 1
        u(x, 0) = 0 ,  u(x, 1) = x      (bordes inferior y superior)
        u(0, y) = 0 ,  u(1, y) = y      (bordes izquierdo y derecho)
        Exacta: u(x, y) = x y

Método: diferencias finitas centrales (plantilla de 5 puntos) + Gauss-Seidel
(y su versión acelerada SOR). Ver metodos.gauss_seidel_poisson.

Dato clave para la defensa: u = xy es un polinomio bilineal, sus derivadas
cuartas son nulas, así que el ERROR DE TRUNCAMIENTO de la plantilla de
5 puntos es CERO: la solución discreta coincide con la exacta en los nodos.
Todo el error que se observe es error de ITERACIÓN (detener Gauss-Seidel
antes de llegar al punto fijo) y de redondeo.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import gauss_seidel_laplace, omega_optimo, ruta_grafica

# ---------------------------- parámetros ------------------------------------
N = 100              # número de intervalos en x y en y
TOL = 1e-8           # tolerancia: máximo cambio entre dos iteraciones
MAX_ITER = 50_000    # número máximo de iteraciones
TOLS_ESTUDIO = [1e-4, 1e-6, 1e-8]   # para mostrar el efecto de la tolerancia


def exacta(x, y):
    """Solución exacta u(x, y) = x y."""
    return x * y


def crear_malla(n):
    """Nodos x_i = i*h, y_j = j*h (i, j = 0..n) y matrices X, Y con
    indexado [i, j] (i -> x, j -> y), igual que u_{i,j} en el libro."""
    x = np.linspace(0.0, 1.0, n + 1)
    y = np.linspace(0.0, 1.0, n + 1)
    X, Y = np.meshgrid(x, y, indexing="ij")
    return x, y, X, Y


def aplicar_condiciones_borde(x, y):
    """Matriz u con las condiciones de borde (Dirichlet) cargadas en las
    filas/columnas extremas y 0 (semilla inicial) en los nodos interiores."""
    u = np.zeros((len(x), len(y)))
    u[:, 0] = 0.0          # u(x, 0) = 0     (borde inferior,  j = 0)
    u[:, -1] = x           # u(x, 1) = x     (borde superior,  j = n)
    u[0, :] = 0.0          # u(0, y) = 0     (borde izquierdo, i = 0)
    u[-1, :] = y           # u(1, y) = y     (borde derecho,   i = n)
    return u


def resolver(n, tol, omega, ordenamiento="rojo-negro"):
    """Arma la malla, aplica los bordes y llama a Gauss-Seidel/SOR."""
    x, y, X, Y = crear_malla(n)
    h = 1.0 / n
    u0 = aplicar_condiciones_borde(x, y)
    u, it, ok, hist, t = gauss_seidel_laplace(u0, h, h, tol=tol,
                                              max_iter=MAX_ITER, omega=omega,
                                              ordenamiento=ordenamiento)
    err = np.abs(u - exacta(X, Y)).max()
    return {"u": u, "X": X, "Y": Y, "it": it, "ok": ok, "hist": hist,
            "t": t, "err": err, "omega": omega, "orden": ordenamiento}


def validar_con_ejercicio_2():
    """Malla 3x3 del Ejercicio 2 (resuelto a mano): 4 incógnitas interiores.
    Sirve para comprobar el programa contra el cálculo manual."""
    r = resolver(3, 1e-12, 1.0, ordenamiento="lexicografico")
    print("\nValidación con la malla 3x3 del Ejercicio 2 (Gauss-Seidel "
          "lexicográfico, tol = 1e-12):")
    print(f"  Iteraciones: {r['it']}")
    for i in (1, 2):
        for j in (1, 2):
            xv, yv = r["X"][i, j], r["Y"][i, j]
            print(f"  u({xv:.4f}, {yv:.4f}) = {r['u'][i, j]:.10f}   "
                  f"exacta = {exacta(xv, yv):.10f}")


def imprimir_comparacion_metodos(resultados):
    print(f"\nMalla {N} x {N}: {(N - 1) ** 2} incógnitas, tol = {TOL:g}")
    cab = (f"{'Método':>28} | {'omega':>7} | {'iteraciones':>11} | "
           f"{'tiempo [s]':>10} | {'error máx |u - xy|':>18}")
    print("-" * len(cab))
    print(cab)
    print("-" * len(cab))
    for nombre, r in resultados.items():
        estado = "" if r["ok"] else " (NO convergió)"
        print(f"{nombre:>28} | {r['omega']:7.4f} | {r['it']:11d} | "
              f"{r['t']:10.3f} | {r['err']:18.4e}{estado}")
    print("-" * len(cab))


def imprimir_tabla_nodos(r, paso=0.2):
    """Tabla tipo 'Tabla 8.1' del libro: valores numéricos y entre
    paréntesis la exacta, en nodos cada 'paso' (filas = y de arriba a abajo)."""
    k = int(round(paso * N))
    u, X, Y = r["u"], r["X"], r["Y"]
    print(f"\nSolución numérica (exacta entre paréntesis) en nodos cada {paso}:")
    print("  y \\ x " + "".join(f"{X[i, 0]:>19.2f}" for i in range(0, N + 1, k)))
    for j in range(N, -1, -k):
        fila = f"  {Y[0, j]:5.2f}"
        for i in range(0, N + 1, k):
            fila += f"  {u[i, j]:8.6f} ({exacta(X[i, j], Y[i, j]):6.4f})"
        print(fila)


def estudiar_tolerancia():
    """Muestra que el criterio 'cambio entre iteraciones < tol' NO garantiza
    error < tol cuando la convergencia es lenta (Gauss-Seidel puro)."""
    print("\nEfecto de la tolerancia (Gauss-Seidel puro, omega = 1):")
    cab = f"{'tol':>8} | {'iteraciones':>11} | {'error máx':>12}"
    print("-" * len(cab))
    print(cab)
    print("-" * len(cab))
    for tol in TOLS_ESTUDIO:
        r = resolver(N, tol, 1.0)
        print(f"{tol:8.0e} | {r['it']:11d} | {r['err']:12.4e}")
    print("-" * len(cab))
    print("  El error real es ~ tol / (1 - rho), con rho ~ 0.999 el radio "
          "espectral de Gauss-Seidel:\n  por eso hay que pedir una tolerancia "
          "bastante menor que el error que se desea.")


def graficar(resultados):
    r = resultados["SOR (rojo-negro)"]
    X, Y, u = r["X"], r["Y"], r["u"]
    fig, ejes = plt.subplots(1, 3, figsize=(17, 4.8))

    # 1) Mapa de color de la solución numérica + isolíneas
    c = ejes[0].contourf(X, Y, u, levels=20, cmap="viridis")
    ejes[0].contour(X, Y, u, levels=10, colors="k", linewidths=0.5)
    fig.colorbar(c, ax=ejes[0], label="u")
    ejes[0].set(xlabel="x", ylabel="y", aspect="equal",
                title=f"Solución numérica (SOR), malla {N}x{N}")

    # 2) Mapa del error absoluto
    c = ejes[1].pcolormesh(X, Y, np.abs(u - exacta(X, Y)), cmap="magma",
                           shading="auto")
    fig.colorbar(c, ax=ejes[1], label="|u - xy|")
    ejes[1].set(xlabel="x", ylabel="y", aspect="equal",
                title="Error absoluto |u - xy|")

    # 3) Historial de convergencia
    for nombre, res in resultados.items():
        ejes[2].semilogy(res["hist"], label=f"{nombre} ({res['it']} it.)")
    ejes[2].axhline(TOL, color="gray", ls=":", label=f"tol = {TOL:g}")
    ejes[2].set(xlabel="iteración", ylabel="máx |u_nuevo - u_viejo|",
                title="Convergencia")
    ejes[2].grid(alpha=0.3, which="both")
    ejes[2].legend(fontsize=8)

    fig.suptitle("Ej. 8.3: Laplace en el cuadrado unitario, exacta u = xy")
    fig.tight_layout()
    ruta = ruta_grafica("cap8_ej03_laplace.png")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 8 - EJERCICIO 3: Laplace, malla 100 x 100, exacta u = xy")
    print("=" * 70)
    validar_con_ejercicio_2()

    w = omega_optimo(N, N)
    print(f"\nomega óptimo de SOR para {N}x{N}: {w:.4f}")
    resultados = {
        "Gauss-Seidel (rojo-negro)": resolver(N, TOL, 1.0),
        "SOR (rojo-negro)": resolver(N, TOL, w),
        "SOR (lexicográfico, libro)": resolver(N, TOL, w, "lexicografico"),
    }
    imprimir_comparacion_metodos(resultados)
    imprimir_tabla_nodos(resultados["SOR (rojo-negro)"])
    estudiar_tolerancia()
    graficar(resultados)


if __name__ == "__main__":
    main()
