"""
Capítulo 8 - Ejercicio 4
========================
Enunciado: Modificar el programa de la Fig. 8.10 para resolver la
distribución de temperatura en la placa dividiéndola en 200 intervalos en x
y 100 en y. Comparar con la solución exacta de la Ec. (8.25).

Problema (Ejemplo 8.1 del libro, Fig. 8.7):
    Placa rectangular de 2 x 1. Ecuación de Laplace (conducción de calor en
    régimen permanente, sin fuentes):
            T_xx + T_yy = 0 ,      0 <= x <= 2 ,  0 <= y <= 1
    Condiciones de borde (Ec. 8.26):
            T(x, 0) = 0             borde inferior
            T(0, y) = 0             borde izquierdo
            T(2, y) = 0             borde derecho
            T(x, 1) = sin(pi x / 2) borde superior
    Solución exacta (Ec. 8.25):
            T(x, y) = sin(pi x / 2) * sinh(pi y / 2) / sinh(pi / 2)

El programa de la Fig. 8.10 usa una malla 8 x 4 (dx = dy = 0.25), semilla
inicial T(i, j) = j * T(i, superior) / 4 (interpolación lineal desde el borde
superior), Gauss-Seidel con tol = 1e-5 y máximo 100 iteraciones. Aquí:
    1) se reproduce ese caso (debe dar 24 iteraciones y la Tabla 8.1),
    2) se cambia a 200 x 100 (dx = dy = 0.01, 199 x 99 = 19 701 incógnitas),
       con Gauss-Seidel puro y con SOR.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metodos import gauss_seidel_laplace, omega_optimo, ruta_grafica

# ---------------------------- parámetros ------------------------------------
LX, LY = 2.0, 1.0        # dimensiones de la placa
NX, NY = 200, 100        # intervalos pedidos en el ejercicio
TOL = 1e-8               # tolerancia (cambio máximo entre iteraciones)
MAX_ITER = 100_000       # máximo de iteraciones

# Caso original del libro (Fig. 8.10)
NX_LIBRO, NY_LIBRO = 8, 4
TOL_LIBRO, MAX_ITER_LIBRO = 1e-5, 100


def exacta(x, y):
    """Ec. (8.25): T = sin(pi x/2) sinh(pi y/2) / sinh(pi/2)."""
    return np.sin(np.pi * x / 2) * np.sinh(np.pi * y / 2) / np.sinh(np.pi / 2)


def crear_malla(nx, ny):
    """Nodos y matrices X, Y con indexado [i, j] (i -> x, j -> y)."""
    x = np.linspace(0.0, LX, nx + 1)
    y = np.linspace(0.0, LY, ny + 1)
    X, Y = np.meshgrid(x, y, indexing="ij")
    return x, y, X, Y


def condiciones_y_semilla(x, y):
    """Carga las condiciones de borde (Ec. 8.26) y la semilla inicial del
    programa de la Fig. 8.10 en los nodos interiores."""
    T = np.zeros((len(x), len(y)))
    T[:, -1] = np.sin(np.pi * x / 2)     # borde superior: T(x,1) = sin(pi x/2)
    # Bordes inferior, izquierdo y derecho quedan en 0.
    # Semilla: interpolación lineal entre T = 0 abajo y el borde superior,
    # T(i, j) = (y_j / LY) * T(i, superior)  (la intención de j*t(i,5)/5 en
    # el libro; tal como está impreso, ese código lee una columna que aún vale
    # 0, así que en la práctica arranca desde cero. Por eso el número de
    # iteraciones puede diferir en 1-2 del "24" del libro; la solución final
    # es la misma porque no depende de la semilla.)
    T[1:-1, 1:-1] = (y[1:-1] / LY)[None, :] * T[1:-1, -1][:, None]
    T[0, :] = 0.0
    T[-1, :] = 0.0                       # sin(pi) ~ 1e-16 -> se fija en 0 exacto
    return T


def resolver(nx, ny, tol, max_iter, omega, ordenamiento="rojo-negro"):
    x, y, X, Y = crear_malla(nx, ny)
    dx, dy = LX / nx, LY / ny
    T0 = condiciones_y_semilla(x, y)
    T, it, ok, hist, t = gauss_seidel_laplace(T0, dx, dy, tol=tol,
                                              max_iter=max_iter, omega=omega,
                                              ordenamiento=ordenamiento)
    err = np.abs(T - exacta(X, Y)).max()
    return {"T": T, "X": X, "Y": Y, "it": it, "ok": ok, "hist": hist,
            "t": t, "err": err, "omega": omega, "nx": nx, "ny": ny}


def imprimir_tabla_81(r):
    """Imprime con el formato de la Tabla 8.1: nodos cada 0.25 (ubicaciones
    de la Fig. 8.9), numérico arriba y exacto entre paréntesis abajo."""
    nx, ny = r["nx"], r["ny"]
    kx, ky = int(round(0.25 / (LX / nx))), int(round(0.25 / (LY / ny)))
    T, X, Y = r["T"], r["X"], r["Y"]
    cols = range(0, nx + 1, kx)
    for j in range(ny, -1, -ky):
        print("  " + "".join(f"{T[i, j]:9.5f}" for i in cols)
              + f"     <- y = {Y[0, j]:.2f}")
        print("  " + "".join(f"({exacta(X[i, j], Y[i, j]):7.5f})" for i in cols))
    print("  " + "".join(f"{'x=' + format(X[i, 0], '.2f'):>9}" for i in cols))


def reproducir_libro():
    """Caso 8 x 4 de la Fig. 8.10 (Gauss-Seidel lexicográfico, como el libro)."""
    r = resolver(NX_LIBRO, NY_LIBRO, TOL_LIBRO, MAX_ITER_LIBRO, 1.0,
                 ordenamiento="lexicografico")
    print(f"\n1) Reproducción del programa de la Fig. 8.10: malla "
          f"{NX_LIBRO} x {NY_LIBRO}, tol = {TOL_LIBRO:g}")
    print(f"   Iteraciones: {r['it']} (el libro reporta 24)   "
          f"error máximo = {r['err']:.4e}")
    imprimir_tabla_81(r)
    return r


def imprimir_comparacion(resultados, r_libro):
    print(f"\n2) Malla {NX} x {NY} (dx = dy = {LX / NX:g}): "
          f"{(NX - 1) * (NY - 1)} incógnitas, tol = {TOL:g}")
    cab = (f"{'Método':>26} | {'omega':>7} | {'iteraciones':>11} | "
           f"{'tiempo [s]':>10} | {'error máx':>11}")
    print("-" * len(cab))
    print(cab)
    print("-" * len(cab))
    for nombre, r in resultados.items():
        estado = "" if r["ok"] else " (NO convergió)"
        print(f"{nombre:>26} | {r['omega']:7.4f} | {r['it']:11d} | "
              f"{r['t']:10.3f} | {r['err']:11.4e}{estado}")
    print("-" * len(cab))
    r = resultados["SOR (rojo-negro)"]
    h_libro, h_fino = LX / NX_LIBRO, LX / NX
    print(f"   Error máximo malla 8x4   : {r_libro['err']:.4e}  (h = {h_libro})")
    print(f"   Error máximo malla 200x100: {r['err']:.4e}  (h = {h_fino})")
    print(f"   Razón de errores = {r_libro['err'] / r['err']:.1f};  "
          f"(h_8x4 / h_200x100)^2 = {(h_libro / h_fino) ** 2:.1f}  -> "
          f"el error de la plantilla de 5 puntos es O(h^2)")


def graficar(resultados):
    r = resultados["SOR (rojo-negro)"]
    X, Y, T = r["X"], r["Y"], r["T"]

    # Figura 1: mapa de temperatura y mapa de error
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8))
    c = ax1.contourf(X, Y, T, levels=20, cmap="inferno")
    ax1.contour(X, Y, T, levels=10, colors="w", linewidths=0.5)
    fig.colorbar(c, ax=ax1, label="T")
    ax1.set(xlabel="x", ylabel="y", aspect="equal",
            title=f"Temperatura numérica (SOR), malla {NX} x {NY}")
    c = ax2.pcolormesh(X, Y, np.abs(T - exacta(X, Y)), cmap="viridis",
                       shading="auto")
    fig.colorbar(c, ax=ax2, label="|T - T exacta|")
    ax2.set(xlabel="x", ylabel="y", aspect="equal",
            title="Error absoluto respecto de la Ec. (8.25)")
    fig.tight_layout()
    ruta = ruta_grafica("cap8_ej04_placa.png")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    print(f"Gráfica guardada en: {ruta}")

    # Figura 2: superficie 3D (equivalente al 'mesh' de la Fig. 8.10) y
    # historial de convergencia
    fig = plt.figure(figsize=(14, 5))
    ax3 = fig.add_subplot(1, 2, 1, projection="3d")
    paso = 4   # se dibuja 1 de cada 4 nodos para que la malla se vea
    ax3.plot_surface(X[::paso, ::paso], Y[::paso, ::paso], T[::paso, ::paso],
                     cmap="inferno", edgecolor="k", linewidth=0.2)
    ax3.set(xlabel="x", ylabel="y", zlabel="T", title="Superficie de temperatura")
    ax4 = fig.add_subplot(1, 2, 2)
    for nombre, res in resultados.items():
        ax4.semilogy(res["hist"], label=f"{nombre} ({res['it']} it.)")
    ax4.axhline(TOL, color="gray", ls=":", label=f"tol = {TOL:g}")
    ax4.set(xlabel="iteración", ylabel="máx |T_nuevo - T_viejo|",
            title="Convergencia de Gauss-Seidel vs SOR")
    ax4.grid(alpha=0.3, which="both")
    ax4.legend()
    fig.tight_layout()
    ruta = ruta_grafica("cap8_ej04_superficie_convergencia.png")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    print(f"Gráfica guardada en: {ruta}")


def main():
    print("=" * 70)
    print("CAPÍTULO 8 - EJERCICIO 4: placa 2 x 1, malla 200 x 100, "
          "exacta Ec. (8.25)")
    print("=" * 70)
    r_libro = reproducir_libro()

    w = omega_optimo(NX, NY)
    print(f"\nomega óptimo de SOR para {NX}x{NY}: {w:.4f}")
    resultados = {
        "Gauss-Seidel (rojo-negro)": resolver(NX, NY, TOL, MAX_ITER, 1.0),
        "SOR (rojo-negro)": resolver(NX, NY, TOL, MAX_ITER, w),
    }
    imprimir_comparacion(resultados, r_libro)
    print("\n3) Malla 200 x 100 en los mismos nodos de la Fig. 8.9 "
          "(numérico / exacto):")
    imprimir_tabla_81(resultados["SOR (rojo-negro)"])
    graficar(resultados)


if __name__ == "__main__":
    main()
