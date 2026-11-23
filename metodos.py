"""
metodos.py
==========
Módulo común con los métodos numéricos usados en la práctica.

    * euler(f, x0, y0, h, xf)            -> Método de Euler (orden 1)
    * rk4(f, x0, y0, h, xf)              -> Runge-Kutta de 4º orden (orden 4)
    * gauss_seidel_poisson(...)          -> Gauss-Seidel / SOR para la ecuación
                                            de Poisson  u_xx + u_yy = f(x, y)
                                            con plantilla de 5 puntos
    * gauss_seidel_laplace(...)          -> caso particular f = 0 (Laplace)
    * omega_optimo(nx, ny)               -> factor de sobre-relajación óptimo
    * utilidades para imprimir tablas y crear la carpeta de gráficas

Referencia: P. Dechaumphai & N. Wansophark, "Numerical Methods in Science and
Engineering", Alpha Science, 2022. Capítulos 7 y 8.
"""

import os
import time

import numpy as np

# Carpeta donde todos los scripts guardan sus gráficas
CARPETA_GRAFICAS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "graficas")


# =============================================================================
#  CAPÍTULO 7: ECUACIONES DIFERENCIALES ORDINARIAS  dy/dx = f(x, y)
# =============================================================================

def _numero_de_pasos(x0, h, xf):
    """Número de pasos n tal que x0 + n*h = xf.

    Se usa round() para evitar errores de redondeo en coma flotante
    (por ejemplo (2 - 1)/0.1 = 9.999999999 en vez de 10).
    """
    n = int(round((xf - x0) / h))
    if n <= 0:
        raise ValueError("El paso h debe ser positivo y menor que (xf - x0).")
    if abs(x0 + n * h - xf) > 1e-9 * max(1.0, abs(xf)):
        raise ValueError("(xf - x0) debe ser múltiplo de h.")
    return n


def euler(f, x0, y0, h, xf):
    """Método de Euler (explícito) para  dy/dx = f(x, y),  y(x0) = y0.

    Idea: la serie de Taylor  y(x+h) = y(x) + h*y'(x) + h^2/2*y''(x) + ...
    se corta después del término lineal:

        y_{i+1} = y_i + h * f(x_i, y_i)

    Es decir, se avanza siguiendo la recta tangente en el punto actual.
    - Error local (por paso)       : O(h^2)
    - Error global (acumulado)     : O(h)    -> si h se divide por 10, el
                                                 error se divide por ~10.

    Devuelve dos arreglos numpy: x (nodos) e y (solución aproximada).
    """
    n = _numero_de_pasos(x0, h, xf)
    x = x0 + h * np.arange(n + 1)          # nodos x_0, x_1, ..., x_n
    y = np.empty(n + 1)
    y[0] = y0                              # condición inicial
    for i in range(n):
        pendiente = f(x[i], y[i])          # pendiente al inicio del intervalo
        y[i + 1] = y[i] + h * pendiente    # avanzar sobre la recta tangente
    return x, y


def rk4(f, x0, y0, h, xf):
    """Método de Runge-Kutta clásico de 4º orden para dy/dx = f(x, y).

    En cada paso se evalúan 4 pendientes dentro del intervalo [x_i, x_i + h]:

        k1 = f(x_i,         y_i)                 pendiente al inicio
        k2 = f(x_i + h/2,   y_i + h/2 * k1)      pendiente en el punto medio (usando k1)
        k3 = f(x_i + h/2,   y_i + h/2 * k2)      pendiente en el punto medio (usando k2)
        k4 = f(x_i + h,     y_i + h   * k3)      pendiente al final

    y se avanza con un promedio ponderado (pesos 1-2-2-1, como la regla de
    Simpson):

        y_{i+1} = y_i + h/6 * (k1 + 2*k2 + 2*k3 + k4)

    Los coeficientes se eligen para que coincida con la serie de Taylor hasta
    el término h^4.
    - Error local  : O(h^5)
    - Error global : O(h^4)  -> si h se divide por 10, el error se divide
                                por ~10 000.

    Devuelve dos arreglos numpy: x (nodos) e y (solución aproximada).
    """
    n = _numero_de_pasos(x0, h, xf)
    x = x0 + h * np.arange(n + 1)
    y = np.empty(n + 1)
    y[0] = y0
    for i in range(n):
        xi, yi = x[i], y[i]
        k1 = f(xi, yi)
        k2 = f(xi + h / 2, yi + h / 2 * k1)
        k3 = f(xi + h / 2, yi + h / 2 * k2)
        k4 = f(xi + h, yi + h * k3)
        y[i + 1] = yi + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return x, y


# =============================================================================
#  CAPÍTULO 8: ECUACIÓN ELÍPTICA (LAPLACE / POISSON) - DIFERENCIAS FINITAS
# =============================================================================

def omega_optimo(nx, ny):
    """Factor de sobre-relajación óptimo para SOR en un rectángulo con
    condiciones de Dirichlet y dx = dy (nx, ny = número de intervalos).

        rho_J   = [cos(pi/nx) + cos(pi/ny)] / 2   (radio espectral de Jacobi)
        w_opt   = 2 / (1 + sqrt(1 - rho_J^2))

    Con w = 1 se recupera Gauss-Seidel puro.
    """
    rho_j = 0.5 * (np.cos(np.pi / nx) + np.cos(np.pi / ny))
    return 2.0 / (1.0 + np.sqrt(1.0 - rho_j ** 2))


def gauss_seidel_poisson(u, dx, dy, f=None, tol=1e-6, max_iter=100_000,
                         omega=1.0, ordenamiento="rojo-negro"):
    """Resuelve  u_xx + u_yy = f(x, y)  con diferencias finitas centrales
    (plantilla de 5 puntos) e iteración de Gauss-Seidel (o SOR si omega > 1).

    Parámetros
    ----------
    u : ndarray de forma (nx+1, ny+1), indexado u[i, j] con i -> x, j -> y
        (igual que T_{i,j} en el libro). Las filas/columnas de los BORDES ya
        deben contener las condiciones de borde (Dirichlet). Los valores
        interiores son la "semilla" inicial de la iteración. Se modifica
        in-place y también se devuelve.
    dx, dy : tamaño de los intervalos en x e y.
    f : ndarray (nx+1, ny+1) con el término fuente, o None para Laplace (f=0).
    tol : tolerancia. Se detiene cuando el máximo cambio |u_nuevo - u_viejo|
          en una iteración completa es menor que tol (criterio del libro).
    max_iter : número máximo de iteraciones.
    omega : factor de relajación. 1 = Gauss-Seidel, 1 < omega < 2 = SOR.
    ordenamiento : "rojo-negro" (vectorizado con numpy, rápido) o
                   "lexicografico" (doble bucle i, j como el programa de la
                   Fig. 8.10 del libro; lento en mallas grandes).

    Plantilla de 5 puntos
    ---------------------
        (u[i-1,j] - 2u[i,j] + u[i+1,j]) / dx^2
      + (u[i,j-1] - 2u[i,j] + u[i,j+1]) / dy^2 = f[i,j]

    despejando u[i,j]:

        u[i,j] = ( (u[i-1,j] + u[i+1,j]) / dx^2
                 + (u[i,j-1] + u[i,j+1]) / dy^2 - f[i,j] ) / (2/dx^2 + 2/dy^2)

    Si dx = dy y f = 0 queda el promedio de los 4 vecinos (Ec. 8.27 del libro):
        u[i,j] = (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1]) / 4

    Devuelve
    --------
    u, iteraciones, convergio (bool), historial (lista con el cambio máximo
    de cada iteración), tiempo (s)
    """
    u = np.asarray(u, dtype=float)
    if f is None:
        f = np.zeros_like(u)

    # Coeficientes de la plantilla (constantes en toda la malla)
    cx = 1.0 / dx ** 2
    cy = 1.0 / dy ** 2
    inv_diag = 1.0 / (2.0 * cx + 2.0 * cy)

    historial = []
    convergio = False
    t0 = time.perf_counter()

    if ordenamiento == "lexicografico":
        nx, ny = u.shape[0] - 1, u.shape[1] - 1
        it = 0
        for it in range(1, max_iter + 1):
            cambio_max = 0.0
            # Recorremos los puntos interiores fila por fila. Cada punto usa
            # inmediatamente los valores ya actualizados de sus vecinos
            # (esa es la diferencia entre Gauss-Seidel y Jacobi).
            for i in range(1, nx):
                for j in range(1, ny):
                    u_gs = (cx * (u[i - 1, j] + u[i + 1, j])
                            + cy * (u[i, j - 1] + u[i, j + 1])
                            - f[i, j]) * inv_diag
                    # SOR: u_nuevo = u_viejo + omega*(u_GS - u_viejo)
                    u_nuevo = u[i, j] + omega * (u_gs - u[i, j])
                    cambio = abs(u_nuevo - u[i, j])
                    if cambio > cambio_max:
                        cambio_max = cambio
                    u[i, j] = u_nuevo
            historial.append(cambio_max)
            if cambio_max < tol:
                convergio = True
                break
        return u, it, convergio, historial, time.perf_counter() - t0

    if ordenamiento != "rojo-negro":
        raise ValueError("ordenamiento debe ser 'rojo-negro' o 'lexicografico'")

    # ---- Ordenamiento rojo-negro (tablero de ajedrez) ----------------------
    # Se pintan los nodos interiores como un tablero: "rojo" si i+j es par y
    # "negro" si i+j es impar. En la plantilla de 5 puntos los 4 vecinos de
    # un nodo rojo son negros y viceversa. Por eso:
    #   1) se actualizan TODOS los rojos a la vez (sólo dependen de negros),
    #   2) luego TODOS los negros, usando los rojos recién calculados.
    # Sigue siendo Gauss-Seidel (cada nodo usa los valores más nuevos de sus
    # vecinos), sólo cambia el orden de recorrido, y permite vectorizar con
    # numpy (cientos de veces más rápido que un doble bucle en Python).
    ii, jj = np.meshgrid(np.arange(1, u.shape[0] - 1),
                         np.arange(1, u.shape[1] - 1), indexing="ij")
    mascaras = [(ii + jj) % 2 == 0, (ii + jj) % 2 == 1]   # [rojos, negros]
    f_int = f[1:-1, 1:-1]
    interior = u[1:-1, 1:-1]          # vista: modificarla modifica u

    it = 0
    for it in range(1, max_iter + 1):
        cambio_max = 0.0
        for mascara in mascaras:
            u_gs = (cx * (u[:-2, 1:-1] + u[2:, 1:-1])
                    + cy * (u[1:-1, :-2] + u[1:-1, 2:])
                    - f_int) * inv_diag
            delta = omega * (u_gs[mascara] - interior[mascara])
            interior[mascara] += delta
            if delta.size:
                cambio_max = max(cambio_max, float(np.abs(delta).max()))
        historial.append(cambio_max)
        if cambio_max < tol:
            convergio = True
            break
    return u, it, convergio, historial, time.perf_counter() - t0


def gauss_seidel_laplace(u, dx, dy, tol=1e-6, max_iter=100_000, omega=1.0,
                         ordenamiento="rojo-negro"):
    """Ecuación de Laplace u_xx + u_yy = 0: caso f = 0 de gauss_seidel_poisson."""
    return gauss_seidel_poisson(u, dx, dy, f=None, tol=tol, max_iter=max_iter,
                                omega=omega, ordenamiento=ordenamiento)


# =============================================================================
#  UTILIDADES (tablas en consola y carpeta de gráficas)
# =============================================================================

def asegurar_carpeta_graficas():
    """Crea la carpeta graficas/ si no existe y devuelve su ruta."""
    os.makedirs(CARPETA_GRAFICAS, exist_ok=True)
    return CARPETA_GRAFICAS


def ruta_grafica(nombre):
    """Ruta completa para guardar una gráfica dentro de graficas/."""
    return os.path.join(asegurar_carpeta_graficas(), nombre)


def error_porcentual(y_num, y_ex):
    """Error relativo verdadero en %:  |y_exacta - y_num| / |y_exacta| * 100.
    Donde la solución exacta vale 0 no se puede dividir: se devuelve NaN."""
    y_num = np.asarray(y_num, dtype=float)
    y_ex = np.asarray(y_ex, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        e = np.abs(y_ex - y_num) / np.abs(y_ex) * 100.0
    e[np.abs(y_ex) < 1e-14] = np.nan
    return e


def indices_tabla(x, paso_tabla):
    """Índices de los nodos que caen en múltiplos de 'paso_tabla' (para no
    imprimir miles de filas cuando h es pequeño)."""
    h = x[1] - x[0]
    k = max(1, int(round(paso_tabla / h)))
    idx = list(range(0, len(x), k))
    if idx[-1] != len(x) - 1:
        idx.append(len(x) - 1)
    return idx


def _fmt_err(e):
    return "         n/d" if np.isnan(e) else f"{e:12.4e}"


def imprimir_tabla(x, soluciones, y_ex, titulo="", paso_tabla=None):
    """Imprime una tabla  x | y numérico | y exacto | error abs | error %.

    soluciones : dict {"Euler": y_euler, "RK4": y_rk4, ...} (mismos nodos x)
    y_ex       : solución exacta en los nodos x
    paso_tabla : si se da, sólo se imprimen los x múltiplos de ese paso
    """
    idx = indices_tabla(x, paso_tabla) if paso_tabla else range(len(x))
    nombres = list(soluciones.keys())

    cab = f"{'x':>8} |"
    for n in nombres:
        cab += f" {('y ' + n):>14} |"
    cab += f" {'y exacta':>14} |"
    for n in nombres:
        cab += f" {('|err| ' + n):>12} | {('err% ' + n):>12} |"
    linea = "-" * len(cab)

    if titulo:
        print("\n" + titulo)
    print(linea)
    print(cab)
    print(linea)
    for k in idx:
        fila = f"{x[k]:8.4f} |"
        for n in nombres:
            fila += f" {soluciones[n][k]:14.8f} |"
        fila += f" {y_ex[k]:14.8f} |"
        for n in nombres:
            e_abs = abs(y_ex[k] - soluciones[n][k])
            e_pct = error_porcentual([soluciones[n][k]], [y_ex[k]])[0]
            fila += f" {e_abs:12.4e} | {_fmt_err(e_pct)} |"
        print(fila)
    print(linea)


def imprimir_convergencia(filas, titulo="Resumen de convergencia"):
    """Tabla resumen del error vs h y ORDEN OBSERVADO del método.

    filas : lista de dicts {"metodo": str, "h": float, "err_max": float,
                            "err_final": float}
            err_max   = máximo error absoluto en todo el intervalo
            err_final = error absoluto en el último punto x = xf

    Orden observado entre dos pasos consecutivos h1 > h2:
        p = ln(E1 / E2) / ln(h1 / h2)
    Debería acercarse a 1 para Euler y a 4 para RK4.
    """
    print("\n" + titulo)
    cab = (f"{'Método':>8} | {'h':>8} | {'err. máx. abs':>14} | "
           f"{'err. en xf':>12} | {'E(h_ant)/E(h)':>13} | {'orden p':>8}")
    print("-" * len(cab))
    print(cab)
    print("-" * len(cab))
    anterior = {}
    for fila in filas:
        m, h, e = fila["metodo"], fila["h"], fila["err_max"]
        razon, orden = "", ""
        if m in anterior:
            h_a, e_a = anterior[m]
            if e > 0 and e_a > 0:
                razon = f"{e_a / e:13.2f}"
                orden = f"{np.log(e_a / e) / np.log(h_a / h):8.2f}"
        print(f"{m:>8} | {h:8.4g} | {e:14.4e} | {fila['err_final']:12.4e} | "
              f"{razon:>13} | {orden:>8}")
        anterior[m] = (h, e)
    print("-" * len(cab))
