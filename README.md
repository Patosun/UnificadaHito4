# Práctica de Métodos Numéricos — Ejercicios computacionales (Cap. 7 y 8)

Libro: P. Dechaumphai & N. Wansophark, *Numerical Methods in Science and Engineering*, Alpha Science, 2022.

- **Capítulo 7 (EDO):** ejercicios 2, 4, 5, 16, 17, 19, 20, 22, más las gráficas de los ejercicios 1, 3, 18 y 21 (los que se resolvieron a mano).
- **Capítulo 8 (EDP elíptica):** ejercicios 3 y 4 (Laplace con diferencias finitas, Gauss-Seidel y SOR).

> **Restricción de la práctica:** en el Cap. 7 solo se usan **Euler** y **Runge-Kutta de 4º orden (RK4)**. Donde el libro pide Heun, Euler modificado, RK2 o RK3 (ejercicios 16, 19 y 20), esos métodos se reemplazaron por Euler y/o RK4.

---

## 1. Cómo ejecutar

### Requisitos
- Python 3 (probado con **Python 3.12.10**)
- `numpy` y `matplotlib` (no hace falta nada más)

### Comandos (PowerShell o CMD, dentro de la carpeta del proyecto)

| Comando | Qué hace |
|---|---|
| `python --version` | Comprueba que Python está instalado y en el PATH. |
| `python -m pip install -r requirements.txt` | Instala las dos librerías necesarias, numpy y matplotlib (solo la primera vez). |
| `python main.py` | Ejecuta **todos** los ejercicios en orden: imprime las tablas en consola y guarda las gráficas en `graficas/`. Al final muestra un resumen (OK/FALLÓ y el tiempo de cada script). Tarda unos 30 s; casi todo ese tiempo es del Cap. 8. |
| `python main.py cap8_ej03 cap8_ej04` | Ejecuta solo los scripts indicados (sin la extensión `.py`). |
| `python cap7_ej16.py` | Ejecuta un único ejercicio. Cada script es independiente. |
| `python main.py > resultados.txt` | Guarda toda la salida de consola en un archivo de texto (útil para el informe). |

> Si `python` no se reconoce, abre una terminal **nueva** (Python se instaló en `%LOCALAPPDATA%\Programs\Python\Python312` y se agregó al PATH del usuario) o usa el lanzador `py main.py`.

### Estructura del proyecto

```
metodos.py              Módulo común: euler, rk4, gauss_seidel_poisson / gauss_seidel_laplace,
                        omega_optimo y utilidades para imprimir tablas.
main.py                 Ejecuta todos los ejercicios en orden.
cap7_ej01_grafica.py    Gráfica del ej. 7.1 (Euler, h = 0.25)
cap7_ej02.py            Euler, h = 0.25 / 0.1 / 0.01
cap7_ej03_grafica.py    Gráfica del ej. 7.3 (Euler, h = 0.2)
cap7_ej04.py            Euler, h = 0.2 / 0.1 / 0.01
cap7_ej05.py            Euler, h = 0.5 / 0.1 / 0.01
cap7_ej16.py            Euler y RK4 (EDO del ej. 12), h = 0.25 / 0.1 / 0.01
cap7_ej17.py            RK4, h = 0.1 / 0.01 / 0.001 / 0.0001
cap7_ej18_grafica.py    Gráfica del ej. 7.18 (RK4, h = 0.25)
cap7_ej19.py            Euler y RK4 (EDO del ej. 18), h = 0.05, con tiempos de cómputo
cap7_ej20.py            RK4, h = 0.1
cap7_ej21_grafica.py    Gráfica del ej. 7.21 (Euler y RK4, h = 0.5)
cap7_ej22.py            Euler y RK4 con elección justificada de h
cap8_ej03.py            Laplace en el cuadrado unitario, malla 100x100, exacta u = xy
cap8_ej04.py            Placa 2x1 (Ejemplo 8.1, Fig. 8.10), malla 200x100, exacta Ec. (8.25)
graficas/               PNG generados por los scripts
```

Todos los scripts siguen la misma estructura: constantes del problema (`X0`, `XF`, `Y0`, pasos), `f(x, y)` (lado derecho de la EDO), `exacta(x)`, funciones para resolver, imprimir y graficar, y un `main()`. Las tablas muestran **x, y numérico, y exacto, error absoluto y error relativo %**. El error % no está definido donde la solución exacta vale 0 y en ese caso aparece `n/d`.

---

## 2. Métodos del Capítulo 7 (resumen para la defensa)

**Euler:** `y_{i+1} = y_i + h·f(x_i, y_i)`. Avanza por la recta tangente. Es la serie de Taylor cortada en el término lineal. El error local es O(h²) y el **global es O(h)**: si h se divide por 10, el error baja ≈10 veces.

**RK4:** se calculan 4 pendientes en el paso
`k1 = f(x_i, y_i)`, `k2 = f(x_i + h/2, y_i + h/2·k1)`, `k3 = f(x_i + h/2, y_i + h/2·k2)`, `k4 = f(x_i + h, y_i + h·k3)`
y se avanza con el promedio ponderado `y_{i+1} = y_i + h/6·(k1 + 2k2 + 2k3 + k4)`. Coincide con Taylor hasta h⁴. El **error global es O(h⁴)**: si h se divide por 10, el error baja ≈10 000 veces. Cuesta 4 evaluaciones de f por paso (Euler cuesta 1).

**Orden observado.** Los scripts calculan `p = ln(E₁/E₂) / ln(h₁/h₂)` entre dos pasos consecutivos. Debe acercarse a 1 en Euler y a 4 en RK4.

---

## 3. Resultados y conclusiones — Capítulo 7

*(Los errores son absolutos, `|y_exacta − y_numérica|`, salvo que se indique %.)*

### Gráficas de los ejercicios hechos a mano
| Ej. | Método | Resultado en el extremo final |
|---|---|---|
| 7.1 | Euler, h = 0.25 | y(2) = 5.2690 vs 5.3863 exacta → error 2.18 % |
| 7.3 | Euler, h = 0.2 | y(3) = 4.5143 vs 5.8741 → error 23.1 % (la solución crece muy rápido y Euler se queda atrás) |
| 7.18 | RK4, h = 0.25 | y(2) = −0.49988 vs −0.5 → error 0.025 % |
| 7.21 | Euler y RK4, h = 0.5 | y(2): Euler 3.5388 (4.66 %), RK4 3.71205 (0.0045 %), exacta 3.71189 |

Archivos: `cap7_ej01_euler.png`, `cap7_ej03_euler.png`, `cap7_ej18_rk4.png`, `cap7_ej21_euler_rk4.png`.

### Ej. 2 — `dy/dx = 1 + y/x`, y(1) = 2, Euler
| h | error máx | orden p |
|---|---|---|
| 0.25 | 1.17e-1 | — |
| 0.1 | 4.88e-2 | 0.96 |
| 0.01 | 4.99e-3 | 0.99 |

**Conclusión:** al dividir h por 10 el error baja ≈10 veces (orden ≈1, como predice la teoría). El error crece a lo largo de x porque se acumula paso a paso.

### Ej. 4 — `dy/dx = 1 + y/x + (y/x)²`, y(1) = 0, Euler
| h | error máx (en x = 3) | orden p |
|---|---|---|
| 0.2 | 1.36 (23 %) | — |
| 0.1 | 0.814 (14 %) | 0.74 |
| 0.01 | 0.100 (1.7 %) | 0.91 |

**Conclusión:** la EDO es no lineal y la solución `x·tan(ln x)` crece cada vez más rápido. El error acumulado es grande y el orden observado se acerca a 1 recién con h pequeño (con h grande todavía no se está en el régimen asintótico). Para este problema Euler necesita h muy pequeño.

### Ej. 5 — `dy/dx = −y + x + 1`, y(0) = 1, Euler
| h | error máx | error en x = 5 | orden p |
|---|---|---|---|
| 0.5 | 1.18e-1 | 5.8e-3 | — |
| 0.1 | 1.92e-2 | 1.6e-3 | 1.13 |
| 0.01 | 1.85e-3 | 1.7e-4 | 1.02 |

**Conclusión:** orden 1 confirmado. Aquí el error **disminuye** hacia el final porque la ecuación es estable: el término `−y` amortigua las perturbaciones, igual que `e^{−x}` decae.

### Ej. 16 — `dy/dx = (1 − y)(4 − y)`, y(0) = 0 (EDO del ej. 12), Euler y RK4
| h | Euler error máx | RK4 error máx |
|---|---|---|
| 0.25 | 4.02e-1 | 4.08e-3 |
| 0.1 | 9.31e-2 | 8.53e-5 (p = 4.2) |
| 0.01 | 7.50e-3 (p = 1.1) | 6.71e-9 (p = 4.1) |

**Conclusiones:**
- RK4 con h = 0.25 ya es 100 veces más preciso que Euler con h = 0.25. Con h = 0.01 la diferencia es de 6 órdenes de magnitud.
- **Detalle llamativo de Euler con h = 0.25:** `y₁ = 0 + 0.25·(1)(4) = 1`, y y = 1 es un **punto de equilibrio** (f = 0), así que la solución se queda "pegada" en 1. Un paso demasiado grande hace saltar a Euler hasta el valor límite.

### Ej. 17 — `dy/dx = y + x²`, y(1) = 1, RK4
| h | error máx | razón | orden p |
|---|---|---|---|
| 0.1 | 9.10e-6 | — | — |
| 0.01 | 9.92e-10 | 9 179 | 3.96 |
| 0.001 | 9.59e-14 | 10 338 | 4.01 |
| 0.0001 | 4.53e-14 | 2.1 | 0.33 |

**Conclusión:** de 0.1 a 0.001 se cumple O(h⁴) (el error baja ×10⁴ por cada ÷10 en h). De 0.001 a 0.0001 el error **deja de bajar**: se llegó al límite del **error de redondeo** de la doble precisión (~16 cifras). Con más pasos se acumulan más errores de redondeo, así que achicar h ya no sirve.

### Ej. 19 — `dy/dx = 1/x² − y² − y/x`, y(1) = −1 (EDO del ej. 18), h = 0.05
| Método | Evaluaciones de f | Tiempo por ejecución* | Error máx | Error máx % |
|---|---|---|---|---|
| Euler | 20 | ≈ 0.012 ms | 2.88e-2 | 5.76 % |
| RK4 | 80 | ≈ 0.039 ms | 2.48e-7 | 5.0e-5 % |

\*Promedio de 2000 ejecuciones. El valor depende del equipo.

**Conclusión:** RK4 tarda ≈3–4 veces más (hace 4 evaluaciones de f por paso en vez de 1), pero su error es ≈10⁵ veces menor. Para igualar esa precisión Euler necesitaría un h unas 10⁵ veces más chico, o sea muchísimo más tiempo. RK4 es claramente más eficiente.

### Ej. 20 — `dy/dx = x²eˣ + 2y/x`, y(1) = 0, RK4 con h = 0.1
Error máximo 1.7e-4 en x = 2 (y = 18.6829 vs 18.6831), y el error relativo es ≤ 0.003 % en todo el intervalo. Con solo 10 pasos RK4 reproduce la exacta con 4–5 cifras significativas.

### Ej. 22 — `dy/dx = x⁴eˣ + 4y/x`, y(1) = 0, Euler y RK4, elección de h
**Criterio:** error relativo máximo < 1 % en todo [1, 2]. Se prueban h = 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001 y para cada método se toma el **mayor h que cumple** (el más barato):

| h | Euler err máx % | RK4 err máx % |
|---|---|---|
| 0.1 | 35.1 | 0.030 ✔ |
| 0.01 | 4.45 | 4.9e-6 ✔ |
| 0.0025 | 1.14 | 2.0e-8 ✔ |
| **0.001** | **0.457 ✔** | 5.2e-10 ✔ |

- RK4 cumple ya con **h = 0.1** (10 pasos).
- Euler recién cumple con **h = 0.001** (1000 pasos).
- La tabla comparativa usa **h = 0.001** (el que exige el método más débil), para comparar ambos en los mismos x. En x = 2: Euler 74.4306 (0.40 %), RK4 74.732388 (1.6e-10 %), exacta 74.732388.

**Justificación:** h se eligió con un criterio de precisión explícito y un estudio de convergencia, no "a ojo". Además se ve la ventaja de RK4: necesita 100 veces menos pasos para la misma precisión.

---

## 4. Capítulo 8 — Explicación conceptual (ecuación elíptica)

### 4.1 Problema
Ecuación de Laplace (o Poisson si hay fuente f):

    u_xx + u_yy = f(x, y)     en un rectángulo, con u conocida en el borde (Dirichlet)

Se divide el dominio en una malla de `nx × ny` intervalos, con `dx = Lx/nx` y `dy = Ly/ny`. Los nodos son `x_i = i·dx` y `y_j = j·dy`, y `u_{i,j} ≈ u(x_i, y_j)`.

### 4.2 Derivación de la plantilla de 5 puntos
Se suman las series de Taylor hacia adelante y hacia atrás en x:

    u(x+dx) = u + dx·u_x + dx²/2·u_xx + dx³/6·u_xxx + dx⁴/24·u_xxxx + ...
    u(x−dx) = u − dx·u_x + dx²/2·u_xx − dx³/6·u_xxx + dx⁴/24·u_xxxx − ...
    -------------------------------------------------------------------- (sumando)
    u(x+dx) + u(x−dx) = 2u + dx²·u_xx + dx⁴/12·u_xxxx + ...

Despejando sale la **diferencia central** de segundo orden:

    u_xx ≈ (u_{i+1,j} − 2u_{i,j} + u_{i−1,j}) / dx²        error O(dx²)

y lo mismo en y. Reemplazando en la ecuación queda la plantilla de 5 puntos:

    (u_{i+1,j} − 2u_{i,j} + u_{i−1,j}) / dx² + (u_{i,j+1} − 2u_{i,j} + u_{i,j−1}) / dy² = f_{i,j}

Si `dx = dy = h` y f = 0, se despeja el nodo central y queda **el promedio de sus 4 vecinos** (Ec. 8.27 del libro):

    u_{i,j} = (u_{i−1,j} + u_{i+1,j} + u_{i,j−1} + u_{i,j+1}) / 4

```
              (i, j+1)
                 |
   (i−1, j) — (i, j) — (i+1, j)        coeficientes:  1 / −4 / 1  (en cruz)
                 |
              (i, j−1)
```

El error de truncamiento depende de las **derivadas cuartas** (`h²/12·(u_xxxx + u_yyyy)`), así que el método es **O(h²)**: si h se divide por 2, el error baja ≈4 veces.

### 4.3 Condiciones de borde
Son de **Dirichlet** (valor conocido). Los nodos del borde **no son incógnitas**: se cargan una sola vez en la matriz `u` (filas/columnas extremas) y no se modifican nunca. Gauss-Seidel solo recorre los nodos interiores `i = 1..nx−1`, `j = 1..ny−1`. Cuando la plantilla de un nodo vecino al borde usa un valor de borde, ese valor conocido funciona como el término independiente del sistema (como el 140 de la Ec. 8.20).
- Ej. 8.3: `u(x,0) = 0`, `u(0,y) = 0`, `u(x,1) = x`, `u(1,y) = y`.
- Ej. 8.4: `T = 0` en los bordes inferior, izquierdo y derecho, y `T(x,1) = sin(πx/2)` arriba.

### 4.4 Gauss-Seidel y criterio de parada
La plantilla aplicada a todos los nodos interiores da un sistema lineal `A·u = b` muy grande (19 701 ecuaciones en el Ej. 8.4) y ralo (cada fila tiene como mucho 5 coeficientes). La matriz es **diagonalmente dominante** (4 en la diagonal y como mucho cuatro −1 por fila), por eso Gauss-Seidel converge. Además ocupa muy poca memoria: no hace falta guardar A, basta la matriz u.

**Algoritmo:**
1. Cargar las condiciones de borde y una semilla inicial en el interior (ceros, o la interpolación lineal del libro).
2. En cada iteración se recorren los nodos interiores y cada uno se reemplaza por el valor de la plantilla, **usando inmediatamente los valores ya actualizados** de los vecinos. Eso es lo que distingue a Gauss-Seidel de Jacobi y lo que lo hace converger más rápido.
3. **Criterio de parada:** se guarda el cambio máximo de la iteración, `máx |u_nuevo − u_viejo|`. Si es **< tol**, se detiene. Si se llega a `max_iter` sin cumplirlo, se informa que no convergió. Los dos son parámetros (`tol`, `max_iter`) de `gauss_seidel_poisson`.

**Advertencia (se ve en el Ej. 8.3):** que el cambio entre iteraciones sea < tol **no** garantiza que el error sea < tol. Con convergencia lenta (radio espectral ρ ≈ 1), el error real es ≈ `tol/(1−ρ)`. Con GS puro en 100×100 y tol = 1e-4, el error real fue ≈ 0.1. Por eso se usa tol = 1e-8.

### 4.5 SOR (sobre-relajación sucesiva)
Gauss-Seidel en mallas finas es lento: necesita O(N²) iteraciones (10 736 en 100×100). SOR "exagera" la corrección de Gauss-Seidel:

    u_nuevo = u_viejo + ω·(u_GS − u_viejo),     1 < ω < 2

- ω = 1 es Gauss-Seidel puro. 1 < ω < 2 es sobre-relajación. Con ω ≥ 2 diverge.
- Para un rectángulo con `dx = dy` el valor óptimo es `ω = 2/(1 + √(1 − ρ_J²))`, con `ρ_J = [cos(π/nx) + cos(π/ny)]/2` (función `omega_optimo`).
- Con ω óptimo hacen falta O(N) iteraciones: **288 en vez de 10 736**.
- SOR converge a **la misma solución** de las ecuaciones de diferencias, solo que llega antes.

### 4.6 Ordenamiento rojo-negro (implementación)
El libro recorre los nodos fila por fila (ordenamiento *lexicográfico*). Ese recorrido está implementado tal cual (`ordenamiento="lexicografico"`), pero en Python puro los dobles bucles son lentos. Por eso también se implementó el ordenamiento **rojo-negro**: los nodos se pintan como un tablero de ajedrez (i+j par = rojo, impar = negro). En la plantilla de 5 puntos los vecinos de un rojo son todos negros, así que se pueden actualizar **todos los rojos a la vez** con numpy y después todos los negros. **Sigue siendo Gauss-Seidel**: cada nodo usa los valores más nuevos de sus vecinos y solo cambia el orden del recorrido. Las tasas de convergencia y la teoría de ω óptimo son las mismas, y en este equipo cada iteración resultó unas 30–35 veces más rápida.

---

## 5. Resultados y conclusiones — Capítulo 8

*(Tiempos medidos en este equipo; cambian según la máquina.)*

### Ej. 3 — Laplace en [0,1]², malla 100×100 (9 801 incógnitas), exacta u = xy, tol = 1e-8

**Validación:** con la malla 3×3 del Ej. 2 (la resuelta a mano) el programa da exactamente u = 1/9, 2/9, 2/9, 4/9 en los nodos interiores, igual a la exacta.

| Método | ω | Iteraciones | Tiempo | Error máx |
|---|---|---|---|---|
| Gauss-Seidel (rojo-negro) | 1.0000 | 10 736 | 2.8 s | 1.0e-5 |
| SOR (rojo-negro) | 1.9391 | **288** | **0.11 s** | 1.1e-7 |
| SOR (lexicográfico, como el libro) | 1.9391 | 362 | 4.5 s | 2.6e-8 |

Efecto de la tolerancia (GS puro): tol 1e-4 → 1 471 it., error 9.6e-2; tol 1e-6 → 6 071 it., error 1.0e-3; tol 1e-8 → 10 736 it., error 1.0e-5.

**Conclusiones:**
- `u = xy` es bilineal y sus derivadas cuartas son nulas, así que **la plantilla de 5 puntos es exacta para esta solución** (error de truncamiento cero). Todo el error que queda es error de **iteración**, es decir, de cortar Gauss-Seidel antes de tiempo. Por eso el error baja con la tolerancia y no con la malla.
- SOR reduce las iteraciones unas 37 veces respecto de Gauss-Seidel.
- Gráfica `cap8_ej03_laplace.png`: mapa de color con isolíneas (hipérbolas xy = cte), mapa de error e historial de convergencia.

### Ej. 4 — Placa 2×1 (Ejemplo 8.1, Fig. 8.10), malla 200×100 (19 701 incógnitas)
Bordes `T = 0` (inferior, izquierdo, derecho), `T = sin(πx/2)` (superior). Exacta, Ec. (8.25): `T = sin(πx/2)·sinh(πy/2)/sinh(π/2)`.

**1) Reproducción del libro (malla 8×4, GS, tol = 1e-5):** 23 iteraciones (el libro informa 24; ver nota) y los valores coinciden con la **Tabla 8.1** en las 3 cifras que publica el libro. Error máximo 2.46e-3.

> *Nota:* el número exacto de iteraciones depende de la semilla y del orden del recorrido. Tal como está impresa, la semilla del libro `t(i,j) = j*t(i,5)/5` lee una columna que todavía vale 0. Con distintas combinaciones de semilla y orden se obtienen entre 23 y 26 iteraciones. La solución final no cambia.

**2) Malla 200×100 (dx = dy = 0.01), tol = 1e-8:**

| Método | ω | Iteraciones | Tiempo | Error máx |
|---|---|---|---|---|
| Gauss-Seidel (rojo-negro) | 1.0000 | 14 540 | 8.7 s | 2.0e-5 |
| SOR (rojo-negro) | 1.9515 | **326** | **0.22 s** | **4.4e-6** |

**3) Verificación del orden O(h²):** el error máximo pasa de 2.46e-3 (h = 0.25) a 4.40e-6 (h = 0.01). La razón es **559**, cerca del valor teórico (0.25/0.01)² = **625**.

**Conclusiones:**
- Con 200×100 la solución numérica coincide con la exacta en 5 decimales en los nodos de la Fig. 8.9.
- El error de discretización que queda (≈4e-6, máximo en el centro de la placa) es O(h²).
- El error de GS puro (2e-5) es mayor que el de SOR con la misma tolerancia. GS converge tan lento que su criterio de parada se cumple cuando todavía falta camino hasta la solución discreta (ver 4.4).
- Gráficas: `cap8_ej04_placa.png` (mapa de temperatura y mapa de error) y `cap8_ej04_superficie_convergencia.png` (superficie 3D, equivalente al `mesh` del libro, y convergencia GS vs SOR).
