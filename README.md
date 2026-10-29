# Detección de Fallas de Sensor y Reconfiguración en un Convertidor DC-DC Boost mediante un Observador de Luenberger con Umbral Dinámico

**Autores:** Edrich Paul Pacheco Camargo (`epachecoc@unsa.edu.pe`) y Cristhian Alvaro Hancco Chura (`chanccoc@unsa.edu.pe`)  
**Institución:** Universidad Nacional de San Agustín de Arequipa (UNSA) - Escuela Profesional de Ingeniería Electrónica

---

## 1. Descripción General
Este repositorio contiene los modelos de Simulink (`.slx`), el código fuente del bloque de umbral dinámico (`Umbral_Dinamico.m`), los scripts de simulación y análisis estadístico (`.m`), la lista de semillas pseudoaleatorias y las tablas de resultados (`.csv`) que respaldan el artículo:
> *"Detección de fallas de sensor y reconfiguración en un convertidor DC-DC Boost mediante un observador de Luenberger con umbral dinámico"* (Réplica y extensión de Li, Zhang & Li, *Sensors* 2018, 18(5):1375).

---

## 2. Entorno de Software y Parámetros de Configuración
* **Versión de MATLAB/Simulink utilizada:** MATLAB R2025b (compatible con R2020a o superior).
* **Configuración del Solver (Simulink y Scripts):** Paso fijo (`Fixed-step`), integrador Runge-Kutta de 4.º orden (`ode4`), tamaño de paso $T_s = 2 \times 10^{-6}\text{ s}$ ($2\text{ }\mu\text{s}$), tiempo de simulación $t \in [0, 0.6]\text{ s}$ (y hasta $1.2\text{ s}$ para pruebas de ruido prolongado).
* **Condiciones iniciales en estado estacionario:** $x(0) = \hat{x}(0) = [7.6\text{ A},\; 380\text{ V}]^T$.
* **Parámetros físicos del convertidor Boost:** $V_{in} = 190\text{ V}$, $L = 5\text{ mH}$, $C_0 = 2.85\text{ mF}$, $R = 0.082\text{ }\Omega$, $R_{load,\text{nom}} = 100\text{ }\Omega$.
* **Ganancia del observador de Luenberger ($H$):**
  $$H = \begin{bmatrix} 36.4677 & -17.0967 \\ 28.5219 & 75.6359 \end{bmatrix}$$
* **Parámetros del Umbral Dinámico Propuesto (`Umbral_Dinamico.m`):**
  * Longitud de ventana móvil: $20\text{ ms} \implies N = 10\,000\text{ muestras}$ (con $N_{\max} = 25\,000$ para soportar hasta $50\text{ ms}$).
  * Factor de sensibilidad estadística: $k = 4$.
  * Cota inferior de seguridad: $J_{th,\min} = 0.02\text{ A}^2+\text{V}^2$ (idéntica al umbral fijo del artículo base).

---

## 3. Estructura de Archivos del Repositorio

### A. Modelos de Simulink (`.slx`)
1. **Modelos de la Réplica Original (Umbral Fijo $J_{th} = 0.02$):**
   * `Boost_Falla_Circuito_Abierto.slx`: Réplica de la Figura 4 del artículo base (falla de circuito abierto en sensor de corriente y de tensión en $t = 0.2\text{ s}$).
   * `Boost_Falla_Desviacion_Ganancia.slx`: Réplica de la Figura 5 del artículo base (desviación de ganancia de $+10\%$ en sensor de corriente y de tensión en $t = 0.2\text{ s}$).
   * `Boost_Falla_Anomalia_Ruido.slx`: Réplica de la Figura 6 del artículo base (anomalía de ruido en sensor de corriente y de tensión en $t = 0.2\text{ s}$).
2. **Modelos del Aporte Propuesto (Con bloque `Umbral_Dinamico`):**
   * `Boost_Umbral_Dinamico.slx`: Modelo con perturbación transitoria de carga del $3.5\%$ ($R_{load}$ de $100\text{ }\Omega$ a $96.6\text{ }\Omega$) iniciada en $t = 0.3\text{ s}$ sin falla de sensor.
   * `Boost_ExpB_Seis_Fallas.slx`: Modelo unificado de evaluación de las 6 fallas en $t = 0.2\text{ s}$ con el bloque `Umbral_Dinamico`.

### B. Scripts de MATLAB (`.m`)
* `Parametros_Convertidor_Boost_Umbral_Dinamico_Aporte.m`: Script de inicialización que carga todos los parámetros físicos, matrices de estado, ganancias y parámetros del umbral dinámico en el *Workspace*.
* `Umbral_Dinamico.m`: Código fuente independiente del Algoritmo 1 (Anexo A) implementado mediante buffer circular y sumas acumuladas en tiempo constante $O(1)$.
* `Generar_ExpB_fig_dinamico.m`: Ejecuta las 6 fallas con el umbral dinámico y genera la figura en panel $2\times 3$ (`fig_dinamico_blanco.png`).
* `Ejecutar_ExpC_Tiempos_Deteccion.m`: Mide con precisión de un paso del solver ($2\text{ }\mu\text{s}$) los tiempos de detección y reconfiguración para umbral fijo y dinámico.
* `Ejecutar_ExpD_Sensibilidad_3x3.m`: Evalúa la perturbación transitoria de carga y construye la matriz de sensibilidad $3\times 3$ (ventanas de $20, 35, 50\text{ ms}$ vs. $k \in \{3, 4, 5\}$), generando `fig_carga_blanco.png`.
* `Ejecutar_Punto1_Ruido_30Corridas.m`: Ejecuta la validación estocástica de Monte Carlo ($30\text{ corridas}$ independientes con `rng(seed, 'twister')`, $\text{seed} \in \{1,\dots,30\}$, $n(t) = 1.35\,\texttt{randn}$) con y sin amplitud mínima de $0.2$.
* `Ejecutar_Punto2_Analisis_Enmascaramiento.m`: Ejecuta el análisis de riesgo de enmascaramiento para desviaciones de ganancia de $\delta \in \{10\%, 5\%, 2\%\}$ en $t = 0.20\text{ s}$ (nominal) y $t = 0.36\text{ s}$ (pico de perturbación de carga), así como para ruido prolongado ($0.2\text{ s}, 0.4\text{ s}, 1.0\text{ s}$).
* `Ejecutar_Punto3_Comparacion_Adaptativos.m`: Compara cuantitativamente el umbral fijo ($0.02$), el umbral dinámico propuesto (ventana móvil $\mu + k\sigma$), el filtro adaptativo EWMA y el estimador robusto Mediana + MAD.
* `Ejecutar_Punto4_Reproducibilidad_y_README.m`: Verifica la equivalencia numérica exacta (Caso de Control) entre los scripts `.m` y el modelo Simulink `.slx`, y genera los archivos auxiliares de documentación.

### C. Tablas de Datos y Semillas Generadas (`.csv`)
* `Lista_Semillas_1_a_30.csv`
* `Tabla_Punto1_Resumen_Estadistico_Ruido.csv` y `Tabla_Punto1_Detalle_30Semillas.csv`
* `Tabla_Punto2AB_Ganancia_10_5_2_Enmascaramiento.csv` y `Tabla_Punto2C_Ruido_Largo_Margenes.csv`
* `Tabla_Punto3_Comparacion_Fijo_Dinamico_EWMA_MAD.csv`

---

## 4. Orden de Ejecución Recomendado
1. Abrir MATLAB y ubicar el *Current Folder* en la carpeta raíz de este repositorio.
2. Ejecutar `Parametros_Convertidor_Boost_Umbral_Dinamico_Aporte.m` para cargar las variables base en el *Workspace*.
3. Para simular cualquiera de los modelos `.slx` en Simulink, abrir el archivo correspondiente y presionar **Run**.
4. Para reproducir las tablas estadísticas y figuras de alta resolución (300 DPI) del manuscrito revisado, ejecutar en orden:
   * `Ejecutar_Punto1_Ruido_30Corridas.m`
   * `Ejecutar_Punto2_Analisis_Enmascaramiento.m`
   * `Ejecutar_Punto3_Comparacion_Adaptativos.m`
   * `Ejecutar_Punto4_Reproducibilidad_y_README.m`

---

## 5. Verificación de Equivalencia Numérica (Caso de Control: Simulink vs. Script)
Dado que tanto los modelos `.slx` de Simulink como los scripts `.m` de evaluación masiva implementan las mismas ecuaciones diferenciales en espacio de estados e integran con el método de Runge-Kutta de 4.º orden (`ode4`) a paso fijo $T_s = 2\text{ }\mu\text{s}$ desde $x(0) = \hat{x}(0) = [7.6\text{ A}, 380\text{ V}]^T$, la evolución temporal de los estados y residuos es bit a bit idéntica:
* En $t = 0.200000\text{ s}$ (muestra $100\,001$), tanto en `Boost_ExpB_Seis_Fallas.slx` como en los scripts `.m`, el estado antes de la falla es $i_L(0.2) = 7.585605\text{ A}$ y $v_C(0.2) = 379.361467\text{ V}$.
* Al inyectar la Falla 1 (circuito abierto en corriente), ambos entregan exactamente $r_1(0.2) = 57.5414\text{ A}^2+\text{V}^2$.
* Al inyectar la Falla 3 (ganancia $+10\%$ en corriente), ambos entregan exactamente $r_1(0.2) = 0.5754\text{ A}^2+\text{V}^2$.
* Al inyectar la Falla 5 (ruido en corriente), ambos entregan exactamente $r_1(0.2) = 0.2323\text{ A}^2+\text{V}^2$ con el generador interno del bloque, y $r_1(0.2) = 0.6012\text{ A}^2+\text{V}^2$ para la Semilla 1 (`rng(1, 'twister')`).
