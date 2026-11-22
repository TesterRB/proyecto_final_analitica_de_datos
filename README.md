# Lluvia en el Valle del Cauca (1981-2003): datos funcionales, regresión funcional y red KAN

## Integrantes

<table>
  <tr>
    <td align="center"><a href="https://github.com/TesterRB"><img src="https://github.com/TesterRB.png?size=200" width="100" alt="TesterRB"/><br/><b>TesterRB</b></a><br/>Juan David Reyes Burbano</td>
    <td align="center"><a href="https://github.com/Majo0209"><img src="https://github.com/Majo0209.png?size=200" width="100" alt="Majo0209"/><br/><b>Majo0209</b></a><br/>María José Suárez Toro</td>
    <td align="center"><a href="https://github.com/Laura25a"><img src="https://github.com/Laura25a.png?size=200" width="100" alt="Laura25a"/><br/><b>Laura25a</b></a><br/>Laura Daniela Astudillo</td>
    <td align="center"><a href="https://github.com/MarcAragon"><img src="https://github.com/MarcAragon.png?size=200" width="100" alt="MarcAragon"/><br/><b>MarcAragon</b></a><br/>Marco Antonio Aragón</td>
    <td align="center"><a href="https://github.com/mariapaulapinillo"><img src="https://github.com/mariapaulapinillo.png?size=200" width="100" alt="mariapaulapinillo"/><br/><b>mariapaulapinillo</b></a><br/>Maria Paula Pinillo</td>
  </tr>
</table>

Proyecto final de Analítica de Datos 2026-2S.

Corrige con estaciones de la CVC la lluvia diaria del producto satelital de la clase. La lluvia de cada estación y cada año se trata como una curva:
- **Ítems 1 y 2:** construcción de las curvas y análisis exploratorio funcional.
- **Ítem 3:** comparación funcional entre la estación y el satélite.
- **Ítem 4:** regresión funcional y una red funcional con cabezas lineal, MLP y KAN.

## Contenido

| Archivo o carpeta | Qué es |
|---|---|
| `Proyecto_final_AnalíticaDeDatos.R` | **El script de R, en un solo archivo.** Contiene los ítems 1, 2, 3 y 4.1, la regresión funcional (concurrente, función a función y rango reducido) y la **evaluación de todos los modelos**, incluidas las redes. En los comentarios del código de Python se le llama «la V2», el nombre de la versión integrada del grupo: no hay otro script de R. |
| `Proyecto_final_AnalíticaDeDatos.py` | **El script de Python, con `main()`.** Contiene la red funcional (codificador funcional, cabeza lineal, MLP o KAN, y decodificador FPCA), la equivalencia de la cabeza lineal con la regresión y la interpretación de la KAN. |
| `Proyecto_final.Rproj` | Proyecto de RStudio: al abrirlo, la carpeta de trabajo es esta. |
| `requirements.txt` | Versiones de Python usadas. |
| `Proyecto final AnalíticaDeDatos.pdf` | El informe. |
| `MANIFIESTO DE USO DE IA - … .docx` | Anexo en Word con el uso de la herramienta de IA (regla 10 y criterio 5 de la rúbrica). |
| `Datos_clase/` | Los datos de la clase, con sus nombres originales: lluvia diaria de las estaciones de la CVC (`P_1930_1979.xlsx`, `P_1980_1990.csv` y `P_1991_2003.csv`), coordenadas (`CVC_Wheather_Stations_Long_Lat.csv` y `.dta`) y catálogo de elevación de las estaciones (`CatalogoEstacionesElevacion.xlsx`), y las péntadas satelitales (`ImagenesSatelitePentadal/`). |
| `Datos_externos/` | Datos que no son de la clase: el límite del Valle (GADM 4.1, `valle_del_cauca.geojson`), el DEM Copernicus GLO-30 (`dem_valle_copernicus30m.tif`) y el índice ONI de la NOAA CPC (`oni.ascii.txt`, la copia usada; el script avisa si encuentra otra versión). |
| `Datos_prueba/` | Fragmento de los datos de la clase del periodo de prueba: `P_2001_2003.csv`, con los registros de 2001-2003 de las estaciones del modelo (las mismas columnas de `P_1991_2003.csv`), y las 216 péntadas satelitales de 2001-2003 en `ImagenesSatelitePentadal/`. Los scripts no lo usan: leen los datos completos de `Datos_clase/`. |
| `Predicciones_R/` | Predicciones de la regresión funcional función a función en validación (1996-2000). |
| `Predicciones_python/` | Predicciones de las tres cabezas de la red, con la semilla elegida por validación. `semillas/` tiene las 9 corridas. |

Las predicciones se entregan ya calculadas (criterio 4.4 de la rúbrica). Al correr los scripts se vuelven a escribir.

**Formato de las predicciones:** `station_id, year, day, value_mm`, con `day` de 1 a 365 y `value_mm` ≥ 0 en mm/día.
- Las curvas viven en el dominio anual [0, 365]. En los años bisiestos, los 366 días se reescalan a ese dominio con t = (d − 0,5)·365/366.
- Cada archivo da la curva en los 365 centros de día.

## Cómo reproducir

### 1. Datos

Los datos de la clase ya vienen en `Datos_clase/` y los externos en `Datos_externos/`. La carpeta queda así:

```
Entrega_final/
├── Datos_clase/
│   ├── P_1930_1979.xlsx
│   ├── P_1980_1990.csv
│   ├── P_1991_2003.csv
│   ├── CVC_Wheather_Stations_Long_Lat.csv
│   ├── CVC_Wheather_Stations_Long_Lat.dta
│   ├── CatalogoEstacionesElevacion.xlsx
│   └── ImagenesSatelitePentadal/        (R usa las péntadas de 1981 a 2003)
├── Datos_externos/
├── Proyecto_final_AnalíticaDeDatos.R
├── Proyecto_final_AnalíticaDeDatos.py
└── ...
```

- R exige todos los archivos al inicio. Si falta alguno, se detiene y lista los que faltan.
- Python solo lee `CVC_Wheather_Stations_Long_Lat.csv`, para los nombres de las estaciones.
- Si el docente tiene una lista de estaciones reservadas, se pone en `Datos_clase/estaciones_reservadas.csv` (columna `station_id`) y R las saca del entrenamiento. Sin ese archivo no se excluye ninguna.

### 2. Preparar los entornos

**R:** R 4.6.1 con readxl 1.5.0, terra 1.9-46, fda 6.3.0, fda.usc 2.2.0 y foreign 0.8-91. El script se detiene al inicio si falta un paquete y lista lo que falta.

**Python:** Python 3.12 con las versiones de `requirements.txt`. torch se instala en su versión para CPU:

```
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install torch==2.15.0 --index-url https://download.pytorch.org/whl/cpu
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Con otras versiones los scripts corren igual, pero las cifras pueden cambiar en los últimos decimales. Las semillas están fijas en los dos scripts.

### 3. Correr en este orden

Los dos scripts usan rutas relativas a esta carpeta. Cada uno se corre completo, uno después del otro:

1. **R:** abrir `Proyecto_final.Rproj` en RStudio y correr `Proyecto_final_AnalíticaDeDatos.R` completo. Desde la consola, en esta carpeta: `Rscript Proyecto_final_AnalíticaDeDatos.R`.
   - Tarda unos 6 minutos.
   - Escribe `Datos_preparados/` (lo que usa Python), `Figuras_proyecto/` y `Predicciones_R/`.
2. **Python:** `python Proyecto_final_AnalíticaDeDatos.py`.
   - Tarda unos 15 minutos.
   - Lee `Datos_preparados/` y escribe `Predicciones_python/` y `Resultados_python/`.
3. **R otra vez:** `Proyecto_final_AnalíticaDeDatos.R` completo, o desde su sección "5. EVALUACIÓN".
   - R lee las predicciones de Python y evalúa todos los modelos con las mismas métricas: L2 relativo, KGE, sesgo del acumulado anual y fecha de los dos picos.
   - Los intervalos son por bootstrap de estaciones, con diferencias pareadas y prueba de Wilcoxon por estación.

Las carpetas de salida se crean solas. Si todavía no existen las predicciones de Python, R corre igual y evalúa solo sus modelos.

Al terminar, además de lo que trae la entrega, la carpeta tiene:

```
├── Datos_preparados/     curvas, tabla de modelado, FPCA y matrices de Gram (de R para Python)
├── Figuras_proyecto/     figuras de R (PDF)
└── Resultados_python/    verificaciones, pesos, tablas y figuras de la red
```

**Librerías:**
- **R:** solo los paquetes de arriba. Las bases, los mínimos cuadrados y la matriz de suavizado del ítem 1 están programados a mano y comprobados contra fda.
- **Python:** la red usa PyTorch y pykan, las librerías permitidas en el ítem 4. numpy, pandas y matplotlib manejan los datos y las figuras. `scipy.linalg` se usa solo para el álgebra lineal del FPCA (Cholesky y valores propios), lo mismo que haría `numpy.linalg`.

## La prueba y los análisis complementarios

- **Partición por años:** entrenamiento 1981-1995, validación 1996-2000 y prueba 2001-2003.
- **Elección con validación:** los modelos, los hiperparámetros y las semillas se eligieron solo con validación.
- **Reglas de la evaluación:** están escritas al inicio de la sección 5 del script de R (`Proyecto_final_AnalíticaDeDatos.R:3343`).

**Estaciones reservadas.** No hay lista de estaciones reservadas.
- El script de R conserva el filtro (un vector vacío que no quita nada).
- Muestra cómo se predice una estación nueva solo con su píxel, su altitud y el ONI (`Predecir_ff_nuevo`, `Proyecto_final_AnalíticaDeDatos.R:3659`), con una comprobación sobre validación.
- Cómo generaliza la red a estaciones que no vio se mide sacando 30 estaciones del entrenamiento: `Resultados_python/validacion_por_estaciones.csv`, generado por `validacion_por_estaciones` (`Proyecto_final_AnalíticaDeDatos.py:1740`).

Seis análisis complementarios no intervienen en la selección de modelos ni cambian las predicciones principales:

1. **Diagnóstico frente a la media de entrenamiento** (script de R, sección 5, después de la evaluación de prueba). Mide si la predicción es solo la climatología: compara el error con el de la media, el apartamiento de la media y el nivel anual.
2. **Análisis de sensibilidad de los datos.** Se corre con `--sensibilidad` y escribe en `Sensibilidad/`.
   - Imputa con estaciones vecinas los pares de entrenamiento incompletos e invalida los picos que las vecinas no respaldan.
   - Las conclusiones principales se mantienen, incluido el empate entre los mejores modelos y entre las tres cabezas. Lo único frágil es el modelo función a función que elige la validación: con los datos de la sensibilidad, en prueba deja de mejorar al satélite.
3. **Experimento de arquitectura.** Se corre con `--arquitectura` y usa solo validación.
   - Compara anchos, profundidad y rejilla de la KAN, y la MLP de una y de dos capas.
   - Ninguna variante mejora a la arquitectura principal más allá de la DE entre semillas. La mejor, la KAN [8, 16, 6], con el doble de parámetros, queda 0,0025 por debajo, frente a una DE de 0,010. Una segunda capa oculta no ayuda.
4. **Interpretación del flujo entre capas de la KAN:** diagrama de la red, efecto de cada nodo en la curva y ablación.
5. **Desfase entre predicción y observación** (script de R, después del diagnóstico frente a la media).
   - Es el adelanto de máxima correlación de cada par.
   - Las medianas son de pocos días, y la dispersión es la misma de la media de entrenamiento: no hay un desfase sistemático.
6. **Colinealidad de las salidas del codificador** (`Resultados_python/colinealidad_z.csv`). Explica por qué los aportes se compensan dentro de las cabezas y por qué se interpretan efectos agregados, no neuronas sueltas.

Para la sensibilidad, el orden es este:

```
Rscript Proyecto_final_AnalíticaDeDatos.R --sensibilidad
python Proyecto_final_AnalíticaDeDatos.py --sensibilidad
Rscript Proyecto_final_AnalíticaDeDatos.R --sensibilidad
```

## Arquitectura de los modelos

Cada figura y tabla del informe indica en su pie el objeto o el archivo que la genera y el script donde está.

### Flujo de la entrega

```
Datos_clase/  +  Datos_externos/
      │
      ▼
┌─ 1. Proyecto_final_AnalíticaDeDatos.R (primera corrida)
│    Ítem 4.1   depuración con 5 reglas, cobertura mínima y trazabilidad
│    Ítem 1     curvas B-spline por mínimos cuadrados (K = 8 estaciones, K = 11 satélite)
│    Ítem 2     análisis exploratorio funcional
│    Ítem 3     comparación estación–píxel
│    Ítem 4     escala raíz, FPCA de las estaciones (M = 6) y regresiones funcionales
└──────────────────────────────────────────────
      │   Datos_preparados/  (curvas, tabla de modelado, FPCA y matrices de Gram)
      ▼
┌─ 2. Proyecto_final_AnalíticaDeDatos.py
│    red funcional con tres cabezas: lineal, MLP y KAN
│    18 configuraciones × 3 semillas por cabeza, elegidas solo con validación
│    equivalencia de la cabeza lineal con la regresión e interpretación de la KAN
└──────────────────────────────────────────────
      │   Predicciones_python/  (Red_lineal.csv, Red_MLP.csv, Red_KAN.csv)
      ▼
┌─ 3. Proyecto_final_AnalíticaDeDatos.R (segunda corrida, sección 5)
│    evalúa todos los modelos con las mismas cuatro métricas:
│    L² relativo, KGE, sesgo del acumulado anual y error en la fecha de los dos picos
│    IC bootstrap por estaciones, Wilcoxon pareado y resultados por fase ENSO
└──────────────────────────────────────────────
```

### Red funcional (Python)

El codificador y el decodificador son los mismos en los tres modelos; solo cambia la cabeza, y esa diferencia es la ablación. La red trabaja en la escala de la raíz, z(t) = √x(t).

```
Curva del satélite x(t)  ──►  16 coeficientes c en la base B-spline de la raíz
      │
      ▼
┌─ CODIFICADOR FUNCIONAL (igual en los tres modelos)
│    6 neuronas funcionales; cada peso w_k(t) en la misma base, con norma L² igual a 1
│    z_k = (1/√365) ∫ w_k(t) x(t) dt = (1/√365) a_kᵀ W_sat c
│    W_sat: matriz de Gram de la base, por cuadratura de Gauss–Legendre
│    se agregan la altitud y el ONI estandarizados  ──►  8 entradas
└──────────────────────────────────────────────
      │
      ▼
┌─ CABEZA (lo único que cambia entre los modelos)
│    Lineal   u = G z + g                                        54 parámetros
│    MLP      8 → 74 neuronas (SiLU) → 6                       1116 parámetros
│    KAN      [8, 8, 6]; cada arista φ(x) = w_b·SiLU(x) + w_s·Σ c_i B_i(x),
│             B-splines cúbicas de rejilla 5                   1120 parámetros (0,36 % más que la MLP)
└──────────────────────────────────────────────
      │   u: 6 puntajes estandarizados
      ▼
┌─ DECODIFICADOR FPCA (igual en los tres modelos, congelado)
│    ẑ(t) = μ(t) + Σ_m √ρ_m u_m ψ_m(t),   m = 1, …, 6   (92,3 % de la varianza)
│    μ y ψ_m: FPCA de las 1679 curvas de estación de entrenamiento
└──────────────────────────────────────────────
      │
      ▼
Curva de la estación en mm/día:  ŷ(t) = máx(ẑ(t), 0)²   (no negativa por construcción)
```

- **Pérdida:** norma L² funcional con la matriz de Gram de la base de las estaciones, más la rugosidad de los pesos: (1/(365 n)) Σ_i (ĉ_i − c_i)ᵀ W_est (ĉ_i − c_i) + (λ/6) Σ_k a_kᵀ Ω a_k.
- **Entrenamiento:** el mismo presupuesto para las tres cabezas. AdamW, hasta 3000 épocas y parada temprana con la pérdida de validación (revisada cada 5 épocas, paciencia de 300). La rejilla es K ∈ {6, 8} neuronas, λ ∈ {0,01; 0,1; 1} y tasa de aprendizaje ∈ {0,001; 0,003; 0,01}, con 3 semillas por configuración.

### Modelos de referencia (R)

Todos predicen la curva de la estación en la escala de la raíz a partir del satélite, la altitud y el ONI, salvo las dos líneas base. La red y las regresiones se evalúan juntas en la sección 5 del script de R.

```
Líneas base
├─ Satélite sin corregir     ŷ(t) = x_sat(t), la curva del píxel
└─ Media de entrenamiento    ŷ(t) = curva media de las estaciones en entrenamiento (no usa el satélite)

Regresiones funcionales (escala raíz: z = √x)
├─ Concurrente ampliado      z(t) = b0(t) + b1(t) z_sat(t) + g1(t) altitud + g2(t) ONI
│                            rugosidad penalizada en todos los coeficientes
├─ Función a función         z(t) = b0(t) + ∫ b1(s, t) z_sat(s) ds + g1(t) altitud + g2(t) ONI
│                            rugosidad penalizada de b1 en s y en t
└─ Rango reducido            función a función con b1(s, t) de rango 6 sobre el FPCA de las estaciones
                             = la red con cabeza lineal (verificado numéricamente, Tabla 12 del informe)
```

## Fuentes de datos

- **Datos de la clase:** lluvia diaria de la CVC, coordenadas, catálogo de elevación y péntadas satelitales.
- **Límite del Valle del Cauca:** GADM 4.1, nivel 1 (https://gadm.org).
- **Altitud de las estaciones sin dato en el catálogo:** Copernicus DEM GLO-30 (https://copernicus-dem-30m.s3.amazonaws.com/).
- **ONI:** NOAA CPC (https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt).
