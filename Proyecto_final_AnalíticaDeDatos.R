#Juan David Reyes Burbano
#María José Suárez Toro
#Marco Antonio Aragón
#Maria Paula Pinillo
#Laura Daniela Astudillo

#Limpiar área de trabajo
rm(list=ls())

#Modo de sensibilidad (Rscript Proyecto_final_AnalíticaDeDatos.R --sensibilidad): regla de picos aislados e imputación de pares de
#entrenamiento incompletos. Escribe en Sensibilidad/ y no interviene en la selección de modelos
Modo_sensibilidad <- "--sensibilidad" %in% commandArgs(trailingOnly = TRUE)     #En RStudio: cambiar a TRUE a mano
Ruta <- function(x) if (Modo_sensibilidad) file.path("Sensibilidad", x) else x
if (Modo_sensibilidad) dir.create("Sensibilidad", showWarnings = FALSE)

#install.packages(c("fda", "fda.usc"))

#Verificación de paquetes y carga de librerías (foreign se usa como foreign::read.dta)
Paquetes   <- c("readxl", "terra", "fda", "fda.usc", "foreign")
Faltan_paq <- Paquetes[!sapply(Paquetes, requireNamespace, quietly = TRUE)]
if (length(Faltan_paq) > 0) stop("Instalar antes de ejecutar: ", paste(Faltan_paq, collapse = ", "))
library(readxl)
library(terra)
library(fda)
library(fda.usc)

#Archivos de entrada. Externos al curso: valle_del_cauca.geojson (GADM 4.1 nivel 1, https://gadm.org),
#dem_valle_copernicus30m.tif (Copernicus DEM GLO-30) y oni.ascii.txt (NOAA CPC)
Archivos_entrada <- c("Datos_clase/P_1930_1979.xlsx",
                      "Datos_clase/P_1980_1990.csv",
                      "Datos_clase/P_1991_2003.csv",
                      "Datos_clase/CVC_Wheather_Stations_Long_Lat.csv",
                      "Datos_clase/CVC_Wheather_Stations_Long_Lat.dta",
                      "Datos_clase/CatalogoEstacionesElevacion.xlsx",
                      "Datos_clase/ImagenesSatelitePentadal",
                      "Datos_externos/valle_del_cauca.geojson",
                      "Datos_externos/dem_valle_copernicus30m.tif",
                      "Datos_externos/oni.ascii.txt")
#ONI: descarga si falta y control de versión (md5 de la copia usada en el informe)
Archivo_oni <- "Datos_externos/oni.ascii.txt"
if (!file.exists(Archivo_oni)) {
  dir.create("Datos_externos", showWarnings = FALSE)
  tryCatch(download.file("https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt", Archivo_oni, quiet = TRUE),
           error = function(e) NULL, warning = function(w) NULL)
  if (file.exists(Archivo_oni) && file.size(Archivo_oni) == 0) file.remove(Archivo_oni)
}
Faltan_arch <- Archivos_entrada[!file.exists(Archivos_entrada)]
if (length(Faltan_arch) > 0) stop("Faltan archivos de entrada:\n  ", paste(Faltan_arch, collapse = "\n  "))
if (unname(tools::md5sum(Archivo_oni)) != "91610417cb5414ca970752d82e865dad")
  warning("El ONI no es la copia usada en el informe (otra versión de la NOAA): las cifras del ENSO pueden variar")

#Carpetas de salida (figuras e insumos para Python)
dir.create(Ruta("Figuras_proyecto"), showWarnings = FALSE)

#Función para mostrar en pantalla una figura ya guardada en PDF (solo en sesión interactiva)
Mostrar <- function(dibujo) {
  if (interactive()) tryCatch(dibujo, error = function(e) message("Panel de gráficos muy pequeño: la figura quedó en el PDF"))
  invisible(NULL)
}

dir.create(Ruta("Datos_preparados"), showWarnings = FALSE)

#Lectura de los tres archivos de lluvia
P_1930_1979 <- as.data.frame(read_excel("Datos_clase/P_1930_1979.xlsx", sheet = "Hoja1"))
P_1980_1990 <- read.csv("Datos_clase/P_1980_1990.csv")
P_1991_2003 <- read.csv("Datos_clase/P_1991_2003.csv")

#Registros, estaciones y años de cada archivo
rbind(Registros  = c(P_1930_1979 = nrow(P_1930_1979), P_1980_1990 = nrow(P_1980_1990), P_1991_2003 = nrow(P_1991_2003)),
      Estaciones = c(length(unique(P_1930_1979$STATION_ID)), length(unique(P_1980_1990$STATION_ID)),
                     length(unique(P_1991_2003$STATION_ID))),
      Desde      = c(min(P_1930_1979$YEAR), min(P_1980_1990$YEAR), min(P_1991_2003$YEAR)),
      Hasta      = c(max(P_1930_1979$YEAR), max(P_1980_1990$YEAR), max(P_1991_2003$YEAR)))

#Columnas extra del archivo histórico: un único valor cada una (Pluviometer, mm)
table(P_1930_1979$SENSORNAME)
table(P_1930_1979$UNITS)

#Mismas cinco columnas en los tres archivos
P_1930_1979 <- P_1930_1979[, c("STATION_ID", "YEAR", "MONTH", "DAY", "VALUE")]

#Unión de los tres periodos
Estaciones <- rbind(P_1930_1979, P_1980_1990, P_1991_2003)

#Fechas del calendario (las inexistentes, como un 30 de febrero, quedan como NA)
Fecha <- as.Date(paste(Estaciones$YEAR, Estaciones$MONTH, Estaciones$DAY, sep = "-"))
#Verificación de la tabla unida: registros, duplicados estación-día, fechas inexistentes, faltantes, estaciones y años
c(Registros = nrow(Estaciones), Duplicados = sum(duplicated(Estaciones[, c("STATION_ID", "YEAR", "MONTH", "DAY")])),
  Fechas_NA = sum(is.na(Fecha)), Faltantes = sum(is.na(Estaciones)), Estaciones = length(unique(Estaciones$STATION_ID)),
  Desde = min(Estaciones$YEAR), Hasta = max(Estaciones$YEAR))

#Código de estación como texto
Estaciones$STATION_ID <- as.character(Estaciones$STATION_ID)

#Días registrados por par estación-año
Dias <- aggregate(VALUE ~ STATION_ID + YEAR, data = Estaciones, FUN = length)
names(Dias)[3] <- "N_dias"

#Días esperados por año (366 en bisiestos)
Bisiesto <- (Dias$YEAR %% 4 == 0 & Dias$YEAR %% 100 != 0) | Dias$YEAR %% 400 == 0
Dias$Dias_esperados <- ifelse(Bisiesto, 366, 365)
Dias$Faltan <- Dias$Dias_esperados - Dias$N_dias

#Rejilla estación x año de 1981-2003 (días registrados por par)
Periodo <- Estaciones[Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 2003, ]
Rejilla <- table(Periodo$STATION_ID, Periodo$YEAR)

#PREPROCESAMIENTO (ítem 4.1)

#Columna de trabajo depurada (VALUE queda como dato crudo)
Estaciones$Lluvia <- Estaciones$VALUE
#Motivo de invalidación de cada registro ("Valido" si ninguna regla aplica)
Estaciones$Marca  <- "Valido"

#Estaciones reservadas por el docente (Datos_clase/estaciones_reservadas.csv, columna station_id). Sin el
#archivo el vector queda vacío y los filtros "%in% Reservadas" no excluyen nada
Reservadas <- character(0)
if (file.exists("Datos_clase/estaciones_reservadas.csv"))
  Reservadas <- read.csv("Datos_clase/estaciones_reservadas.csv", colClasses = "character")$station_id
c(Reservadas = length(Reservadas), En_los_datos = sum(Reservadas %in% Estaciones$STATION_ID))

#Coordenadas de las estaciones (archivo del curso, códigos de 10 dígitos, Latin-1) y límite del Valle del Cauca
Coord_est <- read.csv("Datos_clase/CVC_Wheather_Stations_Long_Lat.csv", fileEncoding = "latin1",
                      colClasses = c(CODIGO = "character"))
#Verificación de coordenadas: un código por fila, 10 dígitos y sin faltantes
c(Filas = nrow(Coord_est), Codigos = length(unique(Coord_est$CODIGO)), Con_10_digitos = sum(nchar(Coord_est$CODIGO) == 10),
  Sin_coord = sum(is.na(Coord_est$Lon) | is.na(Coord_est$lat)))
#Verificación frente a la versión .dta del curso: mismas estaciones y coordenadas
Coord_dta <- foreign::read.dta("Datos_clase/CVC_Wheather_Stations_Long_Lat.dta")
Pos_dta   <- match(Coord_est$CODIGO, sprintf("%.0f", Coord_dta$station_id))
c(Faltan_en_dta = sum(is.na(Pos_dta)), Dif_lon = max(abs(Coord_est$Lon - Coord_dta$lon[Pos_dta])),
  Dif_lat = max(abs(Coord_est$lat - Coord_dta$lat[Pos_dta])))

#Límite departamental (GADM 4.1 nivel 1)
Valle <- vect("Datos_externos/valle_del_cauca.geojson")

#Estaciones con coordenadas dentro del polígono
Pts_est <- vect(data.frame(lon = Coord_est$Lon, lat = Coord_est$lat), geom = c("lon", "lat"), crs = "EPSG:4326")
Coord_est$Dentro_Valle <- is.related(Pts_est, Valle, relation = "intersects")

#Estaciones sin coordenadas: se ubican por el municipio del catálogo de elevación (hoja 1, se reutiliza en la altitud).
#Municipios del Valle: los de estaciones con coordenadas dentro del polígono
Elevacion <- as.data.frame(read_excel("Datos_clase/CatalogoEstacionesElevacion.xlsx", sheet = "Hoja1"))
Elevacion$Municipio <- toupper(trimws(Elevacion$municipio))
Con_coord_cat <- Elevacion$codigo %in% Coord_est$CODIGO
Dentro_cat    <- Elevacion$codigo %in% Coord_est$CODIGO[Coord_est$Dentro_Valle]
Mun_valle <- unique(Elevacion$Municipio[Dentro_cat])
Mun_fuera <- unique(Elevacion$Municipio[Con_coord_cat & !Dentro_cat])
#Validación: ningún municipio con estaciones a ambos lados del límite
c(Con_coord_y_catalogo = sum(Con_coord_cat), Municipios_dentro = length(Mun_valle), Municipios_fuera = length(Mun_fuera),
  Municipios_en_ambos = length(intersect(Mun_valle, Mun_fuera)))
#Estaciones sin coordenadas de municipios del Valle: entran a los ítems 1 y 2, no a los ítems 3 y 4 (sin píxel)
Sin_coord_valle <- Elevacion$codigo[!Con_coord_cat & Elevacion$Municipio %in% Mun_valle]
Sin_coord_valle <- intersect(Sin_coord_valle, unique(Estaciones$STATION_ID))

#Estaciones del Valle: dentro del polígono o, sin coordenadas, con municipio del Valle
Estaciones$Dentro_Valle <- Estaciones$STATION_ID %in% c(Coord_est$CODIGO[Coord_est$Dentro_Valle], Sin_coord_valle)

#Función para contar estaciones, pares y días válidos de 1981-2003 (tabla de trazabilidad)
Contar_validos <- function(Paso) {
  Validos <- Estaciones[!is.na(Estaciones$Lluvia) & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 2003, ]
  Pares   <- unique(Validos[, c("STATION_ID", "YEAR")])
  data.frame(Paso = Paso,
             Estaciones = length(unique(Pares$STATION_ID)),
             Pares = nrow(Pares),
             Dias = nrow(Validos))
}

#Días posibles en 1981-2003 (ambos extremos incluidos)
Dias_posibles <- as.numeric(as.Date("2003-12-31") - as.Date("1981-01-01")) + 1

Trazabilidad <- data.frame(Paso = "Rejilla estación x año 1981-2003 (posibles)",
                           Estaciones = nrow(Rejilla), Pares = length(Rejilla),
                           Dias = Dias_posibles * nrow(Rejilla))
Trazabilidad <- rbind(Trazabilidad, Contar_validos("Pares con al menos un registro"))

#Filtro geográfico: estaciones fuera del Valle del Cauca
Fuera_valle <- !Estaciones$Dentro_Valle
Estaciones$Lluvia[Fuera_valle] <- NA
Estaciones$Marca[Fuera_valle]  <- "Fuera_Valle"
c(Estaciones_sin_coordenadas = length(setdiff(unique(Estaciones$STATION_ID), Coord_est$CODIGO)),
  Sin_coordenadas_municipio_Valle = length(Sin_coord_valle),
  Estaciones_fuera_poligono = length(unique(Estaciones$STATION_ID[Estaciones$STATION_ID %in% Coord_est$CODIGO[!Coord_est$Dentro_Valle]])),
  Estaciones_fuera = length(unique(Estaciones$STATION_ID[Fuera_valle])),
  Registros_fuera = sum(Fuera_valle))
Trazabilidad <- rbind(Trazabilidad, Contar_validos("Dentro del Valle del Cauca (polígono GADM 4.1 o municipio del catálogo)"))

#Regla 1: valores imposibles (negativos o mayores que el récord mundial en 24 h, 1825 mm, OMM)
Limite_fisico <- 1825
Imposible <- Estaciones$Marca == "Valido" & (Estaciones$VALUE > Limite_fisico | Estaciones$VALUE < 0)
sum(Imposible)
#Valores imposibles en todos los registros: 859375 mm de la estación 2611600101, ya excluida por el filtro geográfico
Estaciones[Estaciones$VALUE > Limite_fisico | Estaciones$VALUE < 0, c("STATION_ID", "YEAR", "MONTH", "DAY", "VALUE", "Marca")]

#Invalidación: el día pasa a faltante
Estaciones$Lluvia[Imposible] <- NA
Estaciones$Marca[Imposible]  <- "Imposible"

Trazabilidad <- rbind(Trazabilidad, Contar_validos("Regla 1: valores imposibles"))

#Regla 2: totales mensuales (mes con un único registro, ubicado en su último día)
#Control de alineación del vector Fecha con Estaciones
all(as.numeric(format(Fecha, "%d")) == Estaciones$DAY)

#Registros del mes de cada fila (estación, año, mes)
Estaciones$N_mes <- ave(Estaciones$VALUE, Estaciones$STATION_ID, Estaciones$YEAR, Estaciones$MONTH, FUN = length)

#Último día del mes (el día siguiente es 1)
Ultimo_dia <- as.numeric(format(Fecha + 1, "%d")) == 1

Total_mensual <- Estaciones$N_mes == 1 & Ultimo_dia & Estaciones$Marca == "Valido"
#Evidencia: meses con un solo registro, cuántos caen en el último día y cuántos se esperarían al azar (30,44 días/mes)
Un_registro <- Estaciones$N_mes == 1
c(Meses_un_registro = sum(Un_registro), En_ultimo_dia = sum(Un_registro & Ultimo_dia),
  Esperado_al_azar = round(sum(Un_registro) / 30.44, 1))
#Periodo de esos casos (ninguno en prueba)
table(cut(Estaciones$YEAR[Un_registro & Ultimo_dia], c(1929, 1980, 1995, 2000, 2003),
          labels = c("1930-1980", "Entrenamiento", "Validación", "Prueba")))

#Totales mensuales guardados aparte (no son lluvia diaria y no se imputan)
Totales_mensuales <- Estaciones[Total_mensual, c("STATION_ID", "YEAR", "MONTH", "VALUE")]
names(Totales_mensuales)[4] <- "Total_mm"

#Invalidación en la serie diaria
Estaciones$Lluvia[Total_mensual] <- NA
Estaciones$Marca[Total_mensual]  <- "Total_mensual"

Trazabilidad <- rbind(Trazabilidad, Contar_validos("Regla 2: totales mensuales"))

#Regla 3: ceros falsos (rachas de ceros que corresponden a días sin medición)
#Parte A: fecha como columna y orden por estación y fecha
Estaciones$Fecha <- Fecha
Estaciones <- Estaciones[order(Estaciones$STATION_ID, Estaciones$Fecha), ]
#Vectores sueltos alineados con el orden anterior: se eliminan
rm(Fecha, Ultimo_dia, Imposible, Total_mensual, Un_registro)

#Parte B: rachas de ceros consecutivos
n <- nrow(Estaciones)
#Día válido con lluvia igual a cero
Cero <- !is.na(Estaciones$Lluvia) & Estaciones$Lluvia == 0
#La fila continúa la racha de la fila anterior si: misma estación, día siguiente y ambas en cero
Continua <- c(FALSE, Estaciones$STATION_ID[-1] == Estaciones$STATION_ID[-n] &
                as.numeric(Estaciones$Fecha[-1] - Estaciones$Fecha[-n]) == 1 &
                Cero[-1] & Cero[-n])
#Una racha empieza en un cero que no continúa la anterior; cada inicio recibe un número
Inicia   <- Cero & !Continua
Racha_id <- cumsum(Inicia)
Racha_id[!Cero] <- NA

#Tabla con una fila por racha
Rachas <- data.frame(STATION_ID = Estaciones$STATION_ID[Inicia],
                     Inicio     = Estaciones$Fecha[Inicia],
                     Largo      = as.vector(table(Racha_id)))
#Punto medio de la racha: define su mes y su año
Rachas$Medio <- Rachas$Inicio + floor((Rachas$Largo - 1) / 2)
Rachas$Mes   <- as.numeric(format(Rachas$Medio, "%m"))
Rachas$Año   <- as.numeric(format(Rachas$Medio, "%Y"))

#Control: los días de las rachas suman todos los ceros
sum(Rachas$Largo) == sum(Cero)

#Parte C: condición 1, racha más larga que el cuantil 99,9 % de las rachas de entrenamiento del mismo mes
Entrena <- Rachas$Año >= 1981 & Rachas$Año <= 1995 & !(Rachas$STATION_ID %in% Reservadas)
Largo_critico <- tapply(Rachas$Largo[Entrena], Rachas$Mes[Entrena], quantile, probs = 0.999, names = FALSE)
round(Largo_critico, 1)

#Solo rachas de 1981-2003 (red densa). Control: mínimo de estaciones con dato válido en un día de 1981-2003
min(table(Estaciones$Fecha[!is.na(Estaciones$Lluvia) & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 2003]))
Rachas$Cond1 <- Rachas$Año >= 1981 & Rachas$Año <= 2003 & Rachas$Largo > Largo_critico[Rachas$Mes]
sum(Rachas$Cond1)

#Parte D: señal regional, fracción diaria de estaciones con lluvia (sin reservadas)
Region <- aggregate(Lluvia ~ Fecha, data = Estaciones[!(Estaciones$STATION_ID %in% Reservadas), ], FUN = function(x) mean(x > 0))
names(Region)[2] <- "Frac_lluvia"

#Control: días de entrenamiento completos (periodos de días consecutivos)
Region_ent <- Region[Region$Fecha >= as.Date("1981-01-01") & Region$Fecha <= as.Date("1995-12-31"), ]
nrow(Region_ent) == as.numeric(as.Date("1995-12-31") - as.Date("1981-01-01")) + 1

#Sumas acumuladas de la fracción diaria (la suma de un periodo es una resta)
Sumas_ent <- c(0, cumsum(Region_ent$Frac_lluvia))

#Función para calcular el percentil de la lluvia regional de una racha frente a los periodos de entrenamiento del mismo largo y mes
Percentil_regional <- function(Inicio, Largo, Mes) {
  #Fracción regional media durante la racha
  F_racha <- mean(Region$Frac_lluvia[Region$Fecha >= Inicio & Region$Fecha < Inicio + Largo])
  #Media de cada periodo de entrenamiento con ese largo (el periodo i empieza en el día i)
  k <- nrow(Region_ent) - Largo + 1
  Medias <- (Sumas_ent[(Largo + 1):(Largo + k)] - Sumas_ent[1:k]) / Largo
  #Solo los periodos cuyo punto medio cae en el mismo mes
  Punto_medio <- Region_ent$Fecha[1:k] + floor((Largo - 1) / 2)
  Referencia  <- Medias[as.numeric(format(Punto_medio, "%m")) == Mes]
  #Proporción de periodos de referencia igual de secos o más secos que la racha
  mean(Referencia <= F_racha)
}

#Caso de control: mayo de 1990, estación 2612710103
Percentil_regional(as.Date("1990-04-30"), 45, 5)

#Parte E: condición 2 y marcado
Rachas$Percentil <- NA
Rachas$Percentil[Rachas$Cond1] <- mapply(Percentil_regional,
                                         Rachas$Inicio[Rachas$Cond1], Rachas$Largo[Rachas$Cond1], Rachas$Mes[Rachas$Cond1])

#Si la racha cae en el 10 % de periodos más secos de su época, la región estaba seca y la racha se conserva
alfa_region <- 0.10
Rachas$Cero_falso <- Rachas$Cond1 & !is.na(Rachas$Percentil) & Rachas$Percentil >= alfa_region

#Invalidar los días de las rachas marcadas (el número de cada racha es su fila en la tabla Rachas)
Filas_falsas <- which(Racha_id %in% which(Rachas$Cero_falso))
c(Rachas_cero_falso = sum(Rachas$Cero_falso), Estaciones = length(unique(Rachas$STATION_ID[Rachas$Cero_falso])),
  Dias = length(Filas_falsas))
Estaciones$Lluvia[Filas_falsas] <- NA
Estaciones$Marca[Filas_falsas]  <- "Cero_falso"

Trazabilidad <- rbind(Trazabilidad, Contar_validos("Regla 3: ceros falsos"))

#Regla 4: posibles totales acumulados. El primer registro con lluvia > 0 después de uno o más días sin registro se
#invalida (no se distingue de la lluvia de un día); los ceros después de un hueco se conservan

#Días sin registro (filas ausentes) inmediatamente antes de cada fila, dentro de la misma estación
n <- nrow(Estaciones)
Ausentes_antes <- c(NA, ifelse(Estaciones$STATION_ID[-1] == Estaciones$STATION_ID[-n],
                               as.numeric(Estaciones$Fecha[-1] - Estaciones$Fecha[-n]) - 1, NA))

Acumulado <- !is.na(Ausentes_antes) & Ausentes_antes >= 1 & Estaciones$Marca == "Valido" & Estaciones$Lluvia > 0
#Evidencia (entrenamiento): proporción de días con lluvia sin hueco previo y tras un día ausente; tras 6 a 10 días
#ausentes, fracción de registros por encima del percentil 99 de los días lluviosos sin hueco previo
Ent4 <- Estaciones$Marca == "Valido" & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 1995 & !is.na(Ausentes_antes) &
        !(Estaciones$STATION_ID %in% Reservadas)
c(Sin_hueco = mean(Estaciones$Lluvia[Ent4 & Ausentes_antes == 0] > 0),
  Tras_1_dia = mean(Estaciones$Lluvia[Ent4 & Ausentes_antes == 1] > 0), N_tras_1_dia = sum(Ent4 & Ausentes_antes == 1))
P99_sin_hueco <- quantile(Estaciones$Lluvia[Ent4 & Ausentes_antes == 0 & Estaciones$Lluvia > 0], 0.99, names = FALSE)
Tras_6_10     <- Ent4 & Ausentes_antes >= 6 & Ausentes_antes <= 10 & Estaciones$Lluvia > 0
c(P99 = P99_sin_hueco, N = sum(Tras_6_10), Mediana = median(Estaciones$Lluvia[Tras_6_10]),
  Sobre_P99 = mean(Estaciones$Lluvia[Tras_6_10] > P99_sin_hueco))

#Guardar aparte, con el periodo que cubre (días ausentes + el día del registro)
Totales_acumulados <- data.frame(STATION_ID = Estaciones$STATION_ID[Acumulado],
                                 Desde      = Estaciones$Fecha[Acumulado] - Ausentes_antes[Acumulado],
                                 Hasta      = Estaciones$Fecha[Acumulado],
                                 Dias       = Ausentes_antes[Acumulado] + 1,
                                 Total_mm   = Estaciones$VALUE[Acumulado])

#Invalidar en la serie diaria
Estaciones$Lluvia[Acumulado] <- NA
Estaciones$Marca[Acumulado]  <- "Total_acumulado"

c(Acumulados = sum(Acumulado),
  En_1981_2003 = sum(Totales_acumulados$Hasta >= as.Date("1981-01-01") & Totales_acumulados$Hasta <= as.Date("2003-12-31")))
Trazabilidad <- rbind(Trazabilidad, Contar_validos("Regla 4: posibles totales acumulados tras hueco"))

#Regla 5: meses anómalos. Un mes con al menos 10 días válidos en el que todos superan el percentil 95 de los días
#lluviosos de su estación (percentil de entrenamiento) se invalida completo
Lluviosos_ent <- !is.na(Estaciones$Lluvia) & Estaciones$Lluvia > 0 & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 1995 &
                 !(Estaciones$STATION_ID %in% Reservadas)
P95_estacion  <- tapply(Estaciones$Lluvia[Lluviosos_ent], Estaciones$STATION_ID[Lluviosos_ent], quantile, probs = 0.95, names = FALSE)
#Días válidos de estaciones con referencia y días que superan su percentil 95
Valido_P95 <- !is.na(Estaciones$Lluvia) & Estaciones$STATION_ID %in% names(P95_estacion)
Sobre_P95  <- Valido_P95 & Estaciones$Lluvia > P95_estacion[Estaciones$STATION_ID]
Sobre_P95[is.na(Sobre_P95)] <- FALSE
#Conteos por mes (estación, año, mes), repetidos en cada fila
Mes_clave     <- paste(Estaciones$STATION_ID, Estaciones$YEAR, Estaciones$MONTH)
N_validos_mes <- ave(as.numeric(Valido_P95), Mes_clave, FUN = sum)
N_sobre_mes   <- ave(as.numeric(Sobre_P95), Mes_clave, FUN = sum)
Mes_anomalo   <- Valido_P95 & N_validos_mes >= 10 & N_sobre_mes == N_validos_mes
sum(Mes_anomalo)
#Control: en las estaciones del Valle ningún mes cumple la condición

#Control: fracciones más altas de días sobre el percentil 95 en meses con al menos 10 días válidos (1981-2003)
En_periodo_P95 <- Valido_P95 & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 2003
Fraccion_mes   <- tapply(Sobre_P95[En_periodo_P95], Mes_clave[En_periodo_P95], mean)
Validos_mes    <- tapply(Valido_P95[En_periodo_P95], Mes_clave[En_periodo_P95], sum)
round(sort(Fraccion_mes[Validos_mes >= 10], decreasing = TRUE)[1:3], 3)

#Invalidar en la serie diaria
Estaciones$Lluvia[Mes_anomalo] <- NA
Estaciones$Marca[Mes_anomalo]  <- "Mes_anomalo"
Trazabilidad <- rbind(Trazabilidad, Contar_validos("Regla 5: meses anómalos"))

#Picos aislados (no se invalidan): valores que más que duplican el récord de entrenamiento de su estación (para un
#valor de entrenamiento, el récord de los demás días). Los pares que los contienen se marcan para el análisis de sensibilidad
Validos_ent <- !is.na(Estaciones$Lluvia) & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 1995 & !(Estaciones$STATION_ID %in% Reservadas)
Record_1 <- tapply(Estaciones$Lluvia[Validos_ent], Estaciones$STATION_ID[Validos_ent], max)
Record_2 <- tapply(Estaciones$Lluvia[Validos_ent], Estaciones$STATION_ID[Validos_ent], function(x) sort(x, decreasing = TRUE)[2])
Es_record <- Validos_ent & Estaciones$Lluvia == Record_1[Estaciones$STATION_ID]
Es_record[is.na(Es_record)] <- FALSE
Ref_pico <- ifelse(Es_record, Record_2[Estaciones$STATION_ID], Record_1[Estaciones$STATION_ID])
#Récord de referencia positivo (estaciones con casi todo en cero no permiten comparar)
Pico <- !is.na(Estaciones$Lluvia) & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 2003 &
        !is.na(Ref_pico) & Ref_pico > 0 & Estaciones$Lluvia > 2 * Ref_pico
Picos <- data.frame(STATION_ID = Estaciones$STATION_ID[Pico], Fecha = Estaciones$Fecha[Pico],
                    Lluvia = Estaciones$Lluvia[Pico], Record_referencia = Ref_pico[Pico])
#Contexto regional: fracción de estaciones con lluvia ese día y su percentil frente al mismo mes en entrenamiento
Picos$Frac_region <- Region$Frac_lluvia[match(Picos$Fecha, Region$Fecha)]
Mes_region_ent <- as.numeric(format(Region_ent$Fecha, "%m"))
Picos$Percentil_region <- mapply(function(f, d) mean(Region_ent$Frac_lluvia[Mes_region_ent == as.numeric(format(d, "%m"))] <= f),
                                 Picos$Frac_region, Picos$Fecha)
Picos

#SENSIBILIDAD (picos): se invalida un pico si al menos 2 estaciones de la misma vertiente (primer dígito del código) a
#<= 20 km tienen dato ese día y el pico supera 4 veces el máximo de ellas. Los picos sin coordenadas se conservan
#Función de distancia en km entre estaciones (aproximación equirrectangular; filas id1, columnas id2)
Dist_km <- function(id1, id2) {
  a <- match(id1, Coord_est$CODIGO); b <- match(id2, Coord_est$CODIGO); lat0 <- mean(Coord_est$lat) * pi / 180
  sqrt(outer(Coord_est$Lon[a], Coord_est$Lon[b], "-")^2 * cos(lat0)^2 + outer(Coord_est$lat[a], Coord_est$lat[b], "-")^2) * 111.32
}
if (Modo_sensibilidad) {
  Picos$Vecinas_20km <- NA_integer_; Picos$Max_vecinas <- NA_real_
  for (i in seq_len(nrow(Picos))) {
    s <- Picos$STATION_ID[i]
    if (!(s %in% Coord_est$CODIGO)) next
    Mismo_dia <- Estaciones[Estaciones$Fecha == Picos$Fecha[i] & !is.na(Estaciones$Lluvia) & Estaciones$STATION_ID != s &
                            Estaciones$STATION_ID %in% Coord_est$CODIGO & substr(Estaciones$STATION_ID, 1, 1) == substr(s, 1, 1), ]
    Cerca <- Mismo_dia$Lluvia[Dist_km(s, Mismo_dia$STATION_ID)[1, ] <= 20]
    Picos$Vecinas_20km[i] <- length(Cerca)
    Picos$Max_vecinas[i]  <- if (length(Cerca) > 0) max(Cerca) else NA
  }
  Picos$Razon <- Picos$Lluvia / Picos$Max_vecinas
  Picos$Invalidar <- !is.na(Picos$Vecinas_20km) & Picos$Vecinas_20km >= 2 & Picos$Lluvia > 4 * Picos$Max_vecinas
  print(Picos[, c("STATION_ID", "Fecha", "Lluvia", "Record_referencia", "Vecinas_20km", "Max_vecinas", "Razon", "Invalidar")])
  Fila_pico <- match(paste(Picos$STATION_ID, Picos$Fecha)[Picos$Invalidar], paste(Estaciones$STATION_ID, Estaciones$Fecha))
  Estaciones$Lluvia[Fila_pico] <- NA
  Estaciones$Marca[Fila_pico]  <- "Pico_sin_respaldo"
  Picos <- Picos[!Picos$Invalidar, ]                    #Los que quedan siguen marcados como en el análisis principal
  Trazabilidad <- rbind(Trazabilidad, Contar_validos("Sensibilidad: picos sin respaldo de las vecinas"))
}

#CONSTRUCCIÓN DE DATOS FUNCIONALES (ítem 1)

#Base B-spline programada a mano adaptada del taller de repaso de datos funcionales

#Función de nodos: K funciones de orden m sobre [a, b] usan K - m + 2 cortes equiespaciados y extremos repetidos
Nodos_bspline <- function(a, b, K, m) {
  Cortes <- seq(a, b, length.out = K - m + 2)
  c(rep(a, m - 1), Cortes, rep(b, m - 1))
}

#Función de la matriz Phi (filas = tiempos, columnas = funciones base)
Base_bspline <- function(t, a, b, K, m = 4) {
  tau <- Nodos_bspline(a, b, K, m)
  L <- length(tau) - 1
  #Orden 1: indicadora de cada intervalo [tau_j, tau_j+1)
  B <- sapply(1:L, function(j) as.numeric(t >= tau[j] & t < tau[j + 1]))
  B <- matrix(B, nrow = length(t))
  #El extremo derecho b se asigna al último intervalo no vacío
  Ultimo <- max(which(tau[1:L] < tau[2:(L + 1)]))
  B[t == b, ] <- 0
  B[t == b, Ultimo] <- 1
  #Recursión: cada B de orden k combina dos B vecinas de orden k-1 (con 0/0 = 0)
  for (k in 2:m) {
    B_nueva <- matrix(0, length(t), L - k + 1)
    for (j in 1:(L - k + 1)) {
      Izq <- if (tau[j + k - 1] > tau[j]) (t - tau[j]) / (tau[j + k - 1] - tau[j]) * B[, j] else 0
      Der <- if (tau[j + k] > tau[j + 1]) (tau[j + k] - t) / (tau[j + k] - tau[j + 1]) * B[, j + 1] else 0
      B_nueva[, j] <- Izq + Der
    }
    B <- B_nueva
  }
  B
}

#Verificación de la base con K = 20 en los centros de los 365 días
t_dias <- (1:365) - 0.5
Phi_prueba <- Base_bspline(t_dias, 0, 365, K = 20)

#Verificación: igual a splines::splineDesign, partición de la unidad y no negatividad
c(Dif_splineDesign = max(abs(Phi_prueba - splines::splineDesign(Nodos_bspline(0, 365, 20, 4), t_dias, ord = 4))),
  Suma_min = min(rowSums(Phi_prueba)), Suma_max = max(rowSums(Phi_prueba)), Minimo = min(Phi_prueba))

#Tiempo: centro de cada día, reescalado a [0, 365] en años bisiestos
Estaciones$Dia_año <- as.numeric(format(Estaciones$Fecha, "%j"))
N_dias <- ifelse((Estaciones$YEAR %% 4 == 0 & Estaciones$YEAR %% 100 != 0) | Estaciones$YEAR %% 400 == 0, 366, 365)
Estaciones$t <- (Estaciones$Dia_año - 0.5) * 365 / N_dias
#Pendiente b de Var = a * media^b de la lluvia diaria (entrenamiento, grupos estación x mes con al menos 100 días):
#sin transformar, raíz y log(1 + P) en mm y en cm. Define la escala de la regla de cobertura
Ent_b   <- !is.na(Estaciones$Lluvia) & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 1995 & !(Estaciones$STATION_ID %in% Reservadas)
Grupo   <- paste(Estaciones$STATION_ID[Ent_b], Estaciones$MONTH[Ent_b])
Media_g <- tapply(Estaciones$Lluvia[Ent_b], Grupo, mean)
Var_g   <- tapply(Estaciones$Lluvia[Ent_b], Grupo, var)
N_g     <- tapply(Estaciones$Lluvia[Ent_b], Grupo, length)
Usar_g  <- Media_g > 0 & Var_g > 0 & N_g >= 100
#Pendiente de log(varianza) frente a log(media de la lluvia en mm/día)
Pendiente_b <- function(Var, Media) coef(lm(log(Var) ~ log(Media)))[[2]]
round(c(Sin_transformar = Pendiente_b(Var_g[Usar_g], Media_g[Usar_g]),
        Raiz   = Pendiente_b(tapply(sqrt(Estaciones$Lluvia[Ent_b]), Grupo, var)[Usar_g], Media_g[Usar_g]),
        Log_mm = Pendiente_b(tapply(log1p(Estaciones$Lluvia[Ent_b]), Grupo, var)[Usar_g], Media_g[Usar_g]),
        Log_cm = Pendiente_b(tapply(log1p(Estaciones$Lluvia[Ent_b] / 10), Grupo, var)[Usar_g], Media_g[Usar_g])), 3)
#Escala de la regla de cobertura: raíz (el logaritmo depende de las unidades por su "+ 1")

#Variable de trabajo de las curvas: lluvia diaria en mm (NA si fue invalidada). La raíz se aplica a la curva ya
#suavizada: suavizar raíz(P) estimaría E[raíz(P)(t)] y dejaría la lluvia de los aguaceros fuera de la curva
Estaciones$y <- Estaciones$Lluvia

#Mínimos cuadrados para un par de ejemplo (estación 2613600105, ACUEDUCTO BUGA, 1985) con K = 20
Par <- Estaciones[Estaciones$STATION_ID == "2613600105" & Estaciones$YEAR == 1985 & !is.na(Estaciones$y), ]
K_ej <- 20
Phi <- Base_bspline(Par$t, 0, 365, K_ej)

#c = (Phi^t Phi)^-1 Phi^t y
c_hat <- solve(t(Phi) %*% Phi) %*% t(Phi) %*% Par$y
#Comprobación con lsfit (R base)
max(abs(c_hat - lsfit(Phi, Par$y, intercept = FALSE)$coef))

#Matriz suavizadora H = Phi (Phi^t Phi)^-1 Phi^t: simétrica, idempotente y con traza K
H <- Phi %*% solve(t(Phi) %*% Phi) %*% t(Phi)
c(Asimetria = max(abs(H - t(H))), Idempotencia = max(abs(H %*% H - H)), Traza = sum(diag(H)))

#Control: la curva en mm conserva la lluvia observada del año; la curva de raíz(P) al cuadrado recupera solo una parte
y_hat     <- as.vector(Phi %*% c_hat)
c_hat_raiz <- solve(t(Phi) %*% Phi) %*% t(Phi) %*% sqrt(Par$y)
c(Observado = sum(Par$y), Curva_en_mm = sum(y_hat), Curva_en_raiz_al_cuadrado = sum((Phi %*% c_hat_raiz)^2))

#Selección de K con los pares completos de entrenamiento (1981-1995), lluvia en mm (comparten los días observados y Phi)
Dias_validos <- aggregate(y ~ STATION_ID + YEAR, data = Estaciones, FUN = length)
names(Dias_validos)[3] <- "N_validos"
Bis <- (Dias_validos$YEAR %% 4 == 0 & Dias_validos$YEAR %% 100 != 0) | Dias_validos$YEAR %% 400 == 0
Dias_validos$Completo <- Dias_validos$N_validos == ifelse(Bis, 366, 365)
Completos_ent <- Dias_validos[Dias_validos$Completo & Dias_validos$YEAR >= 1981 & Dias_validos$YEAR <= 1995 &
                              !(Dias_validos$STATION_ID %in% Reservadas), ]

#Matrices de datos (filas = días, columnas = pares): una para años de 365 días y otra para bisiestos
Sel <- Estaciones[paste(Estaciones$STATION_ID, Estaciones$YEAR) %in% paste(Completos_ent$STATION_ID, Completos_ent$YEAR), ]
Sel <- Sel[order(Sel$STATION_ID, Sel$YEAR, Sel$t), ]
Bis_sel <- (Sel$YEAR %% 4 == 0 & Sel$YEAR %% 100 != 0) | Sel$YEAR %% 400 == 0
Y365 <- matrix(Sel$y[!Bis_sel], nrow = 365)
Y366 <- matrix(Sel$y[Bis_sel], nrow = 366)
c(Pares_completos = nrow(Completos_ent), Años_365 = ncol(Y365), Años_366 = ncol(Y366),
  Faltantes = sum(is.na(Y365)) + sum(is.na(Y366)))
t365 <- (1:365) - 0.5
t366 <- ((1:366) - 0.5) * 365 / 366

#Función de ECM, VC (atajo con la diagonal de H) y VCG de cada curva para un K dado
Criterios_K <- function(t, Y, K) {
  Phi  <- Base_bspline(t, 0, 365, K)
  A    <- solve(t(Phi) %*% Phi) %*% t(Phi)   #Lleva los datos a los coeficientes
  s_ii <- rowSums(Phi * t(A))                #Diagonal de H = Phi A, sin construir H completa
  E    <- Y - Phi %*% (A %*% Y)              #Residuos de todas las curvas a la vez
  n    <- nrow(Phi)
  cbind(ECM = colMeans(E^2),
        VC  = colMeans((E / (1 - s_ii))^2),
        VCG = n * colSums(E^2) / (n - K)^2)
}

Ks <- 5:40
Criterios_curva <- lapply(Ks, function(K) rbind(Criterios_K(t365, Y365, K), Criterios_K(t366, Y366, K)))
Tabla_K <- data.frame(K = Ks, t(sapply(Criterios_curva, colMeans)))
round(Tabla_K, 5)
K_VC  <- Ks[which.min(Tabla_K$VC)]
K_VCG <- Ks[which.min(Tabla_K$VCG)]
c(K_VC = K_VC, K_VCG = K_VCG)

#Figura: K vs ECM, VC y VCG
pdf(Ruta("Figuras_proyecto/Fig_Seleccion_K_estaciones.pdf"), width = 10, height = 3.8)
par(mfrow = c(1, 3))
plot(Ks, Tabla_K$ECM, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "ECM", main = "ECM de ajuste")
plot(Ks, Tabla_K$VC, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "VC", main = "Validación cruzada")
abline(v = K_VC, lty = 3, col = "red")
plot(Ks, Tabla_K$VCG, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "VCG", main = "Validación cruzada generalizada")
abline(v = K_VCG, lty = 3, col = "red")
dev.off()

#Diferencia pareada de la VC de cada curva frente al K de mínima VC (K equivalentes: |z| < 2)
VC_curva <- sapply(Criterios_curva, function(x) x[, "VC"])
Dif <- VC_curva - VC_curva[, which(Ks == K_VC)]
Comparacion <- data.frame(K = Ks, Dif_media = colMeans(Dif), EE = apply(Dif, 2, sd) / sqrt(nrow(Dif)))
Comparacion$z <- Comparacion$Dif_media / Comparacion$EE
round(Comparacion[Comparacion$K <= 20, ], 5)
K_equivalentes <- Ks[abs(Comparacion$z) < 2 | is.nan(Comparacion$z)]
K_equivalentes

#K de las estaciones: mínimo de la VC, que coincide con el mínimo de la VCG
K_estaciones <- K_VC
c(K_estaciones = K_estaciones, Coincide_con_VCG = K_VC == K_VCG)

#Control de coherencia (no es criterio de selección): proporción de curvas con exactamente dos picos
t_fino <- seq(0, 365, length.out = 731)
Maximos_locales <- function(K) {
  Phi <- Base_bspline(t365, 0, 365, K)
  C   <- solve(t(Phi) %*% Phi) %*% t(Phi) %*% Y365
  X   <- Base_bspline(t_fino, 0, 365, K) %*% C
  apply(X, 2, function(x) { d <- diff(x); sum(d[-1] < 0 & d[-length(d)] >= 0) })
}
round(sapply(K_equivalentes, function(K) mean(Maximos_locales(K) == 2)), 4)
#Calibración de la cobertura por simulación (pares completos de entrenamiento de 365 días, K = K_estaciones): un hueco
#es aceptable si la curva sin esos días queda dentro de la banda de ±2 EE de la curva completa. Escala raíz
R365 <- sqrt(Y365)
R366 <- sqrt(Y366)
t_fino <- seq(0, 365, length.out = 731)
Phi_d  <- Base_bspline(t365, 0, 365, K_estaciones)    #Base en los días observados
Phi_f  <- Base_bspline(t_fino, 0, 365, K_estaciones)  #Base en una rejilla fina para comparar curvas
Inv_d  <- solve(t(Phi_d) %*% Phi_d)                  #(Phi' Phi)^-1 en los 365 días
C_cob  <- Inv_d %*% t(Phi_d) %*% R365
X_cob  <- Phi_f %*% C_cob
sigma2_cob <- colSums((R365 - Phi_d %*% C_cob)^2) / (365 - K_estaciones)
#EE(t) = raíz( sigma2 * phi(t)' (Phi'Phi)^-1 phi(t) ), una columna por curva
EE_cob <- sqrt(outer(rowSums((Phi_f %*% Inv_d) * Phi_f), sigma2_cob))

#Función de la proporción de curvas que siguen dentro de la banda si faltan los días indicados
Dentro_banda <- function(Faltan) {
  P <- Phi_d[-Faltan, ]
  #Si una base queda sin datos, su coeficiente no se puede estimar: la curva no se puede construir
  if (min(colSums(P)) < 1e-8) return(rep(FALSE, ncol(R365)))
  C <- solve(t(P) %*% P) %*% t(P) %*% R365[-Faltan, ]
  apply(abs(Phi_f %*% C - X_cob) / EE_cob, 2, max) <= 2
}

#Huecos interiores (empiezan entre los días 31 y 335, cada 5 días) y bordes (inicio y fin del año)
Largos <- c(1:6, 10, 15:25, 30, 45, 60)
Sim_interior <- sapply(Largos, function(g) mean(sapply(seq(31, 335 - g, by = 5), function(i) mean(Dentro_banda(i:(i + g - 1))))))
Sim_borde    <- sapply(Largos, function(g) mean(c(mean(Dentro_banda(1:g)), mean(Dentro_banda((366 - g):365)))))
round(data.frame(Largo = Largos, Interior = Sim_interior, Borde = Sim_borde), 3)

#Mayor hueco interior y mayor borde con al menos 95 % de curvas dentro de la banda (descriptivo; la regla es la prueba por patrón)
Hueco_max <- max(Largos[Sim_interior >= 0.95])
Borde_max <- max(Largos[Sim_borde >= 0.95])
c(Hueco_max = Hueco_max, Borde_max = Borde_max)

#Figura de la calibración
pdf(Ruta("Figuras_proyecto/Fig_Calibracion_cobertura.pdf"), width = 6, height = 5)
par(mfrow = c(1, 1))
plot(Largos, Sim_interior, type = "b", pch = 19, ylim = c(0, 1), xlab = "Días sin dato",
     ylab = "Proporción de curvas dentro de la banda", main = "Calibración de la cobertura")
lines(Largos, Sim_borde, type = "b", pch = 17, col = "red")
abline(h = 0.95, lty = 2)
legend("bottomleft", c("Hueco interior", "Borde (inicio o fin)"), col = c("black", "red"), pch = c(19, 17), lty = 1, bty = "n")
dev.off()

#Cobertura por patrón real de faltantes: a las curvas completas de entrenamiento se les quitan los días faltantes del
#par y se reajustan; el par se acepta si al menos el 95 % sigue dentro de su banda de ±2 EE

#Referencia para años bisiestos (la de 365 días ya está en Phi_d, X_cob y EE_cob)
Phi_d366      <- Base_bspline(t366, 0, 365, K_estaciones)
Inv_d366      <- solve(t(Phi_d366) %*% Phi_d366)
C_cob366      <- Inv_d366 %*% t(Phi_d366) %*% R366
X_cob366      <- Phi_f %*% C_cob366
sigma2_cob366 <- colSums((R366 - Phi_d366 %*% C_cob366)^2) / (366 - K_estaciones)
EE_cob366     <- sqrt(outer(rowSums((Phi_f %*% Inv_d366) * Phi_f), sigma2_cob366))

#Función de la proporción de curvas de referencia dentro de su banda si faltan los días indicados (N = días del año)
Prop_patron <- function(Faltan, N) {
  if (length(Faltan) == 0) return(1)
  if (N == 365) {
    P <- Phi_d[-Faltan, , drop = FALSE];    Y <- R365[-Faltan, ]; X <- X_cob;    EE <- EE_cob
  } else {
    P <- Phi_d366[-Faltan, , drop = FALSE]; Y <- R366[-Faltan, ]; X <- X_cob366; EE <- EE_cob366
  }
  #Si alguna base queda sin datos (o quedan muy pocos días), la curva no se puede construir
  if (nrow(P) < K_estaciones || min(colSums(P)) < 1e-8) return(0)
  A <- tryCatch(solve(t(P) %*% P), error = function(e) NULL)
  if (is.null(A)) return(0)
  C <- A %*% t(P) %*% Y
  mean(apply(abs(Phi_f %*% C - X) / EE, 2, max) <= 2)
}

#Días válidos de cada par de 1981-2003, borde y hueco más largo
V <- Estaciones[!is.na(Estaciones$y) & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 2003, c("STATION_ID", "YEAR", "Dia_año")]
V <- V[order(V$STATION_ID, V$YEAR, V$Dia_año), ]
Clave_par <- factor(paste(V$STATION_ID, V$YEAR), levels = unique(paste(V$STATION_ID, V$YEAR)))
Cobertura <- unique(V[, c("STATION_ID", "YEAR")])
N_año <- ifelse((Cobertura$YEAR %% 4 == 0 & Cobertura$YEAR %% 100 != 0) | Cobertura$YEAR %% 400 == 0, 366, 365)
Dias_par <- split(V$Dia_año, Clave_par)
Cobertura$Validos <- sapply(Dias_par, length)
Cobertura$Borde   <- mapply(function(d, N) max(d[1] - 1, N - d[length(d)]), Dias_par, N_año)
Cobertura$Hueco   <- sapply(Dias_par, function(d) if (length(d) > 1) max(diff(d) - 1) else 0)
#Prueba por patrón (tarda unos 30 segundos)
Cobertura$Prop    <- mapply(function(d, N) Prop_patron(setdiff(1:N, d), N), Dias_par, N_año)
Cobertura$Pasa    <- Cobertura$Prop >= 0.95
table(Cobertura$Pasa)
#Comparación con la regla de un solo hueco (hueco <= Hueco_max y borde <= Borde_max)
table(Regla_hueco_borde = Cobertura$Hueco <= Hueco_max & Cobertura$Borde <= Borde_max, Prueba_patron = Cobertura$Pasa)
#Caso de muchos huecos pequeños que juntos dañan la curva
Cobertura[Cobertura$STATION_ID == "2622100301" & Cobertura$YEAR == 1994, ]
table(cut(Cobertura$YEAR[Cobertura$Pasa], c(1980, 1995, 2000, 2003), labels = c("Entrenamiento", "Validación", "Prueba")))

#SENSIBILIDAD (imputación): pares de entrenamiento que no pasaron la cobertura, con a lo sumo 120 días sin dato y al
#menos 2 vecinas con dato en cada día faltante. Vecinas: misma vertiente, <= 20 km y |diferencia de altitud| <= 400 m
#(catálogo o DEM). Valor imputado: media de las vecinas ponderada por 1 / distancia^2, cada una escalada por la razón
#de medias diarias de entrenamiento (estación / vecina). Validación y prueba no se tocan
if (Modo_sensibilidad) {
  Alt_imp <- data.frame(id = Coord_est$CODIGO[Coord_est$Dentro_Valle], stringsAsFactors = FALSE)
  Alt_imp$alt <- Elevacion$alt[match(Alt_imp$id, Elevacion$codigo)]
  Sin_alt <- which(is.na(Alt_imp$alt))
  Pts_imp <- vect(Coord_est[match(Alt_imp$id[Sin_alt], Coord_est$CODIGO), c("Lon", "lat")], geom = c("Lon", "lat"), crs = "EPSG:4326")
  Alt_imp$alt[Sin_alt] <- extract(rast("Datos_externos/dem_valle_copernicus30m.tif"), Pts_imp, method = "bilinear")[, 2]
  Alt_imp <- Alt_imp[!is.na(Alt_imp$alt), ]
  D_imp <- pmax(Dist_km(Alt_imp$id, Alt_imp$id), 0.1)       #(0,1 km como mínimo: dos estaciones en el mismo punto)
  Vecinas_de <- function(s) {
    i <- match(s, Alt_imp$id)
    setdiff(Alt_imp$id[D_imp[i, ] <= 20 & abs(Alt_imp$alt - Alt_imp$alt[i]) <= 400 & substr(Alt_imp$id, 1, 1) == substr(s, 1, 1)], s)
  }
  #Lluvia diaria válida de entrenamiento en una matriz estación x día
  Fechas_ent <- seq(as.Date("1981-01-01"), as.Date("1995-12-31"), by = "day")
  Año_ent    <- as.numeric(format(Fechas_ent, "%Y"))
  E_imp <- Estaciones[Estaciones$STATION_ID %in% Alt_imp$id & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 1995 & !is.na(Estaciones$y), ]
  Mat_imp <- matrix(NA_real_, nrow(Alt_imp), length(Fechas_ent), dimnames = list(Alt_imp$id, NULL))
  Mat_imp[cbind(match(E_imp$STATION_ID, Alt_imp$id), match(E_imp$Fecha, Fechas_ent))] <- E_imp$y
  Media_imp <- rowMeans(Mat_imp, na.rm = TRUE)
  #Candidatos: pares de entrenamiento rechazados, de estaciones con coordenadas y altitud, con <= 120 días sin dato
  Cand <- Cobertura[!Cobertura$Pasa & Cobertura$YEAR <= 1995 & Cobertura$STATION_ID %in% Alt_imp$id &
                    !(Cobertura$STATION_ID %in% Reservadas), ]
  Cand$Faltan <- ifelse((Cand$YEAR %% 4 == 0 & Cand$YEAR %% 100 != 0) | Cand$YEAR %% 400 == 0, 366, 365) - Cand$Validos
  print(c(Rechazados_entrenamiento = sum(!Cobertura$Pasa & Cobertura$YEAR <= 1995), Con_coordenadas_y_altitud = nrow(Cand),
          Hasta_120_dias = sum(Cand$Faltan <= 120)))
  Cand <- Cand[Cand$Faltan <= 120, ]
  Imputados <- list()
  for (r in seq_len(nrow(Cand))) {
    s <- Cand$STATION_ID[r]
    v <- Vecinas_de(s); v <- v[is.finite(Media_imp[v]) & Media_imp[v] > 0]
    if (length(v) < 2) next
    Cols  <- which(Año_ent == Cand$YEAR[r])
    Falta <- Cols[is.na(Mat_imp[s, Cols])]
    Yv <- Mat_imp[v, Falta, drop = FALSE] * (Media_imp[s] / Media_imp[v])
    if (any(colSums(!is.na(Yv)) < 2)) next                #Algún día sin 2 vecinas con dato: el par no se completa
    w <- 1 / D_imp[match(s, Alt_imp$id), match(v, Alt_imp$id)]^2
    Imputados[[length(Imputados) + 1]] <- data.frame(STATION_ID = s, Fecha = Fechas_ent[Falta],
                                                     y = colSums(Yv * w, na.rm = TRUE) / colSums(!is.na(Yv) * w))
  }
  Imputados <- do.call(rbind, Imputados)
  Pares_imp <- unique(data.frame(STATION_ID = Imputados$STATION_ID, YEAR = as.numeric(format(Imputados$Fecha, "%Y"))))
  print(c(Pares_imputados = nrow(Pares_imp), Estaciones = length(unique(Pares_imp$STATION_ID)), Dias_imputados = nrow(Imputados)))
  #Los pares completados pasan la cobertura: sus días válidos no se quitan en el paso siguiente
  Cobertura$Imputado <- paste(Cobertura$STATION_ID, Cobertura$YEAR) %in% paste(Pares_imp$STATION_ID, Pares_imp$YEAR)
  Cobertura$Pasa[Cobertura$Imputado] <- TRUE
  print(table(cut(Cobertura$YEAR[Cobertura$Pasa], c(1980, 1995, 2000, 2003), labels = c("Entrenamiento", "Validación", "Prueba"))))
}

#Los días válidos de los pares que no pasan quedan fuera
No_pasa <-paste(Estaciones$STATION_ID, Estaciones$YEAR) %in% paste(Cobertura$STATION_ID, Cobertura$YEAR)[!Cobertura$Pasa]
Fuera <- No_pasa & !is.na(Estaciones$Lluvia)
Estaciones$Lluvia[Fuera] <- NA
Estaciones$y[Fuera]      <- NA
Estaciones$Marca[Fuera]  <- "Sin_cobertura"
table(Estaciones$Marca)
Trazabilidad <- rbind(Trazabilidad, Contar_validos("Cobertura por patrón real (>= 95 % dentro de la banda)"))
#SENSIBILIDAD: escritura de los días imputados. Con fila (invalidada por una regla) se reemplaza el valor; sin fila
#se agrega una con VALUE = NA y Marca = "Imputado"
if (Modo_sensibilidad) {
  Fila_imp <- match(paste(Imputados$STATION_ID, Imputados$Fecha), paste(Estaciones$STATION_ID, Estaciones$Fecha))
  Con_fila <- !is.na(Fila_imp)
  Estaciones$Lluvia[Fila_imp[Con_fila]] <- Imputados$y[Con_fila]
  Estaciones$y[Fila_imp[Con_fila]]      <- Imputados$y[Con_fila]
  Estaciones$Marca[Fila_imp[Con_fila]]  <- "Imputado"
  Nuevas <- Imputados[!Con_fila, ]
  Filas_nuevas <- Estaciones[rep(1, nrow(Nuevas)), ]
  Filas_nuevas[] <- lapply(Filas_nuevas, function(x) { x[] <- NA; x })
  Filas_nuevas$STATION_ID <- Nuevas$STATION_ID; Filas_nuevas$Fecha <- Nuevas$Fecha
  Filas_nuevas$YEAR  <- as.numeric(format(Nuevas$Fecha, "%Y")); Filas_nuevas$MONTH <- as.numeric(format(Nuevas$Fecha, "%m"))
  Filas_nuevas$DAY   <- as.numeric(format(Nuevas$Fecha, "%d")); Filas_nuevas$Dia_año <- as.numeric(format(Nuevas$Fecha, "%j"))
  Filas_nuevas$t     <- (Filas_nuevas$Dia_año - 0.5) * 365 /
                        ifelse((Filas_nuevas$YEAR %% 4 == 0 & Filas_nuevas$YEAR %% 100 != 0) | Filas_nuevas$YEAR %% 400 == 0, 366, 365)
  Filas_nuevas$Lluvia <- Nuevas$y; Filas_nuevas$y <- Nuevas$y; Filas_nuevas$Marca <- "Imputado"; Filas_nuevas$Dentro_Valle <- TRUE
  Estaciones <- rbind(Estaciones, Filas_nuevas)
  Estaciones <- Estaciones[order(Estaciones$STATION_ID, Estaciones$Fecha), ]
  print(c(Dias_con_fila = sum(Con_fila), Dias_sin_fila = nrow(Nuevas)))
  print(table(Estaciones$Marca))
  Trazabilidad <- rbind(Trazabilidad, Contar_validos("Sensibilidad: imputación por vecinas (pares de entrenamiento)"))
}

#Decisiones sobre agregar o eliminar datos (ítem 4.1)
#1. No se imputan faltantes: mínimos cuadrados ajusta cada curva con sus días válidos y la prueba por patrón asegura que,
#   con los faltantes de cada par aceptado, el 95 % de las curvas queda dentro de su banda. La imputación por vecinas
#   se evalúa en el análisis de sensibilidad. Los totales mensuales y acumulados quedan en Totales_mensuales y
#   Totales_acumulados
#2. No se agregan registros de 1930-1980: los pares de entrenamiento bastan para estimar la covarianza de los
#   coeficientes, la regla 3 necesita una red densa (tabla de abajo) y esos años no tienen satélite
#3. Picos aislados: se conservan y se marcan (Pares_curvas$Pico); la regla por vecinas se aplica solo en el análisis de
#   sensibilidad
#4. Satélite: los valores recortados en 200 mm se conservan (son un mínimo) y se marcan en Curvas_sat$Techo
#5. Filtro geográfico: estaciones con coordenadas dentro del polígono GADM 4.1 o, sin coordenadas, con municipio del
#   Valle en el catálogo de elevación
#6. Las estaciones ubicadas por municipio entran a los ítems 1 y 2, no a los ítems 3 y 4 (sin píxel)
#7. Altitud: a las estaciones con coordenadas sin altitud en el catálogo se les asigna la del DEM Copernicus GLO-30,
#   solo entre 500 y 2000 m y con pendiente menor de 15 grados (Tabla_dem); las demás no entran a la red del ítem 4
#8. Satélite: solo los píxeles con centro dentro del Valle (Pixeles_valle)
#Densidad de la red por década: mediana de estaciones del Valle con registro por día (la red que usa la regla 3)
Por_dia <- table(Estaciones$Fecha[Estaciones$Dentro_Valle])
tapply(as.vector(Por_dia), 10 * floor(as.numeric(substr(names(Por_dia), 1, 4)) / 10), median)
#Días marcados que quedaron dentro de los pares aceptados
Aceptado <- paste(Estaciones$STATION_ID, Estaciones$YEAR) %in% paste(Cobertura$STATION_ID, Cobertura$YEAR)[Cobertura$Pasa]
table(Estaciones$Marca[Aceptado])

#Curvas de estación con K = K_estaciones (mínimos cuadrados por par, lluvia diaria en mm)
Datos_curvas <- Estaciones[!is.na(Estaciones$y) & Estaciones$YEAR >= 1981 & Estaciones$YEAR <= 2003, c("STATION_ID", "YEAR", "t", "y")]
Datos_curvas <- Datos_curvas[order(Datos_curvas$STATION_ID, Datos_curvas$YEAR, Datos_curvas$t), ]
Clave_curva  <- factor(paste(Datos_curvas$STATION_ID, Datos_curvas$YEAR),
                       levels = unique(paste(Datos_curvas$STATION_ID, Datos_curvas$YEAR)))
Pares_curvas <- unique(Datos_curvas[, c("STATION_ID", "YEAR")])
rownames(Pares_curvas) <- NULL
#Marca para el análisis de sensibilidad: el par contiene un pico aislado conservado
Pares_curvas$Pico <- paste(Pares_curvas$STATION_ID, Pares_curvas$YEAR) %in% paste(Picos$STATION_ID, format(Picos$Fecha, "%Y"))
#Picos en pares aceptados (los demás pertenecen a pares que no pasaron la cobertura)
Picos$Par_aceptado <- paste(Picos$STATION_ID, format(Picos$Fecha, "%Y")) %in% paste(Pares_curvas$STATION_ID, Pares_curvas$YEAR)
Picos[, c("STATION_ID", "Fecha", "Lluvia", "Record_referencia", "Par_aceptado")]

#Exportación para Python: pares (t, y) de cada curva de estación (y en mm, t en [0, 365]), pares aceptados y base
write.csv(data.frame(station_id = Datos_curvas$STATION_ID,
                     year       = Datos_curvas$YEAR,
                     t          = Datos_curvas$t,
                     y          = Datos_curvas$y),
          Ruta("Datos_preparados/estaciones_curvas_largo.csv"), row.names = FALSE)
#ubicacion = "municipio_catalogo": estación sin coordenadas (curvas sí; emparejamiento y red no)
#reservada = TRUE: estación reservada; no entra a ningún ajuste
write.csv(data.frame(station_id = Pares_curvas$STATION_ID,
                     year       = Pares_curvas$YEAR,
                     bloque     = cut(Pares_curvas$YEAR, c(1980, 1995, 2000, 2003), labels = c("entrenamiento", "validacion", "prueba")),
                     reservada  = Pares_curvas$STATION_ID %in% Reservadas,
                     pico       = Pares_curvas$Pico,
                     ubicacion  = ifelse(Pares_curvas$STATION_ID %in% Sin_coord_valle, "municipio_catalogo", "coordenadas")),
          Ruta("Datos_preparados/estaciones_pares.csv"), row.names = FALSE)
write.csv(data.frame(nudo = Nodos_bspline(0, 365, K_estaciones, 4)),
          Ruta("Datos_preparados/estaciones_nudos.csv"), row.names = FALSE)
write.csv(data.frame(fuente = "estaciones",
                     K = K_estaciones,
                     orden = 4,
                     dominio_inf = 0,
                     dominio_sup = 365,
                     transformacion = "ninguna (la raiz se aplica a la curva: estaciones_raiz_coef.csv)",
                     unidad_respuesta = "mm/dia",
                     periodo = "1981-2003",
                     n_curvas = nrow(Pares_curvas),
                     n_observaciones = nrow(Datos_curvas),
                     archivo_datos = "estaciones_curvas_largo.csv",
                     archivo_pares = "estaciones_pares.csv",
                     archivo_nudos = "estaciones_nudos.csv"),
          Ruta("Datos_preparados/estaciones_config.csv"), row.names = FALSE)
cat("Estaciones exportadas:", nrow(Datos_curvas), "observaciones en", nrow(Pares_curvas), "curvas\n")

#Función de ajuste de una curva con sus días: c = (Phi'Phi)^-1 Phi'y y varianza residual
Ajustar_curva <- function(t, y, K) {
  Phi   <- Base_bspline(t, 0, 365, K)
  c_hat <- solve(t(Phi) %*% Phi) %*% t(Phi) %*% y
  e     <- y - Phi %*% c_hat
  c(c_hat, sum(e^2) / (length(y) - K))
}
Ajustes <- mapply(Ajustar_curva, split(Datos_curvas$t, Clave_curva), split(Datos_curvas$y, Clave_curva),
                  MoreArgs = list(K = K_estaciones))

#Matriz de coeficientes: K filas y una columna por par estación-año
Coef_est <- Ajustes[1:K_estaciones, ]
Pares_curvas$Sigma2 <- Ajustes[K_estaciones + 1, ]

#Verificación: los coeficientes ajustados por par coinciden con los calculados en bloque (pares completos de 365 días)
C_comp  <- Inv_d %*% t(Phi_d) %*% Y365
Comp365 <- unique(Sel[!Bis_sel, c("STATION_ID", "YEAR")])
Pos <- match(paste(Comp365$STATION_ID, Comp365$YEAR), paste(Pares_curvas$STATION_ID, Pares_curvas$YEAR))
c(Sin_par = sum(is.na(Pos)), Dif_max = max(abs(Coef_est[, Pos] - C_comp)))

#Covarianza residual (pares completos de entrenamiento de 365 días)
E_comp  <- Y365 - Phi_d %*% C_comp      #Residuos: 365 días x 1418 curvas
Sigma_e <- var(t(E_comp))               #Covarianza entre días, estimada a través de las curvas (365 x 365)

#Estructura 1: varianza de cada día, media por mes
Mes_dia <- as.numeric(format(as.Date("1981-01-01") + 0:364, "%m"))
round(tapply(diag(Sigma_e), Mes_dia, mean), 2)
#Estructura 2: correlación media entre residuos separados h días
Cor_e <- cov2cor(Sigma_e)
round(sapply(1:10, function(h) mean(Cor_e[cbind(1:(365 - h), (1 + h):365)])), 3)

#Figura: superficie de covarianza residual y su diagonal
pdf(Ruta("Figuras_proyecto/Fig_Covarianza_residual.pdf"), width = 10, height = 4.8)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
contour(t365, t365, Sigma_e, nlevels = 10, labcex = 0.7, xlab = "s (día del año)", ylab = "t (día del año)",
        main = "Covarianza residual")
plot(t365, diag(Sigma_e), type = "l", xlab = "Día del año", ylab = expression(hat(sigma)^2*(t)),
     main = "Varianza residual a lo largo del año")
dev.off()

#Control: la covarianza residual proyectada, H Sigma_e H', es cero (no sirve para la banda)
H_comp <- Phi_d %*% Inv_d %*% t(Phi_d)
max(abs(H_comp %*% Sigma_e %*% t(H_comp)))

#Propiedades de la matriz suavizadora con K_estaciones: simétrica, idempotente y traza K
c(Asimetria = max(abs(H_comp - t(H_comp))), Idempotencia = max(abs(H_comp %*% H_comp - H_comp)), Traza = sum(diag(H_comp)))

#Varianza del error: con errores independientes Var(e_t) = sigma2 (1 - h_tt), así que sigma2(t) ≈ Var(e_t) / (1 - h_tt)
h_tt      <- diag(H_comp)
Var_error <- diag(Sigma_e) / (1 - h_tt)
#Forma estacional de la varianza (media 1): cada curva la escala con su propio sigma2
Perfil <- Var_error / mean(Var_error)
round(tapply(Perfil, Mes_dia, mean), 2)

#Bandas para el par de ejemplo (2613600105, ACUEDUCTO BUGA, 1985)
i_ej   <- which(Pares_curvas$STATION_ID == "2613600105" & Pares_curvas$YEAR == 1985)
Par_ej <- Datos_curvas[Clave_curva == "2613600105 1985", ]
Phi_ej <- Base_bspline(Par_ej$t, 0, 365, K_estaciones)
B_ej   <- Phi_f %*% solve(t(Phi_ej) %*% Phi_ej) %*% t(Phi_ej)   #Lleva los datos a la curva en la rejilla fina
x_ej   <- as.vector(Phi_f %*% Coef_est[, i_ej])
s2_ej  <- Pares_curvas$Sigma2[i_ej]

#Varianza constante: Var = sigma2 * sum_j B_tj^2
EE_clasico <- sqrt(s2_ej * rowSums(B_ej^2))
#Varianza estacional: Var = sum_j B_tj^2 * sigma2 * Perfil(t_j)
EE_heter   <- sqrt(as.vector(B_ej^2 %*% (s2_ej * Perfil[ceiling(Par_ej$t)])))
Mes_fino <- as.numeric(format(as.Date("1981-01-01") + pmin(floor(t_fino), 364), "%m"))
round(tapply(EE_heter / EE_clasico, Mes_fino, mean), 2)

#Figura (eje vertical cortado en 3 veces el borde superior de la banda; los días por encima se cuentan en el subtítulo)
Tope_ej <- 3 * max(x_ej + 2 * EE_heter)
pdf(Ruta("Figuras_proyecto/Fig_Bandas_ejemplo.pdf"), width = 8, height = 5)
par(mfrow = c(1, 1), mar = c(4.2, 4.2, 2.6, 1))
plot(Par_ej$t, Par_ej$y, pch = 19, cex = 0.5, col = "grey50", xlab = "Día del año", ylim = c(min(x_ej - 2 * EE_heter), Tope_ej),
     ylab = "Lluvia (mm/día)", main = paste0("Estación 2613600105, 1985: curva y bandas punto a punto (K = ", K_estaciones, ")"))
mtext(paste0(sum(Par_ej$y > Tope_ej), " días con más de ", round(Tope_ej), " mm quedan por encima del eje"), side = 3, line = 0.2, cex = 0.75)
abline(h = 0, lty = 3)
lines(t_fino, x_ej, lwd = 2, col = "blue")
lines(t_fino, x_ej + 2 * EE_clasico, lty = 2, col = "blue"); lines(t_fino, x_ej - 2 * EE_clasico, lty = 2, col = "blue")
lines(t_fino, x_ej + 2 * EE_heter, lty = 3, lwd = 2, col = "red"); lines(t_fino, x_ej - 2 * EE_heter, lty = 3, lwd = 2, col = "red")
legend("topright", c("Observado", "Curva estimada", "Banda ±2 EE (varianza constante)", "Banda ±2 EE (varianza estacional)"),
       col = c("grey50", "blue", "blue", "red"), pch = c(19, NA, NA, NA), lty = c(NA, 1, 2, 3), lwd = c(NA, 2, 1, 2), bty = "n", cex = 0.8)
dev.off()

#Proporción de los datos diarios dentro de la banda
mean(Par_ej$y >= approx(t_fino, x_ej - 2 * EE_heter, Par_ej$t)$y & Par_ej$y <= approx(t_fino, x_ej + 2 * EE_heter, Par_ej$t)$y)

#Verificación frente a fda (par de ejemplo): misma base, mismos coeficientes, misma VCG y gl = K
Base_fda <- create.bspline.basis(c(0, 365), nbasis = K_estaciones, norder = 4)
max(abs(eval.basis(t365, Base_fda) - Phi_d))
Ajuste_fda <- smooth.basis(Par_ej$t, Par_ej$y, Base_fda)
max(abs(Ajuste_fda$fd$coefs - Coef_est[, i_ej]))
n_ej <- nrow(Par_ej)
c(VCG_fda = Ajuste_fda$gcv, VCG_a_mano = n_ej * sum((Par_ej$y - Phi_ej %*% Coef_est[, i_ej])^2) / (n_ej - K_estaciones)^2)
Ajuste_fda$df
#Verificación: VC con el atajo de la diagonal de H igual a dejar cada día por fuera y reajustar
H_ej       <- Phi_ej %*% solve(t(Phi_ej) %*% Phi_ej) %*% t(Phi_ej)
VC_atajo   <- mean(((Par_ej$y - H_ej %*% Par_ej$y) / (1 - diag(H_ej)))^2)
Pred_fuera <- numeric(nrow(Par_ej))
for (j in 1:nrow(Par_ej)) {
  c_j <- solve(t(Phi_ej[-j, ]) %*% Phi_ej[-j, ]) %*% t(Phi_ej[-j, ]) %*% Par_ej$y[-j]
  Pred_fuera[j] <- Phi_ej[j, ] %*% c_j
}
c(VC_atajo = VC_atajo, VC_uno_por_fuera = mean((Par_ej$y - Pred_fuera)^2))

#Curvas de estación con tramos negativos (mínimos cuadrados no impone signo; pasa donde la lluvia media es casi 0)
X_dias  <- Phi_d %*% Coef_est                         #Curvas en los 365 centros de día (365 x 3204)
Neg_dia <- X_dias < 0
c(Curvas = sum(colSums(Neg_dia) > 0), Porc_curvas = round(100 * mean(colSums(Neg_dia) > 0), 1),
  Porc_puntos = round(100 * mean(Neg_dia), 2), Mediana = round(median(X_dias[Neg_dia]), 3), Minimo = round(min(X_dias), 2))
#Porcentaje de puntos negativos por mes y en los 15 días de cada borde del año
round(100 * tapply(rowMeans(Neg_dia), Mes_dia, mean), 2)
round(100 * sum(Neg_dia[c(1:15, 351:365), ]) / sum(Neg_dia), 1)
#Lluvia observada en los días en que la curva de su par es negativa
Dia_obs_c <- pmin(ceiling(Datos_curvas$t), 365)
En_neg    <- Neg_dia[cbind(Dia_obs_c, as.integer(Clave_curva))]
round(c(Dias = sum(En_neg), Lluvia_media_mm = mean(Datos_curvas$y[En_neg]), Dias_con_lluvia_pct = 100 * mean(Datos_curvas$y[En_neg] > 0),
        Lluvia_media_resto = mean(Datos_curvas$y[!En_neg]), Pct_de_la_lluvia = 100 * sum(Datos_curvas$y[En_neg]) / sum(Datos_curvas$y)), 3)
#Control: ¿la banda punto a punto (varianza estacional) contiene el 0 en cada punto negativo?
T_curvas <- split(Datos_curvas$t, Clave_curva)        #Días observados de cada par, en el orden de Pares_curvas
Signif   <- unlist(lapply(which(colSums(Neg_dia) > 0), function(i) {
  Ph <- Base_bspline(T_curvas[[i]], 0, 365, K_estaciones)
  d  <- which(Neg_dia[, i])
  B  <- Phi_d[d, , drop = FALSE] %*% solve(t(Ph) %*% Ph) %*% t(Ph)
  X_dias[d, i] + 2 * sqrt(as.vector(B^2 %*% (Pares_curvas$Sigma2[i] * Perfil[ceiling(T_curvas[[i]])]))) < 0
}))
c(Puntos_negativos = length(Signif), Significativos = sum(Signif))
#Lluvia que agrega recortar la curva en 0 (% del total de la curva en el año)
Agrega_recorte <- 100 * (colSums(pmax(X_dias, 0)) / colSums(X_dias) - 1)
round(c(Total = 100 * (sum(pmax(X_dias, 0)) / sum(X_dias) - 1), quantile(Agrega_recorte[colSums(Neg_dia) > 0], c(0.5, 0.95, 1))), 3)
#Decisión: los tramos negativos no se distinguen de 0. La curva se recorta en 0 solo al sacar la raíz (ítem 4) y en
#la lluvia reconstruida; los ítems 2 y 3 usan la curva sin recortar, que conserva el total

#DATOS SATELITALES (péntadas): lectura con terra
Carpeta_sat  <- "Datos_clase/ImagenesSatelitePentadal"
Años_sat     <- 1981:2003
Archivos_sat <- file.path(Carpeta_sat, sprintf("ppt%d%02d.tif", rep(Años_sat, each = 72), rep(1:72, length(Años_sat))))
all(file.exists(Archivos_sat))
Sat <- rast(Archivos_sat)
#Píxeles del Valle: centro dentro del departamento (filtro en memoria; los .tif no se modifican)
Pixeles_valle <- which(!is.na(values(mask(Sat[[1]], Valle))))
Valores_sat <- values(Sat)[Pixeles_valle, , drop = FALSE]    #Filas = píxeles del Valle, columnas = péntadas
c(Pixeles = nrow(Valores_sat), Pentadas = ncol(Valores_sat), Faltantes = sum(is.na(Valores_sat)),
  Minimo = min(Valores_sat), Maximo = max(Valores_sat))

#Trazabilidad de la selección espacial del satélite
Trazabilidad_sat <- data.frame(Paso = c("Malla satelital completa", "Dentro del Valle del Cauca"),
                               Pixeles = c(ncell(Sat), length(Pixeles_valle)),
                               Pares_pixel_año = c(ncell(Sat) * length(Años_sat),
                                                   length(Pixeles_valle) * length(Años_sat)),
                               Porcentaje = round(100 * c(1, length(Pixeles_valle) / ncell(Sat)), 1))
Trazabilidad_sat

#Techo del satélite en 200 mm: valores altos, celdas en 200, píxeles afectados y su rango de longitud
table(Valores_sat[Valores_sat >= 195])
Pixeles_200 <- rowSums(Valores_sat == 200) > 0
c(Valores_200 = sum(Valores_sat == 200), Pixeles_200 = sum(Pixeles_200))
range(xyFromCell(Sat, Pixeles_valle[which(Pixeles_200)])[, "x"])

#Estructura de las péntadas: 6 por mes; las 5 primeras de 5 días y la sexta con el resto del mes (3 a 6 días)
Pent <- data.frame(Año = rep(Años_sat, each = 72), Pentada = rep(1:72, length(Años_sat)))
Pent$Mes     <- ceiling(Pent$Pentada / 6)
Pent$j       <- Pent$Pentada - 6 * (Pent$Mes - 1)        #Posición dentro del mes (1 a 6)
Fin_mes      <- as.numeric(format(as.Date(paste(ifelse(Pent$Mes == 12, Pent$Año + 1, Pent$Año),
                                                ifelse(Pent$Mes == 12, 1, Pent$Mes + 1), 1, sep = "-")) - 1, "%d"))
Pent$Dia_ini <- 5 * (Pent$j - 1) + 1
Pent$Dia_fin <- ifelse(Pent$j < 6, 5 * Pent$j, Fin_mes)
Pent$Dias    <- Pent$Dia_fin - Pent$Dia_ini + 1

#Tiempo: centro del intervalo que cubre la péntada, en día del año reescalado a [0, 365]
Doy_ini   <- as.numeric(format(as.Date(paste(Pent$Año, Pent$Mes, Pent$Dia_ini, sep = "-")), "%j"))
Doy_fin   <- as.numeric(format(as.Date(paste(Pent$Año, Pent$Mes, Pent$Dia_fin, sep = "-")), "%j"))
N_año_sat <- ifelse((Pent$Año %% 4 == 0 & Pent$Año %% 100 != 0) | Pent$Año %% 400 == 0, 366, 365)
Pent$t    <- ((Doy_ini - 1 + Doy_fin) / 2) * 365 / N_año_sat
#Verificación: las 72 péntadas de cada año cubren todos sus días
all(tapply(Pent$Dias, Pent$Año, sum) == tapply(N_año_sat, Pent$Año, unique))

#Variable de trabajo: tasa diaria media de cada péntada en mm/día (la raíz se aplica después a la curva)
Y_sat <- t(Valores_sat) / Pent$Dias    #Filas = péntadas, columnas = píxeles

#Pendiente b de la tasa por péntada (grupos píxel x mes, entrenamiento)
Tasa_ent <- t(Valores_sat[, Pent$Año <= 1995]) / Pent$Dias[Pent$Año <= 1995]     #Péntadas x píxeles
Mes_ent  <- Pent$Mes[Pent$Año <= 1995]
Media_s  <- apply(Tasa_ent, 2, tapply, Mes_ent, mean)                           #12 meses x píxeles del Valle
Var_s    <- apply(Tasa_ent, 2, tapply, Mes_ent, var)
Usar_s   <- Media_s > 0 & Var_s > 0
round(c(Sin_transformar = Pendiente_b(Var_s[Usar_s], Media_s[Usar_s]),
        Raiz   = Pendiente_b(apply(sqrt(Tasa_ent), 2, tapply, Mes_ent, var)[Usar_s], Media_s[Usar_s]),
        Log_mm = Pendiente_b(apply(log1p(Tasa_ent), 2, tapply, Mes_ent, var)[Usar_s], Media_s[Usar_s])), 3)

#Selección de K para el satélite: entrenamiento (1981-1995), píxeles del Valle, tasa en mm/día. Una curva por píxel y
#año (72 péntadas), con una matriz para años normales y otra para bisiestos
Armar_sat <- function(Años) do.call(cbind, lapply(Años, function(a) Y_sat[Pent$Año == a, ]))
Años_ent  <- 1981:1995
Bis_ent   <- (Años_ent %% 4 == 0 & Años_ent %% 100 != 0) | Años_ent %% 400 == 0
Y72  <- Armar_sat(Años_ent[!Bis_ent])
Y72b <- Armar_sat(Años_ent[Bis_ent])
t72  <- Pent$t[Pent$Año == 1981]
t72b <- Pent$t[Pent$Año == 1984]

Criterios_sat <- lapply(Ks, function(K) rbind(Criterios_K(t72, Y72, K), Criterios_K(t72b, Y72b, K)))
Tabla_K_sat   <- data.frame(K = Ks, t(sapply(Criterios_sat, colMeans)))
round(Tabla_K_sat, 5)
K_VC_sat  <- Ks[which.min(Tabla_K_sat$VC)]
K_VCG_sat <- Ks[which.min(Tabla_K_sat$VCG)]
c(K_VC_sat = K_VC_sat, K_VCG_sat = K_VCG_sat)

#Figura: K vs ECM, VC y VCG para el satélite (VC se grafica hasta K = 25 porque luego explota)
pdf(Ruta("Figuras_proyecto/Fig_Seleccion_K_satelite.pdf"), width = 10, height = 3.8)
par(mfrow = c(1, 3))
plot(Ks, Tabla_K_sat$ECM, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "ECM", main = "ECM de ajuste (satélite)")
plot(Ks[Ks <= 25], Tabla_K_sat$VC[Ks <= 25], type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "VC", main = "Validación cruzada (satélite)")
abline(v = K_VC_sat, lty = 3, col = "red")
plot(Ks, Tabla_K_sat$VCG, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "VCG", main = "VCG (satélite)")
abline(v = K_VCG_sat, lty = 3, col = "red")
dev.off()

#Comparación pareada frente al K de mínima VC (mismo criterio que en estaciones)
VC_sat  <- sapply(Criterios_sat, function(x) x[, "VC"])
Dif_sat <- VC_sat - VC_sat[, which(Ks == K_VC_sat)]
z_sat   <- colMeans(Dif_sat) / (apply(Dif_sat, 2, sd) / sqrt(nrow(Dif_sat)))
round(z_sat[Ks <= 14], 2)
K_equivalentes_sat <- Ks[abs(z_sat) < 2 | is.nan(z_sat)]
K_equivalentes_sat

#K elegido para el satélite: el que minimiza la VC, igual que en las estaciones (y también la VCG)
K_satelite <- K_VC_sat
c(K_satelite = K_satelite, Coincide_con_VCG = K_VC_sat == K_VCG_sat)

#Control de coherencia (no es criterio de selección): proporción de curvas con exactamente dos picos
Prop_dos_picos <- function(t, Y, K) {
  Phi <- Base_bspline(t, 0, 365, K)
  X   <- Base_bspline(t_fino, 0, 365, K) %*% solve(t(Phi) %*% Phi) %*% t(Phi) %*% Y
  mean(apply(X, 2, function(x) { d <- diff(x); sum(d[-1] < 0 & d[-length(d)] >= 0) }) == 2)
}
round(sapply(K_equivalentes_sat, function(K) Prop_dos_picos(t72, Y72, K)), 4)

#Exportación para Python: Y_sat (filas = péntadas, columnas = píxeles del Valle, tasa en mm/día), péntadas, píxeles y base
Y_sat_export <- as.data.frame(Y_sat)
names(Y_sat_export) <- paste0("pixel_", Pixeles_valle)
write.csv(Y_sat_export, Ruta("Datos_preparados/satelite_y.csv"), row.names = FALSE)
write.csv(data.frame(year = Pent$Año,
                     pentada = Pent$Pentada,
                     mes = Pent$Mes,
                     pos_mes = Pent$j,
                     dia_ini = Pent$Dia_ini,
                     dia_fin = Pent$Dia_fin,
                     dias = Pent$Dias,
                     t = Pent$t),
          Ruta("Datos_preparados/satelite_pentadas.csv"), row.names = FALSE)
Coord_export <- xyFromCell(Sat, Pixeles_valle)
write.csv(data.frame(pixel = Pixeles_valle,
                     lon = Coord_export[, "x"],
                     lat = Coord_export[, "y"]),
          Ruta("Datos_preparados/satelite_pixeles.csv"), row.names = FALSE)
write.csv(data.frame(nudo = Nodos_bspline(0, 365, K_satelite, 4)),
          Ruta("Datos_preparados/satelite_nudos.csv"), row.names = FALSE)
write.csv(data.frame(fuente = "satelite",
                     K = K_satelite,
                     orden = 4,
                     dominio_inf = 0,
                     dominio_sup = 365,
                     transformacion = "ninguna (la raiz se aplica a la curva: satelite_raiz_coef.csv)",
                     unidad_respuesta = "mm/dia",
                     periodo = "1981-2003",
                     n_pixeles = length(Pixeles_valle),
                     n_pentadas = nrow(Pent),
                     n_curvas = length(Pixeles_valle) * length(Años_sat),
                     archivo_datos = "satelite_y.csv",
                     archivo_tiempo = "satelite_pentadas.csv",
                     archivo_pixeles = "satelite_pixeles.csv",
                     archivo_nudos = "satelite_nudos.csv"),
          Ruta("Datos_preparados/satelite_config.csv"), row.names = FALSE)
cat("Satélite exportado:", nrow(Y_sat), "péntadas x", ncol(Y_sat), "píxeles\n")

#Curvas del satélite con K_satelite: una por píxel y año (cada año comparte Phi)
Curvas_sat <- data.frame(Pixel = rep(Pixeles_valle, length(Años_sat)), Año = rep(Años_sat, each = length(Pixeles_valle)))
Coord_pixel <- xyFromCell(Sat, 1:ncell(Sat))
Curvas_sat$Lon <- Coord_pixel[Curvas_sat$Pixel, "x"]
Curvas_sat$Lat <- Coord_pixel[Curvas_sat$Pixel, "y"]
#Columnas creadas antes del ciclo (asignar por partes a una columna inexistente recicla valores)
Curvas_sat$Sigma2 <- NA
Curvas_sat$Techo  <- NA
Coef_sat <- matrix(NA, K_satelite, nrow(Curvas_sat))
for (a in Años_sat) {
  f     <- Pent$Año == a
  Phi_a <- Base_bspline(Pent$t[f], 0, 365, K_satelite)
  A_a   <- solve(t(Phi_a) %*% Phi_a) %*% t(Phi_a)
  C_a   <- A_a %*% Y_sat[f, ]
  Col   <- which(Curvas_sat$Año == a)
  Coef_sat[, Col] <- C_a
  Curvas_sat$Sigma2[Col] <- colSums((Y_sat[f, ] - Phi_a %*% C_a)^2) / (72 - K_satelite)
  #Marca de techo: al menos una péntada en 200 mm ese año
  Curvas_sat$Techo[Col]  <- rowSums(Valores_sat[, f] == 200) > 0
}
c(Curvas = ncol(Coef_sat), Faltantes = sum(is.na(Coef_sat)), Con_techo = sum(Curvas_sat$Techo))

#Propiedades de la matriz suavizadora con el K seleccionado (año de 365 días)
Phi_72 <- Base_bspline(t72, 0, 365, K_satelite)
H_sat  <- Phi_72 %*% solve(t(Phi_72) %*% Phi_72) %*% t(Phi_72)
c(Asimetria = max(abs(H_sat - t(H_sat))), Idempotencia = max(abs(H_sat %*% H_sat - H_sat)), Traza = sum(diag(H_sat)))

#Covarianza residual del satélite (años de entrenamiento de 365 días, píxeles del Valle del Cauca)
E_sat     <- Y72 - H_sat %*% Y72
Sigma_sat <- var(t(E_sat))
Mes_pent  <- ceiling((1:72) / 6)
round(tapply(diag(Sigma_sat), Mes_pent, mean), 3)
Cor_e_sat <- cov2cor(Sigma_sat)
round(sapply(1:6, function(h) mean(Cor_e_sat[cbind(1:(72 - h), (1 + h):72)])), 3)
#Correlación inducida solo por la proyección, Var(e) = sigma2 (I - H), en el satélite y en las estaciones
round(sapply(1:6, function(h) mean(cov2cor(diag(72) - H_sat)[cbind(1:(72 - h), (1 + h):72)])), 3)
round(sapply(1:5, function(h) mean(cov2cor(diag(365) - H_comp)[cbind(1:(365 - h), (1 + h):365)])), 4)

#Control: la covarianza residual proyectada da cero, como en las estaciones; perfil estacional de la varianza del error
max(abs(H_sat %*% Sigma_sat %*% t(H_sat)))
Var_error_sat <- diag(Sigma_sat) / (1 - diag(H_sat))
Perfil_sat    <- Var_error_sat / mean(Var_error_sat)
round(tapply(Perfil_sat, Mes_pent, mean), 2)
round(range(diag(H_sat)), 3)

#Figura
pdf(Ruta("Figuras_proyecto/Fig_Covarianza_residual_satelite.pdf"), width = 10, height = 4.8)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
contour(t72, t72, Sigma_sat, nlevels = 10, labcex = 0.7, xlab = "s (día del año)", ylab = "t (día del año)",
        main = "Covarianza residual (satélite)")
plot(t72, diag(Sigma_sat), type = "b", pch = 19, cex = 0.6, xlab = "Día del año", ylab = expression(hat(sigma)^2*(t)),
     main = "Varianza residual a lo largo del año (satélite)")
dev.off()

#Bandas para el píxel que contiene a Cali (76,53° O; 3,45° N), año 1985
Pixel_ej  <- cellFromXY(Sat, cbind(-76.53, 3.45))
i_sat     <- which(Curvas_sat$Pixel == Pixel_ej & Curvas_sat$Año == 1985)
i_px      <- match(Pixel_ej, Pixeles_valle)    #Posición del píxel de Cali dentro de los píxeles del Valle
y_px      <- Y_sat[Pent$Año == 1985, i_px]
Phi_f_sat <- Base_bspline(t_fino, 0, 365, K_satelite)
B_sat     <- Phi_f_sat %*% solve(t(Phi_72) %*% Phi_72) %*% t(Phi_72)
x_px      <- as.vector(Phi_f_sat %*% Coef_sat[, i_sat])
s2_px     <- Curvas_sat$Sigma2[i_sat]
EE_clasico_sat <- sqrt(s2_px * rowSums(B_sat^2))
EE_heter_sat   <- sqrt(as.vector(B_sat^2 %*% (s2_px * Perfil_sat)))
round(tapply(EE_heter_sat / EE_clasico_sat, Mes_fino, mean), 2)

pdf(Ruta("Figuras_proyecto/Fig_Bandas_ejemplo_satelite.pdf"), width = 8, height = 5)
par(mfrow = c(1, 1), mar = c(4.2, 4.2, 2, 1))
plot(t72, y_px, pch = 19, cex = 0.7, col = "grey50", xlab = "Día del año", ylab = "Lluvia (mm/día)",
     main = paste0("Satélite, píxel de Cali, 1985: curva y bandas punto a punto (K = ", K_satelite, ")"))
abline(h = 0, lty = 3)
lines(t_fino, x_px, lwd = 2, col = "blue")
lines(t_fino, x_px + 2 * EE_clasico_sat, lty = 2, col = "blue"); lines(t_fino, x_px - 2 * EE_clasico_sat, lty = 2, col = "blue")
lines(t_fino, x_px + 2 * EE_heter_sat, lty = 3, lwd = 2, col = "red"); lines(t_fino, x_px - 2 * EE_heter_sat, lty = 3, lwd = 2, col = "red")
legend("topright", c("Observado", "Curva estimada", "Banda ±2 EE (varianza constante)", "Banda ±2 EE (varianza estacional)"),
       col = c("grey50", "blue", "blue", "red"), pch = c(19, NA, NA, NA), lty = c(NA, 1, 2, 3), lwd = c(NA, 2, 1, 2), bty = "n", cex = 0.8)
dev.off()
mean(y_px >= approx(t_fino, x_px - 2 * EE_heter_sat, t72)$y & y_px <= approx(t_fino, x_px + 2 * EE_heter_sat, t72)$y)

#Verificación contra la librería fda, píxel de ejemplo: misma base, mismos coeficientes y misma VCG
Base_fda_sat <- create.bspline.basis(c(0, 365), nbasis = K_satelite, norder = 4)
max(abs(eval.basis(t72, Base_fda_sat) - Phi_72))
Ajuste_fda_sat <- smooth.basis(t72, y_px, Base_fda_sat)
max(abs(Ajuste_fda_sat$fd$coefs - Coef_sat[, i_sat]))
c(VCG_fda = Ajuste_fda_sat$gcv, VCG_a_mano = 72 * sum((y_px - Phi_72 %*% Coef_sat[, i_sat])^2) / (72 - K_satelite)^2)
#CURVAS EN mm Y TRANSFORMACIÓN PARA MODELAR (ítems 1 y 4)
#1. Control: cada curva conserva la lluvia de su par (las B-splines suman 1 y los residuos suman 0)
N_obs   <- as.vector(table(Clave_curva))                       #Días observados de cada par (orden de Pares_curvas)
Phi_obs <- Base_bspline(Datos_curvas$t, 0, 365, K_estaciones)  #Base en los días observados de todos los pares
Idx_par <- as.integer(Clave_curva)                             #Par al que pertenece cada fila de Datos_curvas
X_obs   <- rowSums(Phi_obs * t(Coef_est)[Idx_par, ])           #Curva de cada par en sus días observados
Lluvia_obs <- as.vector(tapply(Datos_curvas$y, Clave_curva, sum))
c(Identidad_estaciones = max(abs(as.vector(tapply(X_obs, Clave_curva, sum)) - Lluvia_obs)), Pares_sin_lluvia = sum(Lluvia_obs == 0))
#Mismo control en el satélite (72 péntadas por año)
max(sapply(Años_sat, function(a) {
  f <- Pent$Año == a
  max(abs(colSums(Base_bspline(Pent$t[f], 0, 365, K_satelite) %*% Coef_sat[, Curvas_sat$Año == a]) - colSums(Y_sat[f, ])))
}))

#Control: con el mismo K, la curva de raíz(P) al cuadrado recupera solo una parte (%) de la lluvia del par, que
#depende de la frecuencia de días con lluvia
Ajustes_raiz  <- mapply(Ajustar_curva, split(Datos_curvas$t, Clave_curva), split(sqrt(Datos_curvas$y), Clave_curva),
                        MoreArgs = list(K = K_estaciones))
Raiz_obs      <- rowSums(Phi_obs * t(Ajustes_raiz[1:K_estaciones, ])[Idx_par, ])
Recupera_raiz <- 100 * as.vector(tapply(pmax(Raiz_obs, 0)^2, Clave_curva, sum)) / Lluvia_obs
Dias_lluvia_par <- 100 * as.vector(tapply(Datos_curvas$y > 0, Clave_curva, mean))
round(quantile(Recupera_raiz, c(0.05, 0.5, 0.95)), 1)
round(tapply(Recupera_raiz, cut(Dias_lluvia_par, c(0, 20, 30, 40, 50, 60, 80, 100)), median), 1)
#Ejemplo: dos pares de entrenamiento con casi la misma curva en raíz y más del doble de lluvia en uno
Ej_raiz <- match(c("2612810102 1987", "2612810101 1994"), levels(Clave_curva))
data.frame(Par = levels(Clave_curva)[Ej_raiz], Media_curva_raiz = round(as.vector(tapply(Raiz_obs, Clave_curva, mean))[Ej_raiz], 3),
           Lluvia_mm_dia = round(Lluvia_obs[Ej_raiz] / N_obs[Ej_raiz], 2), Dias_lluvia_pct = round(Dias_lluvia_par[Ej_raiz]))

#2. Transformación para modelar (ítem 4): z(t) = raíz(max(x(t), 0)), con regreso exacto x = z^2. Pendiente de la
#   varianza entre años de las curvas de cada estación frente a su nivel (estaciones con al menos 5 curvas)
X_fino_est <- Phi_f %*% Coef_est                                #Curvas de estación en la rejilla fina (mm/día)
X_fino_sat <- Phi_f_sat %*% Coef_sat                            #Curvas del satélite en la rejilla fina (mm/día)
Est_5 <- names(which(table(Pares_curvas$STATION_ID) >= 5))
Pendiente_curvas <- function(Z) {
  Nivel <- sapply(Est_5, function(s) mean(X_fino_est[, Pares_curvas$STATION_ID == s]))
  Var_a <- sapply(Est_5, function(s) mean(apply(Z[, Pares_curvas$STATION_ID == s, drop = FALSE], 1, var)))
  Pendiente_b(Var_a, Nivel)
}
round(c(Estaciones = length(Est_5), En_mm = Pendiente_curvas(X_fino_est), Raiz = Pendiente_curvas(sqrt(pmax(X_fino_est, 0))),
        Log_mm = Pendiente_curvas(log1p(pmax(X_fino_est, 0))), Log_cm = Pendiente_curvas(log1p(pmax(X_fino_est, 0) / 10))), 3)
#Escala elegida: raíz (no depende de las unidades, regreso exacto a mm y misma escala de la regla de cobertura)

#3. La raíz de una curva no es una B-spline cúbica: se representa por mínimos cuadrados en K_z bases sobre t_fino.
#   K_z: menor número de bases con mediana del error L2 relativo del regreso menor que 1 % en entrenamiento
#Funciones de representación de la raíz y de su error relativo de regreso (%)
Representar_raiz <- function(X, Kz) {
  P_z <- Base_bspline(t_fino, 0, 365, Kz)
  solve(t(P_z) %*% P_z, t(P_z) %*% sqrt(pmax(X, 0)))
}
Error_raiz <- function(X, Kz) {
  Regreso <- pmax(Base_bspline(t_fino, 0, 365, Kz) %*% Representar_raiz(X, Kz), 0)^2
  e <- sqrt(colSums((Regreso - pmax(X, 0))^2) / colSums(pmax(X, 0)^2))
  100 * c(Mediana = median(e), P95 = quantile(e, 0.95, names = FALSE), Maximo = max(e))
}
Ks_z    <- 8:25
Ent_par <- Pares_curvas$YEAR <= 1995 & !(Pares_curvas$STATION_ID %in% Reservadas)
Tabla_z_est <- data.frame(K_z = Ks_z, t(sapply(Ks_z, function(Kz) Error_raiz(X_fino_est[, Ent_par], Kz))))
Tabla_z_sat <- data.frame(K_z = Ks_z, t(sapply(Ks_z, function(Kz) Error_raiz(X_fino_sat[, Curvas_sat$Año <= 1995], Kz))))
round(Tabla_z_est, 3); round(Tabla_z_sat, 3)
K_z_est <- min(Ks_z[Tabla_z_est$Mediana < 1])
K_z_sat <- min(Ks_z[Tabla_z_sat$Mediana < 1])
c(K_z_estaciones = K_z_est, K_z_satelite = K_z_sat)

#Coeficientes de la raíz de todas las curvas (las de validación y prueba se representan con la misma regla)
Coef_z_est <- Representar_raiz(X_fino_est, K_z_est)
Coef_z_sat <- Representar_raiz(X_fino_sat, K_z_sat)
#Porcentaje de valores recortados: curva en mm negativa (antes de la raíz) y raíz reproyectada negativa (en el regreso)
round(c(Est_curva_mm = 100 * mean(X_fino_est < 0), Est_raiz = 100 * mean(Base_bspline(t_fino, 0, 365, K_z_est) %*% Coef_z_est < 0),
        Sat_curva_mm = 100 * mean(X_fino_sat < 0), Sat_raiz = 100 * mean(Base_bspline(t_fino, 0, 365, K_z_sat) %*% Coef_z_sat < 0)), 3)

#Exportación para Python: coeficientes de la raíz de las curvas (la red modela z; la lluvia es max(z, 0)^2)
Tabla_coef <- function(Ids, C) { D <- data.frame(Ids, t(C)); names(D)[-(1:ncol(Ids))] <- paste0("c", 1:nrow(C)); D }
write.csv(Tabla_coef(data.frame(station_id = Pares_curvas$STATION_ID, year = Pares_curvas$YEAR), Coef_z_est),
          Ruta("Datos_preparados/estaciones_raiz_coef.csv"), row.names = FALSE)
write.csv(Tabla_coef(data.frame(pixel = Curvas_sat$Pixel, year = Curvas_sat$Año), Coef_z_sat),
          Ruta("Datos_preparados/satelite_raiz_coef.csv"), row.names = FALSE)
write.csv(data.frame(fuente = c("estaciones", "satelite"), K_z = c(K_z_est, K_z_sat), orden = 4, dominio_inf = 0, dominio_sup = 365,
                     curva = "z(t) = raiz(max(x(t), 0)), x = curva en mm/dia", regreso = "lluvia = max(z, 0)^2",
                     archivo_coef = c("estaciones_raiz_coef.csv", "satelite_raiz_coef.csv")),
          Ruta("Datos_preparados/raiz_config.csv"), row.names = FALSE)
#UBICACIÓN DE LAS ESTACIONES Y PÍXEL SATELITAL

#Ubicación de cada estación con curvas y el píxel del satélite que la contiene (NA si queda fuera de la malla)
Ubicacion <- Coord_est[, c("CODIGO", "STATION_NAME", "CATEGORY", "Lon", "lat")]
names(Ubicacion) <- c("STATION_ID", "Nombre", "Categoria", "Lon", "Lat")
Estaciones_curvas <- merge(data.frame(STATION_ID = unique(Pares_curvas$STATION_ID)), Ubicacion, by = "STATION_ID", all.x = TRUE)
Con_xy <- !is.na(Estaciones_curvas$Lon)
Estaciones_curvas$Pixel <- NA
Estaciones_curvas$Pixel[Con_xy] <- cellFromXY(Sat, cbind(Estaciones_curvas$Lon[Con_xy], Estaciones_curvas$Lat[Con_xy]))
#Control: las estaciones con coordenadas tienen píxel del Valle; las ubicadas por municipio no tienen píxel
c(Con_coordenadas = sum(Con_xy), Con_pixel_del_Valle = sum(Estaciones_curvas$Pixel %in% Pixeles_valle),
  Ubicadas_por_municipio = sum(Estaciones_curvas$STATION_ID %in% Sin_coord_valle))

#Regla de asignación (ítem 3): cada estación se empareja con el píxel que la contiene (cellFromXY); las que comparten
#píxel reciben la misma curva satelital
Est_pixel <- table(Estaciones_curvas$Pixel[Estaciones_curvas$Pixel %in% Pixeles_valle])
c(Estaciones_en_malla = sum(Est_pixel), Pixeles = length(Est_pixel), Pixeles_compartidos = sum(Est_pixel > 1),
  Estaciones_en_pixel_compartido = sum(Est_pixel[Est_pixel > 1]), Maximo_por_pixel = max(Est_pixel))
Estaciones_curvas$Comparte_pixel <- Estaciones_curvas$Pixel %in% as.numeric(names(Est_pixel)[Est_pixel > 1])
#Lista de las estaciones que comparten píxel
Compartidas <- Estaciones_curvas[Estaciones_curvas$Comparte_pixel, c("Pixel", "STATION_ID", "Nombre", "Lon", "Lat")]
Compartidas[order(Compartidas$Pixel), ]

#ALTITUD DE LAS ESTACIONES (catálogo de elevación del curso, hoja 1, en metros)
c(Filas = nrow(Elevacion), Codigos = length(unique(Elevacion$codigo)), Sin_altitud = sum(is.na(Elevacion$alt)))
Estaciones_curvas$Altitud_catalogo <- Elevacion$alt[match(Estaciones_curvas$STATION_ID, Elevacion$codigo)]
Estaciones_curvas$Altitud <- Estaciones_curvas$Altitud_catalogo
Estaciones_curvas$Altitud_fuente <- ifelse(!is.na(Estaciones_curvas$Altitud), "Catalogo del curso", "Sin fuente")

#RELLENO DE ALTITUD CON EL DEM COPERNICUS GLO-30 (1 arcsec, ~30 m; tiles N03-N05 x W078-W076 recortados al Valle)
#Altitud de cada estación: mediana de la ventana de 3x3 celdas (el valor bilineal queda como diagnóstico)
Dem_valle <- rast("Datos_externos/dem_valle_copernicus30m.tif")
#Solo las estaciones con coordenadas: las ubicadas por municipio no tienen un punto donde leer el DEM
Fil_xy  <- which(Con_xy)
Pts_alt <- vect(Estaciones_curvas[Fil_xy, c("Lon", "Lat")], geom = c("Lon", "Lat"), crs = "EPSG:4326")
Estaciones_curvas$Altitud_DEM_bilineal <- NA_real_
Estaciones_curvas$Altitud_DEM_bilineal[Fil_xy] <- extract(Dem_valle, Pts_alt, method = "bilinear")[, 2]
e_alt <- extract(Dem_valle, Pts_alt, cells = TRUE)
rc_alt <- rowColFromCell(Dem_valle, e_alt$cell)
nrv_alt <- nrow(Dem_valle); ncv_alt <- ncol(Dem_valle)
Estaciones_curvas$Altitud_DEM_mediana3 <- NA_real_
Estaciones_curvas$Pendiente_DEM <- NA_real_
for (k in seq_along(Fil_xy)) {
  i <- Fil_xy[k]
  r <- rc_alt[k, 1]; cc <- rc_alt[k, 2]
  rr <- max(1, r - 1):min(nrv_alt, r + 1)
  co <- max(1, cc - 1):min(ncv_alt, cc + 1)
  vals <- as.numeric(unlist(Dem_valle[cellFromRowColCombine(Dem_valle, rr, co)]))
  vals <- vals[!is.na(vals)]
  Estaciones_curvas$Altitud_DEM_mediana3[i] <- median(vals)
  if (length(vals) == 9) {
    z <- matrix(vals, nrow = 3, byrow = TRUE)
    dx <- res(Dem_valle)[1] * 111320 * cos(Estaciones_curvas$Lat[i] * pi / 180)
    dy <- res(Dem_valle)[2] * 110574
    dzdx <- ((z[1,3] + 2*z[2,3] + z[3,3]) - (z[1,1] + 2*z[2,1] + z[3,1])) / (8 * dx)
    dzdy <- ((z[3,1] + 2*z[3,2] + z[3,3]) - (z[1,1] + 2*z[1,2] + z[1,3])) / (8 * dy)
    Estaciones_curvas$Pendiente_DEM[i] <- atan(sqrt(dzdx^2 + dzdy^2)) * 180 / pi
  }
}
#El método operativo es la mediana 3x3 de la ventana alrededor de la estación
Estaciones_curvas$Altitud_DEM <- Estaciones_curvas$Altitud_DEM_mediana3
c(Cobertura_DEM = sum(!is.na(Estaciones_curvas$Altitud_DEM)),
  Sin_DEM = sum(is.na(Estaciones_curvas$Altitud_DEM)))

#Validación del DEM: estaciones con altitud de catálogo y con coordenadas (error = DEM - catálogo, en metros)
Val_dem <- !is.na(Estaciones_curvas$Altitud_catalogo) & !is.na(Estaciones_curvas$Altitud_DEM)
Err_dem <- Estaciones_curvas$Altitud_DEM[Val_dem] - Estaciones_curvas$Altitud_catalogo[Val_dem]
Err_bil <- Estaciones_curvas$Altitud_DEM_bilineal[Val_dem] - Estaciones_curvas$Altitud_catalogo[Val_dem]
Resumen_error <- function(e) c(N = length(e), Sesgo = mean(e), MAE = mean(abs(e)), RMSE = sqrt(mean(e^2)))
round(rbind(Mediana_3x3 = Resumen_error(Err_dem), Bilineal = Resumen_error(Err_bil)), 1)
#Error por estrato: banda de altitud según el DEM (lo único que se conoce en una estación sin catálogo) y pendiente
Banda_dem <- cut(Estaciones_curvas$Altitud_DEM[Val_dem], c(-Inf, 500, 2000, Inf), right = FALSE,
                 labels = c("<500", "500-2000", ">=2000"))
Pend_dem  <- ifelse(Estaciones_curvas$Pendiente_DEM[Val_dem] < 15, "<15", ">=15")
Tabla_dem <- do.call(rbind, lapply(split(Err_dem, list(Banda_dem, Pend_dem), drop = TRUE), Resumen_error))
round(Tabla_dem, 1)
#Umbrales del relleno según Tabla_dem: sobre 2000 m el DEM se aleja del catálogo en cientos de metros, bajo 500 m solo
#hay 3 estaciones de catálogo y entre 500 y 2000 m la pendiente >= 15 grados casi duplica el error. Los umbrales se
#eligieron con estas mismas estaciones, así que el error del estrato elegido es optimista
Error_relleno      <- Resumen_error(Err_dem[Banda_dem == "500-2000" & Pend_dem == "<15"])
Error_sensibilidad <- Resumen_error(Err_dem[Banda_dem == "500-2000" & Pend_dem == ">=15"])
round(rbind(Relleno = Error_relleno, Sensibilidad_empinadas = Error_sensibilidad), 1)

#Criterio principal: 500-2000 m y pendiente < 15 grados
rellenar_dem <- is.na(Estaciones_curvas$Altitud_catalogo) &
                !is.na(Estaciones_curvas$Altitud_DEM) &
                Estaciones_curvas$Altitud_DEM >= 500 & Estaciones_curvas$Altitud_DEM < 2000 &
                !is.na(Estaciones_curvas$Pendiente_DEM) & Estaciones_curvas$Pendiente_DEM < 15
Estaciones_curvas$Altitud_error_MAE <- NA_real_
Estaciones_curvas$Altitud_error_RMSE <- NA_real_
Estaciones_curvas$Altitud[rellenar_dem] <- Estaciones_curvas$Altitud_DEM[rellenar_dem]
Estaciones_curvas$Altitud_fuente[rellenar_dem] <- "DEM Copernicus GLO-30"
#Error esperado de la altitud rellenada: el del estrato elegido en la validación
Estaciones_curvas$Altitud_error_MAE[rellenar_dem] <- Error_relleno[["MAE"]]
Estaciones_curvas$Altitud_error_RMSE[rellenar_dem] <- Error_relleno[["RMSE"]]

#Sensibilidad: 500-2000 m sin filtro de pendiente (error de validación: Error_sensibilidad)
rellenar_dem_sens <- is.na(Estaciones_curvas$Altitud_catalogo) &
                     !is.na(Estaciones_curvas$Altitud_DEM) &
                     Estaciones_curvas$Altitud_DEM >= 500 & Estaciones_curvas$Altitud_DEM < 2000 &
                     !is.na(Estaciones_curvas$Pendiente_DEM) & Estaciones_curvas$Pendiente_DEM >= 15
Estaciones_curvas$Altitud_sensibilidad <- Estaciones_curvas$Altitud
Estaciones_curvas$Altitud_fuente_sensibilidad <- Estaciones_curvas$Altitud_fuente
Estaciones_curvas$Altitud_sensibilidad[rellenar_dem_sens] <- Estaciones_curvas$Altitud_DEM[rellenar_dem_sens]
Estaciones_curvas$Altitud_fuente_sensibilidad[rellenar_dem_sens] <- "DEM Copernicus GLO-30 (sensibilidad)"
c(Relleno_principal = sum(rellenar_dem), Relleno_sensibilidad = sum(rellenar_dem_sens))

#Estaciones con curva: con coordenadas, con altitud de catálogo, con relleno DEM y sin altitud
c(Con_coord = sum(Con_xy),
  Altitud_catalogo = sum(!is.na(Estaciones_curvas$Altitud_catalogo)),
  Altitud_DEM = sum(Estaciones_curvas$Altitud_fuente == "DEM Copernicus GLO-30"),
  Altitud_total = sum(!is.na(Estaciones_curvas$Altitud)),
  Sin_altitud = sum(is.na(Estaciones_curvas$Altitud)),
  Altitud_sensibilidad_total = sum(!is.na(Estaciones_curvas$Altitud_sensibilidad)))

#Trazabilidad para los ítems 3 y 4: coordenadas dentro de la malla del Valle y, para la red, altitud
Dias_par_curva <- as.vector(table(Clave_curva))
Con_coord <- Pares_curvas$STATION_ID %in% Estaciones_curvas$STATION_ID[Con_xy]
En_malla  <- Pares_curvas$STATION_ID %in% Estaciones_curvas$STATION_ID[Estaciones_curvas$Pixel %in% Pixeles_valle]
Con_alt_cat  <- En_malla & Pares_curvas$STATION_ID %in% Estaciones_curvas$STATION_ID[!is.na(Estaciones_curvas$Altitud_catalogo)]
Con_alt      <- En_malla & Pares_curvas$STATION_ID %in% Estaciones_curvas$STATION_ID[!is.na(Estaciones_curvas$Altitud)]
Con_alt_sens <- En_malla & Pares_curvas$STATION_ID %in% Estaciones_curvas$STATION_ID[!is.na(Estaciones_curvas$Altitud_sensibilidad)]
Con_sin_alt      <- En_malla & !Con_alt
Con_sin_alt_sens <- En_malla & !Con_alt_sens
Fila_traza <- function(Paso, Sel) data.frame(Paso = Paso, Estaciones = length(unique(Pares_curvas$STATION_ID[Sel])),
                                             Pares = sum(Sel), Dias = sum(Dias_par_curva[Sel]))
Trazabilidad <- rbind(Trazabilidad,
                      Fila_traza("Con coordenadas (archivo del curso)", Con_coord),
                      Fila_traza("Dentro de la malla satelital del Valle", En_malla),
                      Fila_traza("Con altitud de catalogo", Con_alt_cat),
                      Fila_traza("Con altitud de catalogo + relleno DEM 500-2000 m y pendiente <15", Con_alt),
                      Fila_traza("Descartadas por falta de altitud (principal)", Con_sin_alt),
                      Fila_traza("Sensibilidad: + DEM 500-2000 m sin filtro de pendiente", Con_alt_sens),
                      Fila_traza("Descartadas por falta de altitud (sensibilidad)", Con_sin_alt_sens))
Trazabilidad
table(cut(Pares_curvas$YEAR[Con_alt], c(1980, 1995, 2000, 2003), labels = c("Entrenamiento", "Validación", "Prueba")))
table(cut(Pares_curvas$YEAR[Con_alt_sens], c(1980, 1995, 2000, 2003), labels = c("Entrenamiento", "Validación", "Prueba")))
#Pares con pico aislado conservado en cada selección (para el análisis de sensibilidad)
c(Con_curva = sum(Pares_curvas$Pico), En_malla = sum(Pares_curvas$Pico & En_malla),
  Con_altitud = sum(Pares_curvas$Pico & Con_alt),
  Con_altitud_sens = sum(Pares_curvas$Pico & Con_alt_sens))

#Exportación de la tabla de altitud (pixel = NA: estación ubicada por municipio)
write.csv(data.frame(
  station_id = Estaciones_curvas$STATION_ID,
  lon = Estaciones_curvas$Lon,
  lat = Estaciones_curvas$Lat,
  pixel = Estaciones_curvas$Pixel,
  comparte_pixel = Estaciones_curvas$Comparte_pixel,
  ubicacion = ifelse(Estaciones_curvas$STATION_ID %in% Sin_coord_valle, "municipio_catalogo", "coordenadas"),
  altitud_catalogo = Estaciones_curvas$Altitud_catalogo,
  altitud_dem_mediana3 = Estaciones_curvas$Altitud_DEM_mediana3,
  pendiente_dem = Estaciones_curvas$Pendiente_DEM,
  altitud_usada = Estaciones_curvas$Altitud,
  altitud_fuente = Estaciones_curvas$Altitud_fuente,
  altitud_error_mae = Estaciones_curvas$Altitud_error_MAE,
  altitud_error_rmse = Estaciones_curvas$Altitud_error_RMSE,
  altitud_sensibilidad = Estaciones_curvas$Altitud_sensibilidad,
  altitud_fuente_sensibilidad = Estaciones_curvas$Altitud_fuente_sensibilidad),
  Ruta("Datos_preparados/estaciones_altitud.csv"), row.names = FALSE)

#ESTACIONES RESERVADAS (ítem 4.4): para predecir una estación hacen falta su píxel del Valle, su altitud y el ONI del
#año. Sin archivo de reservadas, la tabla sale vacía
n_res <- length(Reservadas)
Pos_res <- match(Reservadas, Coord_est$CODIGO)
Diag_reservadas <- data.frame(STATION_ID = Reservadas,
                              En_los_datos = Reservadas %in% Estaciones$STATION_ID,
                              Del_Valle = Reservadas %in% c(Coord_est$CODIGO[Coord_est$Dentro_Valle], Sin_coord_valle),
                              Con_coord = !is.na(Pos_res),
                              Pixel = rep(NA_real_, n_res),
                              Altitud_catalogo = Elevacion$alt[match(Reservadas, Elevacion$codigo)],
                              Altitud_usada = Estaciones_curvas$Altitud[match(Reservadas, Estaciones_curvas$STATION_ID)],
                              Pares_con_curva = as.vector(table(factor(Pares_curvas$STATION_ID, levels = Reservadas))))
if (any(Diag_reservadas$Con_coord))
  Diag_reservadas$Pixel[Diag_reservadas$Con_coord] <- cellFromXY(Sat, cbind(Coord_est$Lon[Pos_res[Diag_reservadas$Con_coord]],
                                                                            Coord_est$lat[Pos_res[Diag_reservadas$Con_coord]]))
Diag_reservadas$Pixel_en_Valle <- Diag_reservadas$Pixel %in% Pixeles_valle
#Se puede predecir si tiene píxel del Valle y altitud (si solo falta la altitud, aplicar el relleno con el DEM)
Diag_reservadas$Se_puede_predecir <- Diag_reservadas$Pixel_en_Valle &
                                     (!is.na(Diag_reservadas$Altitud_usada) | !is.na(Diag_reservadas$Altitud_catalogo))
Diag_reservadas

#Ambas fuentes sobre el mismo dominio anual: estación 2613600105 (ACUEDUCTO BUGA) y su píxel, 1985 (eje cortado)
Px_ej    <- Estaciones_curvas$Pixel[Estaciones_curvas$STATION_ID == "2613600105"]
i_px_ej  <- match(Px_ej, Pixeles_valle)    #Posición del píxel de la estación dentro de los píxeles del Valle
i_sat_ej <- which(Curvas_sat$Pixel == Px_ej & Curvas_sat$Año == 1985)
x_est <- as.vector(Phi_f %*% Coef_est[, i_ej])
x_sat <- as.vector(Phi_f_sat %*% Coef_sat[, i_sat_ej])
pdf(Ruta("Figuras_proyecto/Fig_Mismo_dominio.pdf"), width = 8, height = 5)
par(mfrow = c(1, 1), mar = c(4.2, 4.2, 2, 1))
plot(Par_ej$t, Par_ej$y, pch = 19, cex = 0.4, col = "grey60", xlim = c(0, 365), xlab = "Día del año (t)",
     ylim = c(min(0, x_est), 3 * max(x_est, x_sat, Y_sat[Pent$Año == 1985, i_px_ej])),
     ylab = "Lluvia (mm/día)", main = "Estación 2613600105 y su píxel satelital, 1985: mismo dominio [0, 365]")
abline(h = 0, lty = 3)
points(t72, Y_sat[Pent$Año == 1985, i_px_ej], pch = 17, cex = 0.7, col = "darkorange")
lines(t_fino, x_est, lwd = 2, col = "blue")
lines(t_fino, x_sat, lwd = 2, col = "red")
legend("topright", c("Estación: datos diarios", "Satélite: datos por péntada",
                     paste0("Curva de la estación (K = ", K_estaciones, ")"),
                     paste0("Curva del satélite (K = ", K_satelite, ")")),
       col = c("grey60", "darkorange", "blue", "red"), pch = c(19, 17, NA, NA), lty = c(NA, NA, 1, 1), lwd = c(NA, NA, 2, 2), bty = "n", cex = 0.8)
dev.off()
#Mapa: lluvia media del satélite en entrenamiento (1981-1995) y estaciones emparejadas, según tengan o no altitud
Lluvia_media_pixel <- rast(Sat[[1]])
#Los píxeles fuera del Valle quedan como NA: la media solo usa las mediciones del departamento
Media_pixel <- rep(NA_real_, ncell(Sat))
Media_pixel[Pixeles_valle] <- rowSums(Valores_sat[, Pent$Año <= 1995]) / sum(Pent$Dias[Pent$Año <= 1995])   #mm/día
values(Lluvia_media_pixel) <- Media_pixel
En_malla_est <- Estaciones_curvas$Pixel %in% Pixeles_valle
Con_alt_cat_est <- En_malla_est & !is.na(Estaciones_curvas$Altitud_catalogo)
Con_alt_dem_est <- En_malla_est & Estaciones_curvas$Altitud_fuente == "DEM Copernicus GLO-30"
Sin_alt_est     <- En_malla_est & is.na(Estaciones_curvas$Altitud)

#Ventana que contiene el raster y el contorno del Valle (el norte del departamento supera el borde del raster)
Lim <- union(ext(Lluvia_media_pixel), ext(Valle))

pdf(Ruta("Figuras_proyecto/Fig_Mapa_estaciones.pdf"), width = 7, height = 7)
par(mfrow = c(1, 1), mar = c(4.2, 4.2, 2, 4))
plot(Lluvia_media_pixel, xlim = c(xmin(Lim), xmax(Lim)), ylim = c(ymin(Lim), ymax(Lim)),
     main = "Lluvia media del satélite 1981-1995 (mm/día) y estaciones", xlab = "Longitud", ylab = "Latitud")
#Halo blanco para que el límite se distinga sobre cualquier color del raster
lines(Valle, col = "white", lwd = 3)
lines(Valle, col = "black", lwd = 1.4)
points(Estaciones_curvas$Lon[Con_alt_cat_est], Estaciones_curvas$Lat[Con_alt_cat_est], pch = 21, bg = "red", cex = 0.9)
points(Estaciones_curvas$Lon[Con_alt_dem_est], Estaciones_curvas$Lat[Con_alt_dem_est], pch = 24, bg = "#e9c46a", cex = 0.85)
points(Estaciones_curvas$Lon[Sin_alt_est], Estaciones_curvas$Lat[Sin_alt_est], pch = 21, bg = "white", cex = 0.9)
#(la leyenda se ubica con coordenadas dentro de la malla para que no quede recortada)
legend(x = -77.53, y = 4.95,
       c("Estación con altitud de catálogo",
         "Estación con altitud DEM (500-2000 m, pendiente <15)",
         "Estación sin altitud",
         "Límite del Valle del Cauca"),
       pch = c(21, 24, 21, NA), pt.bg = c("red", "#e9c46a", "white", NA),
       lty = c(NA, NA, NA, 1), lwd = c(NA, NA, NA, 1.4),
       col = "black", bg = "white", cex = 0.72)
dev.off()

#ÍNDICE ONI (NOAA CPC): anomalía de la temperatura del mar en la región Niño 3.4, promedio móvil de 3 meses
ONI <- read.table("Datos_externos/oni.ascii.txt", header = TRUE)
ONI_periodo <- ONI[ONI$YR >= 1981 & ONI$YR <= 2003, ]
#Control: 12 trimestres por año en los 23 años
table(table(ONI_periodo$YR))

#Índice del año: promedio de los 12 trimestres del año calendario (DEF a NDE)
ONI_año <- aggregate(ANOM ~ YR, data = ONI_periodo, FUN = mean)
names(ONI_año) <- c("Año", "ONI")
#Fase con el umbral de la NOAA: El Niño si ONI >= +0,5; La Niña si ONI <= -0,5; neutro en otro caso
ONI_año$Fase <- ifelse(ONI_año$ONI >= 0.5, "Niño", ifelse(ONI_año$ONI <= -0.5, "Niña", "Neutro"))
#Años de transición: al menos 3 trimestres de El Niño y 3 de La Niña en el mismo año. La media anual los deja
#como neutros; se marcan para repetir las comparaciones por fase sin ellos (análisis de sensibilidad)
ONI_año$Trim_niño  <- as.vector(tapply(ONI_periodo$ANOM >= 0.5, ONI_periodo$YR, sum))
ONI_año$Trim_niña  <- as.vector(tapply(ONI_periodo$ANOM <= -0.5, ONI_periodo$YR, sum))
ONI_año$Transicion <- ONI_año$Trim_niño >= 3 & ONI_año$Trim_niña >= 3
ONI_año
table(ONI_año$Fase, cut(ONI_año$Año, c(1980, 1995, 2000, 2003), labels = c("Entrenamiento", "Validación", "Prueba")))
#Exportación para Python: el ONI del año entra al vector z de la red
write.csv(data.frame(year = ONI_año$Año, oni = ONI_año$ONI, fase = ONI_año$Fase, transicion = ONI_año$Transicion),
          Ruta("Datos_preparados/oni_anual.csv"), row.names = FALSE)

#ANÁLISIS EXPLORATORIO FUNCIONAL (ítem 2), en el orden del script de clase: media, varianza, covarianza y correlación,
#derivadas, profundidad FM y boxplot funcional MBD. Todas las curvas de 1981-2003 (estaciones y píxeles del Valle), en
#mm/día sin recortar y evaluadas en t_fino

X_est <- Phi_f %*% Coef_est
X_sat <- Phi_f_sat %*% Coef_sat
fd_est <- fd(coef = Coef_est, basisobj = Base_fda)
fd_sat <- fd(coef = Coef_sat, basisobj = Base_fda_sat)
c(Curvas_estaciones = ncol(X_est), Estaciones = length(unique(Pares_curvas$STATION_ID)),
  Curvas_satelite = ncol(X_sat), Pixeles = length(Pixeles_valle), Años = length(Años_sat))

#Media funcional paso a paso y verificación con mean.fd
Media_est_paso <- as.vector(Phi_f %*% rowMeans(Coef_est))
Media_est_lib  <- as.vector(eval.fd(t_fino, mean.fd(fd_est)))
dif_media_est  <- max(abs(Media_est_paso - Media_est_lib))
Media_sat_paso <- as.vector(Phi_f_sat %*% rowMeans(Coef_sat))
Media_sat_lib  <- as.vector(eval.fd(t_fino, mean.fd(fd_sat)))
dif_media_sat  <- max(abs(Media_sat_paso - Media_sat_lib))
c(Dif_media_estaciones = dif_media_est, Dif_media_satelite = dif_media_sat)

#Varianza y desviación estándar funcional paso a paso y verificación con var.fd
Var_est_paso <- apply((X_est - Media_est_paso)^2, 1, sum) / (ncol(X_est) - 1)
SD_est_paso  <- sqrt(Var_est_paso)
Var_est_lib  <- diag(eval.bifd(t_fino, t_fino, var.fd(fd_est)))
dif_var_est  <- max(abs(Var_est_paso - Var_est_lib))
Var_sat_paso <- apply((X_sat - Media_sat_paso)^2, 1, sum) / (ncol(X_sat) - 1)
SD_sat_paso  <- sqrt(Var_sat_paso)
Var_sat_lib  <- diag(eval.bifd(t_fino, t_fino, var.fd(fd_sat)))
dif_var_sat  <- max(abs(Var_sat_paso - Var_sat_lib))
c(Dif_var_estaciones = dif_var_est, Dif_var_satelite = dif_var_sat)
CV_est <- SD_est_paso / Media_est_paso
CV_sat <- SD_sat_paso / Media_sat_paso

#Covarianza y correlación funcional de cada fuente
#Covarianza paso a paso: Cov(s, t) = phi(s)' Sigma_c phi(t)
CoefCentrados_est <- Coef_est - rowMeans(Coef_est)
SigmaC_est        <- CoefCentrados_est %*% t(CoefCentrados_est) / (ncol(Coef_est) - 1)
Cov_est_paso      <- Phi_f %*% SigmaC_est %*% t(Phi_f)
CoefCentrados_sat <- Coef_sat - rowMeans(Coef_sat)
SigmaC_sat        <- CoefCentrados_sat %*% t(CoefCentrados_sat) / (ncol(Coef_sat) - 1)
Cov_sat_paso      <- Phi_f_sat %*% SigmaC_sat %*% t(Phi_f_sat)

#Covarianza con var.fd() y correlación con cor.fd()
Cov_est <- eval.bifd(t_fino, t_fino, var.fd(fd_est))
Cor_est <- cor.fd(t_fino, fd_est)
Cov_sat <- eval.bifd(t_fino, t_fino, var.fd(fd_sat))
Cor_sat <- cor.fd(t_fino, fd_sat)
c(Dif_cov_est = max(abs(Cov_est_paso - Cov_est)),
  Dif_cov_sat = max(abs(Cov_sat_paso - Cov_sat)))

#Control: la diagonal de la covarianza es la varianza funcional
c(Dif_diag_cov_est = max(abs(diag(Cov_est) - Var_est_paso)),
  Dif_diag_cov_sat = max(abs(diag(Cov_sat) - Var_sat_paso)))

#Derivadas de la media con la derivada exacta de la base (deriv.fd reproyecta la derivada y la deforma)
DPhi_f      <- eval.basis(t_fino, Base_fda, 1)
D2Phi_f     <- eval.basis(t_fino, Base_fda, 2)
DPhi_f_sat  <- eval.basis(t_fino, Base_fda_sat, 1)
D2Phi_f_sat <- eval.basis(t_fino, Base_fda_sat, 2)
DMedia_est  <- as.vector(DPhi_f %*% rowMeans(Coef_est))
D2Media_est <- as.vector(D2Phi_f %*% rowMeans(Coef_est))
DMedia_sat  <- as.vector(DPhi_f_sat %*% rowMeans(Coef_sat))
D2Media_sat <- as.vector(D2Phi_f_sat %*% rowMeans(Coef_sat))

#Control: la media de las derivadas es la derivada de la media (eval.fd con Lfdobj = 1)
c(Dif_D1_est = max(abs(rowMeans(DPhi_f %*% Coef_est) - eval.fd(t_fino, mean.fd(fd_est), 1))),
  Dif_D1_sat = max(abs(rowMeans(DPhi_f_sat %*% Coef_sat) - eval.fd(t_fino, mean.fd(fd_sat), 1))))

#Control sin librerías: diferencias finitas de la media en la malla fina
h_fino <- t_fino[2] - t_fino[1]
t_medio <- t_fino[-1] - h_fino / 2
c(Dif_finita_est = max(abs(diff(Media_est_paso) / h_fino - approx(t_fino, DMedia_est, t_medio)$y)),
  Dif_finita_sat = max(abs(diff(Media_sat_paso) / h_fino - approx(t_fino, DMedia_sat, t_medio)$y)))

#Función de fecha: el día d del año ocupa el intervalo (d - 1, d] de t
Meses <- c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")
fecha_de <- function(t) {
  Dia <- as.Date("1981-01-01") + pmax(ceiling(t), 1) - 1
  ifelse(is.na(t), "—", paste0(format(Dia, "%d"), "-", Meses[as.integer(format(Dia, "%m"))]))
}

#Función de temporadas de lluvia a partir de las derivadas de la media. Picos: los dos máximos locales más altos.
#Valles: mínimo antes del primer pico, entre los picos y después del segundo (o el borde del año). Inicio: máximo de
#D1 entre el valle anterior y el pico; fin: mínimo de D1 entre el pico y el valle siguiente (D2 = 0 en ambos)
Temporadas <- function(Media, D1) {
  n       <- length(Media)
  Maximos <- which(diff(sign(D1)) == -2) + 1
  Minimos <- which(diff(sign(D1)) == 2) + 1
  Picos   <- sort(Maximos[order(Media[Maximos], decreasing = TRUE)][1:2])
  Valles  <- c(max(c(1, Minimos[Minimos < Picos[1]])),
               Picos[1] - 1 + which.min(Media[Picos[1]:Picos[2]]),
               min(c(n, Minimos[Minimos > Picos[2]])))
  Inicio  <- c(Valles[1] - 1 + which.max(D1[Valles[1]:Picos[1]]),
               Valles[2] - 1 + which.max(D1[Valles[2]:Picos[2]]))
  Fin     <- c(Picos[1] - 1 + which.min(D1[Picos[1]:Valles[2]]),
               Picos[2] - 1 + which.min(D1[Picos[2]:Valles[3]]))
  list(Picos = Picos, Valles = Valles, Inicio = Inicio, Fin = Fin)
}
Temp_est <- Temporadas(Media_est_paso, DMedia_est)
Temp_sat <- Temporadas(Media_sat_paso, DMedia_sat)

#Tabla de temporadas (En_borde = el inicio o el fin cae en los primeros o últimos 15 días del año)
t_inicio <- t_fino[c(Temp_est$Inicio, Temp_sat$Inicio)]
t_fin    <- t_fino[c(Temp_est$Fin, Temp_sat$Fin)]
Temporadas_eda <- data.frame(
  Fuente    = rep(c("Estaciones", "Satélite"), each = 2),
  Temporada = rep(c("Primera", "Segunda"), 2),
  Inicio    = fecha_de(t_inicio),
  Pico      = fecha_de(t_fino[c(Temp_est$Picos, Temp_sat$Picos)]),
  Fin       = fecha_de(t_fin),
  Dias      = t_fin - t_inicio,
  En_borde  = t_inicio < 15 | t_fin > 350)
Temporadas_eda
#Entre fuentes se comparan las fechas, no la duración (depende de K). La base no es periódica: los extremos del año
#tienen más incertidumbre (En_borde)

#Valles de la media (centro de las temporadas secas)
Valle_fecha <- function(v) ifelse(v == 1 | v == length(t_fino), "borde del año", fecha_de(t_fino[v]))
data.frame(Fuente          = c("Estaciones", "Satélite"),
           Valle_previo    = c(Valle_fecha(Temp_est$Valles[1]), Valle_fecha(Temp_sat$Valles[1])),
           Valle_entre     = c(Valle_fecha(Temp_est$Valles[2]), Valle_fecha(Temp_sat$Valles[2])),
           Valle_posterior = c(Valle_fecha(Temp_est$Valles[3]), Valle_fecha(Temp_sat$Valles[3])))

#Figura: media, primera y segunda derivada con las temporadas marcadas
Figura_temporadas <- function(Media, D1, D2, Tm, Fuente) {
  plot(t_fino, Media, type = "l", lwd = 3, col = "slategray", xlim = c(0, 365),
       xlab = "Día del año", ylab = "Lluvia (mm/día)", main = paste0(Fuente, ": media"))
  rect(t_fino[Tm$Inicio], par("usr")[3], t_fino[Tm$Fin], par("usr")[4], col = adjustcolor("skyblue", 0.3), border = NA)
  lines(t_fino, Media, lwd = 3, col = "slategray")
  abline(v = t_fino[Tm$Picos], col = "red", lty = 2)
  legend("topleft", c("Media", "Temporada de lluvias", "Pico"), col = c("slategray", "skyblue", "red"),
         lty = c(1, NA, 2), lwd = c(3, NA, 1), pch = c(NA, 15, NA), bty = "n", cex = 0.75)
  plot(t_fino, D1, type = "l", lwd = 3, col = "slategray", xlim = c(0, 365),
       xlab = "Día del año", ylab = "D1 (mm/día por día)", main = paste0(Fuente, ": primera derivada"))
  abline(h = 0, lty = 3)
  abline(v = t_fino[c(Tm$Inicio, Tm$Fin)], col = "blue", lty = 2)
  abline(v = t_fino[Tm$Picos], col = "red", lty = 2)
  legend("bottomleft", c("D1", "Inicio y fin", "Pico"), col = c("slategray", "blue", "red"),
         lty = c(1, 2, 2), lwd = c(3, 1, 1), bty = "n", cex = 0.75)
  plot(t_fino, D2, type = "l", lwd = 3, col = "purple", xlim = c(0, 365),
       xlab = "Día del año", ylab = expression("D2 (mm/día por día"^2*")"), main = paste0(Fuente, ": segunda derivada"))
  abline(h = 0, lty = 3)
  abline(v = t_fino[c(Tm$Inicio, Tm$Fin)], col = "blue", lty = 2)
  legend("bottomleft", c("D2", "Inicio y fin (D2 = 0)"), col = c("purple", "blue"),
         lty = c(1, 2), lwd = c(3, 1), bty = "n", cex = 0.75)
}
pdf(Ruta("Figuras_proyecto/Fig_EDA_derivadas2.pdf"), width = 11, height = 10)
par(mfcol = c(3, 2), mar = c(4.1, 4.1, 2, 1))
Figura_temporadas(Media_est_paso, DMedia_est, D2Media_est, Temp_est, "Estaciones")
Figura_temporadas(Media_sat_paso, DMedia_sat, D2Media_sat, Temp_sat, "Satélite")
dev.off()
Mostrar({
  par(mfcol = c(3, 2), mar = c(4.1, 4.1, 2, 1))
  Figura_temporadas(Media_est_paso, DMedia_est, D2Media_est, Temp_est, "Estaciones")
  Figura_temporadas(Media_sat_paso, DMedia_sat, D2Media_sat, Temp_sat, "Satélite")
})

#Lluvia media en mm/día: las curvas ya están en mm (Mm_*_fino: curvas sin recortar, que conservan el total)
Mm_est_fino <- X_est
Mm_sat_fino <- X_sat
MediaP_est  <- Media_est_paso
MediaP_sat  <- Media_sat_paso
#Control: la curva media nunca es negativa
round(c(Min_media_est = min(MediaP_est), Min_media_sat = min(MediaP_sat)), 3)
#Verificación con la lluvia observada (promedio de los 365 centros de día), anual y por mes
Dia_centro <- match(t365, t_fino)
Obs_dia    <- tapply(Datos_curvas$y, pmin(ceiling(Datos_curvas$t), 365), mean)
round(c(Observada_est = mean(Obs_dia), Media_est = mean(MediaP_est[Dia_centro]),
        Observada_sat = mean(Y_sat), Media_sat = mean(MediaP_sat[Dia_centro])), 3)
round(rbind(Observada_est = tapply(as.vector(Obs_dia), Mes_dia, mean),
            Media_est     = tapply(MediaP_est[Dia_centro], Mes_dia, mean),
            Observada_sat = tapply(rowMeans(Y_sat), Pent$Mes, mean),
            Media_sat     = tapply(MediaP_sat[Dia_centro], Mes_dia, mean)), 2)

#Figuras EDA por fuente (se dibujan 250 curvas al azar; los cálculos usan todas)
set.seed(20261111)
Muestra_est <- sample(ncol(X_est), 250)
Muestra_sat <- sample(ncol(X_sat), 250)
Figura_fuente <- function(X, Muestra, Media, DE, CV, Cov, Cor, MediaP, Tm, Fuente) {
  par(mfrow = c(2, 3), mar = c(4.1, 4.1, 2, 1))
  matplot(t_fino, X[, Muestra], type = "l", lty = 1, col = "grey88",
          xlab = "Día del año", ylab = "Lluvia (mm/día)", main = paste0(Fuente, ": media ± 2 DE"))
  abline(h = 0, lty = 3)
  lines(t_fino, Media, lwd = 3, col = "blue")
  lines(t_fino, Media + 2 * DE, lwd = 2, lty = 2, col = "red")
  lines(t_fino, Media - 2 * DE, lwd = 2, lty = 2, col = "red")
  plot(t_fino, DE, type = "l", lwd = 3, col = "darkgreen",
       xlab = "Día del año", ylab = "DE funcional", main = paste0(Fuente, ": desviación estándar"))
  plot(t_fino, CV, type = "l", lwd = 3, col = "darkorange",
       xlab = "Día del año", ylab = "CV funcional", main = paste0(Fuente, ": coeficiente de variación"))
  contour(t_fino, t_fino, Cov, nlevels = 10, labcex = 0.6,
          xlab = "s (día del año)", ylab = "t (día del año)", main = paste0(Fuente, ": covarianza"))
  image(t_fino, t_fino, Cor, zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"),
        xlab = "s (día del año)", ylab = "t (día del año)", main = paste0(Fuente, ": correlación"))
  contour(t_fino, t_fino, Cor, levels = seq(-1, 1, 0.25), labcex = 0.5, add = TRUE)
  plot(t_fino, MediaP, type = "l", lwd = 3, col = "slategray", ylim = c(0, max(MediaP)),
       xlab = "Día del año", ylab = "Lluvia (mm/día)", main = paste0(Fuente, ": lluvia media y temporadas"))
  rect(t_fino[Tm$Inicio], par("usr")[3], t_fino[Tm$Fin], par("usr")[4], col = adjustcolor("skyblue", 0.3), border = NA)
  lines(t_fino, MediaP, lwd = 3, col = "slategray")
  abline(v = t_fino[Tm$Picos], col = "red", lty = 2)
}
pdf(Ruta("Figuras_proyecto/Fig_EDA_estaciones.pdf"), width = 12, height = 7.5)
Figura_fuente(X_est, Muestra_est, Media_est_paso, SD_est_paso, CV_est, Cov_est, Cor_est, MediaP_est, Temp_est, "Estaciones")
dev.off()
Mostrar(Figura_fuente(X_est, Muestra_est, Media_est_paso, SD_est_paso, CV_est, Cov_est, Cor_est, MediaP_est, Temp_est, "Estaciones"))
pdf(Ruta("Figuras_proyecto/Fig_EDA_satelite.pdf"), width = 12, height = 7.5)
Figura_fuente(X_sat, Muestra_sat, Media_sat_paso, SD_sat_paso, CV_sat, Cov_sat, Cor_sat, MediaP_sat, Temp_sat, "Satélite")
dev.off()
Mostrar(Figura_fuente(X_sat, Muestra_sat, Media_sat_paso, SD_sat_paso, CV_sat, Cov_sat, Cor_sat, MediaP_sat, Temp_sat, "Satélite"))

#ORIGEN DE LA VARIABILIDAD: Var(t) = varianza entre sitios (orografía) + varianza dentro de cada sitio (año a año, ENSO)
Sitio_est <- match(Pares_curvas$STATION_ID, unique(Pares_curvas$STATION_ID))
Sitio_sat <- match(Curvas_sat$Pixel, Pixeles_valle)
Media_sitio_est <- sapply(split(seq_along(Sitio_est), Sitio_est), function(j) rowMeans(X_est[, j, drop = FALSE]))
Media_sitio_sat <- sapply(split(seq_along(Sitio_sat), Sitio_sat), function(j) rowMeans(X_sat[, j, drop = FALSE]))
Var_entre_est  <- rowSums((Media_sitio_est[, Sitio_est] - Media_est_paso)^2) / (ncol(X_est) - 1)
Var_dentro_est <- rowSums((X_est - Media_sitio_est[, Sitio_est])^2) / (ncol(X_est) - 1)
Var_entre_sat  <- rowSums((Media_sitio_sat[, Sitio_sat] - Media_sat_paso)^2) / (ncol(X_sat) - 1)
Var_dentro_sat <- rowSums((X_sat - Media_sitio_sat[, Sitio_sat])^2) / (ncol(X_sat) - 1)
#Control: las dos partes suman la varianza funcional
c(Dif_est = max(abs(Var_entre_est + Var_dentro_est - Var_est_paso)),
  Dif_sat = max(abs(Var_entre_sat + Var_dentro_sat - Var_sat_paso)))
#Fracción de la varianza que está entre sitios: promedio en el año, mínimo y máximo
round(rbind(Estaciones = c(Media = mean(Var_entre_est / Var_est_paso), Min = min(Var_entre_est / Var_est_paso),
                           Max = max(Var_entre_est / Var_est_paso)),
            Satelite   = c(Media = mean(Var_entre_sat / Var_sat_paso), Min = min(Var_entre_sat / Var_sat_paso),
                           Max = max(Var_entre_sat / Var_sat_paso))), 3)

#LLUVIA POR FASE ENSO. Nivel: lluvia media del año (mm/día). Anomalía: diferencia relativa (%) frente a la media de su
#propio sitio (quita el efecto de qué estaciones hay en cada año)
Nivel_est <- colMeans(Mm_est_fino[Dia_centro, ])
Nivel_sat <- colMeans(Mm_sat_fino[Dia_centro, ])
Anom_est  <- 100 * (Nivel_est / ave(Nivel_est, Pares_curvas$STATION_ID) - 1)
Anom_sat  <- 100 * (Nivel_sat / ave(Nivel_sat, Curvas_sat$Pixel) - 1)
Fase_est  <- ONI_año$Fase[match(Pares_curvas$YEAR, ONI_año$Año)]
Fase_sat  <- ONI_año$Fase[match(Curvas_sat$Año, ONI_año$Año)]
Trans_est <- ONI_año$Transicion[match(Pares_curvas$YEAR, ONI_año$Año)]
Trans_sat <- ONI_año$Transicion[match(Curvas_sat$Año, ONI_año$Año)]
#Sensibilidad: sin los años de transición (1983 y 1998), que la media anual deja como neutros
Tabla_enso <- rbind(Curvas_est              = table(Fase_est),
                    Nivel_est_mm_dia        = tapply(Nivel_est, Fase_est, mean),
                    Anomalia_est_pct        = tapply(Anom_est, Fase_est, mean),
                    Anomalia_est_sin_trans  = tapply(Anom_est[!Trans_est], Fase_est[!Trans_est], mean),
                    Curvas_sat              = table(Fase_sat),
                    Nivel_sat_mm_dia        = tapply(Nivel_sat, Fase_sat, mean),
                    Anomalia_sat_pct        = tapply(Anom_sat, Fase_sat, mean),
                    Anomalia_sat_sin_trans  = tapply(Anom_sat[!Trans_sat], Fase_sat[!Trans_sat], mean))
round(Tabla_enso, 2)
#Anomalía funcional por fase (mm/día): promedio de cada curva en mm menos la curva media de su sitio
Fases <- c("Niña", "Neutro", "Niño")
MediaP_sitio_est <- sapply(split(seq_along(Sitio_est), Sitio_est), function(j) rowMeans(Mm_est_fino[, j, drop = FALSE]))
MediaP_sitio_sat <- sapply(split(seq_along(Sitio_sat), Sitio_sat), function(j) rowMeans(Mm_sat_fino[, j, drop = FALSE]))
Anom_fase_est <- sapply(Fases, function(f) rowMeans(Mm_est_fino[, Fase_est == f] - MediaP_sitio_est[, Sitio_est[Fase_est == f]]))
Anom_fase_sat <- sapply(Fases, function(f) rowMeans(Mm_sat_fino[, Fase_sat == f] - MediaP_sitio_sat[, Sitio_sat[Fase_sat == f]]))

#LLUVIA SEGÚN LA ALTITUD VALIDADA. Bandas: < 500 m (Pacífico), 500-1200 m (fondo del valle), 1200-2000 m (vertientes)
#y >= 2000 m (alta montaña); las estaciones sin altitud validada quedan fuera
Altitud_par <- Estaciones_curvas$Altitud[match(Pares_curvas$STATION_ID, Estaciones_curvas$STATION_ID)]
Banda_alt   <- cut(Altitud_par, c(-Inf, 500, 1200, 2000, Inf), right = FALSE,
                   labels = c("<500 m", "500-1200 m", "1200-2000 m", ">=2000 m"))
Tabla_altitud <- rbind(Estaciones        = tapply(Pares_curvas$STATION_ID, Banda_alt, function(x) length(unique(x))),
                       Curvas            = table(Banda_alt),
                       Nivel_mm_dia      = tapply(Nivel_est, Banda_alt, mean),
                       Anomalia_Niña_pct = tapply(Anom_est[Fase_est == "Niña"], Banda_alt[Fase_est == "Niña"], mean),
                       Anomalia_Niño_pct = tapply(Anom_est[Fase_est == "Niño"], Banda_alt[Fase_est == "Niño"], mean))
round(Tabla_altitud, 2)
c(Curvas_sin_altitud = sum(is.na(Banda_alt)),
  Estaciones_sin_altitud = length(unique(Pares_curvas$STATION_ID[is.na(Banda_alt)])))
MediaP_alt <- sapply(levels(Banda_alt), function(b) rowMeans(Mm_est_fino[, which(Banda_alt == b), drop = FALSE]))

#Mapa de la lluvia media del satélite por píxel (mm/día, todos los años), con las estaciones
Lluvia_pixel <- tapply(Nivel_sat, Curvas_sat$Pixel, mean)
Mapa_lluvia  <- rast(Sat[[1]])
Mapa_lluvia[as.numeric(names(Lluvia_pixel))] <- as.vector(Lluvia_pixel)
Mapa_lluvia  <- crop(Mapa_lluvia, Valle)
round(c(Min = min(Lluvia_pixel), Mediana = median(Lluvia_pixel), Max = max(Lluvia_pixel)), 2)

#Figura: anomalía por fase ENSO de cada fuente, lluvia media por banda de altitud y mapa de lluvia del satélite
Figura_enso_altitud <- function() {
  par(mfrow = c(2, 2), mar = c(4.1, 4.1, 2, 1))
  Col_fase <- c("blue", "grey40", "red")
  matplot(t_fino, Anom_fase_est, type = "l", lty = 1, lwd = 3, col = Col_fase,
          xlab = "Día del año", ylab = "Anomalía (mm/día)", main = "Estaciones: anomalía por fase ENSO")
  abline(h = 0, lty = 3)
  legend("topleft", Fases, col = Col_fase, lwd = 3, bty = "n", cex = 0.8)
  matplot(t_fino, Anom_fase_sat, type = "l", lty = 1, lwd = 3, col = Col_fase,
          xlab = "Día del año", ylab = "Anomalía (mm/día)", main = "Satélite: anomalía por fase ENSO")
  abline(h = 0, lty = 3)
  legend("topleft", Fases, col = Col_fase, lwd = 3, bty = "n", cex = 0.8)
  matplot(t_fino, MediaP_alt, type = "l", lty = 1, lwd = 3, log = "y", col = c("darkgreen", "orange", "brown", "purple"),
          xlab = "Día del año", ylab = "mm/día (escala log)", main = "Estaciones: lluvia media por altitud")
  legend("topleft", levels(Banda_alt), col = c("darkgreen", "orange", "brown", "purple"), lwd = 3, bty = "n", cex = 0.8)
  plot(Mapa_lluvia, col = hcl.colors(30, "YlGnBu", rev = TRUE), main = "Satélite: lluvia media (mm/día)", mar = c(3, 3, 2, 4))
  lines(Valle)
  points(Estaciones_curvas$Lon, Estaciones_curvas$Lat, pch = 21, bg = "white", cex = 0.6)
}
pdf(Ruta("Figuras_proyecto/Fig_EDA_ENSO_altitud.pdf"), width = 11, height = 8.5)
Figura_enso_altitud()
dev.off()
Mostrar(Figura_enso_altitud())

#ANÁLISIS POR AÑO DEL PERIODO DE ENTRENAMIENTO (1981-1995): el detalle año por año se limita a entrenamiento para no
#examinar uno por uno los años de validación y prueba
Años_ent_eda <- 1981:1995
Ent_est <- Pares_curvas$YEAR %in% Años_ent_eda
Ent_sat <- Curvas_sat$Año %in% Años_ent_eda
c(Curvas_est = sum(Ent_est), Estaciones = length(unique(Pares_curvas$STATION_ID[Ent_est])), Curvas_sat = sum(Ent_sat))

#Medias funcionales por año como en el ANOVA funcional de clase: medias = X Z diag(1/n_año)
Año_est   <- factor(Pares_curvas$YEAR[Ent_est], levels = Años_ent_eda)
Año_sat   <- factor(Curvas_sat$Año[Ent_sat], levels = Años_ent_eda)
Z_año_est <- model.matrix(~ Año_est - 1)
Z_año_sat <- model.matrix(~ Año_sat - 1)
n_año_est <- colSums(Z_año_est)
n_año_sat <- colSums(Z_año_sat)
MediaP_año_est <- Mm_est_fino[, Ent_est] %*% Z_año_est %*% diag(1 / n_año_est)
MediaP_año_sat <- Mm_sat_fino[, Ent_sat] %*% Z_año_sat %*% diag(1 / n_año_sat)
colnames(MediaP_año_est) <- colnames(MediaP_año_sat) <- Años_ent_eda
#Control: promedio directo de las curvas de 1985
c(Dif_est = max(abs(MediaP_año_est[, "1985"] - rowMeans(Mm_est_fino[, Pares_curvas$YEAR == 1985]))),
  Dif_sat = max(abs(MediaP_año_sat[, "1985"] - rowMeans(Mm_sat_fino[, Curvas_sat$Año == 1985]))))

#Anomalía de cada curva (%) frente a la media de entrenamiento de su propio sitio
Nivel_ent_est <- Nivel_est[Ent_est]
Nivel_ent_sat <- Nivel_sat[Ent_sat]
Anom_ent_est  <- 100 * (Nivel_ent_est / ave(Nivel_ent_est, Pares_curvas$STATION_ID[Ent_est]) - 1)
Anom_ent_sat  <- 100 * (Nivel_ent_sat / ave(Nivel_ent_sat, Curvas_sat$Pixel[Ent_sat]) - 1)

#Valle que separa las dos mitades del año (curva media de entrenamiento, regla de temporadas)
Temp_ent_est <- Temporadas(as.vector(Phi_f %*% rowMeans(Coef_est[, Ent_est])),
                           as.vector(DPhi_f %*% rowMeans(Coef_est[, Ent_est])))
Temp_ent_sat <- Temporadas(as.vector(Phi_f_sat %*% rowMeans(Coef_sat[, Ent_sat])),
                           as.vector(DPhi_f_sat %*% rowMeans(Coef_sat[, Ent_sat])))
c(Valle_est = fecha_de(t_fino[Temp_ent_est$Valles[2]]), Valle_sat = fecha_de(t_fino[Temp_ent_sat$Valles[2]]))

#Función de fechas de los dos picos de cada año: máximo de la media del año en cada mitad, en los días 15 a 350 (NA si
#queda pegado al límite de su mitad). El inicio y el fin solo se calculan para la media de todo el periodo
Interior <- which(t_fino >= 15 & t_fino <= 350)
Picos_año <- function(MediaP, Valle) {
  Mitad_1 <- Interior[Interior < Valle]
  Mitad_2 <- Interior[Interior > Valle]
  t(apply(MediaP, 2, function(m) {
    p1 <- Mitad_1[which.max(m[Mitad_1])]
    p2 <- Mitad_2[which.max(m[Mitad_2])]
    #(NA_integer_ y no NA: con un NA lógico, t_fino[NA] devuelve un vector de 731 valores y la tabla no se arma)
    if (p1 %in% range(Mitad_1)) p1 <- NA_integer_
    if (p2 %in% range(Mitad_2)) p2 <- NA_integer_
    c(Pico_1 = t_fino[p1], Pico_1_mm = m[p1], Pico_2 = t_fino[p2], Pico_2_mm = m[p2])
  }))
}
Picos_año_est <- Picos_año(MediaP_año_est, Temp_ent_est$Valles[2])
Picos_año_sat <- Picos_año(MediaP_año_sat, Temp_ent_sat$Valles[2])

#Tabla por año: estaciones con curva, lluvia media (mm/día), anomalía media (%) y fechas de los dos picos
Fase_año <- ONI_año$Fase[match(Años_ent_eda, ONI_año$Año)]
Tabla_año <- data.frame(Año        = Años_ent_eda,
                        Fase       = Fase_año,
                        ONI        = round(ONI_año$ONI[match(Años_ent_eda, ONI_año$Año)], 2),
                        Estaciones = n_año_est,
                        Lluvia_est = round(as.vector(tapply(Nivel_ent_est, Año_est, mean)), 2),
                        Anom_est   = round(as.vector(tapply(Anom_ent_est, Año_est, mean)), 1),
                        Pico_1_est = fecha_de(Picos_año_est[, "Pico_1"]),
                        Pico_2_est = fecha_de(Picos_año_est[, "Pico_2"]),
                        Lluvia_sat = round(as.vector(tapply(Nivel_ent_sat, Año_sat, mean)), 2),
                        Anom_sat   = round(as.vector(tapply(Anom_ent_sat, Año_sat, mean)), 1),
                        Pico_1_sat = fecha_de(Picos_año_sat[, "Pico_1"]),
                        Pico_2_sat = fecha_de(Picos_año_sat[, "Pico_2"]))
print(Tabla_año, row.names = FALSE)
#Años sin pico detectado en alguna mitad
c(Sin_pico_est = sum(is.na(Picos_año_est[, c("Pico_1", "Pico_2")])),
  Sin_pico_sat = sum(is.na(Picos_año_sat[, c("Pico_1", "Pico_2")])))
#Variación de las fechas de pico entre años (desviación estándar, días)
round(rbind(Estaciones = apply(Picos_año_est[, c("Pico_1", "Pico_2")], 2, sd, na.rm = TRUE),
            Satelite   = apply(Picos_año_sat[, c("Pico_1", "Pico_2")], 2, sd, na.rm = TRUE)), 1)
#Anomalía media y día medio de cada pico por fase ENSO
Por_fase <- rbind(Anom_est = tapply(Tabla_año$Anom_est, Fase_año, mean),
                  Anom_sat = tapply(Tabla_año$Anom_sat, Fase_año, mean),
                  sapply(split(as.data.frame(Picos_año_est), Fase_año), function(d) colMeans(d[, c("Pico_1", "Pico_2")], na.rm = TRUE)),
                  sapply(split(as.data.frame(Picos_año_sat), Fase_año), function(d) colMeans(d[, c("Pico_1", "Pico_2")], na.rm = TRUE)))
rownames(Por_fase) <- c("Anomalia_est_pct", "Anomalia_sat_pct", "Estaciones_pico_1", "Estaciones_pico_2",
                        "Satelite_pico_1", "Satelite_pico_2")
round(Por_fase, 1)
#Con 4 a 7 años por fase no se sacan conclusiones sobre la fecha de los picos por fase

#Figura: anomalía de cada año por fuente y fechas de los dos picos, coloreadas por fase ENSO
Col_año <- c(Niña = "blue", Neutro = "grey50", Niño = "red")[Fase_año]
Figura_por_año <- function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
  barplot(rbind(Tabla_año$Anom_est, Tabla_año$Anom_sat), beside = TRUE, names.arg = Años_ent_eda, las = 2, cex.names = 0.7,
          col = as.vector(rbind(Col_año, adjustcolor(Col_año, 0.35))), border = NA,
          ylim = 1.15 * range(c(Tabla_año$Anom_est, Tabla_año$Anom_sat)),
          ylab = "Anomalía (%)", main = "Anomalía anual de la lluvia (1981-1995)")
  abline(h = 0)
  legend("topright", c("Niña", "Neutro", "Niño"), fill = c("blue", "grey50", "red"), border = NA, bty = "n", cex = 0.8)
  legend("bottomleft", c("Estaciones (oscuro)", "Satélite (claro)"), fill = c("grey30", "grey85"), border = NA, bty = "n", cex = 0.8)
  plot(Años_ent_eda, Picos_año_est[, "Pico_1"], pch = 19, col = Col_año, ylim = c(60, 345),
       xlab = "Año", ylab = "Día del año", main = "Fechas de los dos picos de lluvia (1981-1995)")
  points(Años_ent_eda, Picos_año_est[, "Pico_2"], pch = 19, col = Col_año)
  points(Años_ent_eda, Picos_año_sat[, "Pico_1"], pch = 1, col = Col_año)
  points(Años_ent_eda, Picos_año_sat[, "Pico_2"], pch = 1, col = Col_año)
  legend("center", c("Estaciones", "Satélite"), pch = c(19, 1), bty = "n", cex = 0.8)
}
pdf(Ruta("Figuras_proyecto/Fig_EDA_por_año.pdf"), width = 11, height = 5)
Figura_por_año()
dev.off()
Mostrar(Figura_por_año())

#PROFUNDIDAD FUNCIONAL Y BOXPLOT FUNCIONAL (todas las curvas, en mm/día sin recortar). FM y MBD dependen solo del
#orden en cada t; la regla de atípicas sí depende de la escala y se aplica en mm
#Profundidad de Fraiman-Muniz con fda.usc
fdata_est  <- fdata(t(X_est), t_fino, range(t_fino))
fdata_sat  <- fdata(t(X_sat), t_fino, range(t_fino))
ProfFM_est <- depth.FM(fdata_est, trim = 0.25)
ProfFM_sat <- depth.FM(fdata_sat, trim = 0.25)
#Verificación a mano: FM_i = promedio en t de 1 - |1/2 - F_t(x_i(t))|; depth.FM() entrega 2 FM - 1
Rango_est <- t(apply(X_est, 1, rank, ties.method = "max"))
Rango_sat <- t(apply(X_sat, 1, rank, ties.method = "max"))
FM_est_mano <- colMeans(1 - abs(0.5 - Rango_est / ncol(X_est)))
FM_sat_mano <- colMeans(1 - abs(0.5 - Rango_sat / ncol(X_sat)))
c(Dif_FM_est = max(abs(ProfFM_est$dep - (2 * FM_est_mano - 1))),
  Dif_FM_sat = max(abs(ProfFM_sat$dep - (2 * FM_sat_mano - 1))))
#MBD a mano: en cada t, la curva de rango r queda en (r - 1)(n - r) + n - 1 de las n(n - 1)/2 bandas
MBD_est_mano <- colMeans((Rango_est - 1) * (ncol(X_est) - Rango_est) + ncol(X_est) - 1) / choose(ncol(X_est), 2)
MBD_sat_mano <- colMeans((Rango_sat - 1) * (ncol(X_sat) - Rango_sat) + ncol(X_sat) - 1) / choose(ncol(X_sat), 2)

#Función de la figura de profundidad: curvas ordenadas por FM con media, mediana FM y media recortada al 25 %, y
#boxplot funcional MBD (fbplot)
Figura_profundidad <- function(X, ProfFM, Media, Fuente) {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
  Orden <- order(ProfFM$dep, decreasing = TRUE)
  Sel   <- Orden[round(seq(1, length(Orden), length.out = 300))]
  matplot(t_fino, X[, rev(Sel)], type = "l", lty = 1, col = rev(gray(seq(0.15, 0.9, length.out = 300))),
          xlab = "Día del año", ylab = "Lluvia (mm/día)",
          main = paste0(Fuente, ": curvas ordenadas por profundidad FM"))
  lines(t_fino, Media, lwd = 3, col = "blue")
  lines(t_fino, as.vector(ProfFM$median$data), lwd = 3, col = "red")
  lines(t_fino, as.vector(ProfFM$mtrim$data), lwd = 2, lty = 2, col = "darkgreen")
  legend("topright", c("Media", "Mediana FM", "Media recortada 25 %"),
         col = c("blue", "red", "darkgreen"), lty = c(1, 1, 2), lwd = c(3, 3, 2), bty = "n", cex = 0.8)
  Box <- fbplot(X, t_fino, method = "MBD", xlim = c(0, 365), ylim = range(X),
                xlab = "Día del año", ylab = "Lluvia (mm/día)",
                main = paste0(Fuente, ": boxplot funcional MBD"))
  legend("topright", c("Mediana", "Región central 50 %", "Bigotes", "Atípicas"),
         col = c("black", "magenta", "blue", "red"), lty = c(1, 1, 1, 2), lwd = c(2, 6, 2, 2), bty = "n", cex = 0.72)
  invisible(Box)
}
pdf(Ruta("Figuras_proyecto/Fig_EDA_profundidad_estaciones.pdf"), width = 11, height = 5.2)
Box_est <- Figura_profundidad(X_est, ProfFM_est, Media_est_paso, "Estaciones")
dev.off()
Mostrar(Figura_profundidad(X_est, ProfFM_est, Media_est_paso, "Estaciones"))
pdf(Ruta("Figuras_proyecto/Fig_EDA_profundidad_satelite.pdf"), width = 11, height = 5.2)
Box_sat <- Figura_profundidad(X_sat, ProfFM_sat, Media_sat_paso, "Satélite")
dev.off()
Mostrar(Figura_profundidad(X_sat, ProfFM_sat, Media_sat_paso, "Satélite"))
#Verificación: la MBD a mano coincide con la de fbplot
c(Dif_MBD_est = max(abs(MBD_est_mano - Box_est$depth)),
  Dif_MBD_sat = max(abs(MBD_sat_mano - Box_sat$depth)))

#Atípicas del boxplot MBD (fuera de la región central del 50 % ampliada 1,5 veces en algún t) y verificación a mano
Atipicas_bp <- function(X, Prof) {
  Centro <- X[, order(Prof, decreasing = TRUE)[1:ceiling(ncol(X) / 2)]]
  Inf_c <- apply(Centro, 1, min); Sup_c <- apply(Centro, 1, max)
  which(colSums(X <= Inf_c - 1.5 * (Sup_c - Inf_c) | X >= Sup_c + 1.5 * (Sup_c - Inf_c)) > 0)
}
Atipicas_est <- Box_est$outpoint
Atipicas_sat <- Box_sat$outpoint
c(Atipicas_estaciones = length(Atipicas_est), Atipicas_satelite = length(Atipicas_sat),
  Regla_a_mano_est = setequal(Atipicas_bp(X_est, Box_est$depth), Atipicas_est),
  Regla_a_mano_sat = setequal(Atipicas_bp(X_sat, Box_sat$depth), Atipicas_sat))

#Curvas atípicas por estación y año (lluvia del año, fase ENSO, profundidad MBD y puesto FM, 1 = menos profunda)
Nombre_est    <- ifelse(is.na(Estaciones_curvas$Nombre), Elevacion$nombre[match(Estaciones_curvas$STATION_ID, Elevacion$codigo)],
                        Estaciones_curvas$Nombre)
Municipio_est <- Elevacion$municipio[match(Estaciones_curvas$STATION_ID, Elevacion$codigo)]
Fila_est      <- match(Pares_curvas$STATION_ID, Estaciones_curvas$STATION_ID)
Tabla_atipicas <- data.frame(STATION_ID = Pares_curvas$STATION_ID, Nombre = Nombre_est[Fila_est], Año = Pares_curvas$YEAR,
                             Fase = Fase_est, Lluvia_mm_año = round(365 * Nivel_est), Prof_MBD = round(Box_est$depth, 4),
                             Rango_FM = rank(ProfFM_est$dep))[Atipicas_est, ]
print(Tabla_atipicas[order(Tabla_atipicas$STATION_ID, Tabla_atipicas$Año), ], row.names = FALSE)
#Resumen por estación de las curvas atípicas
Atip_estacion <- table(Pares_curvas$STATION_ID[Atipicas_est])
Fila_atip     <- match(names(Atip_estacion), Estaciones_curvas$STATION_ID)
Dias_lluvia   <- tapply(Datos_curvas$y > 0, Datos_curvas$STATION_ID, mean)
print(data.frame(STATION_ID      = names(Atip_estacion),
                 Nombre          = Nombre_est[Fila_atip],
                 Municipio       = Municipio_est[Fila_atip],
                 Años_atipicos   = as.vector(Atip_estacion),
                 Años_con_curva  = as.vector(table(Pares_curvas$STATION_ID)[names(Atip_estacion)]),
                 Lluvia_mm_año   = round(365 * as.vector(tapply(Nivel_est, Pares_curvas$STATION_ID, mean)[names(Atip_estacion)])),
                 Dias_lluvia_pct = round(100 * as.vector(Dias_lluvia[names(Atip_estacion)])),
                 Altitud         = Estaciones_curvas$Altitud[Fila_atip],
                 Lon             = round(Estaciones_curvas$Lon[Fila_atip], 3)), row.names = FALSE)
#Referencia: lluvia media anual y porcentaje de días con lluvia de la estación típica (mediana de las estaciones)
round(c(Lluvia_mm_año = median(365 * tapply(Nivel_est, Pares_curvas$STATION_ID, mean)),
        Dias_lluvia_pct = median(100 * Dias_lluvia)))
#Atípicas por fase ENSO, frente al total de curvas de cada fase
rbind(Atipicas = table(factor(Fase_est[Atipicas_est], Fases)), Curvas = table(factor(Fase_est, Fases)))
#Control con FM: cuántas de las curvas menos profundas según FM son atípicas MBD
c(Menos_profundas_FM = length(Atipicas_est),
  De_ellas_atipicas_MBD = sum(order(ProfFM_est$dep)[seq_along(Atipicas_est)] %in% Atipicas_est))
#Función de dirección de las atípicas: por encima o por debajo de la región central ampliada
Direccion_atipicas <- function(X, Prof, Atip) {
  Centro <- X[, order(Prof, decreasing = TRUE)[1:ceiling(ncol(X) / 2)]]
  Inf_c <- apply(Centro, 1, min); Sup_c <- apply(Centro, 1, max); R <- Sup_c - Inf_c
  c(Por_encima = sum(colSums(X[, Atip, drop = FALSE] >= Sup_c + 1.5 * R) > 0),
    Por_debajo = sum(colSums(X[, Atip, drop = FALSE] <= Inf_c - 1.5 * R) > 0))
}
rbind(Estaciones = Direccion_atipicas(X_est, Box_est$depth, Atipicas_est),
      Satelite   = Direccion_atipicas(X_sat, Box_sat$depth, Atipicas_sat))

#Satélite: curvas atípicas por píxel (años atípicos de los 23, coordenadas y lluvia media anual del píxel)
Atip_pixel <- table(Curvas_sat$Pixel[Atipicas_sat])
Fila_px    <- match(as.numeric(names(Atip_pixel)), Curvas_sat$Pixel)
print(data.frame(Pixel         = names(Atip_pixel),
                 Lon           = round(Curvas_sat$Lon[Fila_px], 3),
                 Lat           = round(Curvas_sat$Lat[Fila_px], 3),
                 Años_atipicos = as.vector(Atip_pixel),
                 Lluvia_mm_año = round(365 * as.vector(tapply(Nivel_sat, Curvas_sat$Pixel, mean)[names(Atip_pixel)]))), row.names = FALSE)
#Referencia: lluvia media anual del píxel típico (mediana) y longitud más oriental de los píxeles con curvas atípicas
round(c(Lluvia_mm_año_mediana = median(365 * tapply(Nivel_sat, Curvas_sat$Pixel, mean)),
        Lon_max_atipicas = max(Curvas_sat$Lon[Atipicas_sat])), 3)
#Atípicas del satélite por fase ENSO, frente al total de curvas de cada fase
rbind(Atipicas = table(factor(Fase_sat[Atipicas_sat], Fases)), Curvas = table(factor(Fase_sat, Fases)))
#Satélite: las 10 curvas menos profundas según FM (FM también marca curvas extremas por debajo)
print(data.frame(Pixel = Curvas_sat$Pixel, Año = Curvas_sat$Año, Lon = round(Curvas_sat$Lon, 3), Lat = round(Curvas_sat$Lat, 3),
                 Fase = Fase_sat, Lluvia_mm_año = round(365 * Nivel_sat),
                 Prof_FM = round(ProfFM_sat$dep, 4))[order(ProfFM_sat$dep)[1:10], ], row.names = FALSE)

#RESUMEN DEL EDA POR FUENTE (lluvia en mm/día)
#Picos: máximo de la lluvia media en cada temporada; valle: mínimo entre las dos
Pico_mm <- function(MediaP, Tm, k) Tm$Inicio[k] - 1 + which.max(MediaP[Tm$Inicio[k]:Tm$Fin[k]])
Valle_mm <- function(MediaP, Tm) Tm$Fin[1] - 1 + which.min(MediaP[Tm$Fin[1]:Tm$Inicio[2]])
Resumen_fuente <- function(MediaP, Tm, Fuente, Curvas, Atipicas) {
  p1 <- Pico_mm(MediaP, Tm, 1); p2 <- Pico_mm(MediaP, Tm, 2); v <- Valle_mm(MediaP, Tm)
  data.frame(Fuente = Fuente, Curvas = Curvas,
             Lluvia_mm_dia = round(mean(MediaP[Dia_centro]), 2), Lluvia_mm_año = round(365 * mean(MediaP[Dia_centro])),
             Pico_1 = fecha_de(t_fino[p1]), Pico_1_mm = round(MediaP[p1], 2),
             Pico_2 = fecha_de(t_fino[p2]), Pico_2_mm = round(MediaP[p2], 2),
             Valle = fecha_de(t_fino[v]), Valle_mm = round(MediaP[v], 2),
             Razon_picos_valle = round(mean(MediaP[c(p1, p2)]) / MediaP[v], 2), Atipicas_MBD = length(Atipicas))
}
rbind(Resumen_fuente(MediaP_est, Temp_est, "Estaciones", ncol(X_est), Atipicas_est),
      Resumen_fuente(MediaP_sat, Temp_sat, "Satélite", ncol(X_sat), Atipicas_sat))
#El satélite promedia todos los píxeles del Valle y las estaciones sus sitios: la comparación en los mismos sitios y
#años es la del ítem 3

#Cuantiles funcionales punto a punto (mm/día)
Probs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
Q_est <- t(apply(X_est, 1, quantile, probs = Probs))
Q_sat <- t(apply(X_sat, 1, quantile, probs = Probs))
colnames(Q_est) <- colnames(Q_sat) <- paste0("P", Probs * 100)
Figura_cuantiles <- function(Q, Media, Fuente) {
  plot(t_fino, Q[, "P50"], type = "l", lwd = 3, col = "black", ylim = range(c(Q[, "P10"], Q[, "P90"])),
       xlab = "Día del año", ylab = "Lluvia (mm/día)", main = paste0(Fuente, ": cuantiles funcionales"))
  polygon(c(t_fino, rev(t_fino)), c(Q[, "P10"], rev(Q[, "P90"])), col = adjustcolor("steelblue", 0.20), border = NA)
  polygon(c(t_fino, rev(t_fino)), c(Q[, "P25"], rev(Q[, "P75"])), col = adjustcolor("steelblue", 0.45), border = NA)
  lines(t_fino, Q[, "P50"], lwd = 3, col = "black")
  lines(t_fino, Media, lwd = 2, lty = 2, col = "darkred")
  legend("topright", c("P10-P90", "P25-P75", "Mediana", "Media"), col = c("steelblue", "steelblue", "black", "darkred"),
         lty = c(1, 1, 1, 2), lwd = c(8, 8, 3, 2), bty = "n", cex = 0.75)
}
pdf(Ruta("Figuras_proyecto/Fig_EDA_cuantiles.pdf"), width = 11, height = 5)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
Figura_cuantiles(Q_est, Media_est_paso, "Estaciones")
Figura_cuantiles(Q_sat, Media_sat_paso, "Satélite")
dev.off()
Mostrar({
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
  Figura_cuantiles(Q_est, Media_est_paso, "Estaciones")
  Figura_cuantiles(Q_sat, Media_sat_paso, "Satélite")
})
#COMPARACIÓN FUNCIONAL POR FUENTE (ítem 3)
#1. EMPAREJAMIENTO ESTACIÓN-PÍXEL: cada estación con coordenadas se compara con el píxel de 0,05° que la contiene
#(cellFromXY) en el mismo año; las estaciones ubicadas por municipio no tienen píxel
Pixel_par   <- Estaciones_curvas$Pixel[match(Pares_curvas$STATION_ID, Estaciones_curvas$STATION_ID)]
Idx_sat_par <- match(paste(Pixel_par, Pares_curvas$YEAR), paste(Curvas_sat$Pixel, Curvas_sat$Año))
Con_par     <- !is.na(Idx_sat_par)
#Verificación del emparejamiento: mismo píxel, mismo año y píxel dentro del Valle
c(Pares_estacion_año = nrow(Pares_curvas), Con_pixel = sum(!is.na(Pixel_par)), Con_par = sum(Con_par),
  Pixeles_distintos = length(unique(Pixel_par[Con_par])))
stopifnot(all(Curvas_sat$Pixel[Idx_sat_par[Con_par]] == Pixel_par[Con_par]),
          all(Curvas_sat$Año[Idx_sat_par[Con_par]] == Pares_curvas$YEAR[Con_par]),
          all(Pixel_par[Con_par] %in% Pixeles_valle))

#Geometría del píxel: tamaño en km a la latitud de cada estación (1° de latitud = 110,57 km; de longitud = 111,32 cos(lat))
Km_lon    <- function(lat) 111.32 * cos(lat * pi / 180)
Fil_px    <- Estaciones_curvas$Pixel %in% Pixeles_valle
Centro_px <- xyFromCell(Sat, Estaciones_curvas$Pixel[Fil_px])
Emparejamiento <- data.frame(Estaciones_curvas[Fil_px, c("STATION_ID", "Nombre", "Lon", "Lat", "Pixel")],
                             Lon_pixel = Centro_px[, "x"], Lat_pixel = Centro_px[, "y"])
Dx_km <- abs(Emparejamiento$Lon - Emparejamiento$Lon_pixel) * Km_lon(Emparejamiento$Lat)
Dy_km <- abs(Emparejamiento$Lat - Emparejamiento$Lat_pixel) * 110.57
Emparejamiento$Dist_centro_km <- sqrt(Dx_km^2 + Dy_km^2)
Emparejamiento$Dist_borde_km  <- pmin(res(Sat)[1] / 2 * Km_lon(Emparejamiento$Lat) - Dx_km, res(Sat)[2] / 2 * 110.57 - Dy_km)
Emparejamiento$Est_en_pixel   <- as.vector(Est_pixel[as.character(Emparejamiento$Pixel)])
Emparejamiento$Años_par       <- as.vector(table(factor(Pares_curvas$STATION_ID[Con_par], levels = Emparejamiento$STATION_ID)))
#Control: distancia al centro <= media diagonal del píxel y distancia al borde >= 0
Media_diag <- sqrt((res(Sat)[1] / 2 * Km_lon(Emparejamiento$Lat))^2 + (res(Sat)[2] / 2 * 110.57)^2)
stopifnot(all(Emparejamiento$Dist_centro_km <= Media_diag), all(Emparejamiento$Dist_borde_km >= 0))
round(c(Lado_lon_km = res(Sat)[1] * Km_lon(mean(Emparejamiento$Lat)), Lado_lat_km = res(Sat)[2] * 110.57,
        Media_diagonal_km = mean(Media_diag), Dist_centro_mediana_km = median(Emparejamiento$Dist_centro_km),
        Dist_centro_maxima_km = max(Emparejamiento$Dist_centro_km)), 2)
#Estaciones a menos de 100 m del borde de su píxel (píxel incierto por la precisión de las coordenadas)
Emparejamiento[Emparejamiento$Dist_borde_km < 0.1, c("STATION_ID", "Nombre", "Pixel", "Dist_borde_km")]

#Pares por partición y curvas satelitales compartidas (en un píxel compartido, la misma curva entra en varios pares)
table(cut(Pares_curvas$YEAR[Con_par], c(1980, 1995, 2000, 2003), labels = c("Entrenamiento", "Validación", "Prueba")))
Px_año <- table(paste(Pixel_par[Con_par], Pares_curvas$YEAR[Con_par]))
c(Pares = sum(Con_par), Curvas_satelitales_distintas = length(Px_año), Usadas_mas_de_una_vez = sum(Px_año > 1),
  Pares_que_las_usan = sum(Px_año[Px_año > 1]))

#Mapa del emparejamiento: píxeles con estación (los compartidos en naranja) y estaciones
Figura_emparejamiento <- function() {
  par(mfrow = c(1, 1), mar = c(4.2, 4.2, 2, 1))
  Cx <- xyFromCell(Sat, as.numeric(names(Est_pixel)))
  h  <- res(Sat) / 2
  plot(Emparejamiento$Lon, Emparejamiento$Lat, type = "n", asp = 1, xlab = "Longitud", ylab = "Latitud",
       main = "Emparejamiento estación - píxel (0,05°)")
  lines(Valle, col = "grey50")
  rect(Cx[, 1] - h[1], Cx[, 2] - h[2], Cx[, 1] + h[1], Cx[, 2] + h[2],
       col = ifelse(Est_pixel > 1, adjustcolor("darkorange", 0.5), NA), border = "grey40")
  points(Emparejamiento$Lon, Emparejamiento$Lat, pch = 19, cex = 0.45,
         col = ifelse(Emparejamiento$Est_en_pixel > 1, "darkred", "black"))
  legend("topleft", c("Píxel con una estación", "Píxel compartido", "Estación"), pch = c(22, 22, 19),
         pt.bg = c(NA, adjustcolor("darkorange", 0.5), NA), pt.cex = c(1.6, 1.6, 0.6), bty = "n", cex = 0.8)
}
pdf(Ruta("Figuras_proyecto/Fig_Emparejamiento.pdf"), width = 7, height = 7)
Figura_emparejamiento()
dev.off()
Mostrar(Figura_emparejamiento())

#2. CURVA DE DIFERENCIA d(t) = ESTACIÓN - SATÉLITE (mm/día, curvas sin recortar, píxel de la estación en el mismo año)
D_par   <- Mm_est_fino[, Con_par] - Mm_sat_fino[, Idx_sat_par[Con_par]]
Est_par <- Pares_curvas$STATION_ID[Con_par]
Año_par <- Pares_curvas$YEAR[Con_par]

#Persistencia del sesgo: correlación entre el sesgo medio anual de cada estación en años impares y en años pares
d_anual_par <- colMeans(D_par[Dia_centro, ])                       #Sesgo anual de cada par (mm/día)
Sesgo_impar <- tapply(d_anual_par[Año_par %% 2 == 1], Est_par[Año_par %% 2 == 1], mean)
Sesgo_par   <- tapply(d_anual_par[Año_par %% 2 == 0], Est_par[Año_par %% 2 == 0], mean)
Ambos <- intersect(names(Sesgo_impar), names(Sesgo_par))
round(c(Estaciones = length(Ambos), Correlacion_impar_par = cor(Sesgo_impar[Ambos], Sesgo_par[Ambos])), 3)
#Unidad de análisis: la estación; d_s(t) = media de sus años, con la matriz indicadora del ANOVA de clase
Z_est <- model.matrix(~ factor(Est_par) - 1)
Por_estacion <- function(M) M %*% Z_est %*% diag(1 / colSums(Z_est))
D_est <- Por_estacion(D_par)
#Control: media directa de los años de la primera estación
max(abs(D_est[, 1] - rowMeans(D_par[, Est_par == sort(unique(Est_par))[1], drop = FALSE])))
d_media <- rowMeans(D_est)
EE_d    <- apply(D_est, 1, sd) / sqrt(ncol(D_est))
#Referencia: misma fórmula con los pares tomados como independientes
d_media_pares <- rowMeans(D_par)
EE_pares      <- apply(D_par, 1, sd) / sqrt(ncol(D_par))
round(c(EE_estacion = mean(EE_d), EE_pares = mean(EE_pares), Razon = mean(EE_d / EE_pares)), 3)

#Niveles en los mismos sitios y años (mm/día): por estación (cada sitio pesa igual) y por par
Niveles <- rbind(Por_estacion = c(Estacion = mean(Por_estacion(Mm_est_fino[Dia_centro, Con_par])),
                                  Satelite = mean(Por_estacion(Mm_sat_fino[Dia_centro, Idx_sat_par[Con_par]])),
                                  d = mean(d_media[Dia_centro])),
                 Por_par = c(mean(Mm_est_fino[Dia_centro, Con_par]), mean(Mm_sat_fino[Dia_centro, Idx_sat_par[Con_par]]),
                             mean(d_media_pares[Dia_centro])))
round(cbind(Niveles, d_mm_año = 365 * Niveles[, "d"]), 3)

#Función de intervalos del año en que la banda d(t) ± 2 EE(t) no contiene el cero (Borde: tramos en los primeros o
#últimos 15 días, que no se interpretan)
Intervalos_sin_cero <- function(dm, ee) {
  Sin_cero  <- dm - 2 * ee > 0 | dm + 2 * ee < 0
  Tramos    <- rle(Sin_cero)
  Fin_tramo <- cumsum(Tramos$lengths)
  Ini_tramo <- Fin_tramo - Tramos$lengths + 1
  Tabla <- data.frame(Inicio = fecha_de(t_fino[Ini_tramo]), Fin = fecha_de(t_fino[Fin_tramo]),
                      Dias = t_fino[Fin_tramo] - t_fino[Ini_tramo],
                      Signo = ifelse(dm[Ini_tramo] > 0, "estación > satélite", "satélite > estación"),
                      d_medio = sapply(seq_along(Ini_tramo), function(k) round(mean(dm[Ini_tramo[k]:Fin_tramo[k]]), 2)),
                      Borde = t_fino[Fin_tramo] < 15 | t_fino[Ini_tramo] > 350)[Tramos$values, ]
  list(Tabla = Tabla, Pct_año = round(100 * mean(Sin_cero), 1),
       Pct_interior = round(100 * mean(Sin_cero[t_fino >= 15 & t_fino <= 350]), 1))
}
Int_est   <- Intervalos_sin_cero(d_media, EE_d)
Int_pares <- Intervalos_sin_cero(d_media_pares, EE_pares)
print(Int_est$Tabla, row.names = FALSE)
c(Estacion_pct_año = Int_est$Pct_año, Estacion_pct_dias_15_350 = Int_est$Pct_interior,
  Pares_pct_año = Int_pares$Pct_año, Pares_pct_dias_15_350 = Int_pares$Pct_interior)

#Figura: diferencia media de cada estación (eje cortado en los percentiles 0,5 y 99,5) y banda con la estación como
#unidad frente a la de pares independientes
Figura_diferencia <- function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
  matplot(t_fino, D_est, type = "l", lty = 1, col = "grey80", ylim = quantile(D_est, c(0.005, 0.995)),
          xlab = "Día del año", ylab = "Estación - satélite (mm/día)", main = "Diferencia media de cada estación")
  lines(t_fino, d_media, lwd = 3, col = "darkred")
  abline(h = 0, lty = 3)
  legend("topright", c(paste0("Estaciones (", ncol(D_est), ")"), "Media"), col = c("grey80", "darkred"), lty = 1,
         lwd = c(1, 3), bty = "n", cex = 0.8)
  plot(t_fino, d_media, type = "n", xlab = "Día del año", ylab = "mm/día", main = "Diferencia media ± 2 EE",
       ylim = range(c(d_media - 2 * EE_d, d_media + 2 * EE_d, d_media_pares - 2 * EE_pares, d_media_pares + 2 * EE_pares)))
  rect(c(0, 350), par("usr")[3], c(15, 365), par("usr")[4], col = "grey92", border = NA)   #Bordes del año
  polygon(c(t_fino, rev(t_fino)), c(d_media - 2 * EE_d, rev(d_media + 2 * EE_d)), col = adjustcolor("darkred", 0.2), border = NA)
  lines(t_fino, d_media, lwd = 3, col = "darkred")
  lines(t_fino, d_media_pares, lwd = 1.5, col = "grey40")
  lines(t_fino, d_media_pares + 2 * EE_pares, lty = 2, col = "grey40")
  lines(t_fino, d_media_pares - 2 * EE_pares, lty = 2, col = "grey40")
  abline(h = 0, lty = 3)
  legend("topright", c("Media por estación", "± 2 EE (estación como unidad)", "Referencia: pares independientes", "Borde del año"),
         col = c("darkred", adjustcolor("darkred", 0.2), "grey40", "grey92"), lty = c(1, NA, 2, NA), lwd = c(3, NA, 1, NA),
         pch = c(NA, 15, NA, 15), pt.cex = 2, bty = "n", cex = 0.75)
}
pdf(Ruta("Figuras_proyecto/Fig_Diferencia_fuentes.pdf"), width = 11, height = 5)
Figura_diferencia()
dev.off()
Mostrar(Figura_diferencia())

#3. VERIFICACIÓN CON LOS DATOS SIN SUAVIZAR: DIFERENCIA MENSUAL EN LOS MISMOS PARES (las curvas tienen K distinto).
#Estación: media de los días observados del mes; satélite: lluvia de las 6 péntadas del mes entre los días del mes
#Mes de cada dato diario: t = (día del año - 0,5) 365 / días del año
Dias_año_dc <- ifelse((Datos_curvas$YEAR %% 4 == 0 & Datos_curvas$YEAR %% 100 != 0) | Datos_curvas$YEAR %% 400 == 0, 366, 365)
Dia_año_dc  <- round(Datos_curvas$t * Dias_año_dc / 365 + 0.5)
Mes_dc      <- as.numeric(format(as.Date(paste0(Datos_curvas$YEAR, "-01-01")) + Dia_año_dc - 1, "%m"))
#Control: días distintos de cada mes en años de 365 y de 366 días
rbind(Año_365 = tapply(Dia_año_dc[Dias_año_dc == 365], Mes_dc[Dias_año_dc == 365], function(d) length(unique(d))),
      Año_366 = tapply(Dia_año_dc[Dias_año_dc == 366], Mes_dc[Dias_año_dc == 366], function(d) length(unique(d))))
Clave_dc <- factor(paste(Datos_curvas$STATION_ID, Datos_curvas$YEAR), levels = paste(Pares_curvas$STATION_ID, Pares_curvas$YEAR))
Est_mes  <- tapply(Datos_curvas$y, list(Clave_dc, Mes_dc), mean)                 #Pares x 12 (mm/día)
#Satélite: lluvia mensual de cada píxel-año (mm/día)
Sat_mes <- t(sapply(seq_len(nrow(Curvas_sat)), function(k) {
  f <- Pent$Año == Curvas_sat$Año[k]
  i <- match(Curvas_sat$Pixel[k], Pixeles_valle)
  tapply(Valores_sat[i, f], Pent$Mes[f], sum) / tapply(Pent$Dias[f], Pent$Mes[f], sum)
}))
#Control (1981): los 12 meses por sus días suman la lluvia del año
Dias_mes_1981 <- tapply(Pent$Dias[Pent$Año == 1981], Pent$Mes[Pent$Año == 1981], sum)
k_1981 <- Curvas_sat$Año == 1981
max(abs(Sat_mes[k_1981, ] %*% Dias_mes_1981 - rowSums(Valores_sat[match(Curvas_sat$Pixel[k_1981], Pixeles_valle), Pent$Año == 1981])))
#Diferencia mensual de cada par y media de los años de cada estación (los pares-mes sin ningún día observado no entran)
D_mes <- t(Est_mes[Con_par, ] - Sat_mes[Idx_sat_par[Con_par], ])               #12 x pares
D_mes_est <- sapply(split(seq_len(ncol(D_mes)), Est_par), function(j) rowMeans(D_mes[, j, drop = FALSE], na.rm = TRUE))
#Lo mismo con las curvas: promedio de d_s(t) de cada estación en los días de cada mes
D_mes_curva <- apply(D_est[Dia_centro, ], 2, function(d) tapply(d, Mes_dia, mean))
#Media y EE entre estaciones, como en la banda (en cada mes, las estaciones con datos en ese mes)
Media_EE  <- function(M) list(m = rowMeans(M, na.rm = TRUE), ee = apply(M, 1, sd, na.rm = TRUE) / sqrt(rowSums(!is.na(M))))
Mes_datos <- Media_EE(D_mes_est)
Mes_curva <- Media_EE(D_mes_curva)
Banda_mes <- function(m, ee) ifelse(m - 2 * ee > 0, "estación > satélite", ifelse(m + 2 * ee < 0, "satélite > estación", "contiene el 0"))
Tabla_mes <- data.frame(Mes = Meses, Datos = round(Mes_datos$m, 2), EE_datos = round(Mes_datos$ee, 2),
                        Banda_datos = Banda_mes(Mes_datos$m, Mes_datos$ee),
                        Curva = round(Mes_curva$m, 2), EE_curva = round(Mes_curva$ee, 2),
                        Banda_curva = Banda_mes(Mes_curva$m, Mes_curva$ee))
Tabla_mes$Coinciden <- Tabla_mes$Banda_datos == Tabla_mes$Banda_curva
Tabla_mes
#Año completo: datos frente a curva (donde los meses no coinciden, manda el dato mensual)

Figura_diferencia_mes <- function() {
  par(mfrow = c(1, 1), mar = c(4.2, 4.2, 2, 1))
  Lim <- range(c(Mes_datos$m - 2 * Mes_datos$ee, Mes_datos$m + 2 * Mes_datos$ee, Mes_curva$m - 2 * Mes_curva$ee, Mes_curva$m + 2 * Mes_curva$ee))
  plot(1:12 - 0.12, Mes_datos$m, pch = 19, col = "black", ylim = Lim, xaxt = "n", xlim = c(0.5, 12.5),
       xlab = "Mes", ylab = "Estación - satélite (mm/día)", main = "Diferencia mensual ± 2 EE (estación como unidad)")
  axis(1, at = 1:12, labels = Meses)
  segments(1:12 - 0.12, Mes_datos$m - 2 * Mes_datos$ee, 1:12 - 0.12, Mes_datos$m + 2 * Mes_datos$ee, col = "black")
  points(1:12 + 0.12, Mes_curva$m, pch = 17, col = "darkred")
  segments(1:12 + 0.12, Mes_curva$m - 2 * Mes_curva$ee, 1:12 + 0.12, Mes_curva$m + 2 * Mes_curva$ee, col = "darkred")
  abline(h = 0, lty = 3)
  legend("bottomleft", c("Datos mensuales", "Curvas d(t)"), col = c("black", "darkred"), pch = c(19, 17), lty = 1, bty = "n", cex = 0.8)
}
pdf(Ruta("Figuras_proyecto/Fig_Diferencia_mensual.pdf"), width = 8, height = 5)
Figura_diferencia_mes()
dev.off()
Mostrar(Figura_diferencia_mes())

#4. COVARIANZA DE CADA FUENTE Y CORRELACIÓN CRUZADA (X = satélite, Y = estación). Como en el ANOVA de clase se separan
#entre sitios (media de los años de cada estación y de su píxel en esos años) y año a año (curva de cada año menos la
#media de su estación)
C_est_par <- Coef_est[, Con_par]
C_sat_par <- Coef_sat[, Idx_sat_par[Con_par]]
C_est_sit <- Por_estacion(C_est_par)
C_sat_sit <- Por_estacion(C_sat_par)
Grupo_est <- match(Est_par, sort(unique(Est_par)))           #Columna de Z_est (estación) de cada par
stopifnot(all(colnames(Z_est) == paste0("factor(Est_par)", sort(unique(Est_par)))))
C_est_anu <- C_est_par - C_est_sit[, Grupo_est]
C_sat_anu <- C_sat_par - C_sat_sit[, Grupo_est]
#Control: curva de los coeficientes medios = media de las curvas
c(Estacion = max(abs(Phi_f %*% C_est_sit - Por_estacion(Mm_est_fino[, Con_par]))),
  Satelite = max(abs(Phi_f_sat %*% C_sat_sit - Por_estacion(Mm_sat_fino[, Idx_sat_par[Con_par]]))))
Lecturas <- list(Entre_sitios = list(X = fd(C_sat_sit, Base_fda_sat), Y = fd(C_est_sit, Base_fda)),
                 Año_a_año    = list(X = fd(C_sat_anu, Base_fda_sat), Y = fd(C_est_anu, Base_fda)))
#Año a año: var.fd divide entre N - 1 y no entre N - número de estaciones; la escala cambia igual en las dos fuentes y
#la correlación no cambia

#Función de covarianza de cada fuente y covarianza y correlación cruzadas (var.fd y cor.fd), con verificación paso a paso
Analisis_cruzado <- function(L) {
  Cov_X  <- eval.bifd(t_fino, t_fino, var.fd(L$X))
  Cov_Y  <- eval.bifd(t_fino, t_fino, var.fd(L$Y))
  Cov_XY <- eval.bifd(t_fino, t_fino, var.fd(L$X, L$Y))
  Cor_XY <- cor.fd(t_fino, L$X, t_fino, L$Y)
  X_val <- eval.fd(t_fino, L$X); Y_val <- eval.fd(t_fino, L$Y)
  Cov_XY_paso <- (X_val - rowMeans(X_val)) %*% t(Y_val - rowMeans(Y_val)) / (ncol(X_val) - 1)
  Cor_XY_paso <- Cov_XY_paso / sqrt(outer(diag(Cov_X), diag(Cov_Y)))
  list(Cov_X = Cov_X, Cov_Y = Cov_Y, Cor_X = cor.fd(t_fino, L$X), Cor_Y = cor.fd(t_fino, L$Y), Cor_XY = Cor_XY,
       Diag = diag(Cor_XY), N = ncol(X_val),
       Comprobacion = c(Cov = max(abs(Cov_XY_paso - Cov_XY)), Cor = max(abs(Cor_XY_paso - Cor_XY))))
}
Cruzado <- lapply(Lecturas, Analisis_cruzado)
sapply(Cruzado, function(A) c(N = A$N, A$Comprobacion))

#Diagonal r_XY(t, t): concordancia día a día entre las dos fuentes
round(sapply(Cruzado, function(A) summary(A$Diag[Dia_centro])), 3)
Diag_mes <- t(sapply(Cruzado, function(A) tapply(A$Diag[Dia_centro], Mes_dia, mean)))
colnames(Diag_mes) <- Meses
round(Diag_mes, 2)
#Control de desfase: mayor correlación cruzada fuera de la diagonal frente a la de la diagonal
round(sapply(Cruzado, function(A) c(Max_fuera_diagonal = max(A$Cor_XY), Max_diagonal = max(A$Diag))), 3)

#Comparación de la covarianza de cada fuente (mm/día)^2: varianza media en el año y razón de desviaciones estándar
Correlacion_a <- function(C, Dias) { h <- round(Dias / diff(t_fino[1:2])); mean(C[cbind(1:(nrow(C) - h), (1 + h):nrow(C))]) }
round(sapply(Cruzado, function(A) c(Var_estacion = mean(diag(A$Cov_Y)[Dia_centro]), Var_satelite = mean(diag(A$Cov_X)[Dia_centro]),
                                    Razon_DE = mean(sqrt(diag(A$Cov_Y) / diag(A$Cov_X))[Dia_centro]),
                                    Cor_30_dias_est = Correlacion_a(A$Cor_Y, 30), Cor_30_dias_sat = Correlacion_a(A$Cor_X, 30))), 3)
#Control de representación: el satélite proyectado a la base de las estaciones (K = K_estaciones)
C_sat8_par <- solve(t(Phi_f) %*% Phi_f, t(Phi_f) %*% (Phi_f_sat %*% C_sat_par))
C_sat8_sit <- Por_estacion(C_sat8_par)
Sat_K8 <- list(Entre_sitios = fd(C_sat8_sit, Base_fda), Año_a_año = fd(C_sat8_par - C_sat8_sit[, Grupo_est], Base_fda))
round(sapply(Sat_K8, function(f) c(Var_satelite_K8 = mean(diag(eval.bifd(t_fino, t_fino, var.fd(f)))[Dia_centro]),
                                   Cor_30_dias_sat_K8 = Correlacion_a(cor.fd(t_fino, f), 30))), 3)

#Las superficies se dibujan cada 2 días (una de cada 4 posiciones de t_fino) para que las figuras no pesen demasiado
Fg <- seq(1, length(t_fino), by = 4)
Figura_cor_cruzada <- function() {
  par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2, 1))
  for (k in 1:2) {
    image(t_fino[Fg], t_fino[Fg], Cruzado[[k]]$Cor_XY[Fg, Fg], zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"),
          xlab = "s: satélite (día)", ylab = "t: estación (día)", main = paste("Correlación cruzada:", gsub("_", " ", names(Cruzado)[k])))
    contour(t_fino[Fg], t_fino[Fg], Cruzado[[k]]$Cor_XY[Fg, Fg], levels = seq(-1, 1, 0.2), labcex = 0.7, add = TRUE)
    abline(0, 1, lty = 3, col = "red")
  }
  plot(t_fino, Cruzado$Entre_sitios$Diag, type = "l", lwd = 3, col = "darkblue", ylim = c(-1, 1),
       xlab = "Día del año", ylab = expression(r[XY](t, t)), main = "Diagonal de la correlación cruzada")
  lines(t_fino, Cruzado$Año_a_año$Diag, lwd = 3, lty = 2, col = "darkorange")
  abline(h = 0, lty = 3)
  legend("bottomleft", c("Entre sitios", "Año a año"), col = c("darkblue", "darkorange"), lty = c(1, 2), lwd = 3, bty = "n")
}
Figura_cov_fuentes <- function() {
  par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2, 1))
  for (k in 1:2) {
    A <- Cruzado[[k]]; Lectura <- gsub("_", " ", names(Cruzado)[k])
    plot(t_fino, sqrt(diag(A$Cov_Y)), type = "l", lwd = 3, col = "slategray", ylim = c(0, max(sqrt(diag(A$Cov_Y)))),
         xlab = "Día del año", ylab = "DE (mm/día)", main = paste("Desviación estándar:", Lectura))
    lines(t_fino, sqrt(diag(A$Cov_X)), lwd = 3, lty = 2, col = "darkturquoise")
    legend("topleft", c("Estaciones", "Satélite"), col = c("slategray", "darkturquoise"), lty = c(1, 2), lwd = 3, bty = "n", cex = 0.8)
    for (Fuente in c("Y", "X")) {
      C <- A[[paste0("Cor_", Fuente)]]
      image(t_fino[Fg], t_fino[Fg], C[Fg, Fg], zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"), xlab = "s (día)", ylab = "t (día)",
            main = paste0(ifelse(Fuente == "Y", "Estaciones", "Satélite"), ": correlación (", Lectura, ")"))
      contour(t_fino[Fg], t_fino[Fg], C[Fg, Fg], levels = seq(-1, 1, 0.2), labcex = 0.6, add = TRUE)
    }
  }
}
pdf(Ruta("Figuras_proyecto/Fig_Correlacion_cruzada.pdf"), width = 14, height = 5)
Figura_cor_cruzada()
dev.off()
Mostrar(Figura_cor_cruzada())
pdf(Ruta("Figuras_proyecto/Fig_Covarianza_fuentes.pdf"), width = 13, height = 8.5)
Figura_cov_fuentes()
dev.off()
Mostrar(Figura_cov_fuentes())

#5. COMPARACIÓN POR ALTITUD, VERTIENTE Y FASE ENSO (estación como unidad). Vertiente: primer dígito del código IDEAM
#(2 = Magdalena-Cauca, 5 = Pacífico), columna auxiliar; la comparación principal es por altitud
Est_sit       <- sort(unique(Est_par))                                 #Mismo orden de las columnas de D_est
Banda_sit     <- Banda_alt[match(Est_sit, Pares_curvas$STATION_ID)]
Vertiente_sit <- ifelse(substr(Est_sit, 1, 1) == "5", "Pacífico", "Cauca")
table(Vertiente = Vertiente_sit, Altitud = Banda_sit, useNA = "ifany")
#Nivel anual de cada fuente y sesgo anual de cada estación (mm/día; media de los años de la estación)
Nivel_est_sit <- colMeans(Por_estacion(Mm_est_fino[Dia_centro, Con_par]))
Nivel_sat_sit <- colMeans(Por_estacion(Mm_sat_fino[Dia_centro, Idx_sat_par[Con_par]]))
d_sit <- colMeans(D_est[Dia_centro, ])
max(abs(d_sit - (Nivel_est_sit - Nivel_sat_sit)))
Resumen_grupo <- function(g) {
  data.frame(Estaciones = as.vector(table(g)), Estacion_mm_dia = tapply(Nivel_est_sit, g, mean),
             Satelite_mm_dia = tapply(Nivel_sat_sit, g, mean), d = tapply(d_sit, g, mean),
             EE = tapply(d_sit, g, function(x) sd(x) / sqrt(length(x))),
             d_pct_del_satelite = 100 * tapply(d_sit, g, mean) / tapply(Nivel_sat_sit, g, mean),
             Pct_satelite_mayor = 100 * tapply(d_sit < 0, g, mean))
}
round(Resumen_grupo(Banda_sit), 2)
round(Resumen_grupo(Vertiente_sit), 2)
round(Resumen_grupo(interaction(Vertiente_sit, Banda_sit, drop = TRUE, sep = " ")), 2)

#d(t) por banda de altitud ± 2 EE e intervalos sin cero (la banda < 500 m, con 3 estaciones, es descriptiva)
Dif_banda <- lapply(levels(Banda_sit), function(b) {
  M <- D_est[, which(Banda_sit == b), drop = FALSE]
  list(m = rowMeans(M), ee = apply(M, 1, sd) / sqrt(ncol(M)), N = ncol(M))
})
names(Dif_banda) <- levels(Banda_sit)
Int_banda <- lapply(Dif_banda, function(B) Intervalos_sin_cero(B$m, B$ee))
lapply(Int_banda, function(I) I$Tabla)
sapply(Int_banda, function(I) c(Pct_año = I$Pct_año, Pct_dias_15_350 = I$Pct_interior))

#Fase ENSO: la fase es un atributo del año, así que la réplica es el año. Anomalía de cada par: diferencia relativa (%)
#de su nivel frente a la media de su estación (y del píxel en esos años); anomalía del año: media de sus estaciones
Anom_rel_est <- 100 * (colMeans(Mm_est_fino[Dia_centro, Con_par]) / Nivel_est_sit[Grupo_est] - 1)
Anom_rel_sat <- 100 * (colMeans(Mm_sat_fino[Dia_centro, Idx_sat_par[Con_par]]) / Nivel_sat_sit[Grupo_est] - 1)
Anom_año <- data.frame(Año = sort(unique(Año_par)), Estacion = as.vector(tapply(Anom_rel_est, Año_par, mean)),
                       Satelite = as.vector(tapply(Anom_rel_sat, Año_par, mean)),
                       Estaciones = as.vector(table(Año_par)))
Anom_año <- cbind(Anom_año, ONI_año[match(Anom_año$Año, ONI_año$Año), c("ONI", "Fase", "Transicion")])
Anom_año[, c("Estacion", "Satelite", "ONI")] <- round(Anom_año[, c("Estacion", "Satelite", "ONI")], 2)
Anom_año
Tabla_fase <- function(Usar) {
  Fase_u <- factor(Anom_año$Fase[Usar], levels = Fases)
  Por_f  <- function(x) rbind(Media = tapply(x[Usar], Fase_u, mean), EE = tapply(x[Usar], Fase_u, function(v) sd(v) / sqrt(length(v))))
  M <- rbind(Por_f(Anom_año$Estacion), Por_f(Anom_año$Satelite), Años = table(Fase_u))
  rownames(M) <- c("Estación media (%)", "Estación EE", "Satélite media (%)", "Satélite EE", "Años")
  round(M, 2)
}
Tabla_fase(rep(TRUE, nrow(Anom_año)))
#Sensibilidad: sin los años de transición (1983 y 1998, neutros por la media anual)
Tabla_fase(!Anom_año$Transicion)
#Concordancia año a año de las anomalías regionales y relación con el ONI (23 años)
round(c(Cor_estacion_satelite = cor(Anom_año$Estacion, Anom_año$Satelite), Cor_estacion_ONI = cor(Anom_año$Estacion, Anom_año$ONI),
        Cor_satelite_ONI = cor(Anom_año$Satelite, Anom_año$ONI)), 3)
#Anomalía funcional (mm/día) por fase: media ± 2 EE entre los años de la fase
Anom_cur_est <- Mm_est_fino[, Con_par] - Por_estacion(Mm_est_fino[, Con_par])[, Grupo_est]
Anom_cur_sat <- Mm_sat_fino[, Idx_sat_par[Con_par]] - Por_estacion(Mm_sat_fino[, Idx_sat_par[Con_par]])[, Grupo_est]
Anom_fase_cur <- function(A) {
  A_año <- sapply(split(seq_along(Año_par), Año_par), function(j) rowMeans(A[, j, drop = FALSE]))   #731 x 23
  R <- lapply(Fases, function(f) {
    M <- A_año[, Anom_año$Fase == f, drop = FALSE]
    list(m = rowMeans(M), ee = apply(M, 1, sd) / sqrt(ncol(M)))
  })
  names(R) <- Fases
  R
}
Anom_fase_est3 <- Anom_fase_cur(Anom_cur_est)
Anom_fase_sat3 <- Anom_fase_cur(Anom_cur_sat)
#Fracción del año en que la banda de la anomalía no contiene el cero, por fuente y fase
round(rbind(Estacion = sapply(Anom_fase_est3, function(A) 100 * mean(A$m - 2 * A$ee > 0 | A$m + 2 * A$ee < 0)),
            Satelite = sapply(Anom_fase_sat3, function(A) 100 * mean(A$m - 2 * A$ee > 0 | A$m + 2 * A$ee < 0))), 1)

Figura_altitud_dif <- function() {
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2, 1))
  for (b in names(Dif_banda)) {
    B <- Dif_banda[[b]]
    plot(t_fino, B$m, type = "n", ylim = range(c(B$m - 2 * B$ee, B$m + 2 * B$ee, 0)), xlab = "Día del año",
         ylab = "Estación - satélite (mm/día)", main = paste0("d(t) ± 2 EE: ", b, " (", B$N, " estaciones)"))
    rect(c(0, 350), par("usr")[3], c(15, 365), par("usr")[4], col = "grey92", border = NA)
    polygon(c(t_fino, rev(t_fino)), c(B$m - 2 * B$ee, rev(B$m + 2 * B$ee)), col = adjustcolor("darkred", 0.2), border = NA)
    lines(t_fino, B$m, lwd = 3, col = "darkred")
    abline(h = 0, lty = 3)
  }
}
Figura_enso_fuentes <- function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
  Col_fase <- c("blue", "grey40", "red")
  Lim <- range(sapply(c(Anom_fase_est3, Anom_fase_sat3), function(A) range(c(A$m - 2 * A$ee, A$m + 2 * A$ee))))
  for (Fuente in c("Estaciones", "Satélite")) {
    A3 <- if (Fuente == "Estaciones") Anom_fase_est3 else Anom_fase_sat3
    plot(t_fino, A3[[1]]$m, type = "n", ylim = Lim, xlab = "Día del año", ylab = "Anomalía ± 2 EE entre años (mm/día)",
         main = paste(Fuente, "por fase ENSO"))
    for (k in 1:3) {
      polygon(c(t_fino, rev(t_fino)), c(A3[[k]]$m - 2 * A3[[k]]$ee, rev(A3[[k]]$m + 2 * A3[[k]]$ee)),
              col = adjustcolor(Col_fase[k], 0.15), border = NA)
      lines(t_fino, A3[[k]]$m, lwd = 3, col = Col_fase[k])
    }
    abline(h = 0, lty = 3)
    legend("topleft", Fases, col = Col_fase, lwd = 3, bty = "n", cex = 0.8)
  }
}
pdf(Ruta("Figuras_proyecto/Fig_Diferencia_altitud.pdf"), width = 11, height = 8)
Figura_altitud_dif()
dev.off()
Mostrar(Figura_altitud_dif())
pdf(Ruta("Figuras_proyecto/Fig_ENSO_fuentes.pdf"), width = 11, height = 5)
Figura_enso_fuentes()
dev.off()
Mostrar(Figura_enso_fuentes())

#6. SENSIBILIDAD: SIN LOS PARES CON EL SATÉLITE EN EL TECHO DE 200 mm O CON PICO AISLADO CONSERVADO (banda de d(t),
#estación como unidad)
Techo_par <- Curvas_sat$Techo[Idx_sat_par[Con_par]]
Pico_par  <- Pares_curvas$Pico[Con_par]
c(Pares_techo = sum(Techo_par), Pares_pico = sum(Pico_par), En_ambos = sum(Techo_par & Pico_par))
#Estaciones y bandas de altitud de los pares con techo
table(Estacion = Est_par[Techo_par], Banda = Banda_sit[match(Est_par[Techo_par], Est_sit)], useNA = "ifany")
Sensibilidad <- function(Usar) {
  M  <- sapply(split(which(Usar), Est_par[Usar]), function(j) rowMeans(D_par[, j, drop = FALSE]))   #731 x estaciones
  m  <- rowMeans(M)
  ee <- apply(M, 1, sd) / sqrt(ncol(M))
  I  <- Intervalos_sin_cero(m, ee)
  Banda_M <- Banda_sit[match(colnames(M), Est_sit)]
  c(Pares = sum(Usar), Estaciones = ncol(M), d_mm_dia = mean(m[Dia_centro]), Pct_año = I$Pct_año,
    Pct_dias_15_350 = I$Pct_interior, Intervalos_interiores = sum(!I$Tabla$Borde),
    d_menos_500m = mean(colMeans(M[Dia_centro, which(Banda_M == "<500 m"), drop = FALSE])))
}
Tabla_sensibilidad <- sapply(list(Todos = rep(TRUE, length(Techo_par)), Sin_techo = !Techo_par, Sin_pico = !Pico_par,
                                  Sin_techo_ni_pico = !Techo_par & !Pico_par), Sensibilidad)
round(Tabla_sensibilidad, 3)
#Control: con todos los pares se reproduce la banda del punto 2
c(d = Tabla_sensibilidad["d_mm_dia", "Todos"] - mean(d_media[Dia_centro]), Pct = Tabla_sensibilidad["Pct_año", "Todos"] - Int_est$Pct_año)

#7. ANOVA FUNCIONAL ENTRE FUENTES (ESTACIÓN VS SATÉLITE): medidas repetidas con la estación como bloque (media de sus
#años), X_ij(t) = mu(t) + beta_i(t) + tau_j(t) + e_ij(t), i = 1..n_s, j = estación, satélite. No se usa el ANOVA no
#pareado: su suma de cuadrados dentro de cada fuente incluye la variabilidad entre sitios. Descomposición, pruebas L2 y
#tipo F con corrección de sesgo (Zhang, 2013), F puntual globalizada, permutación y curva crítica según
#Script_anova_funcional.R, con este diseño
n_s       <- ncol(D_est)
Est_anova <- Por_estacion(Mm_est_fino[Dia_centro, Con_par])
Sat_anova <- Por_estacion(Mm_sat_fino[Dia_centro, Idx_sat_par[Con_par]])
D_s       <- D_est[Dia_centro, ]
t_eval    <- t_fino[Dia_centro]
M_anova   <- length(t_eval)
#Pesos del trapecio sobre la malla de 365 centros de día: la integral del año es sum(w_anova * .)
w_anova <- rep(1, M_anova); w_anova[c(1, M_anova)] <- 0.5
#Control: la diferencia de las curvas medias por estación es D_s
max(abs(D_s - (Est_anova - Sat_anova)))

#Descomposición puntual del ANOVA de dos vías con una observación por celda (la interacción es el residual)
d_bar        <- rowMeans(D_s)                                  #Diferencia media entre estaciones (mm/día)
gran_anova   <- rowMeans(cbind(Est_anova, Sat_anova))
media_bloque <- (Est_anova + Sat_anova) / 2
SS_fuente_t  <- n_s * d_bar^2 / 2
SS_bloque_t  <- 2 * rowSums((media_bloque - gran_anova)^2)
SS_resid_t   <- 0.5 * rowSums((D_s - d_bar)^2)
SS_total_t   <- rowSums((Est_anova - gran_anova)^2) + rowSums((Sat_anova - gran_anova)^2)
max(abs(SS_total_t - (SS_fuente_t + SS_bloque_t + SS_resid_t)))
#F(t) del efecto fuente: con dos fuentes equivale al cuadrado del estadístico t pareado, (d_bar(t)/EE(t))^2
F_t   <- (SS_fuente_t / (2 - 1)) / (SS_resid_t / (n_s - 1))
p_F_t <- pf(F_t, 1, n_s - 1, lower.tail = FALSE)
c(Estaciones = n_s, Dias = M_anova, Pct_dias_p_nominal_005 = round(100 * mean(p_F_t < 0.05), 1))
#Verificación: F(t) = (d_bar(t) / EE(t))^2
EE_D_t <- apply(D_s, 1, sd) / sqrt(n_s)
max(abs(F_t - (d_bar / EE_D_t)^2))
#Control: días en que la F puntual coincide con la banda del punto 2
c(t_critico = qt(0.975, n_s - 1),
  Pct_dias_coinciden_con_banda = 100 * mean((p_F_t < 0.05) == (d_media - 2 * EE_d > 0 | d_media + 2 * EE_d < 0)[Dia_centro]))

#Sumas de cuadrados integradas en el año y tabla de ANOVA funcional. Las unidades son (mm/día)^2 x día = mm^2/día
SS_fuente_n <- sum(w_anova * SS_fuente_t)
SS_bloque_n <- sum(w_anova * SS_bloque_t)
SS_resid_n  <- sum(w_anova * SS_resid_t)
SS_total_n  <- sum(w_anova * SS_total_t)
#Escala de las curvas de diferencia: SSR_n = 2 SS_fuente_n y SSE_n = 2 SS_resid_n (el cociente F no cambia)
SSR_n <- n_s * sum(w_anova * d_bar^2)
SSE_n <- 2 * SS_resid_n
max(abs(c(SSR_n - 2 * SS_fuente_n, SSE_n - 2 * SS_resid_n)))
F_n <- (SSR_n / (2 - 1)) / (SSE_n / (n_s - 1))
#Fracción de la variabilidad integrada que aporta cada término
round(c(Fuente = SS_fuente_n / SS_total_n, Bloque = SS_bloque_n / SS_total_n,
        Residual = SS_resid_n / SS_total_n), 4)

#Covarianza residual de las curvas de diferencia y corrección de sesgo de las trazas (Zhang, 2013)
Resid_D   <- D_s - d_bar
G_comb    <- Resid_D %*% t(Resid_D) / (n_s - 1)
DG        <- G_comb * w_anova                                  #diag(w_anova) %*% G
tr_gamma  <- sum(diag(DG))
tr_gamma2 <- sum(DG * t(DG))
N_gl      <- n_s - 1
tr_gamma2_c <- N_gl^2 / ((N_gl - 1) * (N_gl + 2)) * (tr_gamma2 - tr_gamma^2 / N_gl)
tr2_gamma_c <- N_gl * (N_gl + 1) / ((N_gl - 1) * (N_gl + 2)) * (tr_gamma^2 - 2 * tr_gamma2 / (N_gl + 1))
beta_hat <- tr_gamma2_c / tr_gamma
d_hat    <- tr2_gamma_c / tr_gamma2_c
c(beta = beta_hat, d = d_hat)
#Porcentaje de la variabilidad integrada de D_s que es nivel anual de cada estación (dirección dominante si d_hat ≈ 1)
Nivel_D <- colSums(w_anova * Resid_D) / sum(w_anova)               #Promedio en el año de cada curva centrada
round(100 * sum(w_anova) * sum(Nivel_D^2) / sum(w_anova * Resid_D^2), 1)
#Prueba L2: SSR_n ~ beta chi2_{(k-1)d} y prueba tipo F: F_n ~ F_{(k-1)d, (n_s-1)d}
p_L2      <- pchisq(SSR_n / beta_hat, df = (2 - 1) * d_hat, lower.tail = FALSE)
p_F_aprox <- pf(F_n, (2 - 1) * d_hat, (n_s - 1) * d_hat, lower.tail = FALSE)

#Estadístico F puntual globalizado: T_GPF = integral de F(t). La distribución nula depende de la correlación rho(s,t)
T_GPF <- sum(w_anova * F_t)
s_t   <- sqrt(diag(G_comb))
Con_var <- s_t > 0
Rho   <- matrix(0, M_anova, M_anova)
Rho[Con_var, Con_var] <- G_comb[Con_var, Con_var] / outer(s_t[Con_var], s_t[Con_var])
DR      <- Rho * w_anova
tr_rho  <- sum(diag(DR))
tr_rho2 <- sum(DR * t(DR))
beta_rho <- tr_rho2 / ((2 - 1) * tr_rho)
d_rho    <- (2 - 1) * tr_rho^2 / tr_rho2
p_GPF_aprox <- pchisq(T_GPF / beta_rho, df = d_rho, lower.tail = FALSE)

#Permutación pareada: bajo H0 se cambia el signo de D_s(t) en cada estación, con las B permutaciones en forma
#matricial (con píxeles compartidos, el valor p puede ser algo liberal)
set.seed(2026)
B_anova <- 5000
Signos    <- matrix(sample(c(-1, 1), n_s * B_anova, replace = TRUE), n_s, B_anova)
Dbar_null <- D_s %*% Signos / n_s                              #365 x B
SSD_t     <- rowSums(D_s^2)
SSR_t_null <- n_s * Dbar_null^2
SSE_t_null <- SSD_t - n_s * Dbar_null^2
F_t_null   <- SSR_t_null / (SSE_t_null / (n_s - 1))
SSR_n_null <- colSums(w_anova * SSR_t_null)
SSE_n_null <- colSums(w_anova * SSE_t_null)
F_n_null   <- SSR_n_null / (SSE_n_null / (n_s - 1))
GPF_null   <- colSums(w_anova * F_t_null)
Fmax_null  <- apply(F_t_null, 2, max)
min(SSE_t_null)
p_perm <- function(nulos, obs) (1 + sum(nulos >= obs - 1e-12)) / (length(nulos) + 1)
F_max       <- max(F_t)
p_F_perm    <- p_perm(F_n_null, F_n)
p_GPF_perm  <- p_perm(GPF_null, T_GPF)
p_Fmax_perm <- p_perm(Fmax_null, F_max)
c(p_F_perm = p_F_perm, p_GPF_perm = p_GPF_perm, p_Fmax_perm = p_Fmax_perm)
#Con B permutaciones el menor valor p posible es 1/(B + 1): un valor igual a ese se reporta como p < 0,001
c(p_minimo_posible = 1 / (B_anova + 1))

#Curva crítica puntual y valor crítico máximo al 5 %: el máximo controla el error de tipo I en todo el año
alfa_anova <- 0.05
c_puntual  <- apply(F_t_null, 1, quantile, probs = 1 - alfa_anova)
c_maximo   <- quantile(Fmax_null, 1 - alfa_anova)
p_puntual  <- (1 + rowSums(F_t_null >= F_t - 1e-12)) / (B_anova + 1)

#Días con diferencia significativa, sin corrección (puntual) y con el criterio del máximo (todo el año)
Intervalos_F <- function(Indicador) {
  r   <- rle(Indicador)
  Fin <- cumsum(r$lengths); Ini <- Fin - r$lengths + 1
  data.frame(Desde = fecha_de(t_eval[Ini[r$values]]), Hasta = fecha_de(t_eval[Fin[r$values]]),
             Dias = r$lengths[r$values])
}
Intervalos_F(p_puntual < alfa_anova)
Intervalos_F(F_t > c_maximo)
round(c(Pct_año_puntual = 100 * mean(p_puntual < alfa_anova),
        Pct_año_maximo  = 100 * mean(F_t > c_maximo)), 1)

#Sensibilidad a los píxeles compartidos: mismo signo para las estaciones de cada píxel (no recoge el efecto común del año)
Pixel_est  <- Estaciones_curvas$Pixel[match(Est_sit, Estaciones_curvas$STATION_ID)]
Pixeles_u  <- sort(unique(Pixel_est))
Signos_px  <- matrix(sample(c(-1, 1), length(Pixeles_u) * B_anova, replace = TRUE), length(Pixeles_u), B_anova)
Dbar_px    <- D_s %*% Signos_px[match(Pixel_est, Pixeles_u), ] / n_s
F_t_px     <- n_s * Dbar_px^2 / ((SSD_t - n_s * Dbar_px^2) / (n_s - 1))
Fmax_px    <- apply(F_t_px, 2, max)
c_maximo_px <- quantile(Fmax_px, 1 - alfa_anova)
c(Pixeles = length(Pixeles_u), c_maximo_estacion = unname(c_maximo), c_maximo_pixel = unname(c_maximo_px),
  p_Fmax_pixel = p_perm(Fmax_px, F_max))
Intervalos_F(F_t > c_maximo_px)

#Sensibilidad a los bordes del año: pruebas globales solo con los días 15-350 (mismas permutaciones)
Interior_a <- t_eval >= 15 & t_eval <= 350
w_int      <- rep(1, sum(Interior_a)); w_int[c(1, sum(Interior_a))] <- 0.5
c(Dia_F_max = t_eval[which.max(F_t)], F_max = F_max, F_max_dias_15_350 = max(F_t[Interior_a]),
  Pct_GPF_en_bordes = 100 * sum((w_anova * F_t)[!Interior_a]) / T_GPF)
F_n_int       <- sum(w_int * n_s * d_bar[Interior_a]^2) / (sum(w_int * rowSums(Resid_D[Interior_a, ]^2)) / (n_s - 1))
F_n_null_int  <- colSums(w_int * SSR_t_null[Interior_a, ]) / (colSums(w_int * SSE_t_null[Interior_a, ]) / (n_s - 1))
GPF_int       <- sum(w_int * F_t[Interior_a])
GPF_null_int  <- colSums(w_int * F_t_null[Interior_a, ])
Fmax_null_int <- apply(F_t_null[Interior_a, ], 2, max)
c_maximo_int  <- quantile(Fmax_null_int, 1 - alfa_anova)
p_int <- c(p_F_perm = p_perm(F_n_null_int, F_n_int), p_GPF_perm = p_perm(GPF_null_int, GPF_int),
           p_Fmax_perm = p_perm(Fmax_null_int, max(F_t[Interior_a])))
round(c(p_int, c_maximo_dias_15_350 = unname(c_maximo_int)), 4)
Intervalos_F(F_t > c_maximo_int & Interior_a)

#Tabla de ANOVA funcional (sumas de cuadrados integradas) y resumen de las pruebas globales
Tabla_ANOVA_fuentes <- data.frame(
  Fuente         = c("Fuente (estación vs satélite)", "Estación (bloque)", "Fuente x estación (residual)", "Total"),
  Suma_Cuadrados = c(SS_fuente_n, SS_bloque_n, SS_resid_n, SS_total_n),
  gl             = c(1, n_s - 1, n_s - 1, 2 * n_s - 1),
  Cuadrado_Medio = c(SS_fuente_n, SS_bloque_n / (n_s - 1), SS_resid_n / (n_s - 1), NA),
  F              = c(F_n, NA, NA, NA),
  p_aprox        = c(p_F_aprox, NA, NA, NA),
  p_perm         = c(p_F_perm, NA, NA, NA))
Tabla_ANOVA_fuentes
Tabla_pruebas_anova <- data.frame(
  Prueba      = c("L2 (chi-cuadrado corregida)", "Tipo F corregida", "F_n por permutación",
                  "GPF aproximada", "GPF por permutación", "F_max por permutación"),
  Estadistico = c(SSR_n, F_n, F_n, T_GPF, T_GPF, F_max),
  gl1         = c((2 - 1) * d_hat, (2 - 1) * d_hat, NA, d_rho, NA, NA),
  gl2         = c(NA, (n_s - 1) * d_hat, NA, NA, NA, NA),
  Valor_p     = c(p_L2, p_F_aprox, p_F_perm, p_GPF_aprox, p_GPF_perm, p_Fmax_perm))
Tabla_pruebas_anova

#Figura: descomposición puntual, F(t) con las curvas críticas y distribución nula de F_n por permutación
Figura_anova_fuentes <- function() {
  par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2, 1))
  matplot(t_eval, cbind(SS_fuente_t, SS_resid_t), type = "l", lty = 1, lwd = 2, col = c("darkred", "steelblue"),
          xlab = "Día del año", ylab = expression("Suma de cuadrados (mm"^2*"/día)"), main = "Descomposición puntual")
  legend("topright", c("Fuente", "Fuente x estación"), col = c("darkred", "steelblue"), lty = 1, lwd = 2,
         bty = "n", cex = 0.8)
  plot(t_eval, F_t, type = "l", lwd = 2, col = "darkred", xlab = "Día del año", ylab = "F(t)",
       main = "F(t) del efecto fuente", ylim = c(0, max(F_t, c_puntual, c_maximo) * 1.1))
  lines(t_eval, c_puntual, lty = 3, lwd = 2, col = "blue")
  abline(h = c_maximo, lty = 2)
  legend("topleft", c("Observado", "Crítico puntual 0,05", "Crítico máximo 0,05"),
         col = c("darkred", "blue", "black"), lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
  hist(F_n_null, breaks = 50, freq = FALSE, col = "grey85", border = "white",
       xlim = c(0, max(F_n, quantile(F_n_null, 0.995))), xlab = expression(F[n]), main = "Permutación de F")
  curve(df(x, (2 - 1) * d_hat, (n_s - 1) * d_hat), add = TRUE, col = "blue", lwd = 2)
  abline(v = F_n, col = "darkred", lwd = 2)
  legend("topright", c("Aproximación F", "Observado"), col = c("blue", "darkred"), lwd = 2, bty = "n", cex = 0.8)
}
pdf(Ruta("Figuras_proyecto/Fig_ANOVA_fuentes.pdf"), width = 14, height = 4.8)
Figura_anova_fuentes()
dev.off()
Mostrar(Figura_anova_fuentes())

#p_Fmax_perm: diferencia en algún día con control en todo el año; p_F_perm: versión integrada; p_GPF_perm: promedio de
#F(t). Los intervalos con F(t) > c_maximo tienen control en todo el año; los de p_puntual < 0,05 no

#RESUMEN DEL ÍTEM 3 (cifras del informe; mm/día salvo que se indique)
Resumen_item3 <- c(
  Estaciones_emparejadas = sum(Est_pixel), Pixeles = length(Est_pixel), Pixeles_compartidos = sum(Est_pixel > 1),
  Estaciones_en_pixel_compartido = sum(Est_pixel[Est_pixel > 1]), Pares_estacion_año = sum(Con_par),
  Cor_sesgo_años_impares_pares = cor(Sesgo_impar[Ambos], Sesgo_par[Ambos]),
  d_anual_por_estacion = mean(d_media[Dia_centro]), d_anual_por_par = mean(d_media_pares[Dia_centro]),
  Pct_año_sin_cero_estacion = Int_est$Pct_año, Pct_dias_15_350_sin_cero_estacion = Int_est$Pct_interior,
  Pct_año_sin_cero_pares_referencia = Int_pares$Pct_año, Meses_que_coinciden_datos_y_curva = sum(Tabla_mes$Coinciden),
  r_diagonal_mediana_entre_sitios = median(Cruzado$Entre_sitios$Diag[Dia_centro]),
  r_diagonal_mediana_año_a_año = median(Cruzado$Año_a_año$Diag[Dia_centro]),
  Razon_DE_est_sat_entre_sitios = mean(sqrt(diag(Cruzado$Entre_sitios$Cov_Y) / diag(Cruzado$Entre_sitios$Cov_X))[Dia_centro]),
  Razon_DE_est_sat_año_a_año = mean(sqrt(diag(Cruzado$Año_a_año$Cov_Y) / diag(Cruzado$Año_a_año$Cov_X))[Dia_centro]),
  setNames(Resumen_grupo(Banda_sit)$d, paste("d_anual", levels(Banda_sit))),
  Anomalia_Niña_estacion_pct = mean(Anom_año$Estacion[Anom_año$Fase == "Niña"]),
  Anomalia_Niño_estacion_pct = mean(Anom_año$Estacion[Anom_año$Fase == "Niño"]),
  Anomalia_Niña_satelite_pct = mean(Anom_año$Satelite[Anom_año$Fase == "Niña"]),
  Anomalia_Niño_satelite_pct = mean(Anom_año$Satelite[Anom_año$Fase == "Niño"]),
  Cor_anomalia_estacion_ONI = cor(Anom_año$Estacion, Anom_año$ONI), Cor_anomalia_satelite_ONI = cor(Anom_año$Satelite, Anom_año$ONI),
  ANOVA_pct_SS_fuente = 100 * SS_fuente_n / SS_total_n, ANOVA_pct_SS_bloque = 100 * SS_bloque_n / SS_total_n,
  ANOVA_p_L2 = p_L2, ANOVA_p_F_aprox = p_F_aprox, ANOVA_p_F_perm = p_F_perm, ANOVA_p_GPF_perm = p_GPF_perm,
  ANOVA_p_Fmax_perm = p_Fmax_perm, ANOVA_p_minimo_posible = 1 / (B_anova + 1),
  ANOVA_Pct_año_criterio_maximo = 100 * mean(F_t > c_maximo), ANOVA_Pct_año_criterio_maximo_pixel = 100 * mean(F_t > c_maximo_px),
  ANOVA_p_F_perm_dias_15_350 = unname(p_int["p_F_perm"]), ANOVA_p_GPF_perm_dias_15_350 = unname(p_int["p_GPF_perm"]),
  ANOVA_p_Fmax_perm_dias_15_350 = unname(p_int["p_Fmax_perm"]))
data.frame(Valor = round(Resumen_item3, 4))
#Intervalos con control en todo el año (criterio del máximo) y, para el informe, los del interior (días 15-350)
Intervalos_F(F_t > c_maximo)
Intervalos_F(F_t > c_maximo_int & Interior_a)

#Lectura del sesgo del satélite (mismos sitios y años):
#1. Momento: el satélite llega tarde a las dos temporadas (la estación tiene más lluvia en ene-mar y sep, el satélite
#   en may y oct), con balance anual casi nulo por estación (puntos 2 y 3)
#2. Nivel: el satélite aplana el gradiente orográfico; sobreestima en el fondo del valle, subestima sobre 2000 m y en
#   la llanura del Pacífico registra la mitad de la lluvia (punto 5)
#3. Escala: el satélite varía menos que las estaciones entre sitios y entre años (punto 4)
#4. Interanual: el satélite reproduce qué lugares llueven más, pero casi no sigue la variación año a año ni el ENSO
#   (puntos 4 y 5)

#MODELO DE REGRESIÓN FUNCIONAL (ítem 4)
#1. TABLA DE MODELADO Y LÍNEA BASE (CURVA SATELITAL SIN CORREGIR)
#Unidad: el par estación-año; entradas: curva del píxel en el mismo año, altitud de la estación y ONI del año, sin años
#anteriores ni parámetros por estación (sirve para estaciones nuevas). Partición: estimación 1981-1995, elección
#1996-2000 y evaluación 2001-2003
Sel_mod <- Con_par & Con_alt & !(Pares_curvas$STATION_ID %in% Reservadas)
Modelo  <- data.frame(STATION_ID = Pares_curvas$STATION_ID[Sel_mod], YEAR = Pares_curvas$YEAR[Sel_mod],
                      Fila_est = which(Sel_mod), Fila_sat = Idx_sat_par[Sel_mod], Pixel = Pixel_par[Sel_mod],
                      Altitud = Altitud_par[Sel_mod])
Modelo$Bloque <- cut(Modelo$YEAR, c(1980, 1995, 2000, 2003), labels = c("Entrenamiento", "Validación", "Prueba"))
Modelo$ONI    <- ONI_año$ONI[match(Modelo$YEAR, ONI_año$Año)]
Modelo$Fase   <- factor(ONI_año$Fase[match(Modelo$YEAR, ONI_año$Año)], levels = Fases)
#Verificación: cada fila tiene la curva de su estación, la curva de su píxel en el mismo año, altitud y ONI
stopifnot(!anyNA(Modelo), all(Curvas_sat$Pixel[Modelo$Fila_sat] == Modelo$Pixel),
          all(Curvas_sat$Año[Modelo$Fila_sat] == Modelo$YEAR), !anyDuplicated(Modelo[, c("STATION_ID", "YEAR")]))
#Control: todos los píxeles del Valle tienen curva en los 23 años
stopifnot(nrow(Curvas_sat) == length(Pixeles_valle) * length(Años_sat), !anyNA(Coef_z_sat))

#Cobertura de la rejilla estación x año
Años_est <- split(Modelo$YEAR, Modelo$STATION_ID)
Cob_mod  <- data.frame(Años = lengths(Años_est), Primero = sapply(Años_est, min), Ultimo = sapply(Años_est, max))
Cob_mod$Faltan_en_medio <- Cob_mod$Ultimo - Cob_mod$Primero + 1 - Cob_mod$Años
round(c(Estaciones = nrow(Cob_mod), Pares = nrow(Modelo), Con_los_23_años = sum(Cob_mod$Años == 23),
        Empiezan_despues_de_1981 = sum(Cob_mod$Primero > 1981), Terminan_antes_de_2003 = sum(Cob_mod$Ultimo < 2003),
        Con_años_faltantes_en_medio = sum(Cob_mod$Faltan_en_medio > 0),
        Años_faltantes_en_medio = sum(Cob_mod$Faltan_en_medio),
        Pct_rejilla_ocupada = 100 * nrow(Modelo) / (length(Años_sat) * nrow(Cob_mod)),
        Pct_pares_con_año_anterior = 100 * mean(paste(Modelo$STATION_ID, Modelo$YEAR - 1) %in%
                                                  paste(Modelo$STATION_ID, Modelo$YEAR))), 1)
quantile(Cob_mod$Años, c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
#(un modelo con el año anterior de la estación perdería los pares sin año anterior y no serviría para estaciones nuevas)

#Pares y estaciones por bloque, y cuántos son de estaciones que no tienen ningún año de entrenamiento
Est_ent <- unique(Modelo$STATION_ID[Modelo$Bloque == "Entrenamiento"])
Tabla_bloques <- t(sapply(split(Modelo, Modelo$Bloque), function(d)
  c(Pares = nrow(d), Estaciones = length(unique(d$STATION_ID)),
    Estaciones_sin_entrenamiento = length(setdiff(d$STATION_ID, Est_ent)),
    Pares_sin_entrenamiento = sum(!(d$STATION_ID %in% Est_ent)))))
Tabla_bloques
table(Modelo$Bloque, Modelo$Fase)

#Curvas en mm/día en los 365 centros de día (365 x N, una columna por fila de Modelo), recortadas en 0
Obs_mm <- pmax(Mm_est_fino[Dia_centro, Modelo$Fila_est], 0)     #Curva de la estación
Sat_mm <- pmax(Mm_sat_fino[Dia_centro, Modelo$Fila_sat], 0)     #Curva del píxel sin corregir: la línea base
#Segunda referencia, sin satélite: curva media de las estaciones en entrenamiento
Ent_mod <- Modelo$Bloque == "Entrenamiento"
Clim_mm <- matrix(rowMeans(Obs_mm[, Ent_mod]), nrow(Obs_mm), ncol(Obs_mm))

#Función de métricas por par (4.4), con la regla del punto medio en los 365 centros de día:
#  L2_rel = raíz(integral (pred - obs)^2 / integral obs^2); KGE = 1 - raíz((r - 1)^2 + (alfa - 1)^2 + (beta - 1)^2)
#  (Gupta et al., 2009), con r = correlación, alfa = DE(pred) / DE(obs) y beta = media(pred) / media(obs);
#  Sesgo_pct = 100 (acumulado pred - acumulado obs) / acumulado obs
Metricas_par <- function(Pred, Obs) {
  r    <- sapply(seq_len(ncol(Obs)), function(i) cor(Pred[, i], Obs[, i]))
  alfa <- apply(Pred, 2, sd) / apply(Obs, 2, sd)
  beta <- colMeans(Pred) / colMeans(Obs)
  data.frame(L2_rel = sqrt(colSums((Pred - Obs)^2) / colSums(Obs^2)),
             KGE = 1 - sqrt((r - 1)^2 + (alfa - 1)^2 + (beta - 1)^2), r = r, alfa = alfa, beta = beta,
             Sesgo_pct = 100 * (colSums(Pred) - colSums(Obs)) / colSums(Obs))
}
Resumen_metricas <- function(Met, Grupo) {
  t(sapply(split(Met, Grupo, drop = TRUE), function(d)
    c(Pares = nrow(d), L2_mediana = median(d$L2_rel), L2_media = mean(d$L2_rel), KGE_mediana = median(d$KGE),
      r_mediana = median(d$r), alfa_mediana = median(d$alfa), beta_mediana = median(d$beta),
      Sesgo_mediana_pct = median(d$Sesgo_pct), Sesgo_abs_mediana_pct = median(abs(d$Sesgo_pct)))))
}
#Control: una predicción perfecta da error 0, KGE 1 y sesgo 0
stopifnot(all(abs(as.matrix(Metricas_par(Obs_mm[, 1:5], Obs_mm[, 1:5])) - rep(c(0, 1, 1, 1, 1, 0), each = 5)) < 1e-12))

Base_sat  <- Metricas_par(Sat_mm, Obs_mm)
Base_clim <- Metricas_par(Clim_mm, Obs_mm)
#La prueba (2001-2003) se resume solo en la evaluación (sección 5)
Sin_prueba <- Modelo$Bloque != "Prueba"
Bloque_sp  <- droplevels(Modelo$Bloque[Sin_prueba])
#Línea base: curva satelital sin corregir
round(Resumen_metricas(Base_sat[Sin_prueba, ], Bloque_sp), 3)
#Referencia sin satélite: curva media de entrenamiento
round(Resumen_metricas(Base_clim[Sin_prueba, ], Bloque_sp), 3)
#Línea base por fase ENSO y por altitud (las bandas del ítem 3)
round(Resumen_metricas(Base_sat[Sin_prueba, ], paste(Bloque_sp, Modelo$Fase[Sin_prueba])), 3)
Banda_mod <- cut(Modelo$Altitud, c(-Inf, 500, 1200, 2000, Inf), right = FALSE,
                 labels = c("<500 m", "500-1200 m", "1200-2000 m", ">=2000 m"))
round(Resumen_metricas(Base_sat[Ent_mod, ], Banda_mod[Ent_mod]), 3)

#Exportación para Python: tabla de pares de la regresión (la red usa los mismos pares)
write.csv(data.frame(station_id = Modelo$STATION_ID, year = Modelo$YEAR,
                     bloque = c("entrenamiento", "validacion", "prueba")[as.integer(Modelo$Bloque)],
                     pixel = Modelo$Pixel, altitud = Modelo$Altitud, oni = Modelo$ONI, fase = Modelo$Fase),
          Ruta("Datos_preparados/modelo_pares.csv"), row.names = FALSE)

#2. MODELO CONCURRENTE (ejercicio 3 de la clase de regresión funcional), en escala raíz:
#z_est_i(t) = b0(t) + b1(t) z_sat_i(t) + e_i(t); la versión ampliada agrega altitud y ONI (constantes en t). Cada
#b_j(t) = theta(t)' c_j se estima con PENSISE = suma_i integral [z_i - zhat_i]^2 dt + lambda suma_j integral [D2 b_j]^2 dt,
#cuyo mínimo resuelve (A + lambda R) c = d (integrales por Simpson en t_fino)

#Curvas en la escala raíz sobre la rejilla fina: 731 x N, una columna por fila de Modelo
Phi_z_est <- Base_bspline(t_fino, 0, 365, K_z_est)
Phi_z_sat <- Base_bspline(t_fino, 0, 365, K_z_sat)
Z_est <- Phi_z_est %*% Coef_z_est[, Modelo$Fila_est]
Z_sat <- Phi_z_sat %*% Coef_z_sat[, Modelo$Fila_sat]
#Pesos de Simpson: h/3 (1, 4, 2, 4, ..., 2, 4, 1), con un número par de intervalos
h_fino <- t_fino[2] - t_fino[1]
w_fino <- h_fino / 3 * c(1, rep(c(4, 2), (length(t_fino) - 3) / 2), 4, 1)
stopifnot(length(w_fino) == length(t_fino), abs(sum(w_fino) - 365) < 1e-9)

#Covariables escalares estandarizadas con entrenamiento (el ONI, con los 15 años de entrenamiento)
Alt_e   <- (Modelo$Altitud - mean(Modelo$Altitud[Ent_mod])) / sd(Modelo$Altitud[Ent_mod])
ONI_ent <- ONI_año$ONI[ONI_año$Año <= 1995]
ONI_e   <- (Modelo$ONI - mean(ONI_ent)) / sd(ONI_ent)
Fila_cte <- function(x) matrix(x, length(t_fino), length(x), byrow = TRUE)     #Covariable escalar como función constante
X_conc   <- list(Intercepto = Fila_cte(rep(1, nrow(Modelo))), Satelite = Z_sat)
X_conc_e <- c(X_conc, list(Altitud = Fila_cte(Alt_e), ONI = Fila_cte(ONI_e)))

#Base de los coeficientes: la misma de la respuesta (K_z_est funciones), con penalización de la segunda derivada
Base_beta  <- create.bspline.basis(c(0, 365), K_z_est, 4)
Theta_beta <- Phi_z_est                                   #theta(t) en la rejilla fina
R_beta     <- eval.penalty(Base_beta, int2Lfd(2))         #integral D2 theta D2 theta'
stopifnot(max(abs(eval.basis(t_fino, Base_beta) - Theta_beta)) < 1e-10)

#Funciones del sistema lineal (A, d) con las curvas de Filas y del ajuste (matriz K x p, columna j = coeficientes de b_j)
Sistema_concurrente <- function(Xs, Filas) {
  p <- length(Xs); K <- ncol(Theta_beta)
  A <- matrix(0, p * K, p * K); d <- numeric(p * K)
  for (j in 1:p) {
    ij <- (j - 1) * K + 1:K
    d[ij] <- crossprod(Theta_beta, w_fino * rowSums(Xs[[j]][, Filas, drop = FALSE] * Z_est[, Filas, drop = FALSE]))
    for (k in j:p) {
      ik <- (k - 1) * K + 1:K
      A[ij, ik] <- crossprod(Theta_beta, (w_fino * rowSums(Xs[[j]][, Filas, drop = FALSE] * Xs[[k]][, Filas, drop = FALSE])) * Theta_beta)
      A[ik, ij] <- t(A[ij, ik])
    }
  }
  list(A = A, d = d)
}
Ajustar_concurrente <- function(Sis, lambda, Nombres) {
  p <- length(Nombres)
  matrix(solve(Sis$A + lambda * kronecker(diag(p), R_beta), Sis$d), ncol(Theta_beta), p, dimnames = list(NULL, Nombres))
}
#Predicción en la rejilla fina: zhat_i(t) = suma_j x_ij(t) b_j(t)
Predecir_concurrente <- function(Xs, Cb) Reduce(`+`, lapply(seq_along(Xs), function(j) Xs[[j]] * as.vector(Theta_beta %*% Cb[, j])))
#Error cuadrático integrado medio por par, en la escala raíz
SISE_z <- function(Zh, Filas) mean(colSums(w_fino * (Z_est[, Filas] - Zh[, Filas])^2))

Sis_conc   <- Sistema_concurrente(X_conc, Ent_mod)
Sis_conc_e <- Sistema_concurrente(X_conc_e, Ent_mod)
#Control de la cuadratura: el bloque del intercepto es N veces la matriz de Gram de la base
max(abs(Sis_conc$A[1:K_z_est, 1:K_z_est] / sum(Ent_mod) - eval.penalty(Base_beta, 0)))

#Selección de lambda con validación (1996-2000); no se deja una curva fuera porque las curvas de un mismo año
#comparten el clima de ese año
Val_mod <- Modelo$Bloque == "Validación"
Lambdas_conc <- 10^seq(-2, 12, by = 1)
Tabla_lambda_conc <- data.frame(log10_lambda = log10(Lambdas_conc), t(sapply(Lambdas_conc, function(l) {
  Z1 <- Predecir_concurrente(X_conc, Ajustar_concurrente(Sis_conc, l, names(X_conc)))
  Z2 <- Predecir_concurrente(X_conc_e, Ajustar_concurrente(Sis_conc_e, l, names(X_conc_e)))
  c(Ent_satelite = SISE_z(Z1, Ent_mod), Val_satelite = SISE_z(Z1, Val_mod),
    Ent_ampliado = SISE_z(Z2, Ent_mod), Val_ampliado = SISE_z(Z2, Val_mod))
})))
round(Tabla_lambda_conc, 2)
Lambda_conc   <- Lambdas_conc[which.min(Tabla_lambda_conc$Val_satelite)]
Lambda_conc_e <- Lambdas_conc[which.min(Tabla_lambda_conc$Val_ampliado)]
c(log10_lambda_satelite = log10(Lambda_conc), log10_lambda_ampliado = log10(Lambda_conc_e))

Coef_conc   <- Ajustar_concurrente(Sis_conc, Lambda_conc, names(X_conc))
Coef_conc_e <- Ajustar_concurrente(Sis_conc_e, Lambda_conc_e, names(X_conc_e))
Zh_conc     <- Predecir_concurrente(X_conc, Coef_conc)
Zh_conc_e   <- Predecir_concurrente(X_conc_e, Coef_conc_e)
#Error integrado en validación (escala raíz) frente a las dos referencias del punto 1
Media_z_ent <- rowMeans(Z_est[, Ent_mod])
round(c(Satelite_sin_corregir = SISE_z(Z_sat, Val_mod), Media_de_entrenamiento = SISE_z(Fila_cte(rep(0, nrow(Modelo))) + Media_z_ent, Val_mod),
        Concurrente = SISE_z(Zh_conc, Val_mod), Concurrente_ampliado = SISE_z(Zh_conc_e, Val_mod)), 1)

#Verificación frente a fRegress (mismo modelo, base y lambda): coeficientes y criterio PENSISE de cada solución con la
#misma cuadratura (la solución exacta no puede dar un criterio mayor)
PENSISE_conc <- function(Cb) sum(colSums(w_fino * (Z_est[, Ent_mod] - Predecir_concurrente(X_conc_e, Cb)[, Ent_mod])^2)) +
  Lambda_conc_e * sum(diag(t(Cb) %*% R_beta %*% Cb))
Par_beta <- fdPar(Base_beta, int2Lfd(2), Lambda_conc_e)
Fit_freg <- fRegress(fd(Coef_z_est[, Modelo$Fila_est[Ent_mod]], Base_beta),
                     list(rep(1, sum(Ent_mod)), fd(Coef_z_sat[, Modelo$Fila_sat[Ent_mod]], create.bspline.basis(c(0, 365), K_z_sat, 4)),
                          Alt_e[Ent_mod], ONI_e[Ent_mod]),
                     list(Par_beta, Par_beta, Par_beta, Par_beta))
Coef_freg <- sapply(Fit_freg$betaestlist, function(b) b$fd$coefs)
Dif_freg  <- apply(abs(Theta_beta %*% (Coef_conc_e - Coef_freg)), 2, max)
signif(rbind(Diferencia_maxima = Dif_freg, Relativa_al_rango = Dif_freg / apply(Theta_beta %*% Coef_conc_e, 2, function(v) diff(range(v)))), 2)
c(PENSISE_a_mano = PENSISE_conc(Coef_conc_e), PENSISE_fRegress = PENSISE_conc(Coef_freg))
stopifnot(PENSISE_conc(Coef_conc_e) <= PENSISE_conc(Coef_freg))

#Función de R2(t): varianza entre curvas que explica el modelo en cada día (referencia: media de entrenamiento; en
#validación puede ser negativa)
R2_t <- function(Zh, Filas) 1 - rowSums((Z_est[, Filas] - Zh[, Filas])^2) / rowSums((Z_est[, Filas] - Media_z_ent)^2)
R2_conc <- cbind(Ent_satelite = R2_t(Zh_conc, Ent_mod), Val_satelite = R2_t(Zh_conc, Val_mod),
                 Ent_ampliado = R2_t(Zh_conc_e, Ent_mod), Val_ampliado = R2_t(Zh_conc_e, Val_mod))[Dia_centro, ]
round(rbind(Minimo = apply(R2_conc, 2, min), Media = colMeans(R2_conc), Maximo = apply(R2_conc, 2, max)), 3)
#R2(t) medio por mes: en qué parte del año el satélite explica mejor la estación
Mes_dia <- as.integer(format(as.Date(t365, origin = "2001-01-01"), "%m"))
round(apply(R2_conc, 2, function(v) tapply(v, Mes_dia, mean)), 2)

#Bandas de los coeficientes por bootstrap de grupos (1000 remuestreos, percentiles 2,5 y 97,5): estaciones para b0, b1 y
#la altitud, años para el ONI. El sistema (A, d) es la suma de los sistemas de cada grupo
Boot_concurrente <- function(Grupo, B = 1000) {
  Sis_g <- lapply(split(which(Ent_mod), Grupo), function(ii) Sistema_concurrente(X_conc_e, ii))
  A_g <- sapply(Sis_g, function(s) as.vector(s$A)); d_g <- sapply(Sis_g, function(s) s$d)
  Pen <- Lambda_conc_e * kronecker(diag(length(X_conc_e)), R_beta)
  Rep <- replicate(B, {
    n <- tabulate(sample(length(Sis_g), replace = TRUE), length(Sis_g))
    as.vector(Theta_beta[Dia_centro, ] %*% matrix(solve(matrix(A_g %*% n, nrow(Pen)) + Pen, d_g %*% n), ncol(Theta_beta)))
  })
  array(apply(Rep, 1, quantile, c(0.025, 0.975)), c(2, length(Dia_centro), length(X_conc_e)), list(c("Inf", "Sup"), NULL, names(X_conc_e)))
}
set.seed(2026)
Banda_conc_est <- Boot_concurrente(Modelo$STATION_ID[Ent_mod])
Banda_conc_año <- Boot_concurrente(Modelo$YEAR[Ent_mod])
Banda_conc <- Banda_conc_est; Banda_conc[, , "ONI"] <- Banda_conc_año[, , "ONI"]
Beta_conc_e <- Theta_beta[Dia_centro, ] %*% Coef_conc_e                 #365 x 4: coeficientes en los centros de día
Sin_cero <- function(Bd) Bd["Inf", , ] > 0 | Bd["Sup", , ] < 0
round(rbind(Minimo = apply(Beta_conc_e, 2, min), Maximo = apply(Beta_conc_e, 2, max),
            Pct_año_banda_sin_cero = 100 * colMeans(Sin_cero(Banda_conc)),
            Pct_remuestreando_estaciones = 100 * colMeans(Sin_cero(Banda_conc_est)),
            Pct_remuestreando_años = 100 * colMeans(Sin_cero(Banda_conc_año))), 3)
#Porcentaje del año en que la banda de b1(t) queda por debajo de 1 (el modelo atenúa la señal del satélite)
round(100 * mean(Banda_conc["Sup", , "Satelite"] < 1), 1)

#Regreso a mm: lluvia = max(zhat, 0)^2. Como E[z^2] = E[z]^2 + Var(z), el cuadrado pierde lluvia en el total; se compara
#con la variante que suma la varianza residual de entrenamiento (Ampliado_mas_varianza)
A_mm <- function(Zh, Sumar_varianza = FALSE) {
  P <- pmax(Zh[Dia_centro, ], 0)^2
  if (Sumar_varianza) P <- P + rowMeans((Z_est[Dia_centro, Ent_mod] - Zh[Dia_centro, Ent_mod])^2)
  P
}
Met_conc   <- Metricas_par(A_mm(Zh_conc), Obs_mm)
Met_conc_e <- Metricas_par(A_mm(Zh_conc_e), Obs_mm)
Fila_val <- function(Met) Resumen_metricas(Met[Val_mod, ], rep("Validación", sum(Val_mod)))[1, ]
Tabla_conc <- rbind(Satelite_sin_corregir = Fila_val(Base_sat), Media_de_entrenamiento = Fila_val(Base_clim),
                    Concurrente = Fila_val(Met_conc), Concurrente_ampliado = Fila_val(Met_conc_e),
                    Ampliado_mas_varianza = Fila_val(Metricas_par(A_mm(Zh_conc_e, TRUE), Obs_mm)),
                    Curva_observada_en_raiz = Fila_val(Metricas_par(A_mm(Z_est), Obs_mm)))
round(Tabla_conc, 3)
#Porcentaje de valores recortados en 0 al volver a mm
round(100 * c(Concurrente = mean(Zh_conc[Dia_centro, ] < 0), Concurrente_ampliado = mean(Zh_conc_e[Dia_centro, ] < 0)), 4)
#Por fase ENSO en validación: satélite sin corregir frente al concurrente ampliado
round(Resumen_metricas(Base_sat[Val_mod, ], Modelo$Fase[Val_mod]), 3)
round(Resumen_metricas(Met_conc_e[Val_mod, ], Modelo$Fase[Val_mod]), 3)

#Figura: coeficientes con su banda, R2(t) y selección de lambda
Inicio_mes <- cumsum(c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))
Eje_meses  <- function() axis(1, at = Inicio_mes + 15, labels = c("E", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"), tick = FALSE)
Figura_concurrente <- function() {
  par(mfrow = c(2, 3), mar = c(3.2, 4.2, 2.4, 1))
  Titulos <- c("Intercepto b0(t)", "Satélite b1(t)", "Altitud g1(t), por desviación estándar", "ONI g2(t), por desviación estándar")
  for (j in 1:4) {
    plot(t365, Beta_conc_e[, j], type = "n", ylim = range(Banda_conc[, , j], 0, if (j == 2) 1), xaxt = "n", xlab = "",
         ylab = "Raíz de mm/día", main = Titulos[j], cex.main = 0.95)
    polygon(c(t365, rev(t365)), c(Banda_conc["Inf", , j], rev(Banda_conc["Sup", , j])), col = "grey85", border = NA)
    lines(t365, Beta_conc_e[, j], lwd = 2)
    abline(h = 0, lty = 3); if (j == 2) abline(h = 1, lty = 2, col = "red")
    abline(v = Inicio_mes, col = "grey92"); Eje_meses()
  }
  matplot(t365, R2_conc, type = "l", lty = c(2, 1, 2, 1), col = c("grey50", "grey50", "black", "black"), lwd = c(1, 2, 1, 2),
          ylim = c(min(R2_conc), max(R2_conc) + 0.14), xaxt = "n", xlab = "", ylab = expression(R^2*(t)), main = "Varianza explicada",
          cex.main = 0.95)
  abline(h = 0, lty = 3); Eje_meses()
  legend("topleft", c("Solo satélite, validación", "Ampliado, validación", "Entrenamiento (trazos)"), col = c("grey50", "black", "black"),
         lty = c(1, 1, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
  matplot(Tabla_lambda_conc$log10_lambda, Tabla_lambda_conc[, c("Val_satelite", "Val_ampliado")], type = "b", pch = 19, lty = 1,
          col = c("grey50", "black"), xlab = "", ylab = "Error integrado en validación", main = expression("Selección de " * lambda), cex.main = 0.95)
  mtext(expression(log[10](lambda)), 1, line = 2.1, cex = 0.7)
  abline(v = log10(c(Lambda_conc, Lambda_conc_e)), lty = 2, col = c("grey50", "black"))
}
pdf(Ruta("Figuras_proyecto/Fig_Concurrente.pdf"), width = 11, height = 6.5)
Figura_concurrente()
dev.off()
Mostrar(Figura_concurrente())

#3. MODELO FUNCIONAL-FUNCIONAL (función a función) CON COVARIABLES ESCALARES
#z_est_i(t) = b0(t) + integral b1(s, t) z_sat_i(s) ds + g1(t) altitud_i + g2(t) ONI_i + e_i(t)
#Con z_sat_i = phi' c_i, z_est_i = theta' d_i y b1(s, t) = phi(s)' B theta(t): dhat_i' = u_i' Th, u_i = (1, altitud_i,
#ONI_i, W_sat c_i). El mínimo de PENSISE (rugosidad de b1 en s y en t; b0, g1 y g2 sin penalizar) resuelve
#[W_est (x) (U'U + lambda_s P_s) + lambda_t R_est (x) P_t] vec(Th) = vec(U'D W_est), P_s = diag(0, R_sat), P_t = diag(0, W_sat)

#Función de la matriz de Gram por Gauss-Legendre con 4 nodos por intervalo (exacta para el producto de dos B-splines cúbicas)
Gram_bspline <- function(K) {
  Cortes <- unique(Nodos_bspline(0, 365, K, 4))
  Nodo <- c(-sqrt(3 / 7 + 2 / 7 * sqrt(6 / 5)), -sqrt(3 / 7 - 2 / 7 * sqrt(6 / 5)), sqrt(3 / 7 - 2 / 7 * sqrt(6 / 5)), sqrt(3 / 7 + 2 / 7 * sqrt(6 / 5)))
  Peso <- c((18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36)
  a <- head(Cortes, -1); b <- Cortes[-1]
  t_q <- as.vector(outer(Nodo, (b - a) / 2) + rep((a + b) / 2, each = 4))    #4 nodos por intervalo, llevados de [-1, 1] a [a, b]
  w_q <- as.vector(outer(Peso, (b - a) / 2))
  Phi_q <- Base_bspline(t_q, 0, 365, K)
  crossprod(Phi_q, w_q * Phi_q)
}
W_est <- Gram_bspline(K_z_est)                       #K_est x K_est
W_sat <- Gram_bspline(K_z_sat)                       #K_sat x K_sat
Base_sat_z <- create.bspline.basis(c(0, 365), K_z_sat, 4)
#Verificación: igual a la matriz exacta de fda y la base suma 1 (suma de todas las entradas = 365)
c(Dif_est = max(abs(W_est - eval.penalty(Base_beta, 0))), Dif_sat = max(abs(W_sat - eval.penalty(Base_sat_z, 0))),
  Suma_est = sum(W_est), Suma_sat = sum(W_sat))
stopifnot(max(abs(W_est - eval.penalty(Base_beta, 0))) < 1e-10, max(abs(W_sat - eval.penalty(Base_sat_z, 0))) < 1e-10)
R_est <- R_beta                                      #Penalización en t: integral D2 theta D2 theta'
R_sat <- eval.penalty(Base_sat_z, int2Lfd(2))        #Penalización en s: integral D2 phi D2 phi'

#Matrices del modelo: una fila por par de Modelo
D_mod <- t(Coef_z_est[, Modelo$Fila_est])            #N x K_est: coeficientes de la curva de la estación (raíz)
C_mod <- t(Coef_z_sat[, Modelo$Fila_sat])            #N x K_sat: coeficientes de la curva del píxel (raíz)
U_ff  <- cbind(Intercepto = 1, Altitud = Alt_e, ONI = ONI_e, C_mod %*% W_sat)
Fil_b1 <- 4:ncol(U_ff)                               #Filas de Th que forman la matriz B de b1(s, t)
P_s <- matrix(0, ncol(U_ff), ncol(U_ff)); P_s[Fil_b1, Fil_b1] <- R_sat
P_t <- matrix(0, ncol(U_ff), ncol(U_ff)); P_t[Fil_b1, Fil_b1] <- W_sat

#Funciones de ajuste del modelo y de la suma de errores cuadráticos integrados con la matriz de Gram
Ajustar_ff <- function(lambda_s, lambda_t, Filas = Ent_mod, U = U_ff, Ps = P_s, Pt = P_t) {
  UtU <- crossprod(U[Filas, ]); UtD <- crossprod(U[Filas, ], D_mod[Filas, ])
  matrix(solve(kronecker(W_est, UtU + lambda_s * Ps) + lambda_t * kronecker(R_est, Pt), as.vector(UtD %*% W_est)),
         ncol(U), ncol(D_mod))
}
SSE_d <- function(Dh, Filas) { E <- D_mod[Filas, , drop = FALSE] - Dh[Filas, , drop = FALSE]; sum((E %*% W_est) * E) }

#Control sin penalizar: número de parámetros y error en entrenamiento y en validación
Dh_sin <- U_ff %*% Ajustar_ff(0, 0)
round(c(Parametros = length(Ajustar_ff(0, 0)), Ent_sin_penalizar = SSE_d(Dh_sin, Ent_mod) / sum(Ent_mod),
        Val_sin_penalizar = SSE_d(Dh_sin, Val_mod) / sum(Val_mod), Val_concurrente_ampliado = SISE_z(Zh_conc_e, Val_mod)), 1)

#Selección de lambda_s y lambda_t con validación
Lambdas_ff <- 10^(0:14)
Rejilla_ff <- expand.grid(lambda_s = Lambdas_ff, lambda_t = Lambdas_ff)
Rejilla_ff[, c("Ent", "Val")] <- t(sapply(seq_len(nrow(Rejilla_ff)), function(g) {
  Dh <- U_ff %*% Ajustar_ff(Rejilla_ff$lambda_s[g], Rejilla_ff$lambda_t[g])
  c(SSE_d(Dh, Ent_mod) / sum(Ent_mod), SSE_d(Dh, Val_mod) / sum(Val_mod))
}))
round(xtabs(Val ~ log10(lambda_s) + log10(lambda_t), Rejilla_ff), 1)
Mejor_ff <- Rejilla_ff[which.min(Rejilla_ff$Val), ]
Mejor_ff

#Contraste: validación cruzada dejando un año fuera dentro de entrenamiento (no interviene en la elección)
Años_ent_mod <- sort(unique(Modelo$YEAR[Ent_mod]))
Rejilla_ff$VC_año <- sapply(seq_len(nrow(Rejilla_ff)), function(g) sum(sapply(Años_ent_mod, function(a) {
  Fuera <- Ent_mod & Modelo$YEAR == a
  SSE_d(U_ff %*% Ajustar_ff(Rejilla_ff$lambda_s[g], Rejilla_ff$lambda_t[g], Ent_mod & !Fuera), Fuera)
})) / sum(Ent_mod))
VC_año_conc_e <- sum(sapply(Años_ent_mod, function(a) {
  Fuera <- Ent_mod & Modelo$YEAR == a
  Zh <- Predecir_concurrente(X_conc_e, Ajustar_concurrente(Sistema_concurrente(X_conc_e, Ent_mod & !Fuera), Lambda_conc_e, names(X_conc_e)))
  sum(colSums(w_fino * (Z_est[, Fuera] - Zh[, Fuera])^2))
})) / sum(Ent_mod)
Mejor_vc_ff <- Rejilla_ff[which.min(Rejilla_ff$VC_año), ]
Contraste_lambda <- rbind(FF_elegido_con_validacion = unlist(Rejilla_ff[which.min(Rejilla_ff$Val), c("lambda_s", "lambda_t", "Val", "VC_año")]),
                          FF_elegido_con_VC_por_año = unlist(Mejor_vc_ff[, c("lambda_s", "lambda_t", "Val", "VC_año")]),
                          Concurrente_ampliado      = c(Lambda_conc_e, NA, SISE_z(Zh_conc_e, Val_mod), VC_año_conc_e))
Contraste_lambda[, 1:2] <- log10(Contraste_lambda[, 1:2])
colnames(Contraste_lambda)[1:2] <- c("log10_lambda_s", "log10_lambda_t")
round(Contraste_lambda, 1)

Lambda_s_ff <- Mejor_ff$lambda_s; Lambda_t_ff <- Mejor_ff$lambda_t
Th_ff <- Ajustar_ff(Lambda_s_ff, Lambda_t_ff)        #19 x 15
Dh_ff <- U_ff %*% Th_ff                              #N x 15: coeficientes predichos
Zh_ff <- Phi_z_est %*% t(Dh_ff)                      #731 x N: curvas predichas en la rejilla fina
#Control: el error con la matriz de Gram coincide con el de Simpson
c(Gram = SSE_d(Dh_ff, Val_mod) / sum(Val_mod), Simpson = SISE_z(Zh_ff, Val_mod))

#Verificación frente a linmod (función de clase), sin altitud ni ONI (linmod no admite covariables escalares), con las
#mismas bases y lambdas
U_lin <- U_ff[, -(2:3)]
Th_lin <- Ajustar_ff(Lambda_s_ff, Lambda_t_ff, U = U_lin, Ps = P_s[-(2:3), -(2:3)], Pt = P_t[-(2:3), -(2:3)])
Fit_linmod <- linmod(fd(t(C_mod[Ent_mod, ]), Base_sat_z), fd(t(D_mod[Ent_mod, ]), Base_beta),
                     list(fdPar(Base_beta, int2Lfd(2), 0),
                          bifdPar(bifd(matrix(0, K_z_sat, K_z_est), Base_sat_z, Base_beta), int2Lfd(2), int2Lfd(2), Lambda_s_ff, Lambda_t_ff)))
Dif_linmod <- c(Coef_b1 = max(abs(Fit_linmod$beta1estbifd$coefs - Th_lin[-1, ])),
                b0 = max(abs(eval.fd(t_fino, Fit_linmod$beta0estfd) - Phi_z_est %*% Th_lin[1, ])),
                Prediccion = max(abs(eval.fd(t_fino, Fit_linmod$yhatfdobj) - Phi_z_est %*% t(U_lin[Ent_mod, ] %*% Th_lin))))
signif(Dif_linmod, 2)
stopifnot(Dif_linmod < 1e-3)

#Resultados en validación: error integrado (raíz), R2(t) y métricas en mm
round(c(Satelite_sin_corregir = SISE_z(Z_sat, Val_mod), Concurrente_ampliado = SISE_z(Zh_conc_e, Val_mod),
        Funcional_funcional = SISE_z(Zh_ff, Val_mod)), 1)
R2_ff <- cbind(Ent = R2_t(Zh_ff, Ent_mod), Val = R2_t(Zh_ff, Val_mod))[Dia_centro, ]
round(rbind(Minimo = apply(R2_ff, 2, min), Media = colMeans(R2_ff), Maximo = apply(R2_ff, 2, max)), 3)
round(cbind(FF = tapply(R2_ff[, "Val"], Mes_dia, mean), Concurrente_ampliado = tapply(R2_conc[, "Val_ampliado"], Mes_dia, mean)), 2)
Met_ff <- Metricas_par(A_mm(Zh_ff), Obs_mm)
Tabla_ff <- rbind(Tabla_conc[c("Satelite_sin_corregir", "Media_de_entrenamiento", "Concurrente_ampliado"), ],
                  Funcional_funcional = Fila_val(Met_ff))
round(Tabla_ff, 3)
round(100 * mean(Zh_ff[Dia_centro, ] < 0), 4)        #Porcentaje de valores recortados en 0 al volver a mm
round(Resumen_metricas(Met_ff[Val_mod, ], Modelo$Fase[Val_mod]), 3)
#Nivel medio del residuo en validación (raíz de mm/día): positivo = el modelo predice menos lluvia de la observada
round(c(Concurrente_ampliado = mean((Z_est - Zh_conc_e)[Dia_centro, Val_mod]), Funcional_funcional = mean((Z_est - Zh_ff)[Dia_centro, Val_mod])), 3)

#Superficie b1(s, t) en los centros de día
Phi_sat_dia <- Phi_z_sat[Dia_centro, ]; Theta_dia <- Phi_z_est[Dia_centro, ]
Beta1_ff <- Phi_sat_dia %*% Th_ff[Fil_b1, ] %*% t(Theta_dia)          #365 x 365: filas = s (satélite), columnas = t (estación)
Desfase  <- abs(outer(t365, t365, "-"))                              #|s - t| en días
#Proporción de la integral de |b1| a menos de 15 y de 30 días de s = t (8 % y 16 % si b1 fuera constante)
round(100 * c(A_15_dias = sum(abs(Beta1_ff)[Desfase <= 15]), A_30_dias = sum(abs(Beta1_ff)[Desfase <= 30])) / sum(abs(Beta1_ff)), 1)
#Efecto total del satélite sobre el día t: integral en s de b1(s, t), comparable con b1(t) del modelo concurrente
Efecto_total_ff <- colSums(Beta1_ff)
round(rbind(FF_integral_en_s = tapply(Efecto_total_ff, Mes_dia, mean), Concurrente_b1 = tapply(Beta_conc_e[, "Satelite"], Mes_dia, mean)), 2)

#Efecto de borde: |b1(s, t)| en las esquinas de enero y diciembre frente al resto, y curvas con ceros al inicio y al
#final del año
Borde_dia <- t365 < 31 | t365 > 334
Esquina   <- outer(Borde_dia, Borde_dia, "&")
round(c(Media_abs_b1_esquinas = mean(abs(Beta1_ff)[Esquina]), Media_abs_b1_resto = mean(abs(Beta1_ff)[!Esquina]),
        Pct_curvas_con_cero_primeros_10_dias = 100 * mean(apply(Obs_mm[1:10, ], 2, min) == 0),
        Pct_curvas_con_cero_ultimos_10_dias = 100 * mean(apply(Obs_mm[356:365, ], 2, min) == 0)), 4)

Figura_ff <- function() {
  par(mfrow = c(2, 2), mar = c(3.4, 4.2, 2.4, 1))
  Sub <- seq(3, 365, by = 5)
  filled_lim <- max(abs(Beta1_ff))
  image(t365[Sub], t365[Sub], Beta1_ff[Sub, Sub], col = hcl.colors(41, "Blue-Red 3"), zlim = c(-filled_lim, filled_lim), xaxt = "n", yaxt = "n",
        xlab = "", ylab = "t: día de la estación", main = expression(hat(b)[1](s, t) * ": rojo positivo, azul negativo"), cex.main = 0.95)
  contour(t365[Sub], t365[Sub], Beta1_ff[Sub, Sub], add = TRUE, nlevels = 8, col = "grey30", labcex = 0.55)
  abline(0, 1, lty = 2); Eje_meses(); axis(2, at = Inicio_mes + 15, labels = c("E", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"), tick = FALSE, las = 1)
  mtext("s: día del satélite", 1, line = 2.1, cex = 0.75)
  plot(t365, Efecto_total_ff, type = "l", lwd = 2, ylim = range(0, 1, Efecto_total_ff, Beta_conc_e[, "Satelite"]), xaxt = "n", xlab = "",
       ylab = "Raíz de mm/día", main = "Efecto total del satélite sobre el día t", cex.main = 0.95)
  lines(t365, Beta_conc_e[, "Satelite"], lwd = 2, col = "grey50", lty = 2)
  abline(h = 1, lty = 3, col = "red"); Eje_meses()
  legend("bottomleft", c("Integral en s de b1(s, t)", "b1(t) del concurrente"), col = c("black", "grey50"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
  Esc_ff <- Theta_dia %*% t(Th_ff[2:3, ])
  matplot(t365, cbind(Esc_ff, Beta_conc_e[, c("Altitud", "ONI")]), type = "l", lty = c(1, 1, 2, 2), lwd = 2, col = c("darkgreen", "purple"),
          xaxt = "n", xlab = "", ylab = "Raíz de mm/día por desviación estándar", main = "Altitud y ONI", cex.main = 0.95)
  abline(h = 0, lty = 3); Eje_meses()
  legend("bottomleft", c("Altitud", "ONI", "Concurrente (trazos)"), col = c("darkgreen", "purple", "black"), lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.8)
  matplot(t365, cbind(R2_ff[, "Val"], R2_conc[, "Val_ampliado"], R2_ff[, "Ent"]), type = "l", lty = c(1, 1, 2), lwd = c(2, 2, 1),
          col = c("black", "grey50", "black"), ylim = c(min(R2_ff, R2_conc), max(R2_ff, R2_conc) + 0.14), xaxt = "n", xlab = "",
          ylab = expression(R^2*(t)), main = "Varianza explicada", cex.main = 0.95)
  abline(h = 0, lty = 3); Eje_meses()
  legend("topleft", c("Función a función, validación", "Concurrente ampliado, validación", "Función a función, entrenamiento"),
         col = c("black", "grey50", "black"), lty = c(1, 1, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
}
pdf(Ruta("Figuras_proyecto/Fig_Funcional_funcional.pdf"), width = 10, height = 8)
Figura_ff()
dev.off()
Mostrar(Figura_ff())

#4. FPCA DE LAS CURVAS DE ESTACIÓN Y REGRESIÓN DE RANGO REDUCIDO (equivale a la cabeza lineal de la red)
#Decodificador de la red: zhat(t) = mu(t) + suma_m s_m psi_m(t), con el FPCA de entrenamiento congelado

#FPCA en la base B-spline con las curvas de entrenamiento de Modelo: S W b = rho b, resuelto como (L S L') q = rho q con
#W = L'L (Cholesky) y b = L^-1 q, que cumple b' W b = 1
Media_d   <- colMeans(D_mod[Ent_mod, ])                                  #Coeficientes de mu(t)
D_centro  <- sweep(D_mod, 2, Media_d)                                    #Curvas centradas con la media de entrenamiento
S_coef    <- crossprod(D_centro[Ent_mod, ]) / sum(Ent_mod)               #Covarianza de los coeficientes (divisor N, como pca.fd)
L_est     <- chol(W_est)                                                 #W_est = L' L
Propio    <- eigen(L_est %*% S_coef %*% t(L_est), symmetric = TRUE)
V_fpca    <- solve(L_est, Propio$vectors)                                #Columna m = coeficientes de psi_m
#Signo fijo para que el resultado sea reproducible: cada función propia con integral positiva
V_fpca    <- sweep(V_fpca, 2, sign(colSums(W_est %*% V_fpca)), "*")
Var_fpca  <- Propio$values
Tabla_fpca <- data.frame(Componente = 1:8, Valor_propio = Var_fpca[1:8], Pct_varianza = 100 * Var_fpca[1:8] / sum(Var_fpca),
                         Pct_acumulado = 100 * cumsum(Var_fpca)[1:8] / sum(Var_fpca))
round(Tabla_fpca, 2)
#Verificación: ortonormalidad en L2 y coincidencia con pca.fd (valores propios y funciones propias, salvo el signo)
Pca_fda <- pca.fd(fd(t(D_mod[Ent_mod, ]), Base_beta), nharm = 6, centerfns = TRUE)
Dif_fpca <- c(Ortonormalidad = max(abs(t(V_fpca) %*% W_est %*% V_fpca - diag(ncol(V_fpca)))),
              Valores_propios = max(abs(Pca_fda$values[1:6] - Var_fpca[1:6])),
              Funciones_propias = max(abs(abs(Phi_z_est %*% Pca_fda$harmonics$coefs) - abs(Phi_z_est %*% V_fpca[, 1:6]))),
              Varianza_explicada = max(abs(Pca_fda$varprop - Var_fpca[1:6] / sum(Var_fpca))))
signif(Dif_fpca, 2)
#(diferencias del orden de la precisión de inprod)
stopifnot(Dif_fpca < 1e-3)

#Puntajes de todas las curvas con el FPCA de entrenamiento: s_im = integral [z_i(t) - mu(t)] psi_m(t) dt = (d_i - mu)' W b_m
Puntajes <- D_centro %*% W_est %*% V_fpca                                #N x 15
#Control: en entrenamiento los puntajes tienen media 0, varianza rho_m y no están correlacionados
max(abs(crossprod(Puntajes[Ent_mod, ]) / sum(Ent_mod) - diag(Var_fpca)))

#Regresión de rango reducido: s_im = u_i' h_m + e_im con las entradas del modelo funcional-funcional, es decir,
#b1(s, t) = suma_m [phi(s)' h_m] psi_m(t). Por la ortonormalidad de psi_m, el mínimo de PENSISE es (U'U + lambda_s P_s) h_m = U' s_m
Ajustar_rr <- function(M, lambda_s, Filas = Ent_mod) solve(crossprod(U_ff[Filas, ]) + lambda_s * P_s, crossprod(U_ff[Filas, ], Puntajes[Filas, 1:M, drop = FALSE]))
Predecir_rr <- function(H) sweep(U_ff %*% H %*% t(V_fpca[, seq_len(ncol(H)), drop = FALSE]), 2, Media_d, "+")   #N x K_est: coeficientes predichos

#Elección de M (4 a 6) y de lambda_s con validación; techo del decodificador: error de reconstruir la curva observada
#con sus propios puntajes
Lambdas_rr <- 10^(0:14)
Val_rr <- sapply(4:6, function(M) sapply(Lambdas_rr, function(l) SSE_d(Predecir_rr(Ajustar_rr(M, l)), Val_mod) / sum(Val_mod)))
dimnames(Val_rr) <- list(log10_lambda_s = log10(Lambdas_rr), M = 4:6)
round(Val_rr, 1)
Techo_rr <- sapply(4:6, function(M) SSE_d(sweep(Puntajes[, 1:M] %*% t(V_fpca[, 1:M]), 2, Media_d, "+"), Val_mod) / sum(Val_mod))
Tabla_M <- data.frame(M = 4:6, Pct_varianza_ent = Tabla_fpca$Pct_acumulado[4:6],
                      Pct_varianza_val = 100 * (1 - Techo_rr / (SSE_d(matrix(Media_d, nrow(D_mod), ncol(D_mod), byrow = TRUE), Val_mod) / sum(Val_mod))),
                      Techo_decodificador_val = Techo_rr, log10_lambda_s = log10(Lambdas_rr[apply(Val_rr, 2, which.min)]),
                      Error_val = apply(Val_rr, 2, min), Error_val_sin_penalizar = sapply(4:6, function(M) SSE_d(Predecir_rr(Ajustar_rr(M, 0)), Val_mod) / sum(Val_mod)))
round(Tabla_M, 2)
M_fpca      <- Tabla_M$M[which.min(Tabla_M$Error_val)]
Lambda_s_rr <- 10^Tabla_M$log10_lambda_s[Tabla_M$M == M_fpca]
c(M = M_fpca, log10_lambda_s = log10(Lambda_s_rr), Pct_varianza = Tabla_fpca$Pct_acumulado[M_fpca])

H_rr  <- Ajustar_rr(M_fpca, Lambda_s_rr)            #19 x M: fila 1 intercepto, 2 altitud, 3 ONI, 4-19 satélite
Dh_rr <- Predecir_rr(H_rr)
Zh_rr <- Phi_z_est %*% t(Dh_rr)
Met_rr <- Metricas_par(A_mm(Zh_rr), Obs_mm)
#Comparación en validación: error integrado (raíz) y métricas en mm
round(c(Satelite_sin_corregir = SISE_z(Z_sat, Val_mod), Concurrente_ampliado = SISE_z(Zh_conc_e, Val_mod),
        Funcional_funcional = SISE_z(Zh_ff, Val_mod), Rango_reducido = SISE_z(Zh_rr, Val_mod)), 1)
Tabla_rr <- rbind(Tabla_ff, Rango_reducido = Fila_val(Met_rr),
                  Techo_decodificador = Fila_val(Metricas_par(A_mm(Phi_z_est %*% t(sweep(Puntajes[, 1:M_fpca] %*% t(V_fpca[, 1:M_fpca]), 2, Media_d, "+"))), Obs_mm)))
round(Tabla_rr, 3)
round(100 * mean(Zh_rr[Dia_centro, ] < 0), 4)       #Porcentaje de valores recortados en 0 al volver a mm

#EQUIVALENCIA CON LA CABEZA LINEAL DE LA RED (rúbrica 4.2)
#Red lineal: z_k = a_k' W_sat c, s = G_z z + G_e (altitud, ONI)' + g_0 y zhat = mu + V s. Sustituyendo queda el modelo
#funcional-funcional con b0 = mu + V g_0, b1(s, t) = suma_k suma_m w_k(s) G_z,mk psi_m(t) y g_j = V G_e,j (rango <= min(K, M)).
#Al revés, con H_sat' = P diag(sigma) Q', a_k = q_k y G_z = P diag(sigma) dan cualquier H de rango M si K >= M: el mínimo
#de la cabeza lineal sin penalización es Ajustar_rr(M, 0)
K_neuronas <- 8
Dvs   <- svd(t(H_rr[Fil_b1, ]))                                          #H_sat' (M x K_sat) = P diag(sigma) Q'
A_red <- rbind(t(Dvs$v), matrix(0, K_neuronas - M_fpca, K_z_sat))        #K x K_sat: fila k = coeficientes a_k de w_k(s)
G_z   <- cbind(Dvs$u %*% diag(Dvs$d), matrix(0, M_fpca, K_neuronas - M_fpca))   #M x K
G_e   <- t(H_rr[2:3, ]); g_0 <- H_rr[1, ]
#(a) Camino de la red, con la matriz de Gram: z = A W c -> s = G_z z + G_e e + g_0 -> coeficientes mu + V s
Z_red <- C_mod %*% W_sat %*% t(A_red)                                    #N x K: salidas de las neuronas funcionales
S_red <- sweep(Z_red %*% t(G_z) + cbind(Alt_e, ONI_e) %*% t(G_e), 2, g_0, "+")
Zh_red <- Phi_z_est %*% (Media_d + V_fpca[, 1:M_fpca] %*% t(S_red))      #731 x N
#(b) Camino de la regresión: se arma b1(s, t) en la rejilla fina y la integral en s se calcula con Simpson
W_k_fino  <- Phi_z_sat %*% t(A_red)                                      #731 x K: w_k(s)
Psi_fino  <- Phi_z_est %*% V_fpca[, 1:M_fpca]                            #731 x M: psi_m(t)
Beta1_red <- W_k_fino %*% t(G_z) %*% t(Psi_fino)                         #731 x 731: b1(s, t), filas = s
Beta0_red <- as.vector(Phi_z_est %*% Media_d + Psi_fino %*% g_0)
Zh_reg <- Beta0_red + t(Beta1_red) %*% (w_fino * Z_sat) + (Psi_fino %*% G_e[, 1]) %*% t(Alt_e) + (Psi_fino %*% G_e[, 2]) %*% t(ONI_e)
Equivalencia <- c(Red_vs_regresion_integral = max(abs(Zh_red - Zh_reg)), Red_vs_minimos_cuadrados = max(abs(Zh_red - Zh_rr)),
                  Rango_de_b1 = sum(svd(Beta1_red[seq(1, 731, 10), seq(1, 731, 10)])$d > 1e-10 * max(abs(Beta1_red))))
signif(Equivalencia, 3)
stopifnot(Equivalencia[1:2] < 1e-4)
#La diferencia entre (a) y (b) es solo el error de la cuadratura de Simpson frente a la matriz de Gram exacta

#Sensibilidad del FPCA: varianza de validación recuperada con las curvas del modelo frente a todas las curvas de
#entrenamiento del Valle
Fpca_alt <- function(Dx) {
  m <- colMeans(Dx); V <- solve(L_est, eigen(L_est %*% (crossprod(sweep(Dx, 2, m)) / nrow(Dx)) %*% t(L_est), symmetric = TRUE)$vectors[, 1:M_fpca])
  Dc <- sweep(D_mod, 2, m); Rec <- sweep(Dc %*% W_est %*% V %*% t(V), 2, m, "+")
  100 * (1 - SSE_d(Rec, Val_mod) / SSE_d(matrix(m, nrow(D_mod), ncol(D_mod), byrow = TRUE), Val_mod))
}
Ent_todas <- Pares_curvas$YEAR <= 1995 & !(Pares_curvas$STATION_ID %in% Reservadas)
round(c(Curvas_modelo = sum(Ent_mod), Pct_val_modelo = Fpca_alt(D_mod[Ent_mod, ]),
        Curvas_todas_las_del_Valle = sum(Ent_todas), Pct_val_todas = Fpca_alt(t(Coef_z_est[, Ent_todas]))), 2)

#Figura: media, funciones propias (mu ± 2 raíz(rho_m) psi_m), varianza explicada y pesos del satélite h_m(s)
Figura_fpca <- function() {
  par(mfrow = c(2, 4), mar = c(3.2, 4, 2.4, 0.8))
  Mu <- as.vector(Theta_dia %*% Media_d)
  for (m in 1:6) {
    Psi <- as.vector(Theta_dia %*% V_fpca[, m]); Am <- 2 * sqrt(Var_fpca[m])
    plot(t365, Mu, type = "l", lwd = 2, ylim = range(Mu + Am * Psi, Mu - Am * Psi), xaxt = "n", xlab = "", ylab = "Raíz de mm/día",
         main = sprintf("Componente %d (%.1f %%)%s", m, Tabla_fpca$Pct_varianza[m], ifelse(m <= M_fpca, "", " - no usada")), cex.main = 0.95)
    lines(t365, Mu + Am * Psi, col = "red", lty = 2); lines(t365, Mu - Am * Psi, col = "blue", lty = 2)
    Eje_meses()
  }
  plot(1:10, 100 * cumsum(Var_fpca)[1:10] / sum(Var_fpca), type = "b", pch = 19, ylim = c(40, 100), xlab = "", ylab = "% de varianza acumulado",
       main = "Varianza explicada", cex.main = 0.95)
  mtext("Número de componentes", 1, line = 2.1, cex = 0.7); abline(v = M_fpca, lty = 2)
  matplot(t365, Phi_sat_dia %*% H_rr[Fil_b1, ], type = "l", lty = 1, lwd = 2, col = hcl.colors(M_fpca, "Dark 3"), xaxt = "n", xlab = "",
          ylab = "Peso del satélite", main = "h_m(s): qué parte del año usa cada componente", cex.main = 0.85)
  abline(h = 0, lty = 3); Eje_meses()
  legend("topleft", paste("Componente", 1:M_fpca), col = hcl.colors(M_fpca, "Dark 3"), lwd = 2, bty = "n", cex = 0.7)
}
pdf(Ruta("Figuras_proyecto/Fig_FPCA_rango_reducido.pdf"), width = 13, height = 6.5)
Figura_fpca()
dev.off()
Mostrar(Figura_fpca())

#Exportación para Python: decodificador, matrices y escalado (fpca_estaciones.csv: fila 0 = mu, filas 1 a 6 = psi_m)
write.csv(data.frame(componente = 0:6, valor_propio = c(NA, Var_fpca[1:6]), pct_varianza = c(NA, Tabla_fpca$Pct_varianza[1:6]),
                     usada = c(TRUE, 1:6 <= M_fpca), setNames(as.data.frame(rbind(Media_d, t(V_fpca[, 1:6]))), paste0("c", 1:K_z_est))),
          Ruta("Datos_preparados/fpca_estaciones.csv"), row.names = FALSE)
write.csv(W_est, Ruta("Datos_preparados/gram_estaciones.csv"), row.names = FALSE)
write.csv(W_sat, Ruta("Datos_preparados/gram_satelite.csv"), row.names = FALSE)
write.csv(R_sat, Ruta("Datos_preparados/rugosidad_satelite.csv"), row.names = FALSE)
#Escalado de las covariables (media y desviación de entrenamiento) y solución de referencia de la cabeza lineal:
#filas = intercepto, altitud, ONI y los K_sat pesos del satélite (se aplican a W_sat c); columnas = puntajes 1 a M
write.csv(data.frame(variable = c("altitud", "oni"), media = c(mean(Modelo$Altitud[Ent_mod]), mean(ONI_ent)),
                     desviacion = c(sd(Modelo$Altitud[Ent_mod]), sd(ONI_ent))), Ruta("Datos_preparados/modelo_escalado.csv"), row.names = FALSE)
write.csv(data.frame(entrada = c("intercepto", "altitud", "oni", paste0("satelite_", 1:K_z_sat)),
                     setNames(as.data.frame(Ajustar_rr(M_fpca, 0)), paste0("sin_penalizar_s", 1:M_fpca)),
                     setNames(as.data.frame(H_rr), paste0("penalizada_s", 1:M_fpca))),
          Ruta("Datos_preparados/regresion_rango_reducido.csv"), row.names = FALSE)

#5. EVALUACIÓN (4.4): CUATRO MÉTRICAS, INTERVALOS POR BOOTSTRAP, WILCOXON PAREADO Y FASE ENSO
#Por par estación-año, sobre las curvas en mm/día de los 365 centros de día. Validación: el bloque con que se eligieron
#lambda y M; prueba (2001-2003): las mismas funciones con Abrir_prueba <- TRUE. Reglas de la evaluación:
#1. Ningún modelo, lambda, M, K, hiperparámetro de las redes, semilla ni umbral del preprocesamiento se elige con la
#   prueba: todos se eligieron con entrenamiento y validación
#2. Se reportan todos los modelos: satélite sin corregir, media de entrenamiento, concurrente, función a función, rango
#   reducido y las tres redes (cabeza lineal, MLP y KAN) con la semilla que eligió su pérdida de validación, más la
#   media y la desviación estándar de las tres semillas de cada red
#3. Regla de decisión: un modelo supera a otro solo si el IC 95 % (bootstrap de estaciones) de la mediana de las
#   diferencias pareadas excluye el cero; se informan también la media de las diferencias y el Wilcoxon por estación
#4. La prueba no tiene años Niña (2001 y 2003 son neutros y 2002 es Niño): para La Niña, la lectura por fase se apoya en
#   la validación
Abrir_prueba <- TRUE

#Cuarta métrica: error en la fecha de los dos picos (definición del EDA: máximo de cada mitad del año en los días 15 a
#350, mitades separadas por el valle de la media de entrenamiento; NA si queda pegado al límite). Error = día predicho -
#día observado (positivo: el modelo llega tarde)
Dia_valle <- ceiling(t_fino[Temp_ent_est$Valles[2]])
Mitad_1_dia <- 15:(Dia_valle - 1); Mitad_2_dia <- (Dia_valle + 1):350
fecha_de(Dia_valle)
Picos_par <- function(Curvas) t(apply(Curvas, 2, function(m) {
  if (anyNA(m)) return(c(Pico_1 = NA, Pico_2 = NA))      #Par sin predicción: los modelos de Python solo traen su bloque
  p1 <- Mitad_1_dia[which.max(m[Mitad_1_dia])]; p2 <- Mitad_2_dia[which.max(m[Mitad_2_dia])]
  c(Pico_1 = if (p1 %in% range(Mitad_1_dia)) NA else p1, Pico_2 = if (p2 %in% range(Mitad_2_dia)) NA else p2)
}))
Picos_obs <- Picos_par(Obs_mm)
#Pares con pico observado en cada mitad (%), fecha mediana y % de pares con el máximo de la primera mitad en el día 15
round(100 * colMeans(!is.na(Picos_obs)), 1); fecha_de(apply(Picos_obs, 2, median, na.rm = TRUE))
round(100 * mean(apply(Obs_mm[Mitad_1_dia, ], 2, which.max) == 1), 1)

#Función de las cuatro métricas de un modelo por par (una fila por fila de Modelo)
Evaluar_modelo <- function(Pred_mm) {
  Met <- Metricas_par(Pred_mm, Obs_mm)
  Pic <- Picos_par(Pred_mm) - Picos_obs
  data.frame(L2_rel = Met$L2_rel, KGE = Met$KGE, Sesgo_pct = Met$Sesgo_pct, Error_pico_1 = Pic[, "Pico_1"], Error_pico_2 = Pic[, "Pico_2"])
}
#Modelos de R; los de Python se agregan si existen sus predicciones en Predicciones_python (un CSV por modelo con el
#formato de la entrega: station_id, year, day, value_mm; day = 1..365)
Pred_modelos <- list(Satelite_sin_corregir = Sat_mm, Media_de_entrenamiento = Clim_mm, Concurrente_ampliado = A_mm(Zh_conc_e),
                     Funcional_funcional = A_mm(Zh_ff), Rango_reducido = A_mm(Zh_rr))
#Funciones de escritura y lectura de predicciones (Pred_mm: matriz 365 x n, una columna por par)
Escribir_predicciones <- function(Pred_mm, Station_id, Year, Archivo) {
  write.csv(data.frame(station_id = rep(Station_id, each = 365), year = rep(Year, each = 365), day = rep(1:365, length(Station_id)),
                       value_mm = round(as.vector(Pred_mm), 4)), Archivo, row.names = FALSE)
}
Leer_predicciones <- function(Archivo) {
  P <- read.csv(Archivo, colClasses = c(station_id = "character"))
  stopifnot(all(c("station_id", "year", "day", "value_mm") %in% names(P)), all(P$day %in% 1:365), all(P$value_mm >= 0),
            !anyDuplicated(P[, c("station_id", "year", "day")]))                     #Sin filas repetidas
  M <- matrix(NA_real_, 365, nrow(Modelo))
  Col <- match(paste(P$station_id, P$year), paste(Modelo$STATION_ID, Modelo$YEAR))
  M[cbind(P$day[!is.na(Col)], Col[!is.na(Col)])] <- P$value_mm[!is.na(Col)]
  M
}
#Predicciones del modelo función a función en validación; control: al leerlas de vuelta dan las mismas curvas (4 decimales)
dir.create(Ruta("Predicciones_R"), showWarnings = FALSE)
Escribir_predicciones(Pred_modelos$Funcional_funcional[, Val_mod], Modelo$STATION_ID[Val_mod], Modelo$YEAR[Val_mod],
                      Ruta("Predicciones_R/regresion_ff_validacion.csv"))
max(abs(Leer_predicciones(Ruta("Predicciones_R/regresion_ff_validacion.csv"))[, Val_mod] - Pred_modelos$Funcional_funcional[, Val_mod]))
Archivos_py <- list.files(Ruta("Predicciones_python"), pattern = "\\.csv$", full.names = TRUE)
for (a in Archivos_py) Pred_modelos[[sub("\\.csv$", "", basename(a))]] <- Leer_predicciones(a)
names(Pred_modelos)

#Función de resumen de un bloque: medianas entre pares de L2, KGE y sesgo; mediana del error absoluto de cada pico en
#los pares con pico en todos los modelos; sesgo agregado 100 (total predicho / total observado - 1). IC 95 % por
#bootstrap de estaciones completas (2000 remuestreos). Comparación pareada: mediana y media de las diferencias por par
#con su IC y Wilcoxon pareado por estación
B_eval <- 2000
Resumir <- function(x) if (all(is.na(x))) NA else median(x, na.rm = TRUE)
Evaluar_bloque <- function(Filas, Semilla = 2026) {
  Filas <- which(Filas)
  Modelos <- names(Pred_modelos)[sapply(Pred_modelos, function(P) !anyNA(P[, Filas]))]       #Solo modelos con todos los pares del bloque
  Ev <- lapply(Pred_modelos[Modelos], function(P) { E <- Evaluar_modelo(P)[Filas, ]; E[, 4:5] <- abs(E[, 4:5]); E })
  Con_pico <- sapply(c("Error_pico_1", "Error_pico_2"), function(k) !Reduce(`|`, lapply(Ev, function(E) is.na(E[[k]]))))
  for (m in Modelos) for (k in colnames(Con_pico)) Ev[[m]][[k]][!Con_pico[, k]] <- NA
  Metr <- colnames(Ev[[1]]); nM <- length(Metr); nMod <- length(Modelos)
  X <- do.call(cbind, lapply(Ev, as.matrix))                                      #pares x (modelos * 5 métricas)
  Tot <- sapply(Pred_modelos[Modelos], function(P) colSums(P[, Filas])); Tot_obs <- colSums(Obs_mm[, Filas])
  Est_id <- Modelo$STATION_ID[Filas]; Fase <- droplevels(Modelo$Fase[Filas])
  Parejas <- t(combn(nMod, 2))[, 2:1, drop = FALSE]                               #(modelo, frente a), todos los pares de modelos
  Col <- function(m, k) (m - 1) * nM + k
  Dp <- do.call(cbind, lapply(seq_len(nrow(Parejas)), function(j)                 #Diferencias por par de L2 y de KGE
    cbind(X[, Col(Parejas[j, 1], 1)] - X[, Col(Parejas[j, 2], 1)], X[, Col(Parejas[j, 1], 2)] - X[, Col(Parejas[j, 2], 2)])))
  Estadisticos <- function(ii) c(apply(X[ii, , drop = FALSE], 2, Resumir), 100 * (colSums(Tot[ii, , drop = FALSE]) / sum(Tot_obs[ii]) - 1),
                                 apply(Dp[ii, , drop = FALSE], 2, median), colMeans(Dp[ii, , drop = FALSE]))
  Grupos <- split(seq_along(Filas), Est_id)
  set.seed(Semilla)
  Rep <- replicate(B_eval, Estadisticos(unlist(Grupos[sample(length(Grupos), replace = TRUE)], use.names = FALSE)))
  Est <- Estadisticos(seq_along(Filas)); IC <- apply(Rep, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
  n1 <- nMod * nM; n2 <- nMod; n3 <- ncol(Dp)
  Tabla <- data.frame(Modelo = c(rep(Modelos, each = nM), Modelos), Metrica = c(rep(Metr, nMod), rep("Sesgo_agregado_pct", nMod)),
                      Valor = Est[1:(n1 + n2)], IC_inf = IC[1, 1:(n1 + n2)], IC_sup = IC[2, 1:(n1 + n2)], row.names = NULL)
  Tabla <- Tabla[order(match(Tabla$Modelo, Modelos), match(Tabla$Metrica, c(Metr[1:3], "Sesgo_agregado_pct", Metr[4:5]))), ]
  Por_est <- lapply(1:2, function(k) sapply(Ev, function(E) tapply(E[[k]], Est_id, mean)))    #estaciones x modelos, L2 y KGE
  Dif <- do.call(rbind, lapply(seq_len(nrow(Parejas)), function(j) do.call(rbind, lapply(1:2, function(k) {
    a <- Parejas[j, 1]; b <- Parejas[j, 2]; i_med <- n1 + n2 + 2 * (j - 1) + k; i_mea <- i_med + n3
    Mejora <- if (k == 1) Por_est[[k]][, a] < Por_est[[k]][, b] else Por_est[[k]][, a] > Por_est[[k]][, b]
    data.frame(Modelo = Modelos[a], Frente_a = Modelos[b], Metrica = Metr[k], Mediana_dif = Est[i_med], IC_inf = IC[1, i_med], IC_sup = IC[2, i_med],
               IC_excluye_cero = IC[1, i_med] > 0 | IC[2, i_med] < 0, Media_dif = Est[i_mea], Media_inf = IC[1, i_mea], Media_sup = IC[2, i_mea],
               p_Wilcoxon = wilcox.test(Por_est[[k]][, a], Por_est[[k]][, b], paired = TRUE)$p.value, Pct_estaciones_mejora = 100 * mean(Mejora), row.names = NULL)
  }))))
  #Por fase ENSO: mediana de L2 con IC, diferencia pareada de L2 frente al satélite sin corregir con IC y medianas de las
  #demás métricas (bootstrap de las estaciones con pares en esa fase)
  Por_fase <- do.call(rbind, lapply(levels(Fase), function(f) {
    idx <- which(Fase == f); Gr <- split(idx, Est_id[idx]); L2 <- X[, Col(1:nMod, 1), drop = FALSE]; D1 <- L2 - L2[, 1]
    Rf <- replicate(B_eval, { ii <- unlist(Gr[sample(length(Gr), replace = TRUE)], use.names = FALSE); c(apply(L2[ii, , drop = FALSE], 2, median), apply(D1[ii, , drop = FALSE], 2, median)) })
    Q <- apply(Rf, 1, quantile, c(0.025, 0.975))
    data.frame(Fase = f, Modelo = Modelos, Pares = length(idx), L2_rel = apply(L2[idx, , drop = FALSE], 2, median), L2_inf = Q[1, 1:nMod], L2_sup = Q[2, 1:nMod],
               Dif_satelite = apply(D1[idx, , drop = FALSE], 2, median), Dif_inf = Q[1, nMod + 1:nMod], Dif_sup = Q[2, nMod + 1:nMod],
               t(sapply(Ev, function(E) apply(as.matrix(E[idx, -1]), 2, Resumir))), row.names = NULL)
  }))
  list(Tabla = Tabla, Diferencias = Dif, Por_fase = Por_fase, Pares = length(Filas), Estaciones = length(Grupos),
       Pares_con_pico_en_todos = colSums(Con_pico), Con_pico_propio_pct = sapply(Pred_modelos[Modelos], function(P) 100 * colMeans(!is.na(Picos_par(P)[Filas, ]))))
}
#Función de la tabla ancha para el informe: "valor [IC inferior; IC superior]"
Tabla_informe <- function(Ev) {
  Txt <- sprintf(ifelse(grepl("pico|Sesgo", Ev$Tabla$Metrica), "%.1f [%.1f; %.1f]", "%.3f [%.3f; %.3f]"), Ev$Tabla$Valor, Ev$Tabla$IC_inf, Ev$Tabla$IC_sup)
  as.data.frame(matrix(Txt, ncol = length(unique(Ev$Tabla$Metrica)), byrow = TRUE, dimnames = list(unique(Ev$Tabla$Modelo), unique(Ev$Tabla$Metrica))))
}
Imprimir_evaluacion <- function(Ev) {
  print(c(Pares = Ev$Pares, Estaciones = Ev$Estaciones, Ev$Pares_con_pico_en_todos))
  print(Tabla_informe(Ev))
  D <- Ev$Diferencias
  print(data.frame(D[, 1:3], round(D[, 4:6], 3), D[7], round(D[, 8:10], 3), p_Wilcoxon = signif(D$p_Wilcoxon, 2), Pct_est_mejora = round(D$Pct_estaciones_mejora, 1)))
  print(data.frame(Ev$Por_fase[, 1:3], round(Ev$Por_fase[, -(1:3)], 3)))
  print(round(Ev$Con_pico_propio_pct, 1))
}

Eval_val <- Evaluar_bloque(Val_mod)
#Salida: tamaño del bloque; métricas; diferencias pareadas (L2: negativa = el primer modelo es mejor; KGE: positiva =
#mejor); resultados por fase ENSO; % de pares con pico propio en cada mitad
Imprimir_evaluacion(Eval_val)
#Lectura (validación):
#1. Las tres regresiones bajan el L2 del satélite sin corregir (IC de la diferencia sin cero). En KGE, con la regla de
#   decisión (IC de la mediana de las diferencias), no se declara que un modelo supere al satélite
#2. Función a función frente a concurrente: la mediana de las diferencias pareadas no los distingue
#3. Sesgo del acumulado: mediana por par cerca de 0, sesgo agregado negativo en los modelos (el regreso a mm,
#   E[z]^2 < E[z^2], y en el FF y el rango reducido un nivel más bajo en 1996-2000; sección 6)
#4. Ningún modelo mejora la fecha del primer pico; la del segundo mejora 2 a 3 días con FF y rango reducido
#5. Redes (semilla elegida por validación): las tres mejoran al satélite en L2 (KAN -0,041 [-0,066; -0,022]) y no se
#   distinguen entre sí en L2. Son cifras optimistas, porque la validación eligió la parada temprana, los
#   hiperparámetros y la semilla (media de las tres semillas en L2: 0,382 ± 0,011 KAN, 0,374 ± 0,006 MLP, 0,377 ± 0,003 lineal)
table(Modelo$YEAR[Modelo$Bloque == "Prueba"], Modelo$Fase[Modelo$Bloque == "Prueba"])

#Función de la tabla de semillas de las redes (rúbrica 4.3): métricas por corrida y media y DE entre semillas por cabeza.
#En la tabla de modelos entra la semilla con menor pérdida de validación, elegida en Python
Tabla_semillas <- function(Filas) {
  Archivos <- list.files(Ruta("Predicciones_python/semillas"), pattern = "\\.csv$", full.names = TRUE)
  if (length(Archivos) == 0) return(NULL)
  Metr <- c("L2_rel", "KGE", "Sesgo_pct", "Sesgo_agregado_pct", "Error_pico_1", "Error_pico_2")
  Por_corrida <- do.call(rbind, lapply(Archivos, function(a) {
    P <- Leer_predicciones(a)
    if (anyNA(P[, Filas])) return(NULL)                                   #Corrida sin todos los pares del bloque
    E <- Evaluar_modelo(P)[Filas, ]
    data.frame(Cabeza = sub("_s[0-9]+\\.csv$", "", basename(a)), Semilla = as.integer(sub(".*_s([0-9]+)\\.csv$", "\\1", basename(a))),
               L2_rel = median(E$L2_rel), KGE = median(E$KGE), Sesgo_pct = median(E$Sesgo_pct),
               Sesgo_agregado_pct = 100 * (sum(P[, Filas]) / sum(Obs_mm[, Filas]) - 1),
               Error_pico_1 = median(abs(E$Error_pico_1), na.rm = TRUE), Error_pico_2 = median(abs(E$Error_pico_2), na.rm = TRUE))
  }))
  Resumen <- do.call(rbind, lapply(split(Por_corrida, Por_corrida$Cabeza), function(d)
    data.frame(Cabeza = d$Cabeza[1], Semillas = nrow(d),
               t(setNames(c(colMeans(d[, Metr]), apply(d[, Metr], 2, sd)), c(paste0(Metr, "_media"), paste0(Metr, "_DE")))))))
  list(Por_corrida = Por_corrida, Resumen = Resumen)
}
Imprimir_semillas <- function(S) if (!is.null(S)) {
  print(data.frame(S$Por_corrida[, 1:2], round(S$Por_corrida[, -(1:2)], 3)), row.names = FALSE)
  print(data.frame(S$Resumen[, 1:2], round(S$Resumen[, -(1:2)], 3)), row.names = FALSE)
}
Sem_val <- Tabla_semillas(Val_mod)
Imprimir_semillas(Sem_val)

if (Abrir_prueba) {
  Eval_prueba <- Evaluar_bloque(Modelo$Bloque == "Prueba")
  Imprimir_evaluacion(Eval_prueba)
  Imprimir_semillas(Tabla_semillas(Modelo$Bloque == "Prueba"))
}
#Lectura (prueba, 2001-2003: 309 pares de 115 estaciones):
#1. Todos los modelos mejoran al satélite sin corregir en L2 y en KGE (IC de la mediana de las diferencias sin cero). El
#   satélite sobrestima el total del bloque en 16,0 % [6,8; 25,6]; las regresiones y la cabeza lineal lo llevan cerca de 0
#2. Ablación (4.3): la KAN no se distingue de la cabeza lineal ni de la MLP, ni en L2 ni en KGE; la MLP supera a la
#   lineal solo en KGE. No hay evidencia de que la curva necesite una corrección no lineal del sesgo en L2
#3. Frente a la regresión funcional no hay diferencia en L2 (KAN - FF +0,007 [-0,012; 0,018])
#4. Semillas: la DE entre semillas en L2 (0,001 a 0,006) es mucho menor que la mejora frente al satélite (~0,05)
#5. Conclusión: corregir el satélite con regresión funcional o con la red funcional mejora el error y el KGE; la KAN
#   alcanza la precisión de la mejor regresión y de la MLP con activaciones interpretables, y entre los mejores modelos
#   (FF, rango reducido y las tres redes) no hay un ganador claro

#DIAGNÓSTICO FRENTE A LA MEDIA DE ENTRENAMIENTO (descriptivo; no interviene en la selección de modelos)
#Con mu(t) = curva media de entrenamiento, por par: L2_media = raíz(integral (pred - obs)^2 / integral (obs - mu)^2)
#(< 1: mejora a la media; 1 - L2_media^2 es la habilidad de Murphy, 1988); Amplitud = raíz(integral (pred - mu)^2 /
#integral (obs - mu)^2) (0: la predicción es la media); r_anomalia = correlación entre pred - mu y obs - mu. Entre pares:
#Mejora_media_pct (% con L2_media < 1), r_totales y DE_totales (acumulados anuales predichos frente a observados).
#Medianas con IC 95 % por bootstrap de estaciones
Metricas_media <- function(P, Filas) {                     #Filas: índices de columna; devuelve pares x 3
  Mu <- Clim_mm[, 1]; Ao <- Obs_mm[, Filas] - Mu; Ap <- P[, Filas] - Mu; Den <- sqrt(colSums(Ao^2))
  cbind(L2_media = sqrt(colSums((P[, Filas] - Obs_mm[, Filas])^2)) / Den, Amplitud = sqrt(colSums(Ap^2)) / Den,
        r_anomalia = sapply(seq_along(Filas), function(i) cor(Ap[, i], Ao[, i])))
}
#Control: la observación da L2_media = 0, Amplitud = 1 y r_anomalia = 1; la media de entrenamiento, L2_media = 1
Control_obs  <- Metricas_media(Obs_mm, which(Val_mod))
Control_clim <- suppressWarnings(Metricas_media(Clim_mm, which(Val_mod)))
stopifnot(max(abs(Control_obs - rep(c(0, 1, 1), each = nrow(Control_obs)))) < 1e-12,
          max(abs(Control_clim[, 1:2] - rep(c(1, 0), each = nrow(Control_clim)))) < 1e-12)
Diagnostico_media <- function(Filas, Semilla = 2026) {
  Filas <- which(Filas)
  Modelos <- setdiff(names(Pred_modelos)[sapply(Pred_modelos, function(P) !anyNA(P[, Filas]))], "Media_de_entrenamiento")
  Par <- lapply(Pred_modelos[Modelos], Metricas_media, Filas = Filas)
  Tot <- sapply(Pred_modelos[Modelos], function(P) colSums(P[, Filas])); Tot_obs <- colSums(Obs_mm[, Filas])
  Estadisticos <- function(ii) unlist(lapply(Modelos, function(m) c(
    apply(Par[[m]][ii, , drop = FALSE], 2, median, na.rm = TRUE), Mejora_media_pct = 100 * mean(Par[[m]][ii, "L2_media"] < 1),
    r_totales = cor(Tot[ii, m], Tot_obs[ii]), DE_totales = sd(Tot[ii, m]) / sd(Tot_obs[ii]))))
  Grupos <- split(seq_along(Filas), Modelo$STATION_ID[Filas])
  set.seed(Semilla)
  Rep <- replicate(B_eval, Estadisticos(unlist(Grupos[sample(length(Grupos), replace = TRUE)], use.names = FALSE)))
  Est <- Estadisticos(seq_along(Filas)); IC <- apply(Rep, 1, quantile, c(0.025, 0.975))
  Nombres <- c("L2_media", "Amplitud", "r_anomalia", "Mejora_media_pct", "r_totales", "DE_totales")
  Txt <- sprintf(ifelse(rep(Nombres, length(Modelos)) == "Mejora_media_pct", "%.1f [%.1f; %.1f]", "%.3f [%.3f; %.3f]"),
                 Est, IC[1, ], IC[2, ])
  as.data.frame(matrix(Txt, ncol = length(Nombres), byrow = TRUE, dimnames = list(Modelos, Nombres)))
}
Diagnostico_media(Val_mod)
if (Abrir_prueba) Diagnostico_media(Modelo$Bloque == "Prueba")
#Lectura (IC de cada modelo, descriptivos):
#1. El satélite sin corregir no mejora a la media de entrenamiento (L2_media 1,02 en validación y 1,01 en prueba);
#   todos los modelos la mejoran (L2_media 0,86 a 0,93), pero poco, porque la media trae la mayor parte de la forma anual
#2. La predicción se aparta de la media menos que la observación (Amplitud 43 % a 62 %); la cabeza lineal es la más
#   encogida hacia la media
#3. La KAN y la MLP siguen el acumulado de cada estación y año (r_totales 0,82 a 0,85, frente a 0,55 a 0,63 del satélite
#   y de los demás modelos); dentro del año, la KAN sigue peor la forma en prueba (r_anomalia 0,22 [0,15; 0,28])

#DESFASE ENTRE PREDICCIÓN Y OBSERVACIÓN (descriptivo; no interviene en la selección de modelos)
#Por par, Adelanto = d (de -60 a 60 días) que hace máxima la correlación entre pred(t) y obs(t + d) en los días en que
#se solapan (d > 0: la predicción va adelantada). Es distinto del error en la fecha de los picos; no se corrige nada
Adelanto_max_cor <- function(P, O, Max = 60) {            #P, O: 365 x n en mm/día; devuelve el adelanto de cada par
  D <- -Max:Max
  sapply(seq_len(ncol(O)), function(i) {
    r <- sapply(D, function(d) if (d >= 0) cor(P[1:(365 - d), i], O[(1 + d):365, i]) else cor(P[(1 - d):365, i], O[1:(365 + d), i]))
    D[which.max(r)]
  })
}
#Control con un corrimiento conocido: adelantada 5 días da +5; atrasada 5 días da -5
O_ej <- Obs_mm[, which(Val_mod)[1]]
stopifnot(Adelanto_max_cor(matrix(c(O_ej[6:365], O_ej[1:5])), matrix(O_ej)) == 5,
          Adelanto_max_cor(matrix(c(O_ej[361:365], O_ej[1:360])), matrix(O_ej)) == -5)
#Función del desfase por bloque y modelo: mediana, cuartiles y % de pares con |d| <= 7 días
Desfase_bloque <- function(Filas) {
  Filas <- which(Filas)
  Modelos <- names(Pred_modelos)[sapply(Pred_modelos, function(P) !anyNA(P[, Filas]))]
  t(sapply(Pred_modelos[Modelos], function(P) {
    d <- Adelanto_max_cor(P[, Filas, drop = FALSE], Obs_mm[, Filas, drop = FALSE])
    c(Mediana = median(d), P25 = unname(quantile(d, 0.25)), P75 = unname(quantile(d, 0.75)),
      Pct_a_7_dias_o_menos = 100 * mean(abs(d) <= 7))
  }))
}
round(Desfase_bloque(Val_mod), 1)
if (Abrir_prueba) round(Desfase_bloque(Modelo$Bloque == "Prueba"), 1)
#Lectura: las medianas del adelanto son de pocos días y la dispersión entre pares es del mismo tamaño que la de la media
#de entrenamiento: no hay un desfase sistemático que corregir

#Lectura del ANÁLISIS DE SENSIBILIDAD DE LOS DATOS ("Rscript Proyecto_final_AnalíticaDeDatos.R --sensibilidad" después de
#"Proyecto_final_AnalíticaDeDatos.py --sensibilidad"; no interviene en la selección de modelos)
#Cambios en los datos: la regla de picos invalida 4 picos y la imputación completa 185 pares de entrenamiento (1679 ->
#1858 pares de entrenamiento; validación y prueba iguales). Con validación se vuelven a elegir K_z (15 -> 17), los
#lambdas del FF ((1e9, 1e10) -> (1, 1e11)), M (6 -> 5) y, en la red, K (6 -> 8) y la semilla de la cabeza lineal
#Se sostiene: todos los modelos mejoran al satélite en L2 en validación y en prueba, salvo el FF en prueba; las tres
#cabezas no se distinguen en L2; entre los mejores modelos no hay ganador; el diagnóstico frente a la media no cambia
#No se sostiene: el FF elegido por validación es frágil (lambda_s = 1, L2 de validación 0,374 -> 0,399 y en prueba deja
#de mejorar al satélite). La imputación no mejora el error

#Estabilidad entre años: cada año de 1981 a 2000 se predice con un modelo que no lo vio (entrenamiento: validación
#cruzada por año con los lambdas elegidos; validación: el modelo final)
Zo_ff <- Zh_ff; Zo_conc_e <- Zh_conc_e
for (a in Años_ent_mod) {
  Fuera <- Ent_mod & Modelo$YEAR == a
  Zo_ff[, Fuera] <- (Phi_z_est %*% t(U_ff %*% Ajustar_ff(Lambda_s_ff, Lambda_t_ff, Ent_mod & !Fuera)))[, Fuera]
  Zo_conc_e[, Fuera] <- Predecir_concurrente(X_conc_e, Ajustar_concurrente(Sistema_concurrente(X_conc_e, Ent_mod & !Fuera), Lambda_conc_e,
                                                                           names(X_conc_e)))[, Fuera]
}
Por_año <- function(Pred, Columna) tapply(Metricas_par(Pred, Obs_mm)[[Columna]][Sin_prueba], Modelo$YEAR[Sin_prueba], median)
Tabla_año_mod <- data.frame(Año = 1981:2000, Fase = ONI_año$Fase[match(1981:2000, ONI_año$Año)], ONI = ONI_año$ONI[match(1981:2000, ONI_año$Año)],
                            Pares = as.vector(table(Modelo$YEAR[Sin_prueba])),
                            Sesgo_satelite = Por_año(Sat_mm, "Sesgo_pct"), Sesgo_concurrente = Por_año(A_mm(Zo_conc_e), "Sesgo_pct"),
                            Sesgo_ff = Por_año(A_mm(Zo_ff), "Sesgo_pct"), L2_satelite = Por_año(Sat_mm, "L2_rel"),
                            L2_concurrente = Por_año(A_mm(Zo_conc_e), "L2_rel"), L2_ff = Por_año(A_mm(Zo_ff), "L2_rel"))
data.frame(Tabla_año_mod[, 1:2], round(Tabla_año_mod[, -(1:2)], 2))
#Resumen entre los 20 años: media y desviación del sesgo anual, y número de años en que cada modelo mejora el L2 del satélite
round(rbind(Sesgo_medio = colMeans(Tabla_año_mod[, 5:7]), Sesgo_DE_entre_años = apply(Tabla_año_mod[, 5:7], 2, sd),
            Sesgo_abs_medio = colMeans(abs(Tabla_año_mod[, 5:7]))), 1)
c(Años_concurrente_mejor_que_satelite = sum(Tabla_año_mod$L2_concurrente < Tabla_año_mod$L2_satelite),
  Años_ff_mejor_que_satelite = sum(Tabla_año_mod$L2_ff < Tabla_año_mod$L2_satelite),
  Años_ff_mejor_que_concurrente = sum(Tabla_año_mod$L2_ff < Tabla_año_mod$L2_concurrente))
#Lectura: fuera de la muestra los dos modelos quitan el sesgo mediano del satélite, pero el sesgo de un año concreto
#sigue siendo grande y el FF oscila más entre años que el concurrente

#Figura de ejemplos en validación con el modelo función a función: dos mejores, dos centrales y dos peores pares
Figura_ejemplos <- function() {
  par(mfrow = c(2, 3), mar = c(3.2, 4.2, 2.6, 1))
  Ind <- which(Val_mod); Ord <- Ind[order(Met_ff$L2_rel[Ind])]; n <- length(Ord)
  Sel <- Ord[c(1, 2, round(n / 2), round(n / 2) + 1, n - 1, n)]
  Tipo <- rep(c("Mejor", "Central", "Peor"), each = 2)
  for (k in seq_along(Sel)) {
    i <- Sel[k]
    matplot(t365, cbind(Obs_mm[, i], Sat_mm[, i], Pred_modelos$Funcional_funcional[, i]), type = "l", lty = c(1, 2, 1), lwd = c(2, 1.5, 2),
            col = c("black", "blue", "red"), xaxt = "n", xlab = "", ylab = "mm/día",
            main = sprintf("%s: %s, %d (%s)\nL2 satélite %.2f, regresión %.2f", Tipo[k], Modelo$STATION_ID[i], Modelo$YEAR[i], Modelo$Fase[i],
                           Base_sat$L2_rel[i], Met_ff$L2_rel[i]), cex.main = 0.85)
    Eje_meses()
    if (k == 1) legend("bottomleft", c("Estación", "Satélite sin corregir", "Regresión función a función"), col = c("black", "blue", "red"),
                       lty = c(1, 2, 1), lwd = 2, bty = "n", cex = 0.75)
  }
}
pdf(Ruta("Figuras_proyecto/Fig_Regresion_ejemplos.pdf"), width = 11, height = 6.5)
Figura_ejemplos()
dev.off()
Mostrar(Figura_ejemplos())

#PREDICCIÓN DE ESTACIONES QUE NO ESTÁN EN LA TABLA DE MODELADO: el modelo función a función solo necesita el píxel de
#la estación, el año y la altitud (u = (1, altitud, ONI, W_sat c)); devuelve 365 x n en mm/día, no negativa
Predecir_ff_nuevo <- function(Pixel, Año, Altitud) {
  Fila <- match(paste(Pixel, Año), paste(Curvas_sat$Pixel, Curvas_sat$Año))
  stopifnot(!anyNA(Fila), !anyNA(Altitud), Año %in% ONI_año$Año)
  U <- cbind(1, (Altitud - mean(Modelo$Altitud[Ent_mod])) / sd(Modelo$Altitud[Ent_mod]),
             (ONI_año$ONI[match(Año, ONI_año$Año)] - mean(ONI_ent)) / sd(ONI_ent), t(Coef_z_sat[, Fila, drop = FALSE]) %*% W_sat)
  pmax(Theta_dia %*% t(U %*% Th_ff), 0)^2
}
#Control: con los pares de validación reproduce las predicciones del modelo
max(abs(Predecir_ff_nuevo(Modelo$Pixel[Val_mod], Modelo$YEAR[Val_mod], Modelo$Altitud[Val_mod]) - Pred_modelos$Funcional_funcional[, Val_mod]))
#Estaciones reservadas que se pueden predecir (Diag_reservadas): sus 23 años en el formato de entrega (sin reservadas
#no se escribe nada)
Res_pred <- Diag_reservadas[Diag_reservadas$Se_puede_predecir, ]
c(Reservadas = nrow(Diag_reservadas), Se_pueden_predecir = nrow(Res_pred))
if (nrow(Res_pred) > 0) {
  Alt_res <- ifelse(is.na(Res_pred$Altitud_usada), Res_pred$Altitud_catalogo, Res_pred$Altitud_usada)
  Rej_res <- expand.grid(i = seq_len(nrow(Res_pred)), Año = Años_sat)
  Escribir_predicciones(Predecir_ff_nuevo(Res_pred$Pixel[Rej_res$i], Rej_res$Año, Alt_res[Rej_res$i]), Res_pred$STATION_ID[Rej_res$i], Rej_res$Año,
                        Ruta("Predicciones_R/reservadas_regresion_ff.csv"))
}

#6. SENSIBILIDAD A LA ESCALA DE MODELADO: concurrente ampliado y función a función con la curva en mm sin transformar,
#con log(1 + curva) y con la raíz (mismas bases, lambda elegido con validación); comparación en mm
X_est_mod <- pmax(Mm_est_fino[, Modelo$Fila_est], 0); X_sat_mod <- pmax(Mm_sat_fino[, Modelo$Fila_sat], 0)
Proyectar <- function(Y, Phi) solve(crossprod(Phi), crossprod(Phi, Y))
Ajustar_en_escala <- function(g, g_inv) {
  De <- t(Proyectar(g(X_est_mod), Phi_z_est)); Cs <- t(Proyectar(g(X_sat_mod), Phi_z_sat))
  Ze <- Phi_z_est %*% t(De); Zs <- Phi_z_sat %*% t(Cs)
  #Funcional-funcional
  U <- cbind(1, Alt_e, ONI_e, Cs %*% W_sat); UtU <- crossprod(U[Ent_mod, ]); UtD <- crossprod(U[Ent_mod, ], De[Ent_mod, ])
  Ajuste <- function(ls, lt) matrix(solve(kronecker(W_est, UtU + ls * P_s) + lt * kronecker(R_est, P_t), as.vector(UtD %*% W_est)), ncol(U), ncol(De))
  Error  <- function(Dh) { E <- De[Val_mod, ] - Dh[Val_mod, ]; sum((E %*% W_est) * E) }
  Rej <- expand.grid(ls = Lambdas_ff, lt = Lambdas_ff)
  m <- Rej[which.min(sapply(seq_len(nrow(Rej)), function(k) Error(U %*% Ajuste(Rej$ls[k], Rej$lt[k])))), ]
  Pred_ff <- g_inv((Phi_z_est %*% t(U %*% Ajuste(m$ls, m$lt)))[Dia_centro, ])
  #Concurrente ampliado
  Xs <- list(X_conc_e[[1]], Zs, X_conc_e[[3]], X_conc_e[[4]]); K <- ncol(Theta_beta); p <- length(Xs)
  A <- matrix(0, p * K, p * K); d <- numeric(p * K)
  for (j in 1:p) {
    ij <- (j - 1) * K + 1:K; d[ij] <- crossprod(Theta_beta, w_fino * rowSums(Xs[[j]][, Ent_mod] * Ze[, Ent_mod]))
    for (k in j:p) { ik <- (k - 1) * K + 1:K; A[ij, ik] <- crossprod(Theta_beta, (w_fino * rowSums(Xs[[j]][, Ent_mod] * Xs[[k]][, Ent_mod])) * Theta_beta); A[ik, ij] <- t(A[ij, ik]) }
  }
  Pred_z <- function(l) Predecir_concurrente(Xs, matrix(solve(A + l * kronecker(diag(p), R_beta), d), K, p))
  l_c <- Lambdas_conc[which.min(sapply(Lambdas_conc, function(l) sum(colSums(w_fino * (Ze[, Val_mod] - Pred_z(l)[, Val_mod])^2))))]
  Pred_c <- g_inv(Pred_z(l_c)[Dia_centro, ])
  Fila <- function(P) { M <- Metricas_par(P, Obs_mm)[Val_mod, ]
    c(L2_mediana = median(M$L2_rel), L2_media = mean(M$L2_rel), KGE_mediana = median(M$KGE), Sesgo_mediana_pct = median(M$Sesgo_pct),
      Sesgo_agregado_pct = 100 * (sum(P[, Val_mod]) / sum(Obs_mm[, Val_mod]) - 1)) }
  rbind(Concurrente_ampliado = c(log10_lambda_s = log10(l_c), log10_lambda_t = NA, Fila(Pred_c)),
        Funcional_funcional = c(log10_lambda_s = log10(m$ls), log10_lambda_t = log10(m$lt), Fila(Pred_ff)))
}
Escalas <- list(Sin_transformar_mm = Ajustar_en_escala(function(x) x, function(z) pmax(z, 0)),
                Log_1_mas_P = Ajustar_en_escala(log1p, function(z) pmax(expm1(z), 0)),
                Raiz = Ajustar_en_escala(sqrt, function(z) pmax(z, 0)^2))
lapply(Escalas, round, 3)
#Control: la fila Raiz reproduce los modelos de los puntos 2 y 3
stopifnot(abs(Escalas$Raiz["Funcional_funcional", "L2_mediana"] - median(Met_ff$L2_rel[Val_mod])) < 1e-6,
          abs(Escalas$Raiz["Concurrente_ampliado", "L2_mediana"] - median(Met_conc_e$L2_rel[Val_mod])) < 1e-6)
#Lectura: en mm sin transformar es peor; la raíz y el logaritmo quedan prácticamente iguales en L2 y KGE

#7. VALIDACIÓN DE SUPUESTOS DEL MODELO DE REFERENCIA: regresión de rango reducido (punto 4), equivalente de la cabeza
#lineal; sus supuestos son los de la familia lineal. No interviene en la selección de modelos. En la escala raíz:
#  z_i(t) = b0(t) + integral b1(s, t) x_i(s) ds + g1(t) altitud_i + g2(t) ONI_i + e_i(t)
#(S1) forma lineal; (S2) errores independientes entre curvas; (S3) misma covarianza de los errores; (S4) normalidad;
#(S5) significancia de b1 por permutación y R2(t); (S6) estabilidad fuera de muestra. S1 a S5 con residuos de entrenamiento
E_rr     <- Z_est - Zh_rr                                   #731 x N: residuos en raíz sobre la rejilla fina
Nivel_rr <- colSums(w_fino * Zh_rr) / 365                   #Nivel medio del año de la curva ajustada (raíz)
Res_rr   <- colSums(w_fino * E_rr) / 365                    #Residuo medio del año de cada par (raíz)
Tam_rr   <- sqrt(colSums(w_fino * E_rr^2) / 365)            #Tamaño del residuo: raíz del error cuadrático medio del año
#Lectura aproximada de un residuo medio r en % de lluvia, con el nivel medio de entrenamiento z: (z + r)^2 / z^2 - 1
Pct_lluvia <- function(r) 100 * ((mean(Nivel_rr[Ent_mod]) + r)^2 / mean(Nivel_rr[Ent_mod])^2 - 1)
#Función de la media de v por grupo con IC 95 % por bootstrap de estaciones
Media_por_grupo <- function(v, Grupo, Filas = Ent_mod, B = 1000, Semilla = 2026) {
  set.seed(Semilla)
  t(sapply(split(which(Filas), Grupo[Filas], drop = TRUE), function(idx) {
    Gr  <- split(idx, Modelo$STATION_ID[idx])
    Rep <- replicate(B, mean(v[unlist(Gr[sample(length(Gr), replace = TRUE)], use.names = FALSE)]))
    c(Pares = length(idx), Estaciones = length(Gr), Media = mean(v[idx]), IC_inf = unname(quantile(Rep, 0.025)),
      IC_sup = unname(quantile(Rep, 0.975)), Pct_lluvia = Pct_lluvia(mean(v[idx])))
  }))
}

#S1. Forma lineal: residuo medio por banda de altitud, fase ENSO y quintil del nivel ajustado (con buena especificación,
#todos los IC contienen el cero)
Banda_alt <- cut(Modelo$Altitud, c(-Inf, 500, 1000, 1500, 2000, Inf), right = FALSE,
                 labels = c("<500", "500-1000", "1000-1500", "1500-2000", ">=2000"))
Quint_aj  <- cut(Nivel_rr, quantile(Nivel_rr[Ent_mod], 0:5 / 5), include.lowest = TRUE)
S1_alt  <- Media_por_grupo(Res_rr, Banda_alt)
S1_fase <- Media_por_grupo(Res_rr, Modelo$Fase)
S1_aj   <- Media_por_grupo(Res_rr, Quint_aj)
lapply(list(Altitud = S1_alt, Fase_ENSO = S1_fase, Quintil_nivel_ajustado = S1_aj), round, 3)

#S2. Independencia: correlación intraclase del residuo medio entre los años de una estación y entre las estaciones de
#un año, CCI = (MSB - MSW) / (MSB + (n0 - 1) MSW)
CCI <- function(v, Grupo) {
  Grupo <- factor(Grupo); a <- anova(lm(v ~ Grupo)); k <- nlevels(Grupo); N <- length(v)
  n0 <- (N - sum(table(Grupo)^2) / N) / (k - 1); MSB <- a$`Mean Sq`[1]; MSW <- a$`Mean Sq`[2]
  c(Grupos = k, CCI = (MSB - MSW) / (MSB + (n0 - 1) * MSW), F = a$`F value`[1], p = a$`Pr(>F)`[1])
}
S2 <- rbind(Entre_años_de_una_estacion = CCI(Res_rr[Ent_mod], Modelo$STATION_ID[Ent_mod]),
            Entre_estaciones_de_un_año = CCI(Res_rr[Ent_mod], Modelo$YEAR[Ent_mod]))
signif(S2, 3)

#S3. Misma covarianza: tamaño del residuo frente al nivel ajustado (Spearman), por quintil y por vertiente
S3_spearman <- cor.test(Nivel_rr[Ent_mod], Tam_rr[Ent_mod], method = "spearman", exact = FALSE)
Vertiente_mod <- ifelse(substr(Modelo$STATION_ID, 1, 1) == "5", "Pacífico", "Cauca")
S3 <- list(Spearman = c(rho = unname(S3_spearman$estimate), p = S3_spearman$p.value),
           Tamaño_por_quintil = tapply(Tam_rr[Ent_mod], Quint_aj[Ent_mod], mean),
           Tamaño_por_vertiente = tapply(Tam_rr[Ent_mod], Vertiente_mod[Ent_mod], mean),
           Pares_por_vertiente = table(Vertiente_mod[Ent_mod]))
lapply(S3, function(x) if (is.numeric(x)) signif(x, 3) else x)

#S4. Normalidad del residuo medio y del primer puntaje de los residuos: asimetría, curtosis en exceso y Shapiro-Wilk
Ec_rr  <- sweep(E_rr[, Ent_mod], 1, rowMeans(E_rr[, Ent_mod])) * sqrt(w_fino)
Sv_rr  <- svd(t(Ec_rr), nu = 1, nv = 1)
Pt1_rr <- Sv_rr$u[, 1] * Sv_rr$d[1]                         #Primer puntaje de los residuos
Pct_pc1_rr <- 100 * Sv_rr$d[1]^2 / sum(Ec_rr^2)             #Varianza de los residuos que recoge ese componente
Asimetria <- function(x) mean((x - mean(x))^3) / mean((x - mean(x))^2)^1.5
Curtosis  <- function(x) mean((x - mean(x))^4) / mean((x - mean(x))^2)^2 - 3
S4 <- rbind(Residuo_medio = c(Asimetria = Asimetria(Res_rr[Ent_mod]), Curtosis_exceso = Curtosis(Res_rr[Ent_mod]),
                              p_Shapiro = shapiro.test(Res_rr[Ent_mod])$p.value),
            Puntaje_1_residuos = c(Asimetria(Pt1_rr), Curtosis(Pt1_rr), shapiro.test(Pt1_rr)$p.value))
signif(S4, 3); round(Pct_pc1_rr, 1)

#S5. Permutación de H0: b1(s, t) = 0 dadas la altitud y el ONI (B = 999, como Fperm.fd): baja del error integrado de
#entrenamiento al agregar el satélite; R2(t) como en clase
U_ent_rr <- U_ff[Ent_mod, ]; Y_ent_rr <- Puntajes[Ent_mod, 1:M_fpca]
SSE_red  <- function(U, P) sum((Y_ent_rr - U %*% solve(crossprod(U) + P, crossprod(U, Y_ent_rr)))^2)
SSE0_rr  <- SSE_red(U_ent_rr[, 1:3], matrix(0, 3, 3))            #Sin satélite: intercepto, altitud y ONI
T_obs_rr <- SSE0_rr - SSE_red(U_ent_rr, Lambda_s_rr * P_s)
set.seed(2026)
T_perm_rr <- replicate(999, { Up <- U_ent_rr; Up[, Fil_b1] <- Up[sample(nrow(Up)), Fil_b1]; SSE0_rr - SSE_red(Up, Lambda_s_rr * P_s) })
R2_t_rr <- 1 - rowSums(E_rr[, Ent_mod]^2) / rowSums((Z_est[, Ent_mod] - rowMeans(Z_est[, Ent_mod]))^2)
S5 <- c(Baja_del_error_pct = 100 * T_obs_rr / SSE0_rr, Max_permutaciones_pct = 100 * max(T_perm_rr) / SSE0_rr,
        p_permutacion = (1 + sum(T_perm_rr >= T_obs_rr)) / 1000, R2_t_media = mean(R2_t_rr), R2_t_min = min(R2_t_rr),
        R2_t_max = max(R2_t_rr))
round(S5, 4)

#S6. Estabilidad fuera de muestra: residuo medio por bloque (IC por bootstrap de estaciones) y por año
S6_bloque <- Media_por_grupo(Res_rr, Modelo$Bloque, Filas = rep(TRUE, nrow(Modelo)))
S6_año    <- tapply(Res_rr, Modelo$YEAR, mean)
round(S6_bloque, 3); round(S6_año, 3)
#DE de los residuos medios por año en entrenamiento, frente a la que tendrían si los años fueran independientes
c(DE_entre_años = sd(S6_año[as.character(Años_ent_mod)]),
  DE_esperada_si_independientes = sd(Res_rr[Ent_mod]) / sqrt(mean(table(Modelo$YEAR[Ent_mod]))))

Figura_supuestos <- function() {
  par(mfrow = c(2, 3), mar = c(4.5, 4.5, 2.5, 1))
  Puntos_ic <- function(S, titulo, etiqueta_x) {
    x <- seq_len(nrow(S))
    plot(x, S[, "Media"], ylim = range(S[, c("IC_inf", "IC_sup")], 0), pch = 19, xaxt = "n", xlab = etiqueta_x,
         ylab = "Residuo medio del año (raíz)", main = titulo, cex.main = 1)
    axis(1, at = x, labels = rownames(S), cex.axis = 0.8); arrows(x, S[, "IC_inf"], x, S[, "IC_sup"], angle = 90, code = 3, length = 0.04)
    abline(h = 0, lty = 2)
  }
  Puntos_ic(S1_alt, "S1. Residuo por altitud", "Altitud (m)")
  Puntos_ic(S1_fase, "S1. Residuo por fase ENSO", "Fase")
  plot(Nivel_rr[Ent_mod], Tam_rr[Ent_mod], pch = 16, cex = 0.4, col = adjustcolor("grey30", 0.4), xlab = "Nivel ajustado (raíz)",
       ylab = "Tamaño del residuo (raíz)", main = sprintf("S3. Tamaño vs nivel (Spearman %.2f)", S3$Spearman["rho"]), cex.main = 1)
  lines(lowess(Nivel_rr[Ent_mod], Tam_rr[Ent_mod]), col = "red", lwd = 2)
  qqnorm(Res_rr[Ent_mod], pch = 16, cex = 0.4, main = "S4. Residuo medio: Q-Q normal", cex.main = 1,
         xlab = "Cuantiles teóricos", ylab = "Cuantiles de la muestra"); qqline(Res_rr[Ent_mod], col = "red", lwd = 2)
  Cor_e <- cor(t(E_rr[Dia_centro, Ent_mod]))
  contour(t365, t365, Cor_e, nlevels = 10, xlab = "Día (s)", ylab = "Día (t)", main = "Correlación de los residuos", cex.main = 1)
  #El valor observado queda muy lejos de las permutaciones: se grafican en su propia escala y el observado se anota
  hist(100 * T_perm_rr / SSE0_rr, breaks = 30, col = "grey75", border = "white", ylab = "Frecuencia",
       xlab = "Baja del error con el satélite permutado (%)", main = "S5. Permutación de H0: b1 = 0", cex.main = 1)
  legend("topright", bty = "n", cex = 0.9, legend = c(sprintf("Observado: %.1f %%", 100 * T_obs_rr / SSE0_rr),
         sprintf("Máximo permutado: %.2f %%", 100 * max(T_perm_rr) / SSE0_rr), sprintf("p = %.3f (999 permutaciones)", S5["p_permutacion"])))
}
pdf(Ruta("Figuras_proyecto/Fig_Supuestos_rango_reducido.pdf"), width = 13, height = 8)
Figura_supuestos()
dev.off()
Mostrar(Figura_supuestos())

#Versiones de R y de los paquetes usados (ítem 6)
sessionInfo()
 