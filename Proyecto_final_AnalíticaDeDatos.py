# Juan David Reyes Burbano
# María José Suárez Toro
# Marco Antonio Aragón
# Maria Paula Pinillo
# Laura Daniela Astudillo

from __future__ import annotations

import argparse
import copy
import dataclasses
import datetime
import math
import platform
import random
import sys
import tempfile
import time
from dataclasses import dataclass, field
from importlib import metadata
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from kan import KAN
from kan.KANLayer import coef2curve
from kan.utils import SYMBOLIC_LIB, fit_params
from scipy.linalg import cholesky, eigh, solve_triangular
from torch import nn


RUTA_PROYECTO = Path(__file__).resolve().parent
RUTA_DATOS = RUTA_PROYECTO / "Datos_preparados"
ORDEN = 4
K_ESTACIONES = 8
K_SATELITE = 11
KZ_ESTACIONES = 15
KZ_SATELITE = 16
N_FINO = 731


def nodos_bspline(a: float, b: float, K: int, m: int = ORDEN) -> np.ndarray:
    """Nodos de R: extremos repetidos m-1 veces y K-m+2 cortes equiespaciados."""
    cortes = np.linspace(a, b, K - m + 2)
    return np.concatenate([np.repeat(a, m - 1), cortes, np.repeat(b, m - 1)])


def _base_por_orden(t: np.ndarray, a: float, b: float, K: int, m: int = ORDEN):
    """Devuelve las bases de orden 1..m con la recursion de Cox-de Boor de R."""
    tau = nodos_bspline(a, b, K, m)
    t = np.asarray(t, dtype=float)
    L = len(tau) - 1
    bases: Dict[int, np.ndarray] = {}

    B = np.zeros((len(t), L))
    for j in range(L):
        B[:, j] = ((t >= tau[j]) & (t < tau[j + 1])).astype(float)
    ultimo = max(j for j in range(L) if tau[j] < tau[j + 1])
    B[t == b, :] = 0.0
    B[t == b, ultimo] = 1.0
    bases[1] = B

    for k in range(2, m + 1):
        Bn = np.zeros((len(t), L - k + 1))
        for j in range(L - k + 1):
            izq = (
                (t - tau[j]) / (tau[j + k - 1] - tau[j]) * B[:, j]
                if tau[j + k - 1] > tau[j]
                else 0.0
            )
            der = (
                (tau[j + k] - t) / (tau[j + k] - tau[j + 1]) * B[:, j + 1]
                if tau[j + k] > tau[j + 1]
                else 0.0
            )
            Bn[:, j] = izq + der
        B = Bn
        bases[k] = B
    return bases, tau


def base_bspline(t, a: float = 0.0, b: float = 365.0, K: int = K_ESTACIONES, m: int = ORDEN):
    """Matriz Phi (filas = tiempos, columnas = K funciones base)."""
    bases, _ = _base_por_orden(t, a, b, K, m)
    return bases[m]


def base_bspline_deriv2(t, a: float = 0.0, b: float = 365.0, K: int = K_ESTACIONES, m: int = ORDEN):
    """Segunda derivada de las K funciones base, por derivacion de la recursion.

    B''_{j,m} = (m-1)(m-2) [ B_{j,m-2} / ((tau_{j+m-1}-tau_j)(tau_{j+m-2}-tau_j))
                            - B_{j+1,m-2} / ((tau_{j+m-1}-tau_j)(tau_{j+m-1}-tau_{j+1}))
                            - B_{j+1,m-2} / ((tau_{j+m}-tau_{j+1})(tau_{j+m-1}-tau_{j+1}))
                            + B_{j+2,m-2} / ((tau_{j+m}-tau_{j+1})(tau_{j+m}-tau_{j+2})) ]
    Se verifica contra la matriz de rugosidad de la V2 (rugosidad_satelite.csv) en cargar_datos_funcionales."""
    bases, tau = _base_por_orden(t, a, b, K, m)
    B2 = bases[m - 2]
    t = np.asarray(t, dtype=float)
    D = np.zeros((len(t), len(tau) - m))

    def div(num, den):
        num = np.asarray(num, dtype=float)
        den = np.asarray(den, dtype=float)
        out = np.zeros(np.broadcast(num, den).shape)
        np.divide(num, den, out=out, where=den != 0)
        return out

    for j in range(len(tau) - m):
        d = lambda i, k: tau[i] - tau[k]
        D[:, j] = (m - 1) * (m - 2) * (
            div(B2[:, j], d(j + m - 1, j) * d(j + m - 2, j))
            - div(B2[:, j + 1], d(j + m - 1, j) * d(j + m - 1, j + 1))
            - div(B2[:, j + 1], d(j + m, j + 1) * d(j + m - 1, j + 1))
            + div(B2[:, j + 2], d(j + m, j + 1) * d(j + m, j + 2))
        )
    return D


def matriz_gram(
    K: int,
    a: float = 0.0,
    b: float = 365.0,
    m: int = ORDEN,
    derivada: int = 0,
    n_cuadratura: int = 8,
) -> np.ndarray:
    """W = integral phi phi' (o de phi'' phi'') por cuadratura de Gauss-Legendre."""
    tau = np.unique(nodos_bspline(a, b, K, m))
    xg, wg = np.polynomial.legendre.leggauss(n_cuadratura)
    W = np.zeros((K, K))
    for lo, hi in zip(tau[:-1], tau[1:]):
        t = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo)
        w = 0.5 * (hi - lo) * wg
        if derivada == 0:
            B = base_bspline(t, a, b, K, m)
        elif derivada == 2:
            B = base_bspline_deriv2(t, a, b, K, m)
        else:
            raise ValueError("Solo se implementa derivada 0 o 2.")
        W += B.T @ (B * w[:, None])
    return 0.5 * (W + W.T)


def representar_raiz(X_fino: np.ndarray, Kz: int, Phi_z: np.ndarray) -> np.ndarray:
    """Coeficientes B-spline de z(t)=raiz(max(x(t),0)) sobre una rejilla fina."""
    Z = np.sqrt(np.maximum(np.asarray(X_fino, dtype=float), 0.0))
    if Z.ndim == 1:
        return np.linalg.solve(Phi_z.T @ Phi_z, Phi_z.T @ Z)
    return np.linalg.solve(Phi_z.T @ Phi_z, Phi_z.T @ Z)


@dataclass
class DatosFuncionales:
    """Insumos funcionales reconstruidos y sus matrices de bases."""

    pares: pd.DataFrame
    curvas_largo: pd.DataFrame
    emparejamiento: pd.DataFrame
    altitud: pd.DataFrame
    oni: pd.DataFrame
    pixeles: pd.DataFrame
    pentadas: pd.DataFrame
    coef_est_mm: np.ndarray
    sigma2_est: np.ndarray
    coef_est_z: np.ndarray
    coef_sat_mm: Dict[int, np.ndarray]
    coef_sat_z: Dict[int, np.ndarray]
    t_fine: np.ndarray
    phi_est_mm_fine: np.ndarray
    phi_sat_mm_fine: np.ndarray
    phi_est_z_fine: np.ndarray
    phi_sat_z_fine: np.ndarray
    W_est_mm: np.ndarray
    W_sat_mm: np.ndarray
    W_est_z: np.ndarray
    W_sat_z: np.ndarray
    R_est_z: np.ndarray
    R_sat_z: np.ndarray
    nudos_est_mm: np.ndarray
    nudos_sat_mm: np.ndarray
    nudos_est_z: np.ndarray
    nudos_sat_z: np.ndarray
    indice_pixel: Dict[int, int] = field(default_factory=dict)
    verificacion: Dict[str, float] = field(default_factory=dict)

    def curva_mm_estacion(self, i: int) -> np.ndarray:
        return self.phi_est_mm_fine @ self.coef_est_mm[i]

    def curva_z_estacion(self, i: int) -> np.ndarray:
        return self.phi_est_z_fine @ self.coef_est_z[i]

    def curva_mm_satelite(self, pixel: int, year: int) -> np.ndarray:
        return self.phi_sat_mm_fine @ self.coef_sat_mm[int(year)][:, self.indice_pixel[int(pixel)]]

    def curva_z_satelite(self, pixel: int, year: int) -> np.ndarray:
        return self.phi_sat_z_fine @ self.coef_sat_z[int(year)][:, self.indice_pixel[int(pixel)]]


def _leer_csv(nombre: str, **kwargs) -> pd.DataFrame:
    return pd.read_csv(RUTA_DATOS / nombre, **kwargs)


def cargar_datos_funcionales(verificar: bool = True) -> DatosFuncionales:
    """Carga los CSV exportados y reconstruye las curvas paso a paso como en R."""
    # K_z se lee de la configuración de la V2 (con --sensibilidad puede cambiar); K se comprueba contra las constantes
    global KZ_ESTACIONES, KZ_SATELITE
    conf_est, conf_sat = _leer_csv("estaciones_config.csv"), _leer_csv("satelite_config.csv")
    if int(conf_est.K.iloc[0]) != K_ESTACIONES or int(conf_sat.K.iloc[0]) != K_SATELITE:
        raise ValueError("La V2 eligió otro K para las curvas en mm: actualizar K_ESTACIONES y K_SATELITE")
    raiz = _leer_csv("raiz_config.csv").set_index("fuente")
    KZ_ESTACIONES, KZ_SATELITE = int(raiz.loc["estaciones", "K_z"]), int(raiz.loc["satelite", "K_z"])
    pares = _leer_csv("estaciones_pares.csv", dtype={"station_id": str})
    curvas = _leer_csv("estaciones_curvas_largo.csv", dtype={"station_id": str})
    curvas = curvas.sort_values(["station_id", "year", "t"], kind="stable").reset_index(drop=True)
    altitud = _leer_csv("estaciones_altitud.csv", dtype={"station_id": str})
    oni = _leer_csv("oni_anual.csv")
    # Tabla de modelado de la V2: pares estación-año con bloque, píxel, altitud, ONI y fase ENSO
    emp = _leer_csv("modelo_pares.csv", dtype={"station_id": str}).rename(
        columns={"altitud": "altitud_usada", "bloque": "bloque_modelo"}
    )
    emp = emp.merge(oni[["year", "transicion"]], on="year", how="left")
    pixeles = _leer_csv("satelite_pixeles.csv")
    pentadas = _leer_csv("satelite_pentadas.csv")
    y_sat = _leer_csv("satelite_y.csv")

    t_fine = np.linspace(0.0, 365.0, N_FINO)
    phi_est_mm_fine = base_bspline(t_fine, 0, 365, K_ESTACIONES)
    phi_sat_mm_fine = base_bspline(t_fine, 0, 365, K_SATELITE)
    phi_est_z_fine = base_bspline(t_fine, 0, 365, KZ_ESTACIONES)
    phi_sat_z_fine = base_bspline(t_fine, 0, 365, KZ_SATELITE)

    W_est_mm = matriz_gram(K_ESTACIONES)
    W_sat_mm = matriz_gram(K_SATELITE)
    W_est_z = matriz_gram(KZ_ESTACIONES)
    W_sat_z = matriz_gram(KZ_SATELITE)
    R_est_z = matriz_gram(KZ_ESTACIONES, derivada=2)
    R_sat_z = matriz_gram(KZ_SATELITE, derivada=2)

    # Curvas de estacion: un minimo cuadratico por par con sus dias validos.
    grupos = curvas.groupby(["station_id", "year"], sort=False)
    claves = list(grupos.groups.keys())
    pos_par = {(s, int(y)): i for i, (s, y) in enumerate(zip(pares.station_id, pares.year))}
    coef_est_mm = np.zeros((len(pares), K_ESTACIONES))
    sigma2_est = np.zeros(len(pares))
    cache_bases: Dict[Tuple[int, ...], np.ndarray] = {}
    for (sid, year), g in grupos:
        clave_t = tuple(np.round(g.t.to_numpy(dtype=float), 10))
        Phi = cache_bases.get(clave_t)
        if Phi is None:
            Phi = base_bspline(g.t.to_numpy(dtype=float), 0, 365, K_ESTACIONES)
            cache_bases[clave_t] = Phi
        y = g.y.to_numpy(dtype=float)
        c = np.linalg.solve(Phi.T @ Phi, Phi.T @ y)
        e = y - Phi @ c
        i = pos_par[(sid, int(year))]
        coef_est_mm[i] = c
        sigma2_est[i] = e @ e / (len(y) - K_ESTACIONES)

    X_est_fine = phi_est_mm_fine @ coef_est_mm.T
    coef_est_z = representar_raiz(X_est_fine, KZ_ESTACIONES, phi_est_z_fine).T
    coef_est_z = np.asarray(coef_est_z)

    # Satélite: un ajuste por año y píxel; la raíz se reexpresa en K_z bases
    pixeles_orden = pixeles.pixel.to_numpy(dtype=int)
    indice_pixel = {int(p): i for i, p in enumerate(pixeles_orden)}
    coef_sat_mm: Dict[int, np.ndarray] = {}
    coef_sat_z: Dict[int, np.ndarray] = {}
    inv_z16 = np.linalg.inv(phi_sat_z_fine.T @ phi_sat_z_fine)
    for year in sorted(pentadas.year.unique()):
        f = pentadas.year.to_numpy() == year
        t_pent = pentadas.t.to_numpy(dtype=float)[f]
        Phi_p = base_bspline(t_pent, 0, 365, K_SATELITE)
        C = np.linalg.solve(Phi_p.T @ Phi_p, Phi_p.T @ y_sat.to_numpy(dtype=float)[f, :])
        coef_sat_mm[int(year)] = C
        X = phi_sat_mm_fine @ C
        coef_sat_z[int(year)] = inv_z16 @ (phi_sat_z_fine.T @ np.sqrt(np.maximum(X, 0.0)))

    verificacion: Dict[str, float] = {}
    if verificar:
        ref_est = _leer_csv("estaciones_raiz_coef.csv", dtype={"station_id": str})
        ref_est = ref_est.set_index(["station_id", "year"]).loc[
            list(zip(pares.station_id, pares.year))
        ]
        ref_C = ref_est[[f"c{i}" for i in range(1, KZ_ESTACIONES + 1)]].to_numpy(dtype=float)
        verificacion["max_abs_coef_z_estaciones"] = float(np.max(np.abs(coef_est_z - ref_C)))

        ref_sat = _leer_csv("satelite_raiz_coef.csv")
        max_sat = 0.0
        for year, C in coef_sat_z.items():
            bloque = ref_sat[ref_sat.year == year].set_index("pixel").loc[pixeles_orden]
            ref = bloque[[f"c{i}" for i in range(1, KZ_SATELITE + 1)]].to_numpy(dtype=float).T
            max_sat = max(max_sat, float(np.max(np.abs(C - ref))))
        verificacion["max_abs_coef_z_satelite"] = max_sat
        if max(verificacion["max_abs_coef_z_estaciones"], max_sat) > 1e-6:
            raise ValueError("Las curvas reconstruidas no coinciden con las de la V2: revisar K, K_z y los CSV")

        verificacion["max_abs_W_est_sim"] = float(np.max(np.abs(W_est_z - W_est_z.T)))
        verificacion["min_eig_W_est"] = float(np.min(np.linalg.eigvalsh(W_est_z)))
        verificacion["min_eig_W_sat"] = float(np.min(np.linalg.eigvalsh(W_sat_z)))
        verificacion["min_eig_R_est"] = float(np.min(np.linalg.eigvalsh(R_est_z)))
        verificacion["min_eig_R_sat"] = float(np.min(np.linalg.eigvalsh(R_sat_z)))
        # Rugosidad R = integral phi'' phi''': igual a la de la V2 (fda) en la base del satélite; en las dos bases, una
        # recta da 0 y la parábola t^2 da 4 * 365 = 1460
        R_v2 = _leer_csv("rugosidad_satelite.csv").to_numpy(dtype=float)
        verificacion["max_abs_R_sat_vs_V2"] = float(np.max(np.abs(R_sat_z - R_v2)))
        for fuente, Kz, R in (("est", KZ_ESTACIONES, R_est_z), ("sat", KZ_SATELITE, R_sat_z)):
            Phi_t = base_bspline(t_fine, 0, 365, Kz)
            recta = np.linalg.lstsq(Phi_t, t_fine, rcond=None)[0]
            parabola = np.linalg.lstsq(Phi_t, t_fine**2, rcond=None)[0]
            verificacion[f"rugosidad_recta_{fuente}"] = float(recta @ R @ recta)
            verificacion[f"rugosidad_t2_{fuente}_menos_1460"] = float(parabola @ R @ parabola - 4 * 365)
        if (verificacion["max_abs_R_sat_vs_V2"] > 1e-12
                or max(abs(verificacion[f"rugosidad_recta_{f}"]) for f in ("est", "sat")) > 1e-8
                or max(abs(verificacion[f"rugosidad_t2_{f}_menos_1460"]) for f in ("est", "sat")) > 1e-6):
            raise ValueError("La matriz de rugosidad no mide la segunda derivada: revisar base_bspline_deriv2")

        # Control: matriz suavizadora de un año completo simétrica e idempotente
        t365 = np.arange(365) + 0.5
        Phi365 = base_bspline(t365, 0, 365, K_ESTACIONES)
        H = Phi365 @ np.linalg.solve(Phi365.T @ Phi365, Phi365.T)
        verificacion["max_abs_H_simetrica"] = float(np.max(np.abs(H - H.T)))
        verificacion["max_abs_H_idempotente"] = float(np.max(np.abs(H @ H - H)))
        verificacion["traza_H"] = float(np.trace(H))

    return DatosFuncionales(
        pares=pares,
        curvas_largo=curvas,
        emparejamiento=emp,
        altitud=altitud,
        oni=oni,
        pixeles=pixeles,
        pentadas=pentadas,
        coef_est_mm=coef_est_mm,
        sigma2_est=sigma2_est,
        coef_est_z=coef_est_z,
        coef_sat_mm=coef_sat_mm,
        coef_sat_z=coef_sat_z,
        t_fine=t_fine,
        phi_est_mm_fine=phi_est_mm_fine,
        phi_sat_mm_fine=phi_sat_mm_fine,
        phi_est_z_fine=phi_est_z_fine,
        phi_sat_z_fine=phi_sat_z_fine,
        W_est_mm=W_est_mm,
        W_sat_mm=W_sat_mm,
        W_est_z=W_est_z,
        W_sat_z=W_sat_z,
        R_est_z=R_est_z,
        R_sat_z=R_sat_z,
        nudos_est_mm=nodos_bspline(0, 365, K_ESTACIONES),
        nudos_sat_mm=nodos_bspline(0, 365, K_SATELITE),
        nudos_est_z=nodos_bspline(0, 365, KZ_ESTACIONES),
        nudos_sat_z=nodos_bspline(0, 365, KZ_SATELITE),
        indice_pixel=indice_pixel,
        verificacion=verificacion,
    )


RUTA_SALIDA = RUTA_PROYECTO / "Resultados_python"
RUTA_PRED = RUTA_PROYECTO / "Predicciones_python"
# Interruptor del bloque de prueba (2001-2003), igual al Abrir_prueba de la V2. La prueba no interviene en el
# entrenamiento, la parada temprana, la búsqueda de hiperparámetros ni la elección de semilla
ABRIR_PRUEBA = True
SEMILLAS = (0, 1, 2)
# M = 6 componentes: lo eligió la V2 por validación entre 4 y 6 (Tabla_M, regresión de rango reducido)
M_FPCA = 6
# Rejilla de hiperparámetros, la misma para las tres cabezas (2 x 3 x 3 = 18 configuraciones, cada una con SEMILLAS)
K_NEURONAS_GRID = (6, 8)          # neuronas funcionales del codificador (el enunciado pide entre 6 y 8)
LAMBDA_GRID = (1e-2, 1e-1, 1.0)   # penalización de rugosidad de las w_k: siempre > 0, la pérdida la exige
LR_GRID = (1e-3, 3e-3, 1e-2)      # tasa de aprendizaje de AdamW
MAX_EPOCHS = 3000
# Parada temprana: validación cada 5 épocas; se para tras 60 revisiones (300 épocas) sin mejora
PATIENCE_CHECKS = 60
CADA_CHECKS = 5
H_KAN = 8
GRID_KAN = 5
K_KAN = 3
RANGO_KAN = [-4.0, 4.0]

torch.set_num_threads(4)
plt.rcParams.update(
    {
        "figure.dpi": 130,
        "savefig.dpi": 150,
        "font.size": 9,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


@dataclass
class DatosItem4:
    df: pd.DataFrame
    Cx: np.ndarray
    Cy: np.ndarray
    alt: np.ndarray
    oni: np.ndarray
    Xobs: np.ndarray
    Xsat: np.ndarray
    idx_train: np.ndarray
    idx_val: np.ndarray
    idx_test: np.ndarray
    idx_reservadas: np.ndarray
    alt_mu: float
    alt_sd: float
    oni_mu: float
    oni_sd: float


@dataclass
class ResultadoEntrenamiento:
    modelo: nn.Module
    historial: List[Tuple[int, float, float]]
    mejor_val: float
    mejor_epoca: int
    epocas: int
    segundos: float


def fijar_semilla(semilla: int) -> None:
    random.seed(semilla)
    np.random.seed(semilla)
    torch.manual_seed(semilla)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def imprimir_tabla(titulo: str, df: pd.DataFrame | None = None) -> None:
    print("\n" + "=" * 82)
    print(titulo)
    print("=" * 82)
    if df is not None:
        print(df.to_string(index=False))


def preparar_datos(D) -> DatosItem4:
    pares = D.pares.copy()
    emp = D.emparejamiento[
        ["station_id", "year", "bloque_modelo", "pixel", "altitud_usada", "oni", "fase", "transicion"]
    ]
    df = pares.merge(emp, on=["station_id", "year"], how="inner")
    df = df[df.pixel.notna() & df.altitud_usada.notna()].reset_index(drop=True)
    # Los mismos pares y la misma partición que la tabla Modelo de R
    assert len(df) == len(emp), "Faltan pares de modelo_pares.csv en estaciones_pares.csv"
    assert (df.bloque == df.bloque_modelo).all(), "La partición no coincide con la de R"
    df["station_id"] = df["station_id"].astype(str)
    df["reservada"] = df["reservada"].fillna(False).astype(bool)

    pos_est = {
        (s, int(y)): i for i, (s, y) in enumerate(zip(pares.station_id, pares.year))
    }
    df["idx_est"] = [pos_est[(s, int(y))] for s, y in zip(df.station_id, df.year)]
    df["idx_pix"] = [D.indice_pixel[int(p)] for p in df.pixel]

    Cx = np.stack(
        [
            D.coef_sat_z[int(y)][:, int(ip)]
            for y, ip in zip(df.year.to_numpy(), df.idx_pix.to_numpy())
        ]
    ).astype(np.float64)
    Cy = D.coef_est_z[df.idx_est.to_numpy()].astype(np.float64)
    Xobs = np.stack([D.curva_mm_estacion(i) for i in df.idx_est]).astype(np.float64)
    Xsat = np.stack(
        [D.curva_mm_satelite(p, y) for p, y in zip(df.pixel, df.year)]
    ).astype(np.float64)
    alt = df.altitud_usada.to_numpy(dtype=np.float64)
    oni = df.oni.to_numpy(dtype=np.float64)

    idx_train = np.where((df.bloque == "entrenamiento") & (~df.reservada))[0]
    idx_val = np.where((df.bloque == "validacion") & (~df.reservada))[0]
    idx_test = np.where((df.bloque == "prueba") & (~df.reservada))[0]
    idx_reservadas = np.where(df.reservada)[0]

    # Estandarizacion de altitud y ONI solo con entrenamiento.
    alt_mu, alt_sd = alt[idx_train].mean(), alt[idx_train].std()
    oni_mu, oni_sd = oni[idx_train].mean(), oni[idx_train].std()

    return DatosItem4(
        df=df,
        Cx=Cx,
        Cy=Cy,
        alt=alt,
        oni=oni,
        Xobs=Xobs,
        Xsat=Xsat,
        idx_train=idx_train,
        idx_val=idx_val,
        idx_test=idx_test,
        idx_reservadas=idx_reservadas,
        alt_mu=float(alt_mu),
        alt_sd=float(alt_sd),
        oni_mu=float(oni_mu),
        oni_sd=float(oni_sd),
    )


def ajustar_fpca(C: np.ndarray, W: np.ndarray, M: int):
    """FPCA con metrica L2: autofunciones b que cumplen b' W b = 1.

    Igual que la V2: covarianza con divisor N (como pca.fd) y cada función propia con
    integral positiva, para que mu, psi y los valores propios coincidan con los de R.
    """
    mu = C.mean(0)
    Cc = C - mu
    Sigma = Cc.T @ Cc / len(C)
    L = cholesky(W, lower=True)
    valores, vectores = eigh(L.T @ Sigma @ L)
    orden = np.argsort(valores)[::-1]
    valores = valores[orden]
    vectores = vectores[:, orden]
    B = solve_triangular(L.T, vectores, lower=False)
    B = B / np.sqrt(np.einsum("lk,lm,mk->k", B, W, B))[None, :]
    B = B * np.sign((W @ B).sum(0))[None, :]
    var_acum = np.cumsum(valores) / valores.sum()
    return mu, B[:, :M], valores, var_acum, B


class EncoderFuncional(nn.Module):
    def __init__(self, n_neuronas: int, W: np.ndarray, R: np.ndarray, n_base: int, semilla: int):
        super().__init__()
        generador = torch.Generator().manual_seed(semilla)
        self.A = nn.Parameter(torch.randn(n_neuronas, n_base, generator=generador) * 0.1)
        self.register_buffer("W", torch.tensor(W, dtype=torch.float32))
        self.register_buffer("R", torch.tensor(R, dtype=torch.float32))

    def pesos_normalizados(self) -> torch.Tensor:
        AW = self.A @ self.W
        norma = torch.sqrt(torch.sum(AW * self.A, dim=1, keepdim=True) + 1e-12)
        return self.A / norma

    def forward(self, cx: torch.Tensor):
        An = self.pesos_normalizados()
        z = (cx @ (self.W @ An.T)) / math.sqrt(365.0)
        penal = torch.mean(torch.sum((An @ self.R) * An, dim=1))
        return z, penal, An


class CabezaKAN(nn.Module):
    # ocultas y grid solo los usa el experimento de arquitectura; por omisión es la KAN [d, H_KAN, M] del proyecto
    def __init__(self, d: int, M: int, h: int, semilla: int, ocultas: Tuple[int, ...] | None = None,
                 grid: int | None = None):
        super().__init__()
        self.kan = KAN(
            width=[d, *(ocultas if ocultas else (h,)), M],
            grid=GRID_KAN if grid is None else grid,
            k=K_KAN,
            grid_range=RANGO_KAN,
            auto_save=False,
            symbolic_enabled=False,
            seed=semilla,
            ckpt_path=str(RUTA_SALIDA / "kan_ckpt"),  # no se escribe nada: auto_save=False
        )
        self.kan.speed()
        # pykan crea siempre una rama simbólica (symbolic_fun) que, con symbolic_enabled=False, no entra al
        # cálculo ni recibe gradiente. Se congela para que el conteo de parámetros entrenables, sum(p.numel())
        # sobre los parámetros con requires_grad, cuente solo los que de verdad se ajustan
        for p in self.kan.symbolic_fun.parameters():
            p.requires_grad_(False)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.kan(x)


def construir_cabeza(tipo: str, d: int, M: int, semilla: int, h_mlp: int | None = None,
                     arquitectura: Dict[str, object] | None = None):
    """arquitectura (solo el experimento de arquitectura): {"ocultas": anchos de las capas ocultas, "grid": rejilla
    de la KAN}. Sin ella, las cabezas del proyecto: MLP con una capa oculta de h_mlp y KAN [d, H_KAN, M]."""
    arquitectura = arquitectura or {}
    if tipo == "lineal":
        cabeza = nn.Linear(d, M)
        nn.init.zeros_(cabeza.weight)
        nn.init.zeros_(cabeza.bias)
        return cabeza
    if tipo == "mlp":
        capas, entrada = [], d
        for h in arquitectura.get("ocultas", (h_mlp,)):
            capas += [nn.Linear(entrada, h), nn.SiLU()]
            entrada = h
        return nn.Sequential(*capas, nn.Linear(entrada, M))
    if tipo == "kan":
        return CabezaKAN(d, M, H_KAN, semilla, **arquitectura)
    raise ValueError(tipo)


def contar_parametros(module: nn.Module) -> int:
    """Parámetros entrenables: sum(p.numel()) sobre los parámetros con requires_grad."""
    return int(sum(p.numel() for p in module.parameters() if p.requires_grad))


def contar_parametros_total(module: nn.Module) -> int:
    """Todos los parámetros registrados, incluidos los fijos de pykan (rejilla, máscaras, rama simbólica)."""
    return int(sum(p.numel() for p in module.parameters()))


def emparejar_parametros(d: int, M: int, semilla: int = 0):
    """Busca el ancho MLP cuyo numero de parametros entrenables queda a menos del 10 % de la KAN."""
    kan = CabezaKAN(d, M, H_KAN, semilla)
    objetivo = contar_parametros(kan)
    del kan
    mejor = None
    for h in range(2, 500):
        mlp = nn.Sequential(nn.Linear(d, h), nn.SiLU(), nn.Linear(h, M))
        n = contar_parametros(mlp)
        rel = abs(n - objetivo) / objetivo
        if mejor is None or rel < mejor[0]:
            mejor = (rel, h, n)
    rel, h, n = mejor
    if rel > 0.10:
        raise RuntimeError("No se pudo igualar el numero de parametros KAN/MLP")
    return h, objetivo, n, rel


class DecodificadorFPCA(nn.Module):
    """Decodificador con μ y ψ congelados: u -> c_hat = μ + B_s u -> z(t)."""

    def __init__(self, mu: np.ndarray, Bs: np.ndarray, phi_est: np.ndarray):
        super().__init__()
        self.register_buffer("mu", torch.tensor(mu, dtype=torch.float32))
        self.register_buffer("basis", torch.tensor(Bs, dtype=torch.float32))
        self.register_buffer("phi", torch.tensor(phi_est, dtype=torch.float32))

    def forward(self, u: torch.Tensor):
        c_hat = self.mu + u @ self.basis.T
        return c_hat, c_hat @ self.phi.T

    def project(self, c: torch.Tensor, W: torch.Tensor) -> torch.Tensor:
        return (c - self.mu) @ W @ self.basis


class RedFuncional(nn.Module):
    def __init__(
        self,
        encoder: EncoderFuncional,
        cabeza: nn.Module,
        decodificador: DecodificadorFPCA,
        alt_mu: float,
        alt_sd: float,
        oni_mu: float,
        oni_sd: float,
    ):
        super().__init__()
        self.encoder = encoder
        self.cabeza = cabeza
        self.decodificador = decodificador
        self.register_buffer("alt_mu", torch.tensor(alt_mu, dtype=torch.float32))
        self.register_buffer("alt_sd", torch.tensor(alt_sd, dtype=torch.float32))
        self.register_buffer("oni_mu", torch.tensor(oni_mu, dtype=torch.float32))
        self.register_buffer("oni_sd", torch.tensor(oni_sd, dtype=torch.float32))

    def forward(self, cx: torch.Tensor, alt: torch.Tensor, oni: torch.Tensor):
        z, penal, An = self.encoder(cx)
        alt_s = (alt - self.alt_mu) / self.alt_sd
        oni_s = (oni - self.oni_mu) / self.oni_sd
        entrada = torch.cat([z, alt_s[:, None], oni_s[:, None]], dim=1)
        u = self.cabeza(entrada)
        c_hat, z_fine = self.decodificador(u)
        return {
            "c_hat": c_hat,
            "u": u,
            "z_fine": z_fine,
            "pen": penal,
            "An": An,
            "entrada": entrada,
        }


def perdida_datos(c_hat: torch.Tensor, c_obs: torch.Tensor, W: torch.Tensor) -> torch.Tensor:
    d = c_hat - c_obs
    return torch.mean(torch.sum((d @ W) * d, dim=1)) / 365.0


def entrenar(
    tipo: str,
    semilla: int,
    K: int,
    lam: float,
    lr: float,
    D,
    datos: DatosItem4,
    mu: np.ndarray,
    Bs: np.ndarray,
    h_mlp: int | None,
    max_epochs: int,
    arquitectura: Dict[str, object] | None = None,
) -> ResultadoEntrenamiento:
    fijar_semilla(semilla)
    encoder = EncoderFuncional(K, D.W_sat_z, D.R_sat_z * 365.0**2, KZ_SATELITE, semilla)
    cabeza = construir_cabeza(tipo, K + 2, Bs.shape[1], semilla, h_mlp, arquitectura)
    decodificador = DecodificadorFPCA(mu, Bs, D.phi_est_z_fine)
    modelo = RedFuncional(
        encoder,
        cabeza,
        decodificador,
        datos.alt_mu,
        datos.alt_sd,
        datos.oni_mu,
        datos.oni_sd,
    )
    optimizador = torch.optim.AdamW(
        [p for p in modelo.parameters() if p.requires_grad], lr=lr, weight_decay=1e-4
    )
    programador = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizador, mode="min", factor=0.5, patience=5, min_lr=1e-5
    )

    cx = torch.tensor(datos.Cx, dtype=torch.float32)
    cy = torch.tensor(datos.Cy, dtype=torch.float32)
    alt = torch.tensor(datos.alt, dtype=torch.float32)
    oni = torch.tensor(datos.oni, dtype=torch.float32)
    W = torch.tensor(D.W_est_z, dtype=torch.float32)
    it = torch.tensor(datos.idx_train)
    iv = torch.tensor(datos.idx_val)

    mejor_val = float("inf")
    mejor_estado = None
    mejor_epoca = -1
    sin_mejora = 0
    historial: List[Tuple[int, float, float]] = []
    inicio = time.time()

    for epoca in range(max_epochs):
        modelo.train()
        optimizador.zero_grad()
        salida = modelo(cx[it], alt[it], oni[it])
        perdida = perdida_datos(salida["c_hat"], cy[it], W) + lam * salida["pen"]
        perdida.backward()
        nn.utils.clip_grad_norm_(modelo.parameters(), max_norm=1.0)
        optimizador.step()

        if epoca % CADA_CHECKS == 0 or epoca == max_epochs - 1:
            modelo.eval()
            with torch.no_grad():
                salida_val = modelo(cx[iv], alt[iv], oni[iv])
                val = float(perdida_datos(salida_val["c_hat"], cy[iv], W))
            programador.step(val)
            historial.append((epoca, float(perdida.detach()), val))
            if val < mejor_val - 1e-9:
                mejor_val = val
                mejor_estado = copy.deepcopy(modelo.state_dict())
                mejor_epoca = epoca
                sin_mejora = 0
            else:
                sin_mejora += 1
                if sin_mejora >= PATIENCE_CHECKS:
                    break

    if mejor_estado is not None:
        modelo.load_state_dict(mejor_estado)
    modelo.eval()
    return ResultadoEntrenamiento(
        modelo=modelo,
        historial=historial,
        mejor_val=mejor_val,
        mejor_epoca=mejor_epoca,
        epocas=epoca + 1,
        segundos=time.time() - inicio,
    )


def predecir(modelo: nn.Module, datos: DatosItem4, idx: np.ndarray):
    modelo.eval()
    with torch.no_grad():
        salida = modelo(
            torch.tensor(datos.Cx[idx], dtype=torch.float32),
            torch.tensor(datos.alt[idx], dtype=torch.float32),
            torch.tensor(datos.oni[idx], dtype=torch.float32),
        )
    return (
        salida["c_hat"].numpy().astype(np.float64),
        salida["z_fine"].numpy().astype(np.float64),
        salida["u"].numpy().astype(np.float64),
    )


def l2_relativo(pred_mm: np.ndarray, obs_mm: np.ndarray) -> np.ndarray:
    """Error L2 relativo por curva, solo para ordenar los ejemplos y peores casos de las figuras.

    Las métricas del informe (L2, KGE, sesgo del acumulado, fecha de picos, IC bootstrap por
    estaciones, Wilcoxon y fase ENSO) las calcula la V2 con las predicciones de este archivo.
    """
    pred = np.maximum(pred_mm, 0.0)
    obs = np.maximum(obs_mm, 0.0)
    return np.sqrt(np.sum((pred - obs) ** 2, axis=1)) / np.maximum(np.sqrt(np.sum(obs**2, axis=1)), 1e-12)


# Centros de día t = 0,5; 1,5; ...; 364,5: son los puntos impares de t_fine y los t365 de la V2
CENTROS_DIA = np.arange(1, N_FINO, 2)
# Nombre de cada cabeza en las tablas de la V2 (nombre del archivo de predicciones)
NOMBRE_R = {"lineal": "Red_lineal", "mlp": "Red_MLP", "kan": "Red_KAN"}


def escribir_predicciones(df: pd.DataFrame, z_fine: np.ndarray, ruta: Path) -> None:
    """Predicciones en el formato que lee la V2 (Leer_predicciones) y que exige la entrega.

    station_id, year, day = 1..365 y value_mm = max(z(day - 0,5), 0)^2: no negativa por construcción.
    """
    mm = np.maximum(z_fine[:, CENTROS_DIA], 0.0) ** 2
    n = len(df)
    pd.DataFrame(
        {
            "station_id": np.repeat(df.station_id.to_numpy(), 365),
            "year": np.repeat(df.year.to_numpy().astype(int), 365),
            "day": np.tile(np.arange(1, 366), n),
            "value_mm": np.round(mm.ravel(), 4),
        }
    ).to_csv(ruta, index=False)


def verificar_fpca_con_r(mu: np.ndarray, B: np.ndarray, valores: np.ndarray, M: int) -> Dict[str, float]:
    """Compara el FPCA de este archivo con el de la V2 (Datos_preparados/fpca_estaciones.csv)."""
    ref = _leer_csv("fpca_estaciones.csv")
    C = ref[[f"c{i}" for i in range(1, KZ_ESTACIONES + 1)]].to_numpy(dtype=float)
    return {
        "max_abs_mu": float(np.max(np.abs(mu - C[0]))),
        "max_abs_psi": float(np.max(np.abs(B[:, :M] - C[1 : M + 1].T))),
        "max_abs_valor_propio": float(
            np.max(np.abs(valores[:M] - ref.valor_propio.to_numpy(dtype=float)[1 : M + 1]))
        ),
    }


def pesos_simpson(t: np.ndarray) -> np.ndarray:
    """Pesos de la regla de Simpson en una rejilla equiespaciada con un número par de intervalos."""
    n = len(t) - 1
    assert n % 2 == 0 and np.allclose(np.diff(t), t[1] - t[0])
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return w * (t[1] - t[0]) / 3.0


def construir_red_desde_r(
    H: np.ndarray, K: int, datos: DatosItem4, D, mu: np.ndarray, Bs: np.ndarray, valores: np.ndarray, M: int
) -> nn.Module:
    """Escribe la regresión de rango reducido de la V2 (H, 19 x M) como una red con cabeza lineal.

    En la V2 los puntajes son s = h0 + h_alt (alt - m_A)/d_A + h_oni (oni - m_O)/d_O + (C W_sat) H_sat, y en
    la red u_m = s_m / sqrt(rho_m). La parte del satélite exige A' G' = sqrt(365) H_sat diag(1/sqrt(rho)) con
    filas a_k de norma 1 en W_sat. Con W_sat = L L' y la descomposición en valores singulares
    L' sqrt(365) H_sat diag(1/sqrt(rho)) = P S Q', basta tomar A' = L'^-1 P (filas ortonormales en W_sat) y
    G = Q S. Si K > M, las neuronas que sobran tienen peso cero en la cabeza.
    """
    esc = _leer_csv("modelo_escalado.csv").set_index("variable")
    m_A, d_A = esc.loc["altitud", "media"], esc.loc["altitud", "desviacion"]
    m_O, d_O = esc.loc["oni", "media"], esc.loc["oni", "desviacion"]
    raiz = np.sqrt(valores[:M])
    h0, h_alt, h_oni, H_sat = H[0], H[1], H[2], H[3:]
    g_alt = h_alt * datos.alt_sd / d_A / raiz
    g_oni = h_oni * datos.oni_sd / d_O / raiz
    g0 = (h0 + h_alt * (datos.alt_mu - m_A) / d_A + h_oni * (datos.oni_mu - m_O) / d_O) / raiz
    L = cholesky(D.W_sat_z, lower=True)
    P, S, Qt = np.linalg.svd(L.T @ (math.sqrt(365.0) * H_sat / raiz[None, :]), full_matrices=False)
    A = np.zeros((K, KZ_SATELITE))
    A[:M] = solve_triangular(L.T, P, lower=False).T
    G = np.zeros((M, K))
    G[:, :M] = Qt.T * S[None, :]
    if K > M:
        extra = np.random.default_rng(0).standard_normal((K - M, KZ_SATELITE))
        A[M:] = extra / np.sqrt(np.einsum("kl,lm,km->k", extra, D.W_sat_z, extra))[:, None]
    encoder = EncoderFuncional(K, D.W_sat_z, D.R_sat_z, KZ_SATELITE, 0)
    cabeza = nn.Linear(K + 2, M)
    modelo = RedFuncional(
        encoder, cabeza, DecodificadorFPCA(mu, Bs[:, :M], D.phi_est_z_fine),
        datos.alt_mu, datos.alt_sd, datos.oni_mu, datos.oni_sd,
    ).double()
    with torch.no_grad():
        modelo.encoder.A.copy_(torch.tensor(A))
        modelo.cabeza.weight.copy_(torch.tensor(np.column_stack([G, g_alt, g_oni])))
        modelo.cabeza.bias.copy_(torch.tensor(g0))
    return modelo


def verificar_equivalencia(
    modelo: nn.Module,
    datos: DatosItem4,
    D,
    mu: np.ndarray,
    B: np.ndarray,
    M: int,
    Bs: np.ndarray,
    valores: np.ndarray,
) -> Tuple[Dict[str, float], np.ndarray]:
    """Equivalencia entre la red con cabeza lineal y la regresión funcional función a función.

    Demostración. En la red, z_k = (1/sqrt(365)) integral w_k(s) x(s) ds, la cabeza lineal da
    u = G z + g_alt alt + g_oni ONI + g_0 y el decodificador zhat(t) = mu(t) + sum_m u_m phi_m(t), con
    phi_m = sqrt(rho_m) psi_m. Sustituyendo z en u y u en zhat:
        zhat(t) = b0(t) + integral beta(s, t) x(s) ds + gamma_alt(t) alt + gamma_oni(t) ONI,
        beta(s, t) = (1/sqrt(365)) sum_k sum_m w_k(s) G_mk phi_m(t),
        b0(t) = mu(t) + sum_m g_0m phi_m(t),   gamma_j(t) = sum_m g_jm phi_m(t),
    que es el modelo funcional-funcional con covariables escalares, con beta de rango <= min(K, M) (la forma
    beta(s, t) = sum_k sum_m w_k(s) B_km psi_m(t) del enunciado). Al revés, toda regresión de rango M con el
    mismo decodificador se escribe como red si K >= M (construir_red_desde_r).

    Verificación numérica (la rúbrica pide diferencia máxima <= 1e-4):
    1. Red entrenada frente a la integral con beta(s, t) evaluada en la rejilla fina y regla de Simpson.
    2. Regresión de rango reducido de la V2 (regresion_rango_reducido.csv) escrita como red: su salida
       frente a la predicción de la regresión mu + U H Psi'.
    3. Pérdida de la cabeza lineal entrenada frente al óptimo cerrado de la regresión (convergencia).
    """
    modelo = copy.deepcopy(modelo).double()
    modelo.eval()
    n = len(datos.df)
    todos = torch.arange(n)
    cx = torch.tensor(datos.Cx, dtype=torch.float64)
    alt = torch.tensor(datos.alt, dtype=torch.float64)
    oni = torch.tensor(datos.oni, dtype=torch.float64)
    with torch.no_grad():
        salida = modelo(cx[todos], alt[todos], oni[todos])
        An = modelo.encoder.pesos_normalizados().numpy()
    z_red = salida["z_fine"].numpy()
    G_full = modelo.cabeza.weight.detach().numpy().astype(np.float64)
    g0 = modelo.cabeza.bias.detach().numpy().astype(np.float64)
    K = An.shape[0]
    G, g_alt, g_oni = G_full[:, :K], G_full[:, K], G_full[:, K + 1]
    alt_s = (datos.alt - datos.alt_mu) / datos.alt_sd
    oni_s = (datos.oni - datos.oni_mu) / datos.oni_sd

    # 1. beta(s, t) en la rejilla fina e integral por Simpson
    t = D.t_fine
    w_s = D.phi_sat_z_fine @ An.T                       # 731 x K: w_k(s)
    phi_t = D.phi_est_z_fine @ Bs                       # 731 x M: sqrt(rho_m) psi_m(t)
    beta = (w_s @ G.T @ phi_t.T) / math.sqrt(365.0)     # 731 (s) x 731 (t)
    x_s = (D.phi_sat_z_fine @ datos.Cx.T).T             # n x 731: curva del satélite en raíz
    integral = (x_s * pesos_simpson(t)[None, :]) @ beta
    z_formula = (
        (D.phi_est_z_fine @ mu)[None, :]
        + (g0[None, :] + np.outer(alt_s, g_alt) + np.outer(oni_s, g_oni)) @ phi_t.T
        + integral
    )
    Theta = Bs @ G @ An @ D.W_sat_z / math.sqrt(365.0)  # beta en las bases: K_est x K_sat

    # 2. Regresión de rango reducido de la V2 escrita como red
    rr = _leer_csv("regresion_rango_reducido.csv")
    esc = _leer_csv("modelo_escalado.csv").set_index("variable")
    U = np.column_stack(
        [
            np.ones(n),
            (datos.alt - esc.loc["altitud", "media"]) / esc.loc["altitud", "desviacion"],
            (datos.oni - esc.loc["oni", "media"]) / esc.loc["oni", "desviacion"],
            datos.Cx @ D.W_sat_z,
        ]
    )
    # El M de la regresión lo elige la V2 con validación; en el análisis principal coincide con M_FPCA. En el de
    # sensibilidad puede ser otro: la regresión se escribe como red con su propio M (hace falta K >= M_R)
    M_R = sum(c.startswith("sin_penalizar_s") for c in rr.columns)
    if M_R != M:
        print(f"Aviso: la regresión de rango reducido de la V2 usa M = {M_R} y la red M = {M}; "
              f"el control 2 escribe la regresión como red con M = {M_R}")
    controles = {}
    for columna in ["sin_penalizar", "penalizada"]:
        H = rr[[f"{columna}_s{m}" for m in range(1, M_R + 1)]].to_numpy(dtype=float)
        c_regresion = mu[None, :] + (U @ H) @ B[:, :M_R].T
        red_r = construir_red_desde_r(H, K, datos, D, mu, Bs, valores, M_R)
        with torch.no_grad():
            c_red_r = red_r(cx, alt, oni)["c_hat"].numpy()
        controles[columna] = (
            float(np.max(np.abs(c_red_r - c_regresion))),
            float(np.max(np.abs((c_red_r - c_regresion) @ D.phi_est_z_fine.T))),
            float(np.mean(np.einsum("ij,jk,ik->i", c_regresion - datos.Cy, D.W_est_z, c_regresion - datos.Cy)[datos.idx_val]) / 365.0),
        )

    # 3. Pérdida de la cabeza lineal entrenada (misma métrica que el entrenamiento)
    c_red = salida["c_hat"].numpy()
    perd = np.einsum("ij,jk,ik->i", c_red - datos.Cy, D.W_est_z, c_red - datos.Cy) / 365.0

    resultado = {
        "max_abs_red_vs_integral_simpson": float(np.max(np.abs(z_red - z_formula))),
        "rango_beta": float(np.linalg.matrix_rank(Theta)),
        "max_abs_coef_regresion_V2_como_red": controles["sin_penalizar"][0],
        "max_abs_curva_regresion_V2_como_red": controles["sin_penalizar"][1],
        "max_abs_coef_regresion_V2_penalizada_como_red": controles["penalizada"][0],
        "perdida_val_cabeza_lineal": float(perd[datos.idx_val].mean()),
        "perdida_val_regresion_V2_sin_penalizar": controles["sin_penalizar"][2],
        "perdida_val_regresion_V2_penalizada": controles["penalizada"][2],
    }
    return resultado, beta


def verificar_kan_sintetica(semilla: int = 0) -> Tuple[float, float]:
    """KAN [2, 1] sobre f(x) = sin(2 pi x1) + x2^2 / 2 en [-1, 1]^2: RMSE en los 1500 puntos de
    ajuste y en 1000 puntos nuevos (la rúbrica pide RMSE <= 1e-2 antes de usar datos reales)."""
    fijar_semilla(semilla)
    kan = KAN(
        width=[2, 1],
        grid=8,
        k=3,
        grid_range=[-1, 1],
        auto_save=False,
        symbolic_enabled=False,
        seed=semilla,
        ckpt_path=str(RUTA_SALIDA / "kan_ckpt"),  # no se escribe nada: auto_save=False
    )
    x = torch.rand(1500, 2) * 2 - 1
    y = torch.sin(2 * math.pi * x[:, 0]) + 0.5 * x[:, 1] ** 2
    opt = torch.optim.Adam(kan.parameters(), lr=1e-2)
    perdida = None
    for _ in range(1500):
        opt.zero_grad()
        perdida = torch.mean((kan(x)[:, 0] - y) ** 2)
        perdida.backward()
        opt.step()
    x_nuevo = torch.rand(1000, 2, generator=torch.Generator().manual_seed(semilla + 1)) * 2 - 1
    y_nuevo = torch.sin(2 * math.pi * x_nuevo[:, 0]) + 0.5 * x_nuevo[:, 1] ** 2
    with torch.no_grad():
        rmse_nuevo = float(torch.sqrt(torch.mean((kan(x_nuevo)[:, 0] - y_nuevo) ** 2)))
        rmse_ajuste = float(torch.sqrt(torch.mean((kan(x)[:, 0] - y) ** 2)))
    return rmse_ajuste, rmse_nuevo


def parametros_w(modelo: nn.Module) -> np.ndarray:
    with torch.no_grad():
        return modelo.encoder.pesos_normalizados().numpy().astype(np.float64)


MESES = ["ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic"]


def fecha_de(t: float) -> str:
    """Día t del dominio [0, 365] como 'dd-mes', igual que fecha_de() de la V2."""
    dia = datetime.date(1981, 1, 1) + datetime.timedelta(days=max(math.ceil(t), 1) - 1)
    return f"{dia.day:02d}-{MESES[dia.month - 1]}"


def temporadas(media: np.ndarray, d1: np.ndarray) -> Dict[str, List[int]]:
    """Regla Temporadas() de la V2 sobre una curva media y su derivada (índices de la rejilla fina).

    Picos: los dos máximos locales más altos. Valles: el mínimo entre los picos y los mínimos locales de
    cada lado. Inicio y fin de cada temporada: máximo y mínimo de la derivada entre valle y pico.
    """
    n = len(media)
    signo = np.sign(d1)
    maximos = np.where(np.diff(signo) == -2)[0] + 1
    minimos = np.where(np.diff(signo) == 2)[0] + 1
    picos = np.sort(maximos[np.argsort(media[maximos])[::-1][:2]])
    v1 = int(max([0] + [m for m in minimos if m < picos[0]]))
    v2 = int(picos[0] + np.argmin(media[picos[0] : picos[1] + 1]))
    v3 = int(min([n - 1] + [m for m in minimos if m > picos[1]]))
    inicio = [int(v1 + np.argmax(d1[v1 : picos[0] + 1])), int(v2 + np.argmax(d1[v2 : picos[1] + 1]))]
    fin = [int(picos[0] + np.argmin(d1[picos[0] : v2 + 1])), int(picos[1] + np.argmin(d1[picos[1] : v3 + 1]))]
    return {"picos": [int(p) for p in picos], "valles": [v1, v2, v3], "inicio": inicio, "fin": fin}


def temporadas_entrenamiento(D) -> Dict[str, object]:
    """Temporadas de la curva media de entrenamiento de las estaciones (mm/día), con la regla del EDA.

    Solo entrenamiento: ninguna figura de interpretación usa validación o prueba para ubicar las temporadas.
    Periodos del régimen bimodal: T1 y T2 = de inicio a fin de cada temporada de lluvias; Seca_media = entre
    las dos (mitad de año); Seca_dic_feb = el resto (fin y comienzo de año).
    """
    filas = ((D.pares.bloque == "entrenamiento") & (~D.pares.reservada)).to_numpy()
    media = D.phi_est_mm_fine @ D.coef_est_mm[filas].mean(0)
    T = temporadas(media, np.gradient(media, D.t_fine))
    t = D.t_fine
    T["periodos"] = {
        "T1": (t[T["inicio"][0]], t[T["fin"][0]]),
        "Seca_media": (t[T["fin"][0]], t[T["inicio"][1]]),
        "T2": (t[T["inicio"][1]], t[T["fin"][1]]),
    }
    T["media"] = media
    T["n_curvas"] = int(filas.sum())
    return T


def mascara_periodos(t: np.ndarray, T: Dict[str, object]) -> Dict[str, np.ndarray]:
    m = {nombre: (t >= a) & (t < b) for nombre, (a, b) in T["periodos"].items()}
    m["Seca_dic_feb"] = ~np.any(np.column_stack(list(m.values())), axis=1)
    return m


def sombrear_temporadas(ax, T: Dict[str, object]) -> None:
    for nombre in ("T1", "T2"):
        a, b = T["periodos"][nombre]
        ax.axvspan(a, b, color="tab:blue", alpha=0.07)


DIAS_BORDE = 15


def resumir_ventanas(modelo: nn.Module, tipo: str, D, T: Dict[str, object], datos) -> List[dict]:
    """Qué parte del año usa cada neurona funcional: DE, entre los pares de entrenamiento, de la contribución de cada
    periodo a z_k, relativa a la DE de z_k (Bordes_15d: los 15 primeros y últimos días, con efecto de borde)."""
    t = D.t_fine
    w = D.phi_sat_z_fine @ parametros_w(modelo).T
    x = D.phi_sat_z_fine @ datos.Cx[datos.idx_train].T
    ps = pesos_simpson(t)
    borde = (t < DIAS_BORDE) | (t > 365.0 - DIAS_BORDE)
    periodos = {nombre: mk & ~borde for nombre, mk in mascara_periodos(t, T).items()}
    periodos["Bordes_15d"] = borde
    filas = []
    for k in range(w.shape[1]):
        z = (ps * w[:, k]) @ x
        fila: Dict[str, object] = {"modelo": tipo, "neurona": k + 1}
        de = {nombre: float(np.std((ps * mk * w[:, k]) @ x) / np.std(z)) for nombre, mk in periodos.items()}
        for nombre, valor in de.items():
            fila[f"DE_rel_{nombre}"] = valor
        interior = {n: v for n, v in de.items() if n != "Bordes_15d"}
        fila["periodo_principal_interior"] = max(interior, key=interior.get)
        fila["corr_z_con_interior"] = float(np.corrcoef(z, (ps * ~borde * w[:, k]) @ x)[0, 1])
        filas.append(fila)
    return filas


def graficar_fpca(D, mu, B, valores, var_acum, M, destino: Path):
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    m = np.arange(1, 11)
    ax[0].bar(m, valores[:10] / valores.sum() * 100, color="#4C72B0", alpha=0.8)
    ax[0].plot(m, var_acum[:10] * 100, "o-", color="#C44E52", label="Acumulada")
    ax[0].axvline(M, ls="--", color="black", lw=1)
    ax[0].set_xlabel("Componente")
    ax[0].set_ylabel("Varianza explicada (%)")
    ax[0].set_title("FPCA de las curvas de estación")
    ax[0].legend(frameon=False)

    t = D.t_fine
    ax[1].plot(t, D.phi_est_z_fine @ mu, color="black", lw=2, label="Media")
    for j in range(min(3, M)):
        ax[1].plot(t, D.phi_est_z_fine @ B[:, j], lw=1.2, label=f"ψ{j + 1}")
    ax[1].set_xlabel("Día del año")
    ax[1].set_ylabel("z(t) = sqrt(lluvia)")
    ax[1].set_title(f"M = {M}: {var_acum[M - 1] * 100:.1f} % de varianza")
    ax[1].legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(destino / "Fig_FPCA.png", bbox_inches="tight")
    plt.close(fig)


def graficar_w(modelo, tipo, D, T: Dict[str, object], destino: Path):
    An = parametros_w(modelo)
    curvas = D.phi_sat_z_fine @ An.T
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    for k in range(An.shape[0]):
        ax.plot(D.t_fine, curvas[:, k], lw=1.2, label=f"w{k + 1}")
    ax.axhline(0, color="black", lw=0.6)
    sombrear_temporadas(ax, T)
    for a, b in [(0, DIAS_BORDE), (365 - DIAS_BORDE, 365)]:
        ax.axvspan(a, b, color="grey", alpha=0.15)
    for p in T["picos"]:
        ax.axvline(D.t_fine[p], color="tab:blue", lw=0.6, ls=":")
    ax.set_xlabel("Día del año (azul: temporadas de lluvia de entrenamiento; gris: 15 días de cada borde, donde la "
                  "base no es periódica y las curvas se leen con cautela)")
    ax.set_ylabel("w_k(t) (norma L2 = 1)")
    ax.set_title(f"Ventanas aprendidas por el codificador — cabeza {tipo}")
    ax.legend(frameon=False, ncol=4, fontsize=7)
    fig.tight_layout()
    fig.savefig(destino / f"Fig_wk_{tipo}.png", bbox_inches="tight")
    plt.close(fig)


def activacion_kan(modelo, i: int, x: np.ndarray, j: int) -> np.ndarray:
    """phi_ij(x) = w_b SiLU(x) + w_s sum_c c_c B_c(x): arista de la entrada i al nodo oculto j (primera capa)."""
    capa = modelo.cabeza.kan.act_fun[0]
    x_t = torch.tensor(np.asarray(x, dtype=np.float32)).unsqueeze(1)
    with torch.no_grad():
        spline = coef2curve(x_t, capa.grid[i : i + 1], capa.coef[i, j].reshape(1, 1, -1), capa.k).reshape(-1)
        base = capa.scale_base[i, j] * F.silu(x_t.reshape(-1))
        return (capa.scale_sp[i, j] * spline + base).numpy().astype(np.float64)


def nombres_entradas(K: int) -> List[str]:
    return [f"z{k + 1}" for k in range(K)] + ["altitud", "ONI"]


def a_unidades(i: int, x: np.ndarray, K: int, datos: DatosItem4) -> np.ndarray:
    """Entrada estandarizada -> unidades físicas (m para la altitud; ONI); las z_k quedan como están."""
    if i == K:
        return x * datos.alt_sd + datos.alt_mu
    if i == K + 1:
        return x * datos.oni_sd + datos.oni_mu
    return x


def analizar_activaciones_kan(modelo, K: int, datos: DatosItem4, destino: Path):
    """Activaciones de la primera capa de la KAN en el rango de entrenamiento: importancia (DE de phi_ij), R^2 de una
    recta y % de pares fuera de la rejilla. Control: la suma de las aristas reproduce la salida de la capa de pykan."""
    capa = modelo.cabeza.kan.act_fun[0]
    it = datos.idx_train
    with torch.no_grad():
        entrada_t = modelo(
            torch.tensor(datos.Cx[it], dtype=torch.float32),
            torch.tensor(datos.alt[it], dtype=torch.float32),
            torch.tensor(datos.oni[it], dtype=torch.float32),
        )["entrada"]
        h_pykan = capa(entrada_t)[0].numpy()
    entrada = entrada_t.numpy().astype(np.float64)
    n_in, n_out = capa.in_dim, capa.out_dim
    phi = np.stack([np.stack([activacion_kan(modelo, i, entrada[:, i], j) for j in range(n_out)], 1) for i in range(n_in)], 1)
    control = float(np.max(np.abs(phi.sum(1) - h_pykan)))
    nombres = nombres_entradas(K)
    filas = []
    for i in range(n_in):
        fuera = (entrada[:, i] < RANGO_KAN[0]) | (entrada[:, i] > RANGO_KAN[1])
        for j in range(n_out):
            y = phi[:, i, j]
            b, a = np.polyfit(entrada[:, i], y, 1)
            sst = np.sum((y - y.mean()) ** 2)
            r2 = 1 - np.sum((y - (a + b * entrada[:, i])) ** 2) / sst if sst > 0 else 1.0
            filas.append({"entrada": nombres[i], "i": i, "nodo": j + 1, "importancia": float(y.std()),
                          "R2_recta": float(r2), "valores_distintos": int(len(np.unique(entrada[:, i]))),
                          "pct_fuera_rejilla": float(100 * fuera.mean())})
    tabla = pd.DataFrame(filas)
    tabla["importancia_rel"] = tabla.importancia / tabla.importancia.max()

    rangos = {i: np.percentile(entrada[:, i], [1, 99]) for i in range(n_in)}
    rangos[K + 1] = np.array([entrada[:, K + 1].min(), entrada[:, K + 1].max()])   # ONI: los 15 años
    n_filas = int(math.ceil(n_in / 2))
    fig, axes = plt.subplots(n_filas, 2, figsize=(10, 2.2 * n_filas))
    for i, ax in enumerate(axes.ravel()[:n_in]):
        x = np.linspace(*rangos[i], 250)
        xu = a_unidades(i, x, K, datos)
        sub = tabla[tabla.i == i]
        principal = int(sub.loc[sub.importancia.idxmax(), "nodo"]) - 1
        for j in range(n_out):
            imp = float(sub.importancia_rel.iloc[j])
            ax.plot(xu, activacion_kan(modelo, i, x, j), lw=0.6 + 2.0 * imp, alpha=0.25 + 0.75 * imp,
                    color="tab:red" if j == principal else "grey")
        p5, p95 = np.percentile(entrada[:, i], [5, 95])
        ax.axvspan(*a_unidades(i, np.array([p5, p95]), K, datos), color="tab:blue", alpha=0.07)
        ax.axhline(0, color="black", lw=0.5)
        if i >= K:   # altitud y ONI: marcas en los valores distintos de entrenamiento
            valores = a_unidades(i, np.unique(entrada[:, i]), K, datos)
            ax.plot(valores, np.full(len(valores), ax.get_ylim()[0]), "|", color="black", ms=9)
        unidad = " (m)" if i == K else ""
        ax.set_title(f"Entrada {nombres[i]}{unidad}: importancia máx. {sub.importancia.max():.3f}", fontsize=8)
    for ax in axes.ravel()[n_in:]:
        ax.axis("off")
    fig.suptitle("Activaciones de la primera capa de la KAN (rojo: arista más importante de cada entrada; sombreado: "
                 "P5-P95; marcas: valores distintos en entrenamiento)", fontsize=9)
    fig.tight_layout()
    fig.savefig(destino / "Fig_activaciones_KAN.png", bbox_inches="tight")
    plt.close(fig)
    return tabla, entrada, rangos, control


CANDIDATOS_SIMBOLICOS = ["x", "x^2", "x^3", "exp", "tanh", "sin", "gaussian", "abs"]


def ajustar_simbolico(x: np.ndarray, y: np.ndarray) -> pd.DataFrame:
    """Ajuste simbólico de pykan: para cada función candidata f, fit_params busca y ~ c f(a x + b) + d.

    El R^2 se calcula aquí, 1 - SSE/SST de la curva ajustada (el de fit_params suma 1e-4 al denominador).
    Se elige la candidata de menor complejidad (la que asigna pykan) con R^2 >= 0,99; si ninguna llega, la
    de mayor R^2.
    """
    xt = torch.tensor(x, dtype=torch.float32)
    yt = torch.tensor(y, dtype=torch.float32)
    sst = float(np.sum((y - y.mean()) ** 2))
    filas = []
    for nombre in CANDIDATOS_SIMBOLICOS:
        f, _, complejidad, _ = SYMBOLIC_LIB[nombre]
        (a, b, c, d), _ = fit_params(xt, yt, f, verbose=False)
        ajuste = (c * f(a * xt + b) + d).numpy().astype(np.float64)
        filas.append({"funcion": nombre, "complejidad": complejidad, "a": float(a), "b": float(b),
                      "c": float(c), "d": float(d), "R2": 1 - float(np.sum((y - ajuste) ** 2)) / sst})
    tabla = pd.DataFrame(filas)
    buenas = tabla[tabla.R2 >= 0.99]
    elegida = buenas.sort_values(["complejidad", "R2"], ascending=[True, False]).index[0] if len(buenas) else tabla.R2.idxmax()
    tabla["elegida"] = tabla.index == elegida
    return tabla


def formula_simbolica(fila) -> str:
    if fila.funcion == "x":     # c (a x + b) + d es una recta: se escribe con su pendiente e intercepto
        return f"{fila.c * fila.a:.3g}·x {fila.c * fila.b + fila.d:+.3g}"
    return f"{fila.c:.3g}·{fila.funcion}({fila.a:.3g}·x {fila.b:+.3g}) {fila.d:+.3g}"


def simbolico_kan(modelo, K: int, datos: DatosItem4, tabla_act: pd.DataFrame, rangos, destino: Path):
    """Ajuste simbólico de la arista principal de altitud, de ONI y de la z_k más importante, más las dos
    siguientes aristas más importantes que no estén ya incluidas."""
    elegidas = []
    for grupo in (["altitud"], ["ONI"], [f"z{k + 1}" for k in range(K)]):
        sub = tabla_act[tabla_act.entrada.isin(grupo)]
        elegidas.append(sub.importancia.idxmax())
    for idx in tabla_act.sort_values("importancia", ascending=False).index:
        if len(elegidas) >= 5:
            break
        if idx not in elegidas:
            elegidas.append(idx)
    resumen, detalle = [], []
    fig, axes = plt.subplots(1, len(elegidas), figsize=(3.3 * len(elegidas), 3.2))
    for ax, idx in zip(np.atleast_1d(axes), elegidas):
        fila = tabla_act.loc[idx]
        i, j = int(fila.i), int(fila.nodo) - 1
        x = np.linspace(*rangos[i], 200)
        y = activacion_kan(modelo, i, x, j)
        tabla = ajustar_simbolico(x, y)
        tabla.insert(0, "arista", f"{fila.entrada}->h{j + 1}")
        detalle.append(tabla)
        mejor = tabla[tabla.elegida].iloc[0]
        f = SYMBOLIC_LIB[mejor.funcion][0]
        ajuste = (mejor.c * f(torch.tensor(mejor.a * x + mejor.b, dtype=torch.float32)) + mejor.d).numpy()
        xu = a_unidades(i, x, K, datos)
        ax.plot(xu, y, color="tab:red", lw=2, label="Activación KAN")
        ax.plot(xu, ajuste, "--", color="black", lw=1.2, label=f"{mejor.funcion}, R² = {mejor.R2:.4f}")
        ax.set_title(f"{fila.entrada} → nodo {j + 1}", fontsize=9)
        ax.set_xlabel(fila.entrada + (" (m)" if i == K else ""))
        ax.legend(frameon=False, fontsize=7)
        resumen.append({"arista": f"{fila.entrada}->h{j + 1}", "importancia": fila.importancia,
                        "R2_recta": fila.R2_recta, "funcion": mejor.funcion,
                        "expresion (x estandarizada)": formula_simbolica(mejor), "R2": mejor.R2})
    fig.suptitle("Ajuste simbólico (pykan fit_params) de activaciones de la primera capa", fontsize=9)
    fig.tight_layout()
    fig.savefig(destino / "Fig_ajuste_simbolico.png", bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(resumen), pd.concat(detalle, ignore_index=True)


NOMBRES_PERIODO = {"T1": "1.ª temporada de lluvias", "Seca_media": "Seca de mitad de año",
                   "T2": "2.ª temporada de lluvias", "Seca_dic_feb": "Seca dic-feb"}


def aportes_kan_por_condicion(modelo, K: int, datos: DatosItem4, D, T: Dict[str, object], destino: Path):
    """Aportes de la KAN por condición del año (Niña frente a Niño; año húmedo frente a seco), con pares de
    entrenamiento: diferencia de medias (a) de cada grupo de entradas sobre cada nodo oculto y (b) de cada nodo sobre
    la curva, por periodo. Descomposición exacta: h_j = s_j sum_i phi_ji(x_i) + b_j y u_m = s_m sum_j phi_mj(h_j) + b_m;
    los controles comprueban que los aportes suman la salida de la red."""
    kan = modelo.cabeza.kan
    assert np.array_equal(kan.input_id.long().numpy(), np.arange(K + 2)), "orden de entradas inesperado en pykan"
    it = datos.idx_train
    modelo.eval()
    with torch.no_grad():
        out = modelo(torch.tensor(datos.Cx[it], dtype=torch.float32), torch.tensor(datos.alt[it], dtype=torch.float32),
                     torch.tensor(datos.oni[it], dtype=torch.float32))
        _, _, post0, _ = kan.act_fun[0](out["entrada"])                         # pares x ocultos x entradas
        s0 = kan.node_scale[0] * kan.subnode_scale[0]
        b0 = kan.node_scale[0] * kan.subnode_bias[0] + kan.node_bias[0]
        e = post0 * s0[None, :, None]                                           # aporte de la entrada i al nodo j
        h = e.sum(2) + b0[None, :]
        _, _, post1, _ = kan.act_fun[1](h)                                      # pares x M x ocultos
        s1 = kan.node_scale[1] * kan.subnode_scale[1]
        b1 = kan.node_scale[1] * kan.subnode_bias[1] + kan.node_bias[1]
        a = post1 * s1[None, :, None]                                           # aporte del nodo j al puntaje m
        dec = modelo.decodificador
        psi = dec.basis.T @ dec.phi.T                                           # M x rejilla: psi_m(t)
        C = torch.einsum("nmj,mt->njt", a, psi)                                 # aporte del nodo j a z(t)
        control_u = float((a.sum(2) + b1[None, :] - out["u"]).abs().max())
        control_z = float((dec.mu @ dec.phi.T + b1 @ psi + C.sum(1) - out["z_fine"]).abs().max())
        z_obs = datos.Cy[it] @ dec.phi.numpy().T.astype(np.float64)
    e, C = e.numpy().astype(np.float64), C.numpy().astype(np.float64)
    n_oc = e.shape[1]

    # Condiciones del año
    df = datos.df.iloc[it]
    tot = np.maximum(datos.Xobs[it], 0.0) @ pesos_simpson(D.t_fine)
    anom = tot / pd.Series(tot).groupby(df.station_id.to_numpy()).transform("mean").to_numpy() - 1
    q1, q2 = np.quantile(anom, [1 / 3, 2 / 3])
    humedad = np.where(anom <= q1, "Seco", np.where(anom >= q2, "Húmedo", "Normal"))
    fase = df.fase.to_numpy()
    años = df.year.to_numpy()
    contrastes = [("ENSO", "Niña - Niño", fase == "Niña", fase == "Niño"),
                  ("Humedad", "Año húmedo - año seco (tercil del acumulado frente a la media de la estación)",
                   humedad == "Húmedo", humedad == "Seco")]

    # (a) entradas -> nodo oculto, por grupo: satélite (z_1..z_K), altitud, ONI y el total del nodo
    grupos_a = ["Satélite (z1-z%d)" % K, "Altitud", "ONI", "Total del nodo"]
    E_g = np.stack([e[:, :, :K].sum(2), e[:, :, K], e[:, :, K + 1], e.sum(2)], axis=2)        # n x ocultos x 4
    # (b) nodo oculto -> curva, promedio de cada periodo; más la suma de los nodos y la curva observada
    masc = mascara_periodos(D.t_fine, T)
    periodos = ["T1", "Seca_media", "T2", "Seca_dic_feb"]
    Cp = np.stack([C[:, :, masc[p]].mean(2) for p in periodos], axis=2)                       # n x ocultos x 4
    Cp = np.concatenate([Cp, Cp.sum(1, keepdims=True),
                         np.stack([z_obs[:, masc[p]].mean(1) for p in periodos], axis=1)[:, None, :]], axis=1)
    nodos_b = [f"H{j + 1}" for j in range(n_oc)] + ["Suma KAN", "Observado"]

    filas, resumen = [], []
    fig, axes = plt.subplots(2, 2, figsize=(15, 8.5))
    colores_a = ["tab:blue", "tab:brown", "tab:purple", "black"]
    colores_b = {"T1": "tab:blue", "Seca_media": "tab:orange", "T2": "tab:cyan", "Seca_dic_feb": "tab:red"}
    for fila_ax, (nombre, etiqueta, A, B) in zip(axes, contrastes):
        dif_a = E_g[A].mean(0) - E_g[B].mean(0)                                                # ocultos x 4
        dif_b = Cp[A].mean(0) - Cp[B].mean(0)                                                  # (ocultos + 2) x 4
        info = dict(contraste=nombre, pares_A=int(A.sum()), pares_B=int(B.sum()),
                    años_A=int(len(np.unique(años[A]))), años_B=int(len(np.unique(años[B]))))
        for j in range(n_oc):
            for g, gn in enumerate(grupos_a):
                filas.append({**info, "capa": "entradas -> nodo", "nodo": f"H{j + 1}", "grupo_o_periodo": gn,
                              "diferencia": float(dif_a[j, g])})
        for j, nn_ in enumerate(nodos_b):
            for p, per in enumerate(periodos):
                filas.append({**info, "capa": "nodo -> curva (raíz de mm/día)", "nodo": nn_,
                              "grupo_o_periodo": NOMBRES_PERIODO[per], "diferencia": float(dif_b[j, p])})
        for p, per in enumerate(periodos):
            suma_abs = float(np.abs(dif_b[:n_oc, p]).sum())
            resumen.append({**info, "periodo": NOMBRES_PERIODO[per], "suma_KAN": float(dif_b[n_oc, p]),
                            "observado": float(dif_b[n_oc + 1, p]), "suma_abs_nodos": suma_abs,
                            "razon_cancelacion": suma_abs / max(abs(float(dif_b[n_oc, p])), 1e-12),
                            "nodo_mayor": nodos_b[int(np.argmax(np.abs(dif_b[:n_oc, p])))]})
        ancho = 0.2
        ax = fila_ax[0]
        x = np.arange(n_oc)
        for g, gn in enumerate(grupos_a):
            ax.bar(x + (g - 1.5) * ancho, dif_a[:, g], ancho, label=gn, color=colores_a[g])
        ax.set_xticks(x, [f"H{j + 1}" for j in range(n_oc)])
        ax.set_ylabel("Cambio del nodo oculto")
        ax.set_title(f"(a) {etiqueta.split(' (')[0]}: aporte de las entradas a cada nodo\n"
                     f"({info['pares_A']} y {info['pares_B']} pares; {info['años_A']} y {info['años_B']} años)", fontsize=9)
        ax2 = fila_ax[1]
        x = np.arange(len(nodos_b))
        for p, per in enumerate(periodos):
            ax2.bar(x + (p - 1.5) * ancho, dif_b[:, p], ancho, label=NOMBRES_PERIODO[per], color=colores_b[per])
        ax2.axvline(n_oc - 0.5, color="grey", lw=0.8, ls=":")
        ax2.set_xticks(x, [n_.replace(" ", "\n") for n_ in nodos_b], fontsize=8)
        ax2.set_ylabel("Cambio de la curva, raíz de mm/día")
        ax2.set_title(f"(b) {etiqueta.split(' (')[0]}: aporte de cada nodo a la curva, por periodo", fontsize=9)
        for ax_ in (ax, ax2):
            ax_.axhline(0, color="black", lw=0.6)
            ax_.legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle("KAN por condición (entrenamiento): diferencia de medias entre las dos condiciones. Los nodos que se "
                 "compensan entre sí no se interpretan solos; la suma es el cambio de la predicción", fontsize=9)
    fig.tight_layout()
    fig.savefig(destino / "Fig_aportes_KAN_condicion.png", bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(filas), pd.DataFrame(resumen), control_u, control_z, (q1, q2)


def interpretar_capas_kan(modelo, K: int, datos: DatosItem4, D, T: Dict[str, object], destino: Path):
    """Flujo de la información por las dos capas de la KAN (no interviene en la selección de modelos): diagrama de
    pykan, efecto de cada nodo oculto en la curva (P10 y P90 frente a P50) y ablación (cada entrada, arista y nodo
    fijo en su media de entrenamiento; subida de la pérdida de validación)."""
    kan = modelo.cabeza.kan
    dec = modelo.decodificador
    nombres = nombres_entradas(K)
    it, iv = datos.idx_train, datos.idx_val
    modelo.eval()

    def entrada(idx):
        return modelo(torch.tensor(datos.Cx[idx], dtype=torch.float32), torch.tensor(datos.alt[idx], dtype=torch.float32),
                      torch.tensor(datos.oni[idx], dtype=torch.float32))["entrada"]

    with torch.no_grad():
        s0 = kan.node_scale[0] * kan.subnode_scale[0]
        b0 = kan.node_scale[0] * kan.subnode_bias[0] + kan.node_bias[0]
        s1 = kan.node_scale[1] * kan.subnode_scale[1]
        b1 = kan.node_scale[1] * kan.subnode_bias[1] + kan.node_bias[1]

        def capa0(x):
            return kan.act_fun[0](x)[2] * s0[None, :, None]           # pares x ocultos x entradas

        def capa1(h):
            return kan.act_fun[1](h)[2] * s1[None, :, None]           # pares x M x ocultos

        W = torch.tensor(D.W_est_z, dtype=torch.float32)
        cy = torch.tensor(datos.Cy[iv], dtype=torch.float32)

        def perdida(u):
            d = dec.mu + u @ dec.basis.T - cy
            return float(((d @ W) * d).sum(1).mean() / 365.0)

        E_ent, E_val = entrada(it), entrada(iv)
        e_ent = capa0(E_ent)
        h_ent = e_ent.sum(2) + b0
        a_ent = capa1(h_ent)
        e_val = capa0(E_val)
        h_val = e_val.sum(2) + b0
        a_val = capa1(h_val)
        base = perdida(a_val.sum(2) + b1)
        control = abs(base - perdida(modelo(torch.tensor(datos.Cx[iv], dtype=torch.float32),
                                            torch.tensor(datos.alt[iv], dtype=torch.float32),
                                            torch.tensor(datos.oni[iv], dtype=torch.float32))["u"]))

        def subida(p):  # % que sube la pérdida de validación
            return 100.0 * (p / base - 1.0)

        n_in, n_oc, M = e_ent.shape[2], e_ent.shape[1], a_ent.shape[1]
        filas = []
        for i in range(n_in):
            x = E_val.clone()
            x[:, i] = E_ent[:, i].mean()
            h = capa0(x).sum(2) + b0
            filas.append({"pieza": "entrada", "capa": 0, "origen": nombres[i], "destino": "todas",
                          "subida_perdida_pct": subida(perdida(capa1(h).sum(2) + b1))})
        for j in range(n_oc):
            for i in range(n_in):
                e = e_val.clone()
                e[:, j, i] = e_ent[:, j, i].mean()
                filas.append({"pieza": "arista", "capa": 1, "origen": nombres[i], "destino": f"H{j + 1}",
                              "subida_perdida_pct": subida(perdida(capa1(e.sum(2) + b0).sum(2) + b1)),
                              "DE_aporte_entrenamiento": float(e_ent[:, j, i].std())})
            h = h_val.clone()
            h[:, j] = h_ent[:, j].mean()
            filas.append({"pieza": "nodo", "capa": 1, "origen": f"H{j + 1}", "destino": "todas",
                          "subida_perdida_pct": subida(perdida(capa1(h).sum(2) + b1))})
            for m in range(M):
                a = a_val.clone()
                a[:, m, j] = a_ent[:, m, j].mean()
                filas.append({"pieza": "arista", "capa": 2, "origen": f"H{j + 1}", "destino": f"s{m + 1}",
                              "subida_perdida_pct": subida(perdida(a.sum(2) + b1)),
                              "DE_aporte_entrenamiento": float(a_ent[:, m, j].std())})
        tabla = pd.DataFrame(filas)

        # Qué hace cada nodo a la curva (raíz de mm/día) y qué entradas lo mueven
        reparto = (e_ent.std(0) / e_ent.std(0).sum(1, keepdims=True)).numpy()      # ocultos x entradas
        psi = (dec.basis.T @ dec.phi.T).numpy()                                      # M x rejilla
        curvas = {}
        for j in range(n_oc):
            q = torch.quantile(h_ent[:, j], torch.tensor([0.1, 0.5, 0.9]))
            H = h_ent.mean(0, keepdim=True).repeat(3, 1)
            H[:, j] = q
            aporte = capa1(H)[:, :, j].numpy() @ psi                                  # 3 x rejilla
            curvas[j] = (aporte[0] - aporte[1], aporte[2] - aporte[1])

    nodos = tabla[tabla.pieza == "nodo"].set_index("origen").subida_perdida_pct
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.2), sharex=True)
    t = D.t_fine
    for j, ax in enumerate(axes.ravel()[:n_oc]):
        bajo, alto = curvas[j]
        sombrear_temporadas(ax, T)
        ax.axhline(0, color="black", lw=0.6)
        ax.plot(t, alto, color="tab:blue", lw=1.6, label="nodo alto (P90)")
        ax.plot(t, bajo, color="tab:red", lw=1.2, ls="--", label="nodo bajo (P10)")
        top = np.argsort(reparto[j])[::-1][:2]
        ax.set_title(f"H{j + 1} (ablación +{nodos[f'H{j + 1}']:.1f} %) ← " +
                     ", ".join(f"{nombres[i]} {100 * reparto[j, i]:.0f} %" for i in top), fontsize=8)
        if j % 4 == 0:
            ax.set_ylabel("Cambio de la curva\n(raíz de mm/día)")
        if j >= 4:
            ax.set_xlabel("Día del año")
    axes.ravel()[0].legend(frameon=False, fontsize=7)
    fig.suptitle("Qué hace cada nodo oculto de la KAN a la curva anual (frente a su valor mediano) y qué entradas lo "
                 "mueven (reparto de la DE de la capa 1; entrenamiento)", fontsize=9)
    fig.tight_layout()
    fig.savefig(destino / "Fig_nodos_KAN.png", bbox_inches="tight")
    plt.close(fig)

    # Ablación: matrices de aristas y barras de entradas y nodos
    m1 = tabla[(tabla.pieza == "arista") & (tabla.capa == 1)].pivot(index="origen", columns="destino",
                                                                      values="subida_perdida_pct").loc[nombres]
    m1 = m1[[f"H{j + 1}" for j in range(n_oc)]]
    m2 = tabla[(tabla.pieza == "arista") & (tabla.capa == 2)].pivot(index="origen", columns="destino",
                                                                      values="subida_perdida_pct")
    m2 = m2.loc[[f"H{j + 1}" for j in range(n_oc)], [f"s{m + 1}" for m in range(M)]]
    ent = tabla[tabla.pieza == "entrada"].set_index("origen").subida_perdida_pct.loc[nombres]
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.4), gridspec_kw={"width_ratios": [1.1, 0.9, 0.7, 0.7]})
    tope = max(float(np.nanmax(m1.to_numpy())), float(np.nanmax(m2.to_numpy())), 1e-6)
    for ax, mat, titulo in [(axes[0], m1, "Capa 1: entrada → nodo oculto"), (axes[1], m2, "Capa 2: nodo → puntaje FPCA")]:
        im = ax.imshow(mat.to_numpy(), cmap="Blues", vmin=0, vmax=tope, aspect="auto")
        ax.set_xticks(range(mat.shape[1]), mat.columns)
        ax.set_yticks(range(mat.shape[0]), mat.index)
        for (r, c), v in np.ndenumerate(mat.to_numpy()):
            ax.text(c, r, f"{v:.1f}", ha="center", va="center", fontsize=6,
                    color="white" if v > 0.6 * tope else "black")
        ax.set_title(titulo, fontsize=9)
    fig.colorbar(im, ax=axes[:2], shrink=0.8, label="Subida de la pérdida de validación (%)")
    axes[2].barh(range(n_in)[::-1], ent.to_numpy(), color="tab:purple")
    axes[2].set_yticks(range(n_in)[::-1], nombres)
    axes[2].set_title("Entrada reemplazada por su media", fontsize=9)
    axes[3].barh(range(n_oc)[::-1], nodos.loc[[f"H{j + 1}" for j in range(n_oc)]].to_numpy(), color="tab:green")
    axes[3].set_yticks(range(n_oc)[::-1], [f"H{j + 1}" for j in range(n_oc)])
    axes[3].set_title("Nodo oculto reemplazado por su media", fontsize=9)
    for ax in axes[2:]:
        ax.set_xlabel("Subida de la pérdida (%)")
    fig.suptitle("Ablación de la KAN: cuánto empeora la validación si una pieza deja de variar (se fija en su media "
                 "de entrenamiento)", fontsize=9)
    fig.savefig(destino / "Fig_ablacion_KAN.png", bbox_inches="tight")
    plt.close(fig)

    # Diagrama de la red con pykan: forma de cada activación e importancia hacia atrás de cada arista
    kan.save_act = True
    try:
        with torch.no_grad():
            kan(E_ent)
        with tempfile.TemporaryDirectory() as carpeta:
            kan.plot(folder=carpeta, beta=3, metric="backward", scale=1.0,
                     in_vars=["alt" if n == "altitud" else n for n in nombres],
                     out_vars=[f"s{m + 1}" for m in range(M)], title=f"KAN [{n_in}, {n_oc}, {M}] (pykan, entrenamiento)")
            plt.savefig(destino / "Fig_red_KAN.png", bbox_inches="tight", dpi=150)
            plt.close("all")
    finally:
        kan.save_act = False
    return tabla, base, control, reparto


def colinealidad_z(mejor_por_modelo, datos: DatosItem4) -> pd.DataFrame:
    """Colinealidad de las salidas del codificador (entrenamiento; descriptivo): mayor |r| entre z_k, |r| medio entre
    z1, z2 y z3, número de condición y, en la cabeza lineal, razón de cancelación de cada puntaje del FPCA."""
    it = datos.idx_train
    filas = []
    for tipo in ("kan", "lineal"):
        modelo = mejor_por_modelo[tipo].modelo
        modelo.eval()
        with torch.no_grad():
            u = modelo(torch.tensor(datos.Cx[it], dtype=torch.float32), torch.tensor(datos.alt[it], dtype=torch.float32),
                       torch.tensor(datos.oni[it], dtype=torch.float32))["entrada"].numpy().astype(np.float64)
        K = u.shape[1] - 2
        fuera = np.abs(np.corrcoef(u[:, :K].T) - np.eye(K))
        i, j = np.unravel_index(np.argmax(fuera), fuera.shape)
        fila = {"cabeza": tipo, "K": K, "max_abs_r": float(fuera[i, j]), "par_max": f"z{i + 1}-z{j + 1}",
                "abs_r_medio_z1_z3": float(np.mean([fuera[0, 1], fuera[0, 2], fuera[1, 2]])),
                "condicion_cor_z": float(np.linalg.cond(np.corrcoef(u[:, :K].T)))}
        if tipo == "lineal":
            W = modelo.cabeza.weight.detach().numpy().astype(np.float64)       # M x (K + 2)
            c = u[:, None, :] * W[None, :, :]                                 # pares x M x entradas
            razon = c.std(0).sum(1) / np.maximum(c.sum(2).std(0), 1e-12)
            fila.update(cancelacion_min=float(razon.min()), cancelacion_mediana=float(np.median(razon)),
                        cancelacion_max=float(razon.max()))
        filas.append(fila)
    return pd.DataFrame(filas)


def _predecir_mm(modelo, cx, alt, oni) -> Tuple[np.ndarray, np.ndarray]:
    with torch.no_grad():
        z = modelo(cx, alt, oni)["z_fine"].numpy().astype(np.float64)
    return z, np.maximum(z, 0.0) ** 2


def contraste_cabezas(resultados, mejor_por_modelo, datos: DatosItem4, D, T: Dict[str, object], destino: Path):
    """Contraste de la KAN y la MLP con la cabeza lineal por el ONI y la altitud (entrenamiento): dependencia parcial
    (% frente a ONI = 0 o a la altitud mediana), coeficiente local d zhat(t) / d v y referencia empírica del cociente
    estación / satélite."""
    it = datos.idx_train
    cx = torch.tensor(datos.Cx[it], dtype=torch.float32)
    alt0 = datos.alt[it]
    oni0 = datos.oni[it]
    ps = pesos_simpson(D.t_fine)
    masc = mascara_periodos(D.t_fine, T)

    def lluvia_media(mm):           # mm/día medios de cada periodo y del año, promediados entre pares
        salida = {"Anual": float(np.mean(mm @ ps) / 365.0)}
        for nombre, mk in masc.items():
            salida[nombre] = float(np.mean(mm[:, mk]))
        return salida

    def correr(modelo, variable, v):
        alt = np.full(len(it), v) if variable == "altitud" else alt0
        oni = np.full(len(it), v) if variable == "ONI" else oni0
        return _predecir_mm(modelo, cx, torch.tensor(alt, dtype=torch.float32), torch.tensor(oni, dtype=torch.float32))

    alt_ref = float(np.median(alt0))
    niveles = {"ONI": [-0.75, 0.0, 0.75], "altitud": [float(np.percentile(alt0, 10)), alt_ref, float(np.percentile(alt0, 90))]}
    rejillas = {
        "ONI": (np.unique(np.r_[np.linspace(oni0.min(), oni0.max(), 21), niveles["ONI"]]), 0.0, 0.1),
        "altitud": (np.unique(np.r_[np.linspace(*np.percentile(alt0, [2, 98]), 21), niveles["altitud"]]), alt_ref, 50.0),
    }
    filas_pd = []
    for variable, (rejilla, ref, _) in rejillas.items():
        for tipo, lista in resultados.items():
            for k, res in enumerate(lista):
                base = lluvia_media(correr(res.modelo, variable, ref)[1])
                for v in rejilla:
                    lm = lluvia_media(correr(res.modelo, variable, v)[1])
                    for periodo, valor in lm.items():
                        filas_pd.append({"variable": variable, "cabeza": tipo, "semilla_n": k, "elegida": res is mejor_por_modelo[tipo],
                                         "valor": float(v), "periodo": periodo, "cambio_pct": 100 * (valor / base[periodo] - 1)})
    tabla_pd = pd.DataFrame(filas_pd)

    # Referencia empírica en entrenamiento: cociente estación / satélite (lluvia anual)
    cociente = (datos.Xobs[it] @ ps) / (datos.Xsat[it] @ ps)
    fase = datos.df.iloc[it].fase.to_numpy()
    años = datos.df.iloc[it].year.to_numpy()
    base_neutro = np.median(cociente[fase == "Neutro"])
    ref_años = pd.DataFrame({"year": np.unique(años)})
    ref_años["valor"] = [float(oni0[años == a][0]) for a in ref_años.year]
    ref_años["fase"] = [fase[años == a][0] for a in ref_años.year]
    ref_años["pares"] = [int((años == a).sum()) for a in ref_años.year]
    ref_años["cambio_pct"] = [100 * (np.median(cociente[años == a]) / base_neutro - 1) for a in ref_años.year]
    ref_oni = pd.DataFrame({"grupo": ["Niña", "Neutro", "Niño"]})
    ref_oni["valor"] = [float(oni0[fase == g].mean()) for g in ref_oni.grupo]
    ref_oni["cambio_pct"] = [100 * (np.median(cociente[fase == g]) / base_neutro - 1) for g in ref_oni.grupo]
    base_alt = np.median(cociente[(alt0 >= 1200) & (alt0 < 2000)])
    deciles = pd.qcut(alt0, 10)
    ref_alt = (pd.DataFrame({"alt": alt0, "c": cociente, "q": deciles}).groupby("q", observed=True)
               .agg(valor=("alt", "median"), pares=("c", "size"), cociente=("c", "median")).reset_index())
    ref_alt["cambio_pct"] = 100 * (ref_alt.cociente / base_alt - 1)
    ref_alt["q"] = ref_alt.q.astype(str)
    bajo_500 = alt0 < 500
    nota_500 = (f"<500 m (fuera del rango graficado): {int(bajo_500.sum())} pares, "
                f"{100 * (np.median(cociente[bajo_500]) / base_alt - 1):+.0f} %") if bajo_500.any() else ""

    # Coeficiente local en el año, con la semilla elegida de cada cabeza
    escala = {"ONI": 1.0, "altitud": 100.0}            # por unidad de ONI y por 100 m
    locales = {}
    for variable, (_, _, delta) in rejillas.items():
        for tipo, res in mejor_por_modelo.items():
            for v in (niveles[variable] if tipo == "kan" else niveles[variable][1:2]):
                z_mas = correr(res.modelo, variable, v + delta)[0]
                z_menos = correr(res.modelo, variable, v - delta)[0]
                locales[(variable, tipo, v)] = escala[variable] * (z_mas - z_menos).mean(0) / (2 * delta)
    # Control: en la cabeza lineal el coeficiente local es exactamente gamma(t) de sus pesos
    lin = mejor_por_modelo["lineal"].modelo
    Kn = lin.encoder.A.shape[0]
    Wc = lin.cabeza.weight.detach().numpy().astype(np.float64)
    phi_t = D.phi_est_z_fine @ lin.decodificador.basis.numpy().astype(np.float64)
    control_gamma = float(np.max(np.abs(locales[("ONI", "lineal", 0.0)] - phi_t @ Wc[:, Kn + 1] / datos.oni_sd)))

    colores = {"lineal": "black", "mlp": "tab:orange", "kan": "tab:green"}
    tonos_kan = ["tab:blue", "tab:green", "tab:red"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.8))
    for fila_ax, variable in zip(axes, ["ONI", "altitud"]):
        ax = fila_ax[0]
        sub = tabla_pd[(tabla_pd.variable == variable) & (tabla_pd.periodo == "Anual")]
        for (tipo, k), g in sub.groupby(["cabeza", "semilla_n"]):
            el = bool(g.elegida.iloc[0])
            ax.plot(g.valor, g.cambio_pct, color=colores[tipo], lw=2.2 if el else 0.8, alpha=1 if el else 0.45,
                    label=tipo if el else None)
        if variable == "ONI":
            ax.plot(ref_años.valor, ref_años.cambio_pct, "o", color="tab:purple", ms=4, label="cada año (datos)")
            for _, r in ref_años.iterrows():
                ax.annotate(str(int(r.year))[2:], (r.valor, r.cambio_pct), fontsize=6, color="tab:purple",
                            xytext=(2, 2), textcoords="offset points")
            ax.plot(ref_oni.valor, ref_oni.cambio_pct, "s", color="tab:purple", ms=7, mfc="none", label="cada fase (datos)")
        else:
            ax.plot(ref_alt.valor, ref_alt.cambio_pct, "s", color="tab:purple", ms=6, label="cada decil (datos)")
            ax.text(0.02, 0.97, nota_500, transform=ax.transAxes, fontsize=7, va="top")
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xlabel(variable + (" (m)" if variable == "altitud" else ""))
        ax.set_ylabel("Cambio de la lluvia anual (%)")
        ax.set_title(f"Dependencia parcial: {variable} (gruesa = semilla elegida; puntos: estación/satélite)", fontsize=9)
        ax.legend(frameon=False, fontsize=7, loc="lower left" if variable == "ONI" else "lower right")
        ax = fila_ax[1]
        for (var, tipo, v), curva in locales.items():
            if var != variable:
                continue
            if tipo == "kan":
                ax.plot(D.t_fine, curva, color=tonos_kan[niveles[variable].index(v)], lw=1.4, label=f"kan en {v:.4g}")
            else:
                ax.plot(D.t_fine, curva, ":" if tipo == "mlp" else "-", color=colores[tipo], lw=1.4,
                        label="lineal: gamma(t)" if tipo == "lineal" else f"mlp en {v:.4g}")
        sombrear_temporadas(ax, T)
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xlabel("Día del año")
        ax.set_ylabel("d zhat(t) / d " + ("ONI" if variable == "ONI" else "(100 m)"))
        ax.set_title(f"Coeficiente local de {variable} a lo largo del año", fontsize=9)
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(destino / "Fig_contraste_cabezas.png", bbox_inches="tight")
    plt.close(fig)
    return tabla_pd, ref_oni, ref_años, ref_alt, nota_500, control_gamma, niveles


def validacion_por_estaciones(config, K, h_mlp, D, datos: DatosItem4, mu, Bs, max_epochs: int,
                              pliegues: int = 3, por_pliegue: int = 30, semillas=(0, 1)) -> pd.DataFrame:
    """Pérdida de validación en estaciones no vistas: en cada pliegue se sacan 30 estaciones de entrenamiento (sorteo
    fijo), se entrena con las demás y se mide la validación de las sacadas. Solo entrenamiento y validación."""
    est = datos.df.station_id.to_numpy()
    orden = np.random.default_rng(2026).permutation(np.unique(est[datos.idx_train]))
    filas = []
    for p in range(pliegues):
        fuera = np.isin(est, orden[p * por_pliegue : (p + 1) * por_pliegue])
        d2 = dataclasses.replace(datos, idx_train=datos.idx_train[~fuera[datos.idx_train]],
                                 idx_val=datos.idx_val[~fuera[datos.idx_val]])
        idx_nuevas = datos.idx_val[fuera[datos.idx_val]]
        for tipo, (lam, lr) in config.items():
            for semilla in semillas:
                res = entrenar(tipo, semilla, K, lam, lr, D, d2, mu, Bs, h_mlp if tipo == "mlp" else None, max_epochs)
                perdidas = {}
                for nombre, idx in [("vistas", d2.idx_val), ("nuevas", idx_nuevas)]:
                    c, _, _ = predecir(res.modelo, datos, idx)
                    dif = c - datos.Cy[idx]
                    perdidas[nombre] = float(np.mean(np.einsum("ij,jk,ik->i", dif, D.W_est_z, dif)) / 365.0)
                filas.append({"pliegue": p, "cabeza": tipo, "semilla": semilla, "pares_nuevas": len(idx_nuevas),
                              "perdida_val_estaciones_vistas": perdidas["vistas"],
                              "perdida_val_estaciones_nuevas": perdidas["nuevas"]})
    return pd.DataFrame(filas)


def nombres_estaciones() -> pd.DataFrame:
    """Nombre y cuenca de cada estación, del catálogo del curso (CVC, Latin-1)."""
    ruta = RUTA_PROYECTO / "Datos_clase" / "CVC_Wheather_Stations_Long_Lat.csv"
    cat = pd.read_csv(ruta, encoding="latin-1", dtype={"CODIGO": str})
    cat = cat.rename(columns={"CODIGO": "station_id", "STATION_NAME": "nombre", "BASIN": "cuenca"})
    return cat[["station_id", "nombre", "cuenca"]].drop_duplicates("station_id")


def graficar_predicciones(df, datos: DatosItem4, preds: Dict[str, np.ndarray], idx: np.ndarray, D, T, destino: Path,
                          n_peores: int = 6, sufijo: str = ""):
    """Ejemplos (mejor y típico de la KAN por fase ENSO) y peores casos (mayor L2 de la KAN) del bloque predicho, con
    la estación, el satélite, las tres cabezas y la media de entrenamiento. Devuelve la tabla de peores casos."""
    obs, sat, t = datos.Xobs[idx], datos.Xsat[idx], D.t_fine
    obs_0, sat_0 = np.maximum(obs, 0.0), np.maximum(sat, 0.0)
    clim = np.maximum(datos.Xobs[datos.idx_train], 0.0).mean(axis=0)
    ps = pesos_simpson(t)
    l2 = {tipo: l2_relativo(p, obs) for tipo, p in preds.items()}
    l2_sat = l2_relativo(sat, obs)
    l2_clim = l2_relativo(np.broadcast_to(clim, obs.shape), obs)
    colores = {"lineal": "black", "mlp": "tab:orange", "kan": "tab:green"}
    catalogo = nombres_estaciones().set_index("station_id")

    def nombre(sid):
        return catalogo.nombre.get(str(sid), "") if str(sid) in catalogo.index else ""

    def panel(ax, j, titulo):
        ax.plot(t, obs_0[j], color="tab:blue", lw=1.8, label="Estación")
        ax.plot(t, sat_0[j], color="grey", lw=1.2, ls="--", label="Satélite")
        ax.plot(t, clim, color="darkred", lw=1.0, ls="-.", label="Media de entrenamiento")
        for tipo, p in preds.items():
            ax.plot(t, p[j], color=colores[tipo], lw=1.0 if tipo != "kan" else 1.6, ls=":" if tipo == "mlp" else "-", label=tipo)
        sombrear_temporadas(ax, T)
        ax.set_title(titulo, fontsize=8)
        ax.set_xlabel("Día del año")
        ax.set_ylabel("mm/día")
        ax.legend(frameon=False, fontsize=6, ncol=3)

    fases = df.iloc[idx].fase.to_numpy()
    ejemplos = []
    for f in ["Niña", "Neutro", "Niño"]:
        en_fase = np.where(fases == f)[0]
        if len(en_fase) == 0:
            continue
        orden = en_fase[np.argsort(l2["kan"][en_fase], kind="stable")]
        ejemplos.append((f, [("Mejor", int(orden[0])), ("Típico (mediana)", int(orden[len(orden) // 2]))]))
    fig, axes = plt.subplots(len(ejemplos), 2, figsize=(14, 3.0 * len(ejemplos)), squeeze=False)
    for fila_ax, (f, casos) in zip(axes, ejemplos):
        for ax, (cual, j) in zip(fila_ax, casos):
            fila = df.iloc[idx[j]]
            panel(ax, j, f"{cual}, {f}: {fila.station_id} {nombre(fila.station_id)} — {int(fila.year)} — "
                         f"L2 KAN {l2['kan'][j]:.3f} (satélite {l2_sat[j]:.3f}; media {l2_clim[j]:.3f})")
    fig.tight_layout()
    fig.savefig(destino / f"Fig_predicciones_ejemplo{sufijo}.png", bbox_inches="tight")
    plt.close(fig)

    peores = np.argsort(l2["kan"])[::-1][:n_peores]
    fig, axes = plt.subplots(int(math.ceil(n_peores / 2)), 2, figsize=(11, 2.9 * math.ceil(n_peores / 2)))
    filas = []
    for ax, j in zip(axes.ravel(), peores):
        fila = df.iloc[idx[j]]
        panel(ax, j, f"{fila.station_id} {nombre(fila.station_id)} — {int(fila.year)} — {fila.fase} — "
                     f"{fila.altitud_usada:.0f} m — L2 KAN {l2['kan'][j]:.2f}")
        sid = str(fila.station_id)
        filas.append({"station_id": sid, "nombre": nombre(sid),
                      "cuenca": catalogo.cuenca.get(sid, "") if sid in catalogo.index else "",
                      "year": int(fila.year), "fase": fila.fase,
                      "vertiente": "Pacífico" if sid.startswith("5") else "Cauca",
                      "altitud": float(fila.altitud_usada), "L2_satelite": float(l2_sat[j]),
                      **{f"L2_{tipo}": float(v[j]) for tipo, v in l2.items()},
                      "mm_anual_estacion": float(obs[j] @ ps), "mm_anual_satelite": float(sat[j] @ ps),
                      "mm_anual_kan": float(preds["kan"][j] @ ps)})
    fig.tight_layout()
    fig.savefig(destino / f"Fig_peores_casos{sufijo}.png", bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(filas)


def cancelacion_por_capa(modelo, tipo: str, datos: DatosItem4) -> List[float]:
    """Por capa de la cabeza: mediana entre los nodos de salida de sum_i DE(aporte_i) / DE(sum_i aporte_i), con pares
    de entrenamiento (1 = los aportes se suman sin anularse; > 1 = se compensan)."""
    it = datos.idx_train
    modelo.eval()
    with torch.no_grad():
        out = modelo(torch.tensor(datos.Cx[it], dtype=torch.float32), torch.tensor(datos.alt[it], dtype=torch.float32),
                     torch.tensor(datos.oni[it], dtype=torch.float32))
        x, razones = out["entrada"], []

        def razon(c):  # c: pares x salidas x entradas
            de_partes = c.std(0).sum(1)
            de_suma = c.sum(2).std(0).clamp_min(1e-12)
            return float(torch.median(de_partes / de_suma))

        if tipo == "kan":
            kan = modelo.cabeza.kan
            for l in range(len(kan.act_fun)):
                _, _, post, _ = kan.act_fun[l](x)
                s = kan.node_scale[l] * kan.subnode_scale[l]
                c = post * s[None, :, None]
                razones.append(razon(c))
                x = c.sum(2) + (kan.node_scale[l] * kan.subnode_bias[l] + kan.node_bias[l])[None, :]
        else:
            capas = [m for m in modelo.cabeza if isinstance(m, nn.Linear)]
            for k, capa in enumerate(capas):
                razones.append(razon(capa.weight[None, :, :] * x[:, None, :]))
                x = capa(x)
                if k < len(capas) - 1:
                    x = F.silu(x)
        if float((x - out["u"]).abs().max()) > 1e-4:
            raise RuntimeError("La descomposición por capas no reproduce la salida de la cabeza")
    return razones


def experimento_arquitectura(semillas: List[int], max_epochs: int) -> None:
    """Experimento de arquitectura de la cabeza (solo validación; no interviene en la selección de modelos).

    Con el codificador (K y lambda de la corrida principal), el FPCA, los datos y la parada temprana fijos, cada
    arquitectura se entrena con las semillas y cada lr de LR_GRID, y se elige el lr de menor pérdida de validación
    media. Se miden la pérdida, el L2 en mm, el nivel, la forma, la captura de los contrastes húmedo - seco y
    Niña - Niño y la cancelación por capa.
    """
    salida = RUTA_SALIDA / "arquitectura"
    salida.mkdir(parents=True, exist_ok=True)
    print("=== Experimento de arquitectura (solo validación) ===")
    archivo_busqueda = RUTA_SALIDA / "busqueda_hiperparametros_resumen.csv"
    if not archivo_busqueda.exists():
        raise FileNotFoundError(f"Falta {archivo_busqueda}: correr antes este archivo sin --arquitectura")
    elegidos = pd.read_csv(archivo_busqueda).query("elegido").set_index("cabeza")
    K = int(elegidos.K_neuronas.iloc[0])
    D = cargar_datos_funcionales(verificar=True)
    datos = preparar_datos(D)
    mu, Bm, valores, _, _ = ajustar_fpca(datos.Cy[datos.idx_train], D.W_est_z, M_FPCA)
    Bs = Bm * np.sqrt(valores[:M_FPCA])[None, :]
    d = K + 2
    h_mlp = emparejar_parametros(d, M_FPCA, 0)[0]

    # Referencias de validación (mm en los centros de día; la media es la de entrenamiento, como en la V2)
    it, iv = datos.idx_train, datos.idx_val
    O_ent = np.maximum(datos.Xobs[it][:, CENTROS_DIA], 0.0)
    O_val = np.maximum(datos.Xobs[iv][:, CENTROS_DIA], 0.0)
    mu_mm = O_ent.mean(0)
    z_obs_val = (datos.Cy[iv] @ D.phi_est_z_fine.T).mean(1)            # media anual de la raíz observada
    est_ent, est_val = datos.df.station_id.to_numpy()[it], datos.df.station_id.to_numpy()[iv]
    media_est = pd.Series(O_ent.sum(1)).groupby(est_ent).mean()
    con_media = np.isin(est_val, media_est.index)
    anom = np.full(len(iv), np.nan)
    anom[con_media] = O_val.sum(1)[con_media] / media_est.loc[est_val[con_media]].to_numpy() - 1
    q1, q2 = np.nanquantile(anom, [1 / 3, 2 / 3])
    fase_val = datos.df.fase.to_numpy()[iv]
    contrastes = {"humedo_seco": (anom >= q2, anom <= q1), "nina_nino": (fase_val == "Niña", fase_val == "Niño")}

    def medir(modelo) -> Dict[str, float]:
        _, z, _ = predecir(modelo, datos, iv)
        P = np.maximum(z[:, CENTROS_DIA], 0.0) ** 2
        Ap, Ao = P - mu_mm, O_val - mu_mm
        r_anom = [np.corrcoef(Ap[i], Ao[i])[0, 1] for i in range(len(iv))]
        tot_p, tot_o = P.sum(1), O_val.sum(1)
        zp = z.mean(1)
        m = {"L2_val_mediana": float(np.median(np.sqrt(((P - O_val) ** 2).sum(1) / (O_val ** 2).sum(1)))),
             "amplitud": float(np.median(np.sqrt((Ap ** 2).sum(1) / (Ao ** 2).sum(1)))),
             "r_anomalia": float(np.nanmedian(r_anom)), "r_totales": float(np.corrcoef(tot_p, tot_o)[0, 1]),
             "DE_totales": float(tot_p.std() / tot_o.std())}
        for nombre, (A, B) in contrastes.items():
            obs = z_obs_val[A].mean() - z_obs_val[B].mean()
            m[f"captura_{nombre}"] = float((zp[A].mean() - zp[B].mean()) / obs) if A.any() and B.any() else float("nan")
        return m

    variantes = [
        (f"KAN [{d}, 8, 6] (principal)", "kan", {}),
        (f"KAN [{d}, 4, 6]", "kan", {"ocultas": (4,)}),
        (f"KAN [{d}, 16, 6]", "kan", {"ocultas": (16,)}),
        (f"KAN [{d}, 8, 8, 6]", "kan", {"ocultas": (8, 8)}),
        (f"KAN [{d}, 4, 4, 6]", "kan", {"ocultas": (4, 4)}),
        (f"KAN [{d}, 8, 6], rejilla 3", "kan", {"grid": 3}),
        (f"KAN [{d}, 8, 6], rejilla 10", "kan", {"grid": 10}),
        (f"MLP [{d}, {h_mlp}, 6] (principal)", "mlp", {}),
        (f"MLP [{d}, 26, 26, 6]", "mlp", {"ocultas": (26, 26)}),
    ]
    filas = []
    for nombre, tipo, arq in variantes:
        lam = float(elegidos.loc[tipo, "lambda_rugosidad"])
        prueba_lr = {}
        for lr in LR_GRID:
            prueba_lr[lr] = [entrenar(tipo, semilla, K, lam, lr, D, datos, mu, Bs, h_mlp, max_epochs, arq)
                             for semilla in semillas]
        lr = min(prueba_lr, key=lambda v: np.mean([r.mejor_val for r in prueba_lr[v]]))
        for semilla, res in zip(semillas, prueba_lr[lr]):
            canc = cancelacion_por_capa(res.modelo, tipo, datos)
            filas.append({"arquitectura": nombre, "cabeza": tipo, "semilla": semilla, "lambda": lam, "lr": lr,
                          "parametros": contar_parametros(res.modelo.cabeza), "val_loss": res.mejor_val,
                          "epocas": res.epocas, "segundos": res.segundos, **medir(res.modelo),
                          "cancelacion_capas": " / ".join(f"{c:.2f}" for c in canc),
                          "cancelacion_max": max(canc)})
    tabla = pd.DataFrame(filas)
    tabla.to_csv(salida / "arquitectura_semillas.csv", index=False)
    cols = ["val_loss", "L2_val_mediana", "amplitud", "r_anomalia", "r_totales", "DE_totales", "captura_humedo_seco",
            "captura_nina_nino", "cancelacion_max"]
    resumen = tabla.groupby("arquitectura", sort=False).agg(
        cabeza=("cabeza", "first"), parametros=("parametros", "first"), lr=("lr", "first"),
        **{f"{c}_media": (c, "mean") for c in cols}, val_loss_DE=("val_loss", "std"), L2_DE=("L2_val_mediana", "std"),
    ).reset_index()
    resumen.to_csv(salida / "arquitectura_resumen.csv", index=False)
    pd.set_option("display.width", 250)
    imprimir_tabla("Arquitecturas (validación; media de las semillas)", resumen.round(4))
    imprimir_tabla("Cancelación por capa (entrenamiento), por semilla",
                   tabla[["arquitectura", "semilla", "cancelacion_capas"]])
    print(f"Terciles del acumulado en validación frente a la media de entrenamiento de la estación: seco <= "
          f"{100 * q1:+.1f} %, húmedo >= {100 * q2:+.1f} %; pares Niña / Niño en validación: "
          f"{int(contrastes['nina_nino'][0].sum())} / {int(contrastes['nina_nino'][1].sum())}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 0.45 * len(variantes) + 2.2), sharey=True)
    y = np.arange(len(variantes))[::-1]
    nombres = [v[0] for v in variantes]
    for ax, (col, titulo) in zip(axes, [("val_loss", "Pérdida de validación (menor es mejor)"),
                                        ("captura_humedo_seco", "Captura del contraste húmedo − seco"),
                                        ("cancelacion_max", "Cancelación de la capa más cancelada")]):
        for k, nombre in enumerate(nombres):
            v = tabla.loc[tabla.arquitectura == nombre, col].to_numpy()
            color = "tab:green" if "KAN" in nombre else "tab:orange"
            ax.plot(v, np.full(len(v), y[k]), "o", color=color, alpha=0.55, ms=5)
            ax.plot([v.mean()], [y[k]], "|", color="black", ms=14, mew=2)
        ax.set_title(titulo, fontsize=9)
        if col == "captura_humedo_seco":
            ax.axvline(1, color="grey", lw=0.8, ls=":")
        if col == "cancelacion_max":
            ax.axvline(1, color="grey", lw=0.8, ls=":")
    axes[0].set_yticks(y, nombres)
    fig.suptitle("Experimento de arquitectura (solo validación; puntos = semillas, raya = media). Exploratorio: no "
                 "cambia el modelo principal", fontsize=9)
    fig.tight_layout()
    fig.savefig(salida / "Fig_arquitectura.png", bbox_inches="tight")
    plt.close(fig)
    print(f"\nResultados en {_relativa(salida)}")


def _relativa(ruta: Path) -> Path:
    """Ruta relativa a la carpeta del proyecto, para que la consola no muestre rutas absolutas del equipo."""
    try:
        return ruta.relative_to(RUTA_PROYECTO)
    except ValueError:
        return ruta


def main() -> None:
    global RUTA_DATOS, RUTA_SALIDA, RUTA_PRED
    # En Windows la consola usa cp1252 y no puede imprimir símbolos como ↔ o ψ
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-epochs", type=int, default=MAX_EPOCHS)
    parser.add_argument("--semillas", type=str, default=",".join(str(s) for s in SEMILLAS))
    parser.add_argument("--rapido", action="store_true")
    parser.add_argument("--sensibilidad", action="store_true",
                        help="análisis de sensibilidad de los datos: carpetas dentro de Sensibilidad/")
    parser.add_argument("--arquitectura", action="store_true",
                        help="solo el experimento de arquitectura de la cabeza (validación; escribe en "
                             "Resultados_python/arquitectura/)")
    args = parser.parse_args()
    if args.sensibilidad:
        base = RUTA_PROYECTO / "Sensibilidad"
        RUTA_DATOS, RUTA_SALIDA, RUTA_PRED = base / "Datos_preparados", base / "Resultados_python", base / "Predicciones_python"
        print(f"ANÁLISIS DE SENSIBILIDAD DE LOS DATOS: carpetas dentro de {_relativa(base)}")
    if args.rapido:
        args.max_epochs = min(args.max_epochs, 250)
        args.semillas = "0"

    semillas = [int(s) for s in args.semillas.split(",") if s.strip()]
    if args.arquitectura:
        experimento_arquitectura(semillas, args.max_epochs)
        return
    salida = RUTA_SALIDA
    figuras = salida / "figuras"
    for ruta in (salida, figuras, RUTA_PRED, RUTA_PRED / "semillas"):
        ruta.mkdir(parents=True, exist_ok=True)

    versiones = {"python": platform.python_version()}
    for paquete in ["torch", "pykan", "numpy", "pandas", "scipy", "matplotlib"]:
        versiones[paquete] = metadata.version(paquete)
    pd.DataFrame([versiones]).to_csv(salida / "versiones.csv", index=False)
    imprimir_tabla("Versiones", pd.DataFrame([versiones]))
    inicio = time.time()
    D = cargar_datos_funcionales(verificar=True)
    datos = preparar_datos(D)
    df = datos.df
    tabla_verificacion = pd.DataFrame(
        [{"control": k, "valor": v} for k, v in D.verificacion.items()]
    )
    tabla_verificacion.to_csv(salida / "verificacion_construccion.csv", index=False)
    imprimir_tabla("Verificación de la construcción funcional", tabla_verificacion)

    print("\nPartición:")
    for bloque, grupo in df.groupby("bloque"):
        print(f"  {bloque}: {len(grupo)} pares, {grupo.station_id.nunique()} estaciones, años {grupo.year.min()}-{grupo.year.max()}")
    # Bloque predicho y entregado a la V2: validación y, con ABRIR_PRUEBA, prueba
    idx_eval = np.concatenate([datos.idx_val, datos.idx_test]) if ABRIR_PRUEBA else datos.idx_val
    assert np.allclose(D.t_fine[CENTROS_DIA], np.arange(365) + 0.5)

    # FPCA con las curvas de estación de entrenamiento de la tabla de modelado (el decodificador de la regresión de
    # rango reducido de la V2)
    C_fpca = datos.Cy[datos.idx_train]
    mu, Bm, valores, var_acum, B = ajustar_fpca(C_fpca, D.W_est_z, M_FPCA)
    Bs = Bm * np.sqrt(valores[:M_FPCA])[None, :]
    ortonormal = float(np.max(np.abs(Bm.T @ D.W_est_z @ Bm - np.eye(M_FPCA))))
    print(f"FPCA: {len(C_fpca)} curvas, M={M_FPCA}, varianza acumulada={var_acum[M_FPCA - 1] * 100:.2f} %, "
          f"error ortonormalidad={ortonormal:.2e}")
    imprimir_tabla(
        "FPCA de este archivo frente al de la V2 (fpca_estaciones.csv)",
        pd.DataFrame([verificar_fpca_con_r(mu, B, valores, M_FPCA)]),
    )
    tabla_fpca = pd.DataFrame(
        {
            "componente": np.arange(1, 11),
            "autovalor": valores[:10],
            "varianza_pct": valores[:10] / valores.sum() * 100,
            "varianza_acumulada_pct": var_acum[:10] * 100,
        }
    )
    imprimir_tabla("Varianza explicada por la FPCA", tabla_fpca.round(4))
    graficar_fpca(D, mu, B, valores, var_acum, M_FPCA, figuras)

    rmse_kan, rmse_kan_nuevo = verificar_kan_sintetica(0)
    print(f"KAN sintética: RMSE en los puntos de ajuste={rmse_kan:.5f}, en puntos nuevos={rmse_kan_nuevo:.5f} "
          f"(exigido <= 1e-2)")
    if max(rmse_kan, rmse_kan_nuevo) > 1e-2:
        raise RuntimeError("La KAN no superó la verificación sintética")

    # Búsqueda de hiperparámetros: la misma rejilla K x lambda x lr para las tres cabezas, cada configuración con las
    # semillas del modelo final y resumida con la media de la pérdida de validación. K (común a las tres cabezas): el de
    # menor promedio, entre cabezas, de la mejor media de cada una. Con ese K, cada cabeza toma, entre las
    # configuraciones a menos del 1 % de su mejor media, la de mayor lambda y, en empate, la de menor media
    parametros_por_K = {K: emparejar_parametros(K + 2, M_FPCA, 0) for K in K_NEURONAS_GRID}
    busqueda = []
    for tipo in ["lineal", "mlp", "kan"]:
        for K in K_NEURONAS_GRID:
            for lam in LAMBDA_GRID:
                for lr in LR_GRID:
                    for semilla in semillas:
                        res = entrenar(
                            tipo, semilla, K, lam, lr, D, datos, mu, Bs,
                            parametros_por_K[K][0] if tipo == "mlp" else None, args.max_epochs,
                        )
                        busqueda.append(
                            {
                                "cabeza": tipo,
                                "K_neuronas": K,
                                "lambda_rugosidad": lam,
                                "lr": lr,
                                "semilla": semilla,
                                "val_loss": res.mejor_val,
                                "mejor_epoca": res.mejor_epoca,
                                "epocas": res.epocas,
                                "segundos": res.segundos,
                            }
                        )
                    vals = [b["val_loss"] for b in busqueda[-len(semillas):]]
                    print(f"  búsqueda {tipo}, K={K}, lambda={lam:g}, lr={lr:g}: val={' / '.join(f'{v:.5f}' for v in vals)} "
                          f"(media {np.mean(vals):.5f})")
    tabla_busqueda = pd.DataFrame(busqueda)
    claves_config = ["cabeza", "K_neuronas", "lambda_rugosidad", "lr"]
    resumen_busqueda = tabla_busqueda.groupby(claves_config, as_index=False).agg(
        semillas=("val_loss", "size"), val_media=("val_loss", "mean"), val_DE=("val_loss", "std"),
        mejor_epoca_media=("mejor_epoca", "mean"), epocas_max=("epocas", "max"))
    mejor_por_K = resumen_busqueda.groupby(["K_neuronas", "cabeza"]).val_media.min().unstack("cabeza")
    mejor_por_K["promedio"] = mejor_por_K.mean(axis=1)
    K = int(mejor_por_K.promedio.idxmin())
    imprimir_tabla(
        "Mejor media de validación (entre semillas) por K y cabeza (se elige el K de menor promedio)",
        mejor_por_K.reset_index().round(5),
    )
    resumen_busqueda["elegido"] = False
    config: Dict[str, Tuple[float, float]] = {}
    for tipo in ["lineal", "mlp", "kan"]:
        sub = resumen_busqueda[(resumen_busqueda.cabeza == tipo) & (resumen_busqueda.K_neuronas == K)]
        candidatos = sub[sub.val_media <= sub.val_media.min() * 1.01]
        elegido = candidatos.sort_values(["lambda_rugosidad", "val_media"], ascending=[False, True]).iloc[0]
        resumen_busqueda.loc[elegido.name, "elegido"] = True
        config[tipo] = (float(elegido.lambda_rugosidad), float(elegido.lr))
    tabla_busqueda = tabla_busqueda.merge(resumen_busqueda[claves_config + ["elegido"]], on=claves_config)
    tabla_busqueda.to_csv(salida / "busqueda_hiperparametros.csv", index=False)
    resumen_busqueda.to_csv(salida / "busqueda_hiperparametros_resumen.csv", index=False)
    imprimir_tabla(
        f"Búsqueda de hiperparámetros (solo validación, {len(resumen_busqueda) // 3} configuraciones por cabeza, "
        f"{len(semillas)} semillas cada una: media y DE de la pérdida de validación)",
        resumen_busqueda.sort_values(["cabeza", "val_media"]).round(5),
    )
    print(f"K elegido = {K}; " + "; ".join(f"{t}: lambda={l:g}, lr={r:g}" for t, (l, r) in config.items()))

    h_mlp, n_kan, n_mlp, rel_param = parametros_por_K[K]
    print(f"Parámetros entrenables de la cabeza: KAN={n_kan}, MLP(h={h_mlp})={n_mlp}, diferencia={rel_param * 100:.2f} %")

    resultados: Dict[str, List[ResultadoEntrenamiento]] = {"lineal": [], "mlp": [], "kan": []}
    for tipo in ["lineal", "mlp", "kan"]:
        lam, lr = config[tipo]
        for semilla in semillas:
            res = entrenar(
                tipo,
                semilla,
                K,
                lam,
                lr,
                D,
                datos,
                mu,
                Bs,
                h_mlp if tipo == "mlp" else None,
                args.max_epochs,
            )
            resultados[tipo].append(res)
    # Control: cada corrida final repite la de la búsqueda (misma configuración y semilla; entrenamiento determinista)
    dif_busqueda = 0.0
    for tipo in resultados:
        lam, lr = config[tipo]
        for semilla, res in zip(semillas, resultados[tipo]):
            fila = tabla_busqueda[(tabla_busqueda.cabeza == tipo) & (tabla_busqueda.K_neuronas == K)
                                  & (tabla_busqueda.lambda_rugosidad == lam) & (tabla_busqueda.lr == lr)
                                  & (tabla_busqueda.semilla == semilla)]
            dif_busqueda = max(dif_busqueda, abs(res.mejor_val - float(fila.val_loss.iloc[0])))
    print(f"Control: corridas finales frente a las mismas corridas de la búsqueda, max |diferencia| = {dif_busqueda:.1e}")

    # Predicciones de cada corrida en Predicciones_python/semillas/ (tabla de semillas de la V2); en la carpeta
    # principal, la semilla de menor pérdida de validación de cada cabeza
    filas_semilla = []
    for tipo, lista in resultados.items():
        for k, res in enumerate(lista):
            _, z_eval, _ = predecir(res.modelo, datos, idx_eval)
            escribir_predicciones(
                df.iloc[idx_eval],
                z_eval,
                RUTA_PRED / "semillas" / f"{NOMBRE_R[tipo]}_s{semillas[k]}.csv",
            )
            filas_semilla.append(
                {
                    "modelo": tipo,
                    "semilla": semillas[k],
                    "val_loss": res.mejor_val,
                    "mejor_epoca": res.mejor_epoca,
                    "epocas": res.epocas,
                    "segundos": res.segundos,
                    "pct_z_negativos": float(np.mean(z_eval[:, CENTROS_DIA] < 0) * 100),
                }
            )
    tabla_semilla = pd.DataFrame(filas_semilla)
    imprimir_tabla("Pérdida de validación por semilla", tabla_semilla.round(5))

    mejor_por_modelo = {
        tipo: min(lista, key=lambda r: r.mejor_val) for tipo, lista in resultados.items()
    }
    for tipo, res in mejor_por_modelo.items():
        _, z_eval, _ = predecir(res.modelo, datos, idx_eval)
        escribir_predicciones(df.iloc[idx_eval], z_eval, RUTA_PRED / f"{NOMBRE_R[tipo]}.csv")
        print(f"  {NOMBRE_R[tipo]}.csv: semilla {semillas[resultados[tipo].index(res)]} (menor pérdida de validación), "
              f"{len(idx_eval)} pares")

    # Recorte y no negatividad en los centros de día del bloque predicho
    filas_clip = []
    for tipo, res in mejor_por_modelo.items():
        _, z_fine, _ = predecir(res.modelo, datos, idx_eval)
        z = z_fine[:, CENTROS_DIA]
        filas_clip.append(
            {
                "modelo": tipo,
                "pct_puntos_z_negativos": float(np.mean(z < 0) * 100),
                "pct_curvas_con_negativos": float(np.mean(np.any(z < 0, axis=1)) * 100),
                "masa_mm_recortada_pct": float(
                    np.sum(np.minimum(z, 0) ** 2) / max(np.sum(z**2), 1e-12) * 100
                ),
                "min_z": float(z.min()),
                "max_z": float(z.max()),
            }
        )
    tabla_clip = pd.DataFrame(filas_clip)
    imprimir_tabla("Recorte y no negatividad (bloque predicho)", tabla_clip.round(6))

    # Equivalencia cabeza lineal - regresión funcional función a función (demostración en el docstring)
    mejor_lineal = mejor_por_modelo["lineal"]
    equivalencia, beta_lineal = verificar_equivalencia(
        mejor_lineal.modelo, datos, D, mu, B, M_FPCA, Bs, valores
    )
    tabla_equivalencia = pd.DataFrame([equivalencia]).T.reset_index()
    tabla_equivalencia.columns = ["control", "valor"]
    tabla_equivalencia.to_csv(salida / "equivalencia_cabeza_lineal.csv", index=False)
    imprimir_tabla("Equivalencia cabeza lineal ↔ regresión funcional", tabla_equivalencia)
    if max(equivalencia["max_abs_red_vs_integral_simpson"], equivalencia["max_abs_curva_regresion_V2_como_red"]) > 1e-4:
        raise RuntimeError("La verificación de la equivalencia no cumple la tolerancia de 1e-4")

    # Parámetros entrenables (los que se igualan entre KAN y MLP) y total registrado
    filas_param = []
    for tipo, lista in resultados.items():
        modelo = lista[0].modelo
        n_cabeza = contar_parametros(modelo.cabeza)
        n_encoder = contar_parametros(modelo.encoder)
        filas_param.append(
            {
                "modelo": tipo,
                "parametros_cabeza_entrenables": n_cabeza,
                "parametros_cabeza_registrados": contar_parametros_total(modelo.cabeza),
                "parametros_encoder": n_encoder,
                "parametros_entrenables_totales": n_cabeza + n_encoder,
                "diferencia_kan_mlp_pct": rel_param * 100,
            }
        )
    tabla_param = pd.DataFrame(filas_param)
    tabla_param.to_csv(salida / "parametros.csv", index=False)
    imprimir_tabla("Conteo de parámetros", tabla_param)

    mejor_val_global = min(mejor_por_modelo, key=lambda t: mejor_por_modelo[t].mejor_val)

    # Pesos finales: funciones w_k(t) y todos los parámetros entrenables
    filas_w = []
    filas_pesos = []
    for tipo, res in mejor_por_modelo.items():
        curvas_w = D.phi_sat_z_fine @ parametros_w(res.modelo).T
        for k in range(curvas_w.shape[1]):
            for t, valor in zip(D.t_fine, curvas_w[:, k]):
                filas_w.append(
                    {"modelo": tipo, "neurona": k + 1, "t": float(t), "w": float(valor)}
                )
        for nombre, parametro in res.modelo.named_parameters():
            forma = "x".join(str(d) for d in parametro.shape) or "escalar"
            for indice, valor in enumerate(parametro.detach().cpu().numpy().ravel()):
                filas_pesos.append(
                    {
                        "modelo": tipo,
                        "parametro": nombre,
                        "forma": forma,
                        "indice": indice,
                        "valor": float(valor),
                    }
                )
    pd.DataFrame(filas_w).to_csv(salida / "pesos_w_finales.csv", index=False)
    pd.DataFrame(filas_pesos).to_csv(salida / "pesos_cabezas_finales.csv", index=False)

    # INTERPRETACIÓN (4.5): semilla elegida de cada cabeza y datos de entrenamiento, salvo los ejemplos y los peores
    # casos (bloque predicho)
    T = temporadas_entrenamiento(D)
    tt = D.t_fine
    print(
        f"\nTemporadas de entrenamiento (curva media de {T['n_curvas']} curvas de estación, regla del EDA): "
        f"T1 {fecha_de(tt[T['inicio'][0]])} → pico {fecha_de(tt[T['picos'][0]])} → {fecha_de(tt[T['fin'][0]])}; "
        f"T2 {fecha_de(tt[T['inicio'][1]])} → pico {fecha_de(tt[T['picos'][1]])} → {fecha_de(tt[T['fin'][1]])}; "
        f"valle {fecha_de(tt[T['valles'][1]])}"
    )
    filas_ventanas = []
    for tipo, res in mejor_por_modelo.items():
        graficar_w(res.modelo, tipo, D, T, figuras)
        filas_ventanas.extend(resumir_ventanas(res.modelo, tipo, D, T, datos))
    tabla_ventanas = pd.DataFrame(filas_ventanas)
    tabla_ventanas.to_csv(salida / "ventanas_wk.csv", index=False)
    imprimir_tabla(
        "Ventanas de las neuronas: DE entre pares de la contribución de cada periodo a z_k, relativa a la DE de z_k",
        tabla_ventanas.round(3),
    )

    modelo_kan = mejor_por_modelo["kan"].modelo
    tabla_act, _, rangos, control_act = analizar_activaciones_kan(modelo_kan, K, datos, figuras)
    tabla_act.to_csv(salida / "activaciones_kan.csv", index=False)
    print(f"\nControl: suma de las aristas frente a la salida de la capa de pykan = {control_act:.2e}")
    principal = tabla_act.loc[tabla_act.groupby("entrada", sort=False).importancia.idxmax()]
    resumen_entradas = (
        tabla_act.groupby("entrada", sort=False)
        .agg(importancia_suma=("importancia", "sum"), R2_recta_min=("R2_recta", "min"),
             valores_distintos=("valores_distintos", "first"), pct_fuera_rejilla=("pct_fuera_rejilla", "first"))
        .reset_index()
        .merge(principal[["entrada", "nodo", "importancia", "R2_recta"]].rename(
            columns={"nodo": "nodo_principal", "importancia": "importancia_principal", "R2_recta": "R2_recta_principal"}),
            on="entrada")
        .sort_values("importancia_suma", ascending=False)
    )
    imprimir_tabla("Activaciones de la KAN por entrada (importancia = DE de la activación en entrenamiento)",
                   resumen_entradas.round(4))
    resumen_simb, detalle_simb = simbolico_kan(modelo_kan, K, datos, tabla_act, rangos, figuras)
    detalle_simb.to_csv(salida / "ajuste_simbolico_candidatos.csv", index=False)
    resumen_simb.to_csv(salida / "ajuste_simbolico.csv", index=False)
    imprimir_tabla("Ajuste simbólico de activaciones de la primera capa", resumen_simb.round(4))

    # Aportes de la KAN por condición (entrenamiento)
    tabla_cond, resumen_cond, control_u, control_z, (q1, q2) = aportes_kan_por_condicion(modelo_kan, K, datos, D, T, figuras)
    tabla_cond.to_csv(salida / "aportes_kan_condicion.csv", index=False)
    print(f"\nControl de la descomposición de la KAN: puntajes {control_u:.2e}; curva z(t) {control_z:.2e}. "
          f"Terciles del acumulado relativo (seco <= {100 * q1:+.1f} %, húmedo >= {100 * q2:+.1f} %)")
    imprimir_tabla("KAN por condición: cambio de la curva en cada periodo (raíz de mm/día), suma de los nodos, cambio "
                   "observado y razón de cancelación = suma |aporte de cada nodo| / |suma|", resumen_cond.round(3))
    imprimir_tabla("KAN por condición: entradas -> nodo oculto (cambio del nodo)",
                   tabla_cond[tabla_cond.capa == "entradas -> nodo"].pivot_table(
                       index=["contraste", "nodo"], columns="grupo_o_periodo", values="diferencia", sort=False)
                   .round(3).reset_index())
    imprimir_tabla("KAN por condición: nodo oculto -> curva (raíz de mm/día)",
                   tabla_cond[tabla_cond.capa != "entradas -> nodo"].pivot_table(
                       index=["contraste", "nodo"], columns="grupo_o_periodo", values="diferencia", sort=False)
                   .round(3).reset_index())
    # Flujo entre capas: diagrama de pykan, efecto de cada nodo en la curva, ablación y cancelación
    tabla_capas, base_capas, control_capas, reparto = interpretar_capas_kan(modelo_kan, K, datos, D, T, figuras)
    tabla_capas.to_csv(salida / "ablacion_kan.csv", index=False)
    canc_capas = cancelacion_por_capa(modelo_kan, "kan", datos)
    print(f"\nAblación de la KAN: pérdida de validación base {base_capas:.5f} (control frente a la red {control_capas:.1e}); "
          f"razón de cancelación por capa (mediana de los nodos, entrenamiento): "
          + " / ".join(f"{c:.2f}" for c in canc_capas))
    # Colinealidad de las z y cancelación en la cabeza lineal (descriptivo)
    tabla_colin = colinealidad_z(mejor_por_modelo, datos)
    tabla_colin.to_csv(salida / "colinealidad_z.csv", index=False)
    imprimir_tabla("Colinealidad de las z_k (entrenamiento) y cancelación en la cabeza lineal", tabla_colin.round(3))
    imprimir_tabla("Ablación: entradas y nodos ocultos (subida de la pérdida de validación, %)",
                   tabla_capas[tabla_capas.pieza != "arista"][["pieza", "origen", "subida_perdida_pct"]].round(2))
    imprimir_tabla("Ablación: las 10 aristas más importantes de cada capa",
                   tabla_capas[tabla_capas.pieza == "arista"].sort_values("subida_perdida_pct", ascending=False)
                   .groupby("capa").head(10).round(3))
    imprimir_tabla("Reparto de la DE de la capa 1 en cada nodo oculto (filas H1..H8; entrenamiento)",
                   pd.DataFrame(reparto, columns=nombres_entradas(K), index=[f"H{j + 1}" for j in range(reparto.shape[0])])
                   .round(3).reset_index())

    tabla_pd, ref_oni, ref_años, ref_alt, nota_500, control_gamma, niveles = contraste_cabezas(
        resultados, mejor_por_modelo, datos, D, T, figuras
    )
    tabla_pd.to_csv(salida / "dependencia_parcial.csv", index=False)
    print(f"\nControl: coeficiente local de la cabeza lineal frente a gamma_ONI(t) de sus pesos = {control_gamma:.2e}")
    for variable in ("ONI", "altitud"):
        sub = tabla_pd[(tabla_pd.variable == variable) & tabla_pd.valor.isin(niveles[variable])]
        elegida = sub[sub.elegida].pivot_table(index=["cabeza", "valor"], columns="periodo", values="cambio_pct").reset_index()
        rango = (
            sub[sub.periodo == "Anual"].groupby(["cabeza", "valor"]).cambio_pct.agg(["min", "max"])
            .rename(columns={"min": "Anual_min_semillas", "max": "Anual_max_semillas"}).reset_index()
        )
        imprimir_tabla(
            f"Dependencia parcial de {variable}: cambio de la lluvia predicha (%) frente a "
            + ("ONI = 0" if variable == "ONI" else f"{niveles['altitud'][1]:.0f} m")
            + " (semilla elegida; rango entre las 3 semillas)",
            elegida.merge(rango, on=["cabeza", "valor"]).round(2),
        )
    imprimir_tabla("Referencia empírica (entrenamiento): cociente estación/satélite por año, frente a los años neutros",
                   ref_años.round(2))
    imprimir_tabla("Referencia empírica por fase ENSO", ref_oni.round(2))
    imprimir_tabla("Referencia empírica por decil de altitud, frente a 1200-2000 m", ref_alt.round(2))
    print(nota_500)

    # Ejemplos y peores casos por bloque: validación y, con ABRIR_PRUEBA, prueba
    bloques_figura = [("validacion", datos.idx_val)] + ([("prueba", datos.idx_test)] if ABRIR_PRUEBA else [])
    for nombre_bloque, idx_bloque in bloques_figura:
        preds = {
            tipo: np.maximum(predecir(res.modelo, datos, idx_bloque)[1], 0.0) ** 2
            for tipo, res in mejor_por_modelo.items()
        }
        tabla_peores = graficar_predicciones(df, datos, preds, idx_bloque, D, T, figuras, sufijo=f"_{nombre_bloque}")
        tabla_peores.to_csv(salida / f"peores_casos_{nombre_bloque}.csv", index=False)
        imprimir_tabla(f"Peores casos de {nombre_bloque} (mayor L2 de la KAN)", tabla_peores.round(3))

    tabla_estaciones = validacion_por_estaciones(config, K, h_mlp, D, datos, mu, Bs, args.max_epochs)
    tabla_estaciones.to_csv(salida / "validacion_por_estaciones.csv", index=False)
    resumen_est = tabla_estaciones.groupby(["pliegue", "cabeza"])[
        ["perdida_val_estaciones_vistas", "perdida_val_estaciones_nuevas"]].mean().unstack("cabeza")
    imprimir_tabla("Estaciones fuera del entrenamiento: pérdida de validación media por pliegue (2 semillas)",
                   resumen_est.round(4).reset_index())
    media_est = tabla_estaciones.groupby("cabeza")[["perdida_val_estaciones_vistas", "perdida_val_estaciones_nuevas"]].mean()
    mejora = (100 * (1 - media_est / media_est.loc["lineal"])).round(1)
    imprimir_tabla("Mejora frente a la cabeza lineal (%), promedio de los pliegues", mejora.reset_index())

    fig, ax = plt.subplots(figsize=(7, 3.6))
    for tipo, lista in resultados.items():
        hist = lista[0].historial
        ax.plot([h[0] for h in hist], [h[2] for h in hist], label=tipo, lw=1.2)
    ax.set_yscale("log")
    ax.set_xlabel("Época")
    ax.set_ylabel("Pérdida de validación")
    ax.set_title("Entrenamiento (semilla 0)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figuras / "Fig_entrenamiento.png", bbox_inches="tight")
    plt.close(fig)

    print("\n=== Configuración final ===")
    print(
        f"FPCA: M={M_FPCA}, {var_acum[M_FPCA - 1] * 100:.2f} % de varianza; "
        f"encoder K={K}; {len(semillas)} semillas; {args.max_epochs} épocas máximas; "
        + "; ".join(f"{t}: lambda={l:g}, lr={r:g}" for t, (l, r) in config.items())
    )
    print(
        f"Modelo propuesto: KAN (semilla elegida por validación). "
        f"Menor val_loss global: {mejor_val_global.upper()} "
        f"({mejor_por_modelo[mejor_val_global].mejor_val:.5f} vs "
        f"{mejor_por_modelo['kan'].mejor_val:.5f} de KAN)."
    )
    print("\nArchivos relevantes:")
    for ruta in [
        salida / "verificacion_construccion.csv",
        salida / "pesos_w_finales.csv",
        salida / "pesos_cabezas_finales.csv",
        *[RUTA_PRED / f"{NOMBRE_R[t]}.csv" for t in NOMBRE_R],
        RUTA_PRED / "semillas",
    ]:
        print(f"  {_relativa(ruta)}")
    print(f"Tiempo total: {time.time() - inicio:.1f}s")


if __name__ == "__main__":
    main()
