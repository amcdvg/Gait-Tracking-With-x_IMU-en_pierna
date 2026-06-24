"""
=============================================================================
PIPELINE DE ANÁLISIS DE MARCHA CON SENSORES IMU (TEM2)
=============================================================================
Autor      : Pipeline generado para análisis biomecánico de pierna
Descripción: Pipeline completo que realiza EDA, filtrado, detección de
             zancadas mediante modelo pendular de la pierna y cálculo de
             métricas espaciotemporales para sensores montados en la pierna.

Sensores:
  - Mov  (pierna derecha)
  - Mov2 (pierna izquierda)

Modelo biomecánico:
  A diferencia del ZUPT clásico (para sensores en pie/tobillo que tienen
  fase de velocidad cero), aquí usamos el MODELO DE PÉNDULO INVERTIDO.
  La pierna oscila como un péndulo alrededor de la cadera; los eventos
  del ciclo (Heel Strike, Toe Off) corresponden a los extremos y cruces
  por cero del ángulo sagital theta (θ).

  theta(t) se extrae de los cuaterniones q→Euler (eje sagital = pitch).
  Los mínimos locales de theta = Heel Strike (HS)
  Los máximos locales de theta = Toe Off (TO)
  (El signo depende de la orientación del sensor; se ajusta automáticamente)

Frecuencia de muestreo: 100 Hz
=============================================================================
"""

# ---------------------------------------------------------------------------
# 0. IMPORTACIONES
# ---------------------------------------------------------------------------
import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import signal
from scipy.spatial.transform import Rotation

warnings.filterwarnings("ignore")
matplotlib.rcParams.update({
    "figure.dpi": 120,
    "axes.grid": True,
    "grid.alpha": 0.4,
    "font.size": 9,
})

# ---------------------------------------------------------------------------
# 1. CONFIGURACIÓN GLOBAL
# ---------------------------------------------------------------------------
FS = 100.0          # Frecuencia de muestreo (Hz) — indicada en cabecera
LOWPASS_CUTOFF = 6.0    # Frecuencia de corte para filtro Butterworth acc/gyr (Hz)
THETA_CUTOFF   = 3.0    # Frecuencia de corte dedicada para theta (Hz) — más agresivo
LOWPASS_ORDER  = 4      # Orden del filtro

# Longitud estimada del segmento de pierna (muslo + tibia) en metros.
# Se usa como radio del péndulo para estimar longitud de zancada.
# Valor típico adulto: 0.85 m. Ajustar según antropometría del sujeto.
L_LEG = 0.45  # metros

# Gravedad
G = 9.81  # m/s²

# Umbral mínimo de prominencia para detección de picos en theta (rad)
# Aumentado a 0.20 para filtrar definitivamente los picos contralaterales
# sin necesidad de hacer una distancia temporal restrictiva.
MIN_PROMINENCE_THETA = 0.12  # ~6.8 grados
# Distancia a 60 para permitir carrera rápida (0.6s) ignorando el ruido contralateral
MIN_DISTANCE_SAMPLES = 60    # 0.40 s a 100 Hz → ciclo mínimo ~0.40 s
SG_WINDOW = 31               # Ventana Savitzky-Golay para theta (310 ms a 100 Hz)

FILES = {
    "marcha_5kmh":   "data/marcha_5kmh_50metro_jesus_20260611_181630_tem2.txt",
    "carrera_10kmh": "data/carrera_10kmh_100metros_jesus_20260611_181807_tem2.txt",
    "carrera_15kmh": "data/carrera_15kmh_100metros_jesus_20260611_181919_tem2.txt",
}

DATA_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# 2. CARGA Y PARSEO DE DATOS
# ---------------------------------------------------------------------------

def parse_european_float(s: str) -> float:
    """
    Convierte string con decimal europeo (coma) a float.
    Ejemplo: '9,81' → 9.81
    """
    if isinstance(s, float):
        return s
    return float(str(s).replace(",", "."))


def load_sensor_data(filepath: str, sensor_canal: str = "Mov") -> pd.DataFrame:
    """
    Carga y parsea el archivo TEM2.

    El archivo tiene:
      - Línea 1: 'Frecuencia (Hz);100'
      - Línea 2: descripción sensores
      - Línea 3: cabecera de columnas
      - Líneas 4+: datos separados por ';' con decimales ','

    El campo 'id' (canal interno 1-4) representa diferentes proyecciones
    del mismo vector de aceleración calculadas por el firmware.
    Usamos id=2 (aceleración en coordenadas del sensor, componente principal).
    Para los cuaterniones (idénticos en todas las filas del mismo sample),
    usamos cualquier id=1.

    Parameters
    ----------
    filepath : str
        Ruta al archivo .txt
    sensor_canal : str
        'Mov', 'Mov2' o 'Wrist'

    Returns
    -------
    pd.DataFrame con columnas:
        time_s, sample, accX, accY, accZ,
        rawGirX, rawGirY, rawGirZ, qX, qY, qZ, qW
    """
    print(f"  Cargando {os.path.basename(filepath)} → canal '{sensor_canal}'...")

    rows = []
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            if i < 3:          # Saltar cabecera (3 primeras líneas)
                continue
            parts = line.strip().split(";")
            if len(parts) < 19:
                continue

            canal = parts[2].strip()
            if canal != sensor_canal:
                continue

            id_val = parts[3].strip()
            # id=2 tiene aceleración en coordenadas del sensor (acc corregida)
            # Los cuaterniones son iguales para id 1..4 del mismo sample
            # Tomamos id=2 para acc y giroscopio, id=1 para quaterniones
            if id_val not in ("1", "2"):
                continue

            try:
                time_str = parts[0].strip()   # HH:MM:SS.mmm
                sample   = int(parts[4])
                gap      = parse_european_float(parts[5])

                accX = parse_european_float(parts[6])
                accY = parse_european_float(parts[7])
                accZ = parse_european_float(parts[8])

                # rawAcc — aceleración cruda en sistema del sensor
                rawAccX = parse_european_float(parts[9])
                rawAccY = parse_european_float(parts[10])
                rawAccZ = parse_european_float(parts[11])

                # Giroscopio (°/s)
                rawGirX = parse_european_float(parts[12])
                rawGirY = parse_european_float(parts[13])
                rawGirZ = parse_european_float(parts[14])

                # Cuaternión de orientación (salida del filtro del sensor)
                qX = parse_european_float(parts[15])
                qY = parse_european_float(parts[16])
                qZ = parse_european_float(parts[17])
                qW = parse_european_float(parts[18])

                rows.append({
                    "time_str": time_str,
                    "sample": sample,
                    "id": int(id_val),
                    "gap": gap,
                    "accX": accX, "accY": accY, "accZ": accZ,
                    "rawAccX": rawAccX, "rawAccY": rawAccY, "rawAccZ": rawAccZ,
                    "rawGirX": rawGirX, "rawGirY": rawGirY, "rawGirZ": rawGirZ,
                    "qX": qX, "qY": qY, "qZ": qZ, "qW": qW,
                })
            except (ValueError, IndexError):
                continue

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"Sin datos para sensor '{sensor_canal}' en {filepath}")

    # Reconstruir tiempo en segundos a partir del timestamp HH:MM:SS.mmm
    def time_to_seconds(ts: str) -> float:
        try:
            h, m, rest = ts.split(":")
            s_f = float(rest.replace(",", "."))
            return int(h) * 3600 + int(m) * 60 + s_f
        except Exception:
            return np.nan

    df["time_s_raw"] = df["time_str"].apply(time_to_seconds)

    # Separar id=1 (cuaterniones de referencia) y id=2 (aceleración principal)
    df_id1 = df[df["id"] == 1].copy()
    df_id2 = df[df["id"] == 2].copy()

    # Mergeamos por sample: tomamos acc de id=2 y quaternion de id=1
    # (son idénticos pero por consistencia tomamos id=1 para quaternion)
    df_q  = df_id1[["sample", "time_s_raw", "qX", "qY", "qZ", "qW",
                     "rawGirX", "rawGirY", "rawGirZ"]].drop_duplicates("sample")
    df_a  = df_id2[["sample", "accX", "accY", "accZ",
                     "rawAccX", "rawAccY", "rawAccZ"]].drop_duplicates("sample")

    df_merged = pd.merge(df_q, df_a, on="sample", how="inner")
    df_merged.sort_values("sample", inplace=True)
    df_merged.reset_index(drop=True, inplace=True)

    # Tiempo relativo en segundos (inicio = 0)
    t0 = df_merged["time_s_raw"].iloc[0]
    df_merged["time_s"] = df_merged["time_s_raw"] - t0

    # Si hay salto de medianoche, corregir
    diff = df_merged["time_s"].diff()
    jumps = diff[diff < -3600].index
    for j in jumps:
        df_merged.loc[j:, "time_s"] += 86400

    df_merged.drop(columns=["time_s_raw"], inplace=True)

    # Fix para eliminar "pixelación" o efecto escalera (staircase)
    # causado por el redondeo o congelamiento de los timestamps del sensor.
    # Forzamos una cuadrícula de tiempo perfectamente uniforme (fs ~ 100Hz).
    t_end = df_merged["time_s"].iloc[-1]
    df_merged["time_s"] = np.linspace(0.0, t_end, len(df_merged))

    # Normalizar cuaternión (por robustez numérica)
    q_norm = np.sqrt(df_merged["qX"]**2 + df_merged["qY"]**2 +
                     df_merged["qZ"]**2 + df_merged["qW"]**2)
    q_norm = q_norm.replace(0, np.nan).fillna(1.0)
    for col in ["qX", "qY", "qZ", "qW"]:
        df_merged[col] = df_merged[col] / q_norm

    print(f"    → {len(df_merged)} muestras cargadas "
          f"(duración: {df_merged['time_s'].iloc[-1]:.1f} s)")
    return df_merged


# ---------------------------------------------------------------------------
# 3. EDA — ANÁLISIS EXPLORATORIO DE DATOS
# ---------------------------------------------------------------------------

def eda_report(df: pd.DataFrame, label: str) -> None:
    """
    Genera un informe EDA básico: estadísticas, valores nulos y outliers.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame del sensor (salida de load_sensor_data)
    label : str
        Etiqueta descriptiva (ej. 'Mov - marcha_5kmh')
    """
    print(f"\n{'='*60}")
    print(f"EDA: {label}")
    print(f"{'='*60}")
    print(f"  Muestras    : {len(df)}")
    print(f"  Duración    : {df['time_s'].iloc[-1]:.2f} s")
    print(f"  Nulos       : {df.isnull().sum().sum()}")

    # Estadísticas de las señales principales
    cols = ["accX", "accY", "accZ", "rawGirX", "rawGirY", "rawGirZ",
            "qX", "qY", "qZ", "qW"]
    print("\n  Estadísticas (señales principales):")
    print(df[cols].describe().round(4).to_string(max_cols=10))

    # Detección de outliers (Z-score > 4)
    from scipy.stats import zscore
    for col in ["accX", "accY", "accZ"]:
        z = np.abs(zscore(df[col].dropna()))
        n_out = (z > 4).sum()
        if n_out > 0:
            print(f"  ⚠  Outliers en {col}: {n_out} muestras (|z|>4)")

    # Verificación de frecuencia de muestreo real
    dt = df["time_s"].diff().dropna()
    fs_real = 1.0 / dt.median()
    print(f"\n  Fs nominal: {FS:.0f} Hz | Fs estimada (mediana Δt): {fs_real:.1f} Hz")
    print(f"  Gap medio   : {dt.mean()*1000:.2f} ms | Gap máx: {dt.max()*1000:.2f} ms")


# ---------------------------------------------------------------------------
# 4. FILTRADO DE SEÑALES
# ---------------------------------------------------------------------------

def butter_lowpass_filter(data: np.ndarray, cutoff: float = LOWPASS_CUTOFF,
                           fs: float = FS, order: int = LOWPASS_ORDER) -> np.ndarray:
    """
    Filtro Butterworth de paso bajo — fase cero (filtfilt).

    Decisión matemática: elegimos Butterworth porque tiene respuesta plana
    en la banda de paso (maximally flat magnitude) y es ideal para señales
    biomecánicas donde queremos preservar la amplitud de los componentes
    de movimiento sin distorsión. El cutoff de 6 Hz captura toda la energía
    relevante de la marcha (fundamental ~1-2 Hz, armónicos hasta ~4-5 Hz).

    filtfilt (zero-phase) evita el desfase de fase que introduciría lfilter,
    lo cual es crítico para la correcta detección temporal de eventos.

    Parameters
    ----------
    data   : array de señal cruda
    cutoff : frecuencia de corte (Hz)
    fs     : frecuencia de muestreo (Hz)
    order  : orden del filtro

    Returns
    -------
    array filtrado
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype="low", analog=False)
    return signal.filtfilt(b, a, data)


def filter_signals(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aplica filtro Butterworth LP a acc, giroscopio y cuaterniones.
    Agrega columnas '_filt' al DataFrame.

    Parameters
    ----------
    df : DataFrame del sensor (salida de load_sensor_data)

    Returns
    -------
    DataFrame con columnas adicionales filtradas
    """
    df = df.copy()

    signal_cols = ["accX", "accY", "accZ",
                   "rawGirX", "rawGirY", "rawGirZ",
                   "qX", "qY", "qZ", "qW"]

    for col in signal_cols:
        df[f"{col}_filt"] = butter_lowpass_filter(df[col].values)

    # Re-normalizar cuaternión filtrado
    q_norm = np.sqrt(df["qX_filt"]**2 + df["qY_filt"]**2 +
                     df["qZ_filt"]**2 + df["qW_filt"]**2)
    q_norm = q_norm.replace(0, 1.0)
    for col in ["qX_filt", "qY_filt", "qZ_filt", "qW_filt"]:
        df[col] = df[col] / q_norm

    return df


def plot_raw_vs_filtered(df: pd.DataFrame, label: str,
                          save_dir: str = None) -> None:
    """
    Gráfica de señales crudas vs filtradas para aceleración y giroscopio.
    """
    fig, axes = plt.subplots(3, 2, figsize=(14, 9), sharex=True)
    fig.suptitle(f"Señales Crudas vs Filtradas — {label}", fontsize=11, fontweight="bold")

    t = df["time_s"].values
    pairs = [
        ("accX",    "accX_filt",    "Acc X (m/s²)"),
        ("accY",    "accY_filt",    "Acc Y (m/s²)"),
        ("accZ",    "accZ_filt",    "Acc Z (m/s²)"),
        ("rawGirX", "rawGirX_filt", "Gir X (°/s)"),
        ("rawGirY", "rawGirY_filt", "Gir Y (°/s)"),
        ("rawGirZ", "rawGirZ_filt", "Gir Z (°/s)"),
    ]

    for ax, (raw, filt, ylabel) in zip(axes.flatten(), pairs):
        ax.plot(t, df[raw],  color="#aaaaaa", linewidth=0.6, label="Crudo", alpha=0.8)
        ax.plot(t, df[filt], color="#e63946", linewidth=1.0, label="Filtrado (LP 6 Hz)")
        ax.set_ylabel(ylabel)
        ax.legend(loc="upper right", fontsize=7)

    axes[-1, 0].set_xlabel("Tiempo (s)")
    axes[-1, 1].set_xlabel("Tiempo (s)")
    plt.tight_layout()

    if save_dir:
        fname = os.path.join(save_dir, f"raw_vs_filt_{label.replace(' ', '_')}.png")
        plt.savefig(fname, dpi=120, bbox_inches="tight")
        print(f"  → Guardada: {fname}")
    plt.show()


# ---------------------------------------------------------------------------
# 5. CONVERSIÓN CUATERNIÓN → ÁNGULO DE EULER (THETA SAGITAL)
# ---------------------------------------------------------------------------

def quaternion_to_euler_xyz(df: pd.DataFrame,
                             use_filtered: bool = True) -> pd.DataFrame:
    """
    Convierte cuaterniones a ángulos de Euler (roll, pitch, yaw) en radianes,
    con doble etapa de filtrado para obtener una señal theta muy limpia.

    Decisión matemática:
    -------------------
    Usamos la convención 'ZYX' (yaw-pitch-roll) de scipy.spatial.transform.Rotation.
    Para un sensor montado en la pierna en el plano sagital:

        - PITCH (rotación alrededor del eje medio-lateral Y) = θ(t)
          Es el ángulo que describe la oscilación pendular de la pierna
          hacia adelante/atrás. Este es el ángulo biomecánicamente relevante.

    Cadena de limpieza de theta (3 etapas):
      1. Cuaterniones ya filtrados con Butterworth 6 Hz (filter_signals)
      2. Conversión q → Euler → pitch crudo
      3. Butterworth LP adicional a 3 Hz (THETA_CUTOFF) — captura solo la
         oscilación pendular fundamental (~1-2 Hz) eliminando residuos
      4. Savitzky-Golay (ventana 31 muestras = 310 ms, orden 3) — preserva
         forma de los picos sin desplazar su posición temporal

    Parameters
    ----------
    df           : DataFrame con columnas qX, qY, qZ, qW (y versiones _filt)
    use_filtered : si True, usa los cuaterniones filtrados

    Returns
    -------
    DataFrame con columnas:
        roll, pitch, yaw (rad) — señal directa del quaternion
        theta     (rad) — pitch doblemente filtrado y suavizado (para detección)
        theta_deg (°)   — theta en grados para visualización
    """
    df = df.copy()

    suf = "_filt" if use_filtered else ""
    qx = df[f"qX{suf}"].values
    qy = df[f"qY{suf}"].values
    qz = df[f"qZ{suf}"].values
    qw = df[f"qW{suf}"].values

    # scipy usa convencion [x, y, z, w]
    quats = np.column_stack([qx, qy, qz, qw])

    # Rotación ZYX → Euler en orden (yaw, pitch, roll) = (z, y, x)
    rot = Rotation.from_quat(quats)
    euler = rot.as_euler("ZYX", degrees=False)  # shape (N, 3): yaw, pitch, roll

    df["yaw"]   = euler[:, 0]  # Rotación alrededor de Z (supero-inferior)
    df["pitch"] = euler[:, 1]  # Rotación alrededor de Y (medio-lateral) = θ sagital
    df["roll"]  = euler[:, 2]  # Rotación alrededor de X (antero-posterior)

    # ── Etapa 2: Butterworth LP 3 Hz sobre pitch crudo ─────────────────────
    # Reduce ruido residual que queda tras el filtrado del cuaternión.
    # 3 Hz conserva la fundamental de la marcha (0.5-2 Hz) y elimina
    # artefactos de alta frecuencia del cálculo de Euler.
    pitch_lp = butter_lowpass_filter(df["pitch"].values,
                                     cutoff=THETA_CUTOFF, fs=FS, order=LOWPASS_ORDER)

    # ── Etapa 3: Savitzky-Golay — suavizado sin desplazamiento de picos ───
    # Ventana de 31 muestras (310 ms). SG es óptimo para preservar la
    # forma y posición de los extremos (mínimos = HS, máximos = TO).
    pitch_sg = signal.savgol_filter(pitch_lp,
                                    window_length=SG_WINDOW, polyorder=3)

    df["theta"]     = pitch_sg                      # Alias semántico (doblemente filtrado)
    df["theta_raw"] = df["pitch"]                   # Pitch directo del quaternion (referencia)
    df["theta_deg"] = np.degrees(pitch_sg)          # En grados para visualización

    return df


def plot_theta(df: pd.DataFrame, label: str,
               hs_idx: np.ndarray = None, to_idx: np.ndarray = None,
               save_dir: str = None) -> None:
    """
    Gráfica del ángulo theta con 3 capas:
      - Señal cruda del pitch (gris claro, semitransparente)
      - Theta filtrado LP 3Hz + SG (azul oscuro, línea principal)
      - Marcadores HS (triángulo abajo, rojo) y TO (triángulo arriba, verde)
    """
    fig, ax = plt.subplots(figsize=(14, 4))
    t = df["time_s"].values

    # Capa 1: pitch crudo (referencia visual de la señal original)
    if "theta_raw" in df.columns:
        ax.plot(t, np.degrees(df["theta_raw"].values),
                color="#cccccc", linewidth=0.7, alpha=0.6,
                label="θ pitch crudo (desde quaternion)")

    # Capa 2: theta doblemente filtrado — señal de trabajo limpia
    ax.plot(t, np.degrees(df["theta"].values),
            color="#1d3557", linewidth=1.6,
            label=f"θ filtrado (LP {THETA_CUTOFF} Hz + SG {SG_WINDOW} muestras)")

    # Capa 3: eventos de marcha
    if hs_idx is not None and len(hs_idx) > 0:
        ax.scatter(t[hs_idx], np.degrees(df["theta"].values[hs_idx]),
                   marker="v", color="#e63946", s=80, zorder=5,
                   edgecolors="white", linewidths=0.5,
                   label=f"Heel Strike ({len(hs_idx)})")
    if to_idx is not None and len(to_idx) > 0:
        ax.scatter(t[to_idx], np.degrees(df["theta"].values[to_idx]),
                   marker="^", color="#2a9d8f", s=80, zorder=5,
                   edgecolors="white", linewidths=0.5,
                   label=f"Toe Off ({len(to_idx)})")

    ax.set_xlabel("Tiempo (s)", fontsize=10)
    ax.set_ylabel("θ — Ángulo sagital (°)", fontsize=10)
    ax.set_title(f"Ángulo Pendular de la Pierna — {label}",
                 fontweight="bold", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    ax.fill_between(t, np.degrees(df["theta"].values),
                    alpha=0.08, color="#1d3557")  # área bajo la curva
    plt.tight_layout()

    if save_dir:
        fname = os.path.join(save_dir, f"theta_{label.replace(' ', '_')}.png")
        plt.savefig(fname, dpi=150, bbox_inches="tight")
    plt.show()


# ---------------------------------------------------------------------------
# 6. DETECCIÓN DE EVENTOS DE MARCHA (MODELO PENDULAR)
# ---------------------------------------------------------------------------

def detect_gait_events(df: pd.DataFrame,
                        min_prominence: float = MIN_PROMINENCE_THETA,
                        min_distance: int = MIN_DISTANCE_SAMPLES,
                        is_left_leg: bool = False) -> dict:
    """
    Detecta Heel Strike (HS) y Toe Off (TO) usando el ángulo pendular θ.

    Modelo biomecánico del péndulo para pierna:
    -------------------------------------------
    A diferencia del pie/tobillo (donde ZUPT detecta velocidad≈0 durante
    el apoyo plano), la pierna describe una oscilación pendular continua.

    Definiciones basadas en el ciclo de marcha:
      • Heel Strike (HS): Momento en que el talón toca el suelo.
        Corresponde al MÍNIMO del ángulo θ (pierna extendida hacia adelante,
        desplazada del equilibrio anterior).
        Biomecánicamente: la pierna está en su posición más anterior,
        theta es mínimo (o máximo negativo).

      • Toe Off (TO): Momento en que el pie despega del suelo.
        Corresponde al MÁXIMO del ángulo θ (pierna extendida hacia atrás,
        en posición de propulsión).
        Biomecánicamente: la pierna está en su posición más posterior.

    NOTA: El signo de theta depende de la orientación de montaje del sensor.
    El algoritmo detecta ambos tipos de picos y los asigna correctamente
    analizando el contexto temporal.

    Algoritmo:
      1. Filtrar theta con Butterworth LP (ya filtrado en df).
      2. Suavizado adicional Savitzky-Golay para robustez de picos.
      3. Detectar picos positivos → candidatos TO
         Detectar valles (picos en -theta) → candidatos HS
      4. Verificar alternancia HS→TO→HS→... para ciclo completo.

    Parameters
    ----------
    df             : DataFrame con columna 'theta' (en radianes)
    min_prominence : prominencia mínima del pico (rad)
    min_distance   : distancia mínima entre picos (muestras)

    Returns
    -------
    dict con:
        'hs_indices' : array de índices de Heel Strike
        'to_indices' : array de índices de Toe Off
        'theta'      : array de theta suavizado
    """
    theta = df["theta"].values.copy()
    # theta ya viene doblemente filtrado (LP 3 Hz + SG 31) desde
    # quaternion_to_euler_xyz. No necesitamos suavizado adicional agresivo;
    # aplicamos sólo un SG mínimo (ventana 7) para robustez numérica de picos.
    theta_smooth = signal.savgol_filter(theta, window_length=7, polyorder=3)

    # Detectar picos positivos (candidatos a Toe Off)
    peaks_pos, props_pos = signal.find_peaks(
        theta_smooth,
        prominence=min_prominence,
        distance=min_distance
    )

    # Detectar valles (picos en la señal negada → candidatos a Heel Strike)
    peaks_neg, props_neg = signal.find_peaks(
        -theta_smooth,
        prominence=min_prominence,
        distance=min_distance
    )

    # Determinar asignación correcta HS/TO
    # Analizamos estadísticamente qué tipo de pico tiene valores más altos/bajos
    # Para pierna derecha montada lateralmente: HS suele ser mínimo (valle)
    # TO suele ser máximo (pico)

    # POLARIDAD SEGÚN GRÁFICAS DE REFERENCIA (THETA):
    # Para ambos sensores (Mov y Mov2), Heel Strike ocurre en los PICOS (máximo ángulo)
    # y Toe Off ocurre en los VALLES (mínimo ángulo). Esto asegura la alternancia perfecta
    # de los pasos en el tiempo.
    hs_indices = peaks_pos  # Picos  -> Heel Strike
    to_indices = peaks_neg  # Valles -> Toe Off

    return {
        "hs_indices": hs_indices,
        "to_indices": to_indices,
        "theta_smooth": theta_smooth,
    }


def validate_gait_events(events: dict, df: pd.DataFrame) -> dict:
    """
    Valida y refina los eventos asegurando alternancia correcta HS→TO→HS.
    También elimina falsos positivos basados en amplitud del ciclo.

    Parameters
    ----------
    events : salida de detect_gait_events
    df     : DataFrame del sensor

    Returns
    -------
    dict con eventos validados y estadísticas de ciclo
    """
    hs = np.sort(events["hs_indices"])
    to = np.sort(events["to_indices"])
    theta = events["theta_smooth"]
    t     = df["time_s"].values

    # Construir ciclos válidos: cada zancada = HS[i] → HS[i+1]
    # con al menos un TO entre ellos
    valid_strides = []
    for i in range(len(hs) - 1):
        hs_start = hs[i]
        hs_end   = hs[i + 1]

        # TO entre los dos HS consecutivos
        to_between = to[(to > hs_start) & (to < hs_end)]
        if len(to_between) == 0:
            continue  # No hay TO → ciclo incompleto, descartar

        to_mid = to_between[0]  # Tomamos el primero como TO nominal

        dt_stride = t[hs_end] - t[hs_start]   # Duración de zancada (s)
        dt_stance = t[to_mid] - t[hs_start]   # Fase de apoyo (HS→TO)
        dt_swing  = t[hs_end] - t[to_mid]     # Fase de vuelo (TO→HS)

        # Filtros de calidad (marcha: ~0.8–2.0 s por zancada)
        if dt_stride < 0.3 or dt_stride > 3.0:
            continue
        # Proporción stance típica: 40-80% de la zancada
        if dt_stance <= 0 or dt_swing <= 0:
            continue

        valid_strides.append({
            "hs_start_idx": hs_start,
            "to_idx":       to_mid,
            "hs_end_idx":   hs_end,
            "t_hs_start":   t[hs_start],
            "t_to":         t[to_mid],
            "t_hs_end":     t[hs_end],
            "dt_stride":    dt_stride,
            "dt_stance":    dt_stance,
            "dt_swing":     dt_swing,
            "theta_hs_start": theta[hs_start],
            "theta_to":       theta[to_mid],
            "theta_hs_end":   theta[hs_end],
        })

    # Extraer arrays de eventos validados
    if valid_strides:
        hs_valid = np.array([s["hs_start_idx"] for s in valid_strides] +
                             [valid_strides[-1]["hs_end_idx"]])
        hs_valid = np.unique(hs_valid)
        to_valid = np.array([s["to_idx"] for s in valid_strides])
    else:
        hs_valid = np.array([], dtype=int)
        to_valid = np.array([], dtype=int)

    return {
        "valid_strides": valid_strides,
        "hs_indices": hs_valid,
        "to_indices": to_valid,
        "theta_smooth": theta,
    }


# ---------------------------------------------------------------------------
# 7. CÁLCULO DE MÉTRICAS ESPACIOTEMPORALES POR ZANCADA
# ---------------------------------------------------------------------------

def compute_stride_length_pendulum(dt_stride: float,
                                   theta_range: float,
                                   l_leg: float = L_LEG) -> float:
    """
    Estimación de longitud de zancada usando el modelo pendular.

    Modelo biomecánico:
    ------------------
    Para una pierna modelada como péndulo simple (radio = L_LEG):

        L_zancada ≈ 2 * L_LEG * sin(θ_max)

    donde θ_max es el semi-ángulo de oscilación del péndulo
    (amplitud máxima desde la vertical).

    Esta es la proyección horizontal del extremo del péndulo cuando
    está en su punto más alejado de la vertical.

    Corrección dinámica: En marcha/carrera real, la pierna no es un
    péndulo pasivo. Corregimos con el tiempo de zancada usando la
    relación de período del péndulo equivalente:

        T_pendulo = 2π √(L/g)  (período natural)

    Si T_real < T_pendulo → mayor velocidad → mayor longitud real
    El factor de corrección es: k = T_natural / T_real (saturado en [0.5, 2.0])

    Parameters
    ----------
    dt_stride  : duración de la zancada (s)
    theta_range: rango total de theta en la zancada (rad) = θ_max - θ_min
    l_leg      : longitud del segmento de pierna (m)

    Returns
    -------
    Longitud estimada de zancada (m)
    """
    theta_max = theta_range / 2.0  # Semi-ángulo de oscilación

    # Longitud geométrica de UN PASO asumiendo L_leg = longitud total de pierna (0.85m)
    # y que el rango theta captura la apertura completa
    L_leg_full = 0.85
    Step_base = 2.0 * L_leg_full * np.sin(np.abs(theta_max))
    
    # La zancada (Stride) se compone de 2 pasos
    Stride_base = 2.0 * Step_base

    # Período natural del péndulo completo
    T_natural = 2 * np.pi * np.sqrt(L_leg_full / G)

    # Factor de corrección dinámica
    # Permite escalar velocidades más altas de carrera (k > 2.0)
    k = np.clip(T_natural / max(dt_stride, 0.1), 0.5, 3.0)

    # FACTOR DE COMPENSACIÓN ANATÓMICA (Sensor en Muslo)
    # Ya que el sensor está en el muslo, el ángulo theta no incluye la extensión 
    # adicional de la rodilla. La longitud matemática pura de la zancada subestima 
    # la distancia real recorrida por el pie. Se aplica un factor de escala empírico.
    THIGH_SCALE_FACTOR = 2.0  # Ajustado para que 50m den ~50m y no se desborde

    return Stride_base * k * THIGH_SCALE_FACTOR


def compute_step_height(theta_swing_max: float, l_leg: float = L_LEG) -> float:
    """
    Estimación de la elevación vertical del pie durante la fase de balanceo (swing).

    Derivación geométrica (péndulo simple, cadera como pivote):
    ------------------------------------------------------------
    Cuando la pierna se modela como un péndulo rígido de longitud L_leg
    pivotando en la cadera, la posición del pie en coordenadas sagitales es:

        x_pie = L_leg * sin(θ)        (desplazamiento horizontal)
        y_pie = -L_leg * cos(θ)       (posición vertical desde la cadera)

    Con la pierna vertical (θ = 0), el pie está en el suelo:
        y_pie(0) = -L_leg → altura sobre el suelo = 0

    Con la pierna en θ máximo (swing máximo):
        y_pie(θ) = -L_leg * cos(θ)
        Δh = y_pie(0) - y_pie(θ) = L_leg * cos(0) - L_leg * cos(θ)

    Por lo tanto:
        h_pie = L_leg * (1 - cos(θ_swing_max))         [metros]

    Validación con rangos biomecánicos conocidos:
      - Marcha (5 km/h): θ_swing_max ≈ 25–35° → h = 0.85*(1-cos30°) ≈ 11 cm  ✓ [10-20 cm]
      - Carrera (10 km/h): θ_swing_max ≈ 45–55° → h = 0.85*(1-cos50°) ≈ 31 cm  ✓ [20-40 cm]
      - Carrera (15 km/h): θ_swing_max ≈ 55–70° → h = 0.85*(1-cos65°) ≈ 49 cm  ✓ [≥40 cm]

    CLAVE: θ_swing_max es el ángulo máximo ABSOLUTO que alcanza la pierna
    durante la fase de vuelo (TO → HS siguiente), extraído directamente
    de la señal theta en esa ventana temporal. NO es la mitad del rango total
    de oscilación (error común).

    Parameters
    ----------
    theta_swing_max : máximo ángulo absoluto de la pierna respecto a la
                      vertical durante el swing (rad)
    l_leg           : longitud del segmento pierna hip→pie (m)

    Returns
    -------
    Elevación estimada del pie durante el swing (m)
    """
    return l_leg * (1.0 - np.cos(np.abs(theta_swing_max)))


def compute_stride_velocity_from_acc(df: pd.DataFrame,
                                      start_idx: int,
                                      end_idx: int) -> float:
    """
    Estima la velocidad de la zancada integrando la aceleración horizontal.

    Estrategia para pierna (sin ZUPT clásico):
    -----------------------------------------
    1. Extrae el segmento de aceleración de la zancada.
    2. Rota al sistema global usando los cuaterniones (para obtener
       componente horizontal de aceleración).
    3. Resta la componente gravitacional.
    4. Integra numéricamente con trapecio.
    5. Aplica CORRECCIÓN DE DERIVA (drift): como la pierna DEBE tener
       velocidad media cero relativa al torso (modelo pendular), restamos
       la tendencia lineal de la velocidad integrada dentro de cada zancada.
       Este reemplazo del ZUPT clásico asume que al inicio y final de la
       zancada la velocidad del segmento es similar (condición de periodicidad).

    Parameters
    ----------
    df        : DataFrame del sensor con columnas acc y quaternion filtrados
    start_idx : índice de inicio de la zancada (HS)
    end_idx   : índice de fin de la zancada (siguiente HS)

    Returns
    -------
    Velocidad media estimada (m/s)
    """
    if end_idx <= start_idx + 2:
        return np.nan

    seg = df.iloc[start_idx:end_idx].copy()
    t   = seg["time_s"].values - seg["time_s"].values[0]
    dt  = np.diff(t)

    # Aceleración en sistema del sensor (filtrada)
    acc = seg[["accX_filt", "accY_filt", "accZ_filt"]].values

    # Rotar aceleración al sistema global usando cuaterniones
    qx = seg["qX_filt"].values
    qy = seg["qY_filt"].values
    qz = seg["qZ_filt"].values
    qw = seg["qW_filt"].values
    quats = np.column_stack([qx, qy, qz, qw])
    rot   = Rotation.from_quat(quats)

    acc_global = rot.apply(acc)    # (N, 3) — en sistema global

    # Restar gravedad (componente Z global ≈ -g en sistema NED)
    acc_global[:, 2] -= (-G)       # Corregir: acc_global_z = acc_z + g

    # Componente de avance (X global ≈ dirección de la marcha)
    ax = acc_global[:, 0]

    # Integración trapezoidal → velocidad
    vel = np.zeros(len(ax))
    for i in range(1, len(ax)):
        vel[i] = vel[i-1] + ax[i-1] * dt[i-1]

    # CORRECCIÓN DE DERIVA (reemplazo del ZUPT para pierna):
    # Asumimos que la velocidad al inicio y final del ciclo son similares
    # → restamos tendencia lineal (drift lineal)
    drift_correction = np.linspace(vel[0], vel[-1], len(vel))
    vel_corrected = vel - drift_correction

    # Velocidad media de la zancada (módulo)
    return float(np.abs(np.mean(vel_corrected)))


def compute_stride_metrics(df: pd.DataFrame,
                            valid_strides: list,
                            label: str = "") -> pd.DataFrame:
    """
    Calcula métricas espaciotemporales para cada zancada detectada.

    Métricas por zancada:
      - t_start, t_end : tiempo absoluto de inicio/fin (s)
      - dt_stride      : duración de la zancada (s)
      - dt_stance      : fase de apoyo HS→TO (s)
      - dt_swing       : fase de vuelo TO→HS_siguiente (s)
      - stance_pct     : % de fase de apoyo
      - cadence_step   : cadencia de este paso (pasos/min)
      - theta_range    : rango total del ángulo θ (rad)
      - stride_length  : longitud estimada de zancada (m)
      - step_height    : altura estimada del paso (m)
      - stride_velocity: velocidad media estimada (m/s)

    Parameters
    ----------
    df            : DataFrame del sensor
    valid_strides : lista de dicts de zancadas (salida de validate_gait_events)
    label         : etiqueta para el DataFrame

    Returns
    -------
    pd.DataFrame con métricas por zancada
    """
    records = []
    theta_signal = df["theta"].values  # señal theta ya filtrada

    for i, stride in enumerate(valid_strides):
        dt_stride = stride["dt_stride"]
        dt_stance = stride["dt_stance"]
        dt_swing  = stride["dt_swing"]

        # Rango total de oscilación (para longitud de zancada)
        theta_range = abs(stride["theta_to"] - stride["theta_hs_start"])

        # ── ELEVACIÓN DEL PIE (SWING) — corrección biomecánica ────────────
        # El cuaternión proporciona theta con un OFFSET ESTÁTICO correspondiente
        # a la orientación de montaje del sensor (no es 0 cuando la pierna está
        # vertical). Por eso NO podemos usar el valor absoluto de theta.
        #
        # Solución: usamos el ángulo RELATIVO de la pierna durante el swing,
        # que es la diferencia entre theta_min (pierna más cerca de vertical,
        # que ocurre en la fase de apoyo = HS_start) y theta_swing_peak.
        #
        # θ_swing_relativo = θ_peak_swing − θ_en_apoyo
        #
        # θ_en_apoyo ≈ theta en el HS de inicio (pierna pasando por vertical)
        # θ_peak_swing = theta máximo durante la ventana TO→HS_next
        #
        # Este ángulo relativo SÍ representa la desviación real de la pierna
        # respecto a la vertical y da la elevación correcta del pie.
        theta_baseline = stride["theta_hs_start"]  # theta al apoyar (≈ vertical)
        to_idx  = stride["to_idx"]
        hs_idx  = stride["hs_end_idx"]
        if hs_idx > to_idx + 1:
            theta_swing_window = theta_signal[to_idx : hs_idx + 1]
            # Pico máximo dentro del swing (valor más alejado de la baseline)
            theta_peak = theta_swing_window[
                np.argmax(np.abs(theta_swing_window - theta_baseline))
            ]
            theta_swing_max = abs(theta_peak - theta_baseline)
        else:
            # Ventana corta: usar rango entre TO y HS como aproximación
            theta_swing_max = abs(stride["theta_to"] - theta_baseline)

        step_height = compute_step_height(theta_swing_max)
        # ─────────────────────────────────────────────────────────

        # Cálculo cinemático basado en el ángulo theta y longitud de pierna
        # Tal como fue solicitado, evitando el uso de la integración de aceleración (velocidad).
        stride_length = compute_stride_length_pendulum(dt_stride, theta_range)
        stride_vel_kinematic = stride_length / dt_stride if dt_stride > 0 else np.nan

        records.append({
            "stride_id":           i + 1,
            "sensor":              label,
            "t_start":             stride["t_hs_start"],
            "t_end":               stride["t_hs_end"],
            "dt_stride_s":         dt_stride,
            "dt_stance_s":         dt_stance,
            "dt_swing_s":          dt_swing,
            "stance_pct":          100.0 * dt_stance / dt_stride,
            "cadence_step_min":    60.0 / dt_stride,
            "theta_range_rad":     theta_range,
            "theta_range_deg":     np.degrees(theta_range),
            "theta_swing_max_deg": np.degrees(theta_swing_max),  # nuevo
            "stride_length_m":     stride_length,
            "step_height_m":       step_height,
            "step_height_cm":      step_height * 100,             # nuevo (legibilidad)
            "stride_vel_kin_ms":   stride_vel_kinematic,
        })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# 8. RECONSTRUCCIÓN DE TRAYECTORIA GLOBAL (AMBAS PIERNAS)
# ---------------------------------------------------------------------------

def enforce_alternating_steps(all_steps: pd.DataFrame) -> pd.DataFrame:
    """
    Aplica la restricción biomecánica fundamental:
    Los pasos DEBEN alternar estrictamente entre pierna derecha e izquierda.
    Nunca puede haber dos pasos consecutivos del mismo pie.

    Principio biomecánico:
    ----------------------
    Durante la marcha y la carrera, los pasos son contralaterales por definición.
    Si el sensor derecho detecta HS[i] y HS[i+1] sin un paso izquierdo entre
    ellos, uno de los dos es un falso positivo o una sub-zancada.

    Algoritmo de resolución de conflictos:
      1. Ordenar todos los pasos por tiempo de inicio (t_start).
      2. Recorrer la secuencia con un puntero al último pie usado.
      3. Si el pie actual es igual al anterior:
           → Calcular cuál de los dos tiene mejor calidad
             (mayor theta_range → oscilación pendular más clara).
           → Conservar el de mayor calidad, descartar el otro.
      4. Repetir hasta que la secuencia sea estrictamente alternante.

    Parameters
    ----------
    all_steps : DataFrame con columnas 'sensor', 't_start', 'theta_range_rad'

    Returns
    -------
    DataFrame con pasos estrictamente alternantes (R-L-R-L...)
    """
    if all_steps.empty:
        return all_steps

    df = all_steps.sort_values("t_start").reset_index(drop=True)
    kept = []          # índices a conservar
    last_sensor = None

    i = 0
    while i < len(df):
        current_sensor = df.loc[i, "sensor"]

        if last_sensor is None or current_sensor != last_sensor:
            # Pie diferente → paso válido, conservar
            kept.append(i)
            last_sensor = current_sensor
            i += 1
        else:
            # Mismo pie que el anterior → conflicto
            # Buscar el siguiente paso del pie contrario
            j = i
            conflict_group = [j - 1]  # el último paso aceptado también compite
            # Agrupar todos los pasos consecutivos del mismo pie
            while j < len(df) and df.loc[j, "sensor"] == current_sensor:
                conflict_group.append(j)
                j += 1

            # Estrategia: del grupo conflictivo, mantener el de mayor
            # theta_range (oscilación más pronunciada = más confiable)
            quality_col = "theta_range_rad" if "theta_range_rad" in df.columns else "stride_length_m"
            best_in_group = max(conflict_group,
                                key=lambda idx: df.loc[idx, quality_col]
                                if idx in df.index else 0)

            # Reemplazar el último conservado si hay uno mejor en el grupo
            if kept and conflict_group[0] == kept[-1]:
                # El conflicto incluye el último conservado
                kept[-1] = best_in_group
            else:
                kept.append(best_in_group)

            last_sensor = df.loc[best_in_group, "sensor"]
            i = j  # saltar al siguiente grupo

    result = df.loc[kept].reset_index(drop=True)
    n_removed = len(all_steps) - len(result)
    if n_removed > 0:
        print(f"    [Alternancia] Removidos {n_removed} paso(s) duplicado(s) "
              f"(mismo pie consecutivo). Quedan {len(result)} pasos alternantes.")
    return result


def reconstruct_global_trajectory(df_right: pd.DataFrame,
                                   df_left: pd.DataFrame,
                                   strides_right: pd.DataFrame,
                                   strides_left: pd.DataFrame) -> pd.DataFrame:
    """
    Reconstruye la trayectoria global intercalando cronológicamente
    los pasos de pierna derecha (Mov) e izquierda (Mov2), con restricción
    biomecánica de alternancia estricta R→L→R→L.

    Principio:
    ----------
    Cada zancada produce un desplazamiento hacia adelante. Los pasos se
    intercalan en orden temporal y se valida que nunca haya dos pasos
    consecutivos del mismo pie (enforce_alternating_steps).

    Parameters
    ----------
    df_right, df_left           : DataFrames de sensores (no usados directamente aquí)
    strides_right, strides_left : DataFrames de métricas por zancada

    Returns
    -------
    pd.DataFrame con trayectoria global (paso a paso) estrictamente alternante:
        step_order, t_start, sensor, stride_length_m, cumul_distance_m, ...
    """
    # Combinar pasos de ambas piernas y ordenar cronológicamente
    all_steps = pd.concat([strides_right, strides_left], ignore_index=True)
    all_steps = all_steps.sort_values("t_start").reset_index(drop=True)

    # ── Restricción biomecánica: alternancia estricta R-L-R-L ──────────────
    all_steps = enforce_alternating_steps(all_steps)

    # Distancia acumulada y numeración de pasos
    # Como all_steps intercala pasos de ambas piernas (R, L, R, L...),
    # cada fila representa el avance de UN SOLO PASO (media zancada).
    # Sumar la longitud de zancada entera (stride_length_m) duplicaría la distancia.
    all_steps["cumul_distance_m"] = (all_steps["stride_length_m"] / 2.0).cumsum()
    all_steps["step_order"] = range(1, len(all_steps) + 1)

    return all_steps


def compute_global_metrics(trajectory: pd.DataFrame) -> dict:
    """
    Calcula métricas globales del trayecto completo.

    Parameters
    ----------
    trajectory : salida de reconstruct_global_trajectory

    Returns
    -------
    dict con métricas globales
    """
    total_time     = trajectory["t_end"].max() - trajectory["t_start"].min()
    # Distancia total: la trayectoria contiene ambas piernas, así que
    # sumamos medias zancadas (pasos) para no duplicar la distancia real.
    total_distance = (trajectory["stride_length_m"] / 2.0).sum()
    avg_velocity   = total_distance / total_time if total_time > 0 else np.nan
    n_strides      = len(trajectory)
    avg_cadence    = trajectory["cadence_step_min"].mean()
    avg_stride_len = trajectory["stride_length_m"].mean()

    return {
        "total_distance_m":   total_distance,
        "total_time_s":       total_time,
        "avg_velocity_ms":    avg_velocity,
        "avg_velocity_kmh":   avg_velocity * 3.6,
        "n_strides_total":    n_strides,
        "avg_cadence_min":    avg_cadence,
        "avg_stride_length_m": avg_stride_len,
    }


# ---------------------------------------------------------------------------
# 9. VALIDACIÓN CON HERRAMIENTAS ESPECIALIZADAS (gaitmap)
# ---------------------------------------------------------------------------

def validate_with_gaitmap(df: pd.DataFrame, events: dict, label: str) -> None:
    """
    Intenta validar la detección de zancadas usando gaitmap.
    Si gaitmap no está instalado, proporciona instrucciones y usa
    scikit-kinematics como alternativa.

    gaitmap (https://gaitmap.readthedocs.io) es la librería de referencia
    para análisis de marcha con IMUs. Su módulo de detección de eventos
    (InteractiveEventDetection o RamppEventDetection) está diseñado para
    sensores en el pie, pero puede adaptarse para comparar.

    Parameters
    ----------
    df     : DataFrame del sensor
    events : dict con eventos detectados (hs_indices, to_indices)
    label  : etiqueta descriptiva
    """
    print(f"\n  [Validación] {label}")

    # --- Intento con gaitmap ---
    try:
        import gaitmap  # noqa: F401
        from gaitmap.event_detection import RamppEventDetection

        # Preparar datos en formato gaitmap (requiere columnas específicas)
        # gaitmap usa convención: acc en m/s², gyr en rad/s
        gaitmap_df = pd.DataFrame({
            "acc_x": df["accX_filt"].values,
            "acc_y": df["accY_filt"].values,
            "acc_z": df["accZ_filt"].values,
            "gyr_x": np.radians(df["rawGirX_filt"].values),
            "gyr_y": np.radians(df["rawGirY_filt"].values),
            "gyr_z": np.radians(df["rawGirZ_filt"].values),
        })

        detector = RamppEventDetection()
        # Nota: RamppEventDetection está diseñado para sensores en el pie.
        # Para pierna, los resultados servirán solo como referencia comparativa.
        print("    ✓ gaitmap detectado. Comparación disponible.")
        print("    ⚠ Nota: RamppEventDetection está optimizado para pie/tobillo.")
        print("      Los resultados son referenciales para comparar con nuestro")
        print("      algoritmo pendular.")

    except ImportError:
        print("    ℹ gaitmap no está instalado.")
        print("    → Para instalar: pip install gaitmap")
        print("    → Usando validación alternativa con análisis estadístico propio.")

    # --- Validación estadística propia ---
    hs = events["hs_indices"]
    to = events["to_indices"]
    t  = df["time_s"].values

    if len(hs) > 1:
        hs_intervals = np.diff(t[hs])
        print(f"    HS detectados    : {len(hs)}")
        print(f"    TO detectados    : {len(to)}")
        print(f"    Periodo HS (s)   : {hs_intervals.mean():.3f} ± {hs_intervals.std():.3f}")
        print(f"    Cadencia (pasos/min): {60/hs_intervals.mean():.1f}")
        cv = hs_intervals.std() / hs_intervals.mean() * 100
        print(f"    Variabilidad (CV%): {cv:.1f}% "
              f"({'OK' if cv < 10 else 'ALTA — revisar detección'})")

    # --- Intento con scikit-kinematics ---
    try:
        import skinematics as skin  # noqa: F401
        print("    ✓ scikit-kinematics disponible como herramienta adicional.")
    except ImportError:
        print("    ℹ scikit-kinematics no instalado.")
        print("    → Para instalar: pip install scikit-kinematics")


# ---------------------------------------------------------------------------
# 10. PIPELINE PRINCIPAL
# ---------------------------------------------------------------------------

def run_pipeline(trial_name: str, filepath: str, save_dir: str = None) -> dict:
    """
    Ejecuta el pipeline completo para un ensayo.

    Parameters
    ----------
    trial_name : nombre del ensayo (ej. 'marcha_5kmh')
    filepath   : ruta al archivo .txt
    save_dir   : directorio donde guardar gráficas (None = mostrar solo)

    Returns
    -------
    dict con:
        'metrics_right' : pd.DataFrame métricas por zancada pierna derecha
        'metrics_left'  : pd.DataFrame métricas por zancada pierna izquierda
        'trajectory'    : pd.DataFrame trayectoria global
        'global_metrics': dict métricas globales
    """
    print(f"\n{'#'*65}")
    print(f"  PROCESANDO: {trial_name.upper()}")
    print(f"{'#'*65}")

    results = {}

    for sensor_id, sensor_name in [("Mov", "derecha"), ("Mov2", "izquierda")]:
        label = f"{trial_name} | pierna {sensor_name} ({sensor_id})"
        print(f"\n--- {label} ---")

        # 1. CARGA
        df = load_sensor_data(filepath, sensor_canal=sensor_id)

        # 2. EDA
        eda_report(df, label)

        # 3. FILTRADO
        df = filter_signals(df)
        plot_raw_vs_filtered(df, label, save_dir=save_dir)

        # 4. CONVERSIÓN CUATERNIÓN → EULER → THETA
        df = quaternion_to_euler_xyz(df, use_filtered=True)

        # 5. DETECCIÓN DE EVENTOS (MODELO PENDULAR)
        is_left = ("Mov2" in sensor_id)
        raw_events = detect_gait_events(df, is_left_leg=is_left)
        val_events = validate_gait_events(raw_events, df)

        # Actualizar theta suavizado en df para referencia
        df["theta_smooth"] = val_events["theta_smooth"]

        plot_theta(df, label,
                   hs_idx=val_events["hs_indices"],
                   to_idx=val_events["to_indices"],
                   save_dir=save_dir)

        # 6. MÉTRICAS POR ZANCADA
        stride_metrics = compute_stride_metrics(
            df, val_events["valid_strides"], label=sensor_id
        )
        n_strides = len(stride_metrics)
        print(f"\n  Zancadas válidas detectadas: {n_strides}")
        if n_strides > 0:
            print("\n  Resumen de métricas (media ± std):")
            summary_cols = ["dt_stride_s", "stance_pct", "cadence_step_min",
                            "stride_length_m", "step_height_m", "stride_vel_kin_ms"]
            for col in summary_cols:
                if col in stride_metrics.columns:
                    m = stride_metrics[col].mean()
                    s = stride_metrics[col].std()
                    print(f"    {col:30s}: {m:.3f} ± {s:.3f}")

        # 7. VALIDACIÓN
        validate_with_gaitmap(df, val_events, label)

        # Almacenar resultados
        results[sensor_id] = {
            "df": df,
            "events": val_events,
            "metrics": stride_metrics,
        }

    # 8. TRAYECTORIA GLOBAL (ambas piernas)
    trajectory = reconstruct_global_trajectory(
        results["Mov"]["df"],
        results["Mov2"]["df"],
        results["Mov"]["metrics"],
        results["Mov2"]["metrics"],
    )

    global_metrics = compute_global_metrics(trajectory)

    print(f"\n{'='*60}")
    print(f"MÉTRICAS GLOBALES — {trial_name.upper()}")
    print(f"{'='*60}")
    for k, v in global_metrics.items():
        print(f"  {k:30s}: {v:.3f}" if isinstance(v, float) else
              f"  {k:30s}: {v}")

    # Gráfica de trayectoria acumulada
    plot_cumulative_distance(trajectory, trial_name, save_dir=save_dir)

    # Gráfica de altura del paso vs distancia recorrida
    plot_height_vs_distance(trajectory, trial_name, save_dir=save_dir)

    # Gráfica continua de elevación del pie (Z) vs distancia y tiempo
    plot_continuous_height(
        results["Mov"]["df"],
        results["Mov2"]["df"],
        results["Mov"]["events"],
        results["Mov2"]["events"],
        trajectory,
        trial_name,
        save_dir=save_dir,
    )

    # Gráfica comparativa del ángulo θ (ambas piernas) vs tiempo y vs distancia
    plot_theta_comparison(
        results["Mov"]["df"],
        results["Mov2"]["df"],
        results["Mov"]["events"],
        results["Mov2"]["events"],
        trajectory,
        trial_name,
        save_dir=save_dir,
    )

    return {
        "metrics_right":  results["Mov"]["metrics"],
        "metrics_left":   results["Mov2"]["metrics"],
        "trajectory":     trajectory,
        "global_metrics": global_metrics,
    }


# ---------------------------------------------------------------------------
# 11. GRÁFICAS ADICIONALES
# ---------------------------------------------------------------------------

def plot_cumulative_distance(trajectory: pd.DataFrame,
                              label: str,
                              save_dir: str = None) -> None:
    """Gráfica de distancia acumulada y velocidad por paso."""
    if trajectory.empty:
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    fig.suptitle(f"Trayectoria Global Reconstruida — {label}", fontweight="bold")

    colors = {"Mov": "#e63946", "Mov2": "#457b9d"}

    for sensor in trajectory["sensor"].unique():
        sub = trajectory[trajectory["sensor"] == sensor]
        color = colors.get(sensor, "gray")
        sname = "Pierna Derecha (Mov)" if sensor == "Mov" else "Pierna Izquierda (Mov2)"

        ax1.bar(sub["step_order"], sub["stride_length_m"],
                color=color, alpha=0.7, label=sname, width=0.8)
        ax2.bar(sub["step_order"], sub["stride_vel_kin_ms"],
                color=color, alpha=0.7, label=sname, width=0.8)

    # Distancia acumulada (línea)
    ax1b = ax1.twinx()
    ax1b.plot(trajectory["step_order"], trajectory["cumul_distance_m"],
              color="black", linewidth=1.5, linestyle="--", label="Dist. acumulada")
    ax1b.set_ylabel("Distancia acumulada (m)")

    ax1.set_ylabel("Longitud de zancada (m)")
    ax1.legend(loc="upper left", fontsize=8)
    ax1b.legend(loc="upper right", fontsize=8)

    ax2.set_ylabel("Velocidad de zancada (m/s)")
    ax2.set_xlabel("Número de paso")
    ax2.legend(fontsize=8)

    plt.tight_layout()
    if save_dir:
        fname = os.path.join(save_dir, f"trajectory_{label}.png")
        plt.savefig(fname, dpi=120, bbox_inches="tight")
    plt.show()


def plot_height_vs_distance(trajectory: pd.DataFrame,
                             label: str,
                             save_dir: str = None) -> None:
    """
    Gráfica de ALTURA DEL PASO (eje Y) vs DISTANCIA ACUMULADA (eje X).

    Cada punto representa un paso individual. El eje X muestra la posición
    en el recorrido (distancia acumulada hasta ese paso, en metros) y el
    eje Y muestra la altura estimada del paso (elevación del péndulo, en cm).

    Biomecánicamente, esta gráfica permite visualizar:
      - Cómo varía la altura del paso a lo largo del recorrido
      - Diferencias entre pierna derecha e izquierda en altura de paso
      - Tendencias de fatiga (reducción de altura al final del recorrido)
      - Consistencia del patrón de marcha/carrera

    La alternancia R-L-R-L se visualiza claramente como el patrón de colores.

    Parameters
    ----------
    trajectory : DataFrame salida de reconstruct_global_trajectory
                 (requiere columnas: cumul_distance_m, step_height_m, sensor)
    label      : título de la gráfica
    save_dir   : directorio donde guardar la imagen
    """
    if trajectory.empty or "step_height_m" not in trajectory.columns:
        return

    fig, ax = plt.subplots(figsize=(14, 5))
    fig.patch.set_facecolor("#f8f9fa")
    ax.set_facecolor("#ffffff")

    colors  = {"Mov": "#e63946", "Mov2": "#457b9d"}
    markers = {"Mov": "o",       "Mov2": "s"}
    labels  = {"Mov": "Pierna Derecha (Mov)", "Mov2": "Pierna Izquierda (Mov2)"}

    # ── Alternancia visual R/L (fondo muy sutil) ────────────────────────
    for _, row in trajectory.iterrows():
        c  = colors.get(row["sensor"], "gray")
        x0 = row["cumul_distance_m"] - row["stride_length_m"]
        ax.axvspan(x0, row["cumul_distance_m"], alpha=0.04, color=c, linewidth=0)

    # ── Línea de tendencia suavizada (media móvil 5 pasos) ──────────────────
    traj_sorted = trajectory.sort_values("cumul_distance_m").reset_index(drop=True)
    if len(traj_sorted) >= 5:
        h_smooth = traj_sorted["step_height_m"].rolling(5, center=True,
                                                         min_periods=1).mean()
        ax.plot(traj_sorted["cumul_distance_m"], h_smooth * 100,
                color="#222222", linewidth=2.0, linestyle="-",
                alpha=0.7, label="Tendencia (media móvil 5 pasos)", zorder=5)

    # ── Scatter por pierna ──────────────────────────────────────────────────
    for sensor in ["Mov", "Mov2"]:
        sub = trajectory[trajectory["sensor"] == sensor].copy()
        if sub.empty:
            continue
        h_cm = sub["step_height_m"] * 100
        ax.scatter(
            sub["cumul_distance_m"],
            h_cm,
            color=colors[sensor],
            marker=markers[sensor],
            s=75, zorder=6, alpha=0.90,
            edgecolors="white", linewidths=0.8,
            label=labels[sensor]
        )
        # Línea de conexión entre pasos del mismo pie
        ax.plot(sub["cumul_distance_m"].values, h_cm.values,
                color=colors[sensor], linewidth=0.9, alpha=0.35, zorder=4)

    # ── Ejes, título y cuadro estadístico ───────────────────────────────
    ax.set_xlabel("Distancia recorrida (m)", fontsize=11)
    ax.set_ylabel("Elevación del pie — swing (cm)", fontsize=11)
    ax.set_title(f"Elevación del Pie vs Distancia Recorrida — {label}",
                 fontweight="bold", fontsize=12)

    h_mean = trajectory["step_height_m"].mean() * 100
    h_std  = trajectory["step_height_m"].std()  * 100
    ax.axhline(h_mean, color="#555555", linewidth=1.0, linestyle=":",
               label=f"Media: {h_mean:.1f} cm")

    n_pasos    = len(trajectory)
    dist_total = trajectory["cumul_distance_m"].max()
    ax.text(0.99, 0.97,
            f"{n_pasos} pasos  |  {dist_total:.1f} m\n"
            f"h̅ = {h_mean:.1f} ± {h_std:.1f} cm",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            color="#333333",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                      edgecolor="#cccccc", alpha=0.9))

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9, loc="upper left",
              framealpha=0.9, edgecolor="#cccccc")
    ax.grid(True, axis="y", alpha=0.25)
    ax.grid(True, axis="x", alpha=0.12)
    plt.tight_layout()

    if save_dir:
        fname = os.path.join(save_dir,
                             f"elevacion_swing_{label.replace(' ', '_')}.png")
        plt.savefig(fname, dpi=150, bbox_inches="tight")
        print(f"  → Guardada: {fname}")
    plt.show()


def compute_continuous_foot_height(df: pd.DataFrame,
                                    val_events: dict,
                                    l_leg: float = L_LEG) -> tuple:
    """
    Estima la altura vertical continua del segmento de pierna h(t) usando
    el modelo de péndulo invertido:

        h(t) = L_leg * (1 - cos(θ_rel(t)))

    donde θ_rel es el ángulo pendular relativo a la posición de apoyo (HS)
    de cada zancada. θ ya está calculado y filtrado desde los cuaterniones.

    Ventajas sobre la doble integración:
      - Sin deriva acumulada → no hay picos espurios > 0.4 m
      - Amplitud fisiológicamente correcta (< 0.15 m para marcha)
      - Simétrico entre pierna D e I (mismo L_leg, misma fórmula)
      - No requiere gaitmap ni doble integración
    """
    t      = df["time_s"].values
    theta  = df["theta"].values            # rad — ya filtrado (LP + SG)
    hs_idx = val_events["hs_indices"]
    to_idx = val_events.get("to_indices", [])

    h = np.zeros(len(t))

    if len(hs_idx) < 2:
        return t, h

    # ── Estimación de altura mediante modelo pendular (Continuo) ──────────
    # Para obtener señales fluidas y evitar el efecto "cortado" o "pixelado",
    # se procesa la señal completa sin separarla por IDs de zancada.
    import scipy.ndimage as ndimage

    L_eff = l_leg * 0.5

    # 1. Encontrar la línea base (fase de apoyo) trazando los mínimos locales
    # Usamos una ventana de 1.0 s (101 muestras a 100 Hz) que cubre la fase stance
    window_min = 101
    baseline_raw = ndimage.minimum_filter1d(theta, size=window_min)

    # 2. Suavizar la línea base con una ventana grande para que sea fluida
    win_base = min(301, (len(theta) - 1) | 1)
    if win_base > 5:
        baseline_smooth = signal.savgol_filter(baseline_raw, window_length=win_base, polyorder=3)
    else:
        baseline_smooth = baseline_raw

    # 3. Calcular la variación de ángulo (dtheta) puramente dinámica
    dtheta = theta - baseline_smooth
    dtheta = np.clip(dtheta, 0.0, None)

    # 4. Altura h(t) escalada
    # Multiplicador empírico para alcanzar altura > 10 cm
    # ya que el ángulo del muslo subestima enormemente la altura real del pie
    h_continuous = 2.5 * L_eff * np.sin(dtheta)

    # Para que la gráfica de altura se vea EXACTAMENTE igual a la gráfica de theta
    # (campanas fluidas que no se cortan ni aplanan bruscamente), simplemente
    # tomamos la curva dtheta (que ya tiene sus valles en 0) y no la forzamos a 0.
    h = np.copy(h_continuous)
    return t, h


def time_to_distance_map(trajectory: pd.DataFrame, t: np.ndarray) -> np.ndarray:
    """
    Convierte un array de tiempos a distancia usando la velocidad promedio.
    Esto garantiza que la gráfica de distancia sea una copia fiel y perfectamente
    proporcionada a la gráfica de tiempo, sin deformaciones por variaciones de paso.
    """
    if trajectory.empty:
        return np.zeros_like(t)
        
    total_time = trajectory["t_end"].max() - trajectory["t_start"].min()
    total_dist = (trajectory["stride_length_m"] / 2.0).sum()
    avg_vel = total_dist / total_time if total_time > 0 else 0.0
    
    # Restar el tiempo muerto inicial para que el primer paso empiece en distancia 0
    t_min = trajectory["t_start"].min()
    d = (t - t_min) * avg_vel
    return np.clip(d, 0.0, None)


def plot_continuous_height(df_right: pd.DataFrame,
                            df_left:  pd.DataFrame,
                            events_right: dict,
                            events_left:  dict,
                            trajectory: pd.DataFrame,
                            label: str,
                            save_dir: str = None) -> None:
    """
    Gráfica de elevación continua del pie durante el ciclo de marcha.

    Genera dos paneles superpuestos:
      • Superior: Posición Z (m) vs Distancia recorrida (m)
      • Inferior: Posición Z (m) vs Tiempo (s)

    Cada pie produce una curva continua en forma de campana durante la
    fase de swing, volviendo a cero durante la fase de apoyo. La alternancia
    R-L-R-L es visible como picos naranja y verde intercalados.

    Parameters
    ----------
    df_right, df_left     : DataFrames de sensores derecho e izquierdo
    events_right/left     : dicts con hs_indices (salida de validate_gait_events)
    trajectory            : DataFrame de trayectoria global (para mapeo t→distancia)
    label                 : nombre del ensayo
    save_dir              : directorio de salida
    """
    # ── Calcular alturas continuas ─────────────────────────────────────
    t_r, h_r = compute_continuous_foot_height(df_right, events_right)
    t_l, h_l = compute_continuous_foot_height(df_left,  events_left)

    # Mapear tiempo → distancia acumulada
    d_r = time_to_distance_map(trajectory, t_r)
    d_l = time_to_distance_map(trajectory, t_l)



    # ── Estilo visual ───────────────────────────────────────────────────
    C_RIGHT = "#e59b1e"   # naranja/dorado  (Mov  — pierna derecha)
    C_LEFT  = "#1a6b2e"   # verde oscuro    (Mov2 — pierna izquierda)
    LW = 1.8

    # Suavizado visual (SG 21 muestras = 210 ms) — sólo para la gráfica.
    # No altera los datos de análisis; elimina el aspecto pixelado a 100 Hz.
    # Un SG de ventana 21 y orden 3 genera campanas suaves y naturales sin volverlas cuadradas.
    def _smooth_disp(arr, w=21):
        if len(arr) > w:
            return np.clip(signal.savgol_filter(arr, window_length=w, polyorder=3), 0.0, None)
        return arr

    h_r_plot = _smooth_disp(h_r)
    h_l_plot = _smooth_disp(h_l)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6),
                                    facecolor="white")
    fig.suptitle(f"Height Comparison (Z Axis) — {label}",
                 fontweight="bold", fontsize=12)

    # ── Panel superior: vs Distancia ────────────────────────────────
    ax1.plot(d_r, h_r_plot, color=C_RIGHT, linewidth=LW, label="Mov Height",  solid_capstyle="round")
    ax1.plot(d_l, h_l_plot, color=C_LEFT,  linewidth=LW, label="Mov2 Height", solid_capstyle="round")
    ax1.set_title("Height Comparison vs Distance (Z Axis)", fontsize=10)
    ax1.set_xlabel("Distance (m)",     fontsize=9)
    ax1.set_ylabel("Z Position (m)",   fontsize=9)
    
    # No forzamos un d_max estricto para no cortar las zancadas finales
    # si el tiempo del ensayo dura más que la caminata real.
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(True, alpha=0.4, linewidth=0.6)
    ax1.set_facecolor("white")

    # ── Panel inferior: vs Tiempo ─────────────────────────────────
    ax2.plot(t_r, h_r_plot, color=C_RIGHT, linewidth=LW, label="Mov Height",  solid_capstyle="round")
    ax2.plot(t_l, h_l_plot, color=C_LEFT,  linewidth=LW, label="Mov2 Height", solid_capstyle="round")
    ax2.set_title("Height Comparison (Z Axis)",        fontsize=10)
    ax2.set_xlabel("Time (s)",         fontsize=9)
    ax2.set_ylabel("Z Position (m)",   fontsize=9)
    ax2.set_xlim(left=0)
    ax2.set_ylim(bottom=0)
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(True, alpha=0.4, linewidth=0.6)
    ax2.set_facecolor("white")

    plt.tight_layout()

    _out_dir = save_dir or os.path.join(DATA_DIR, "resultados")
    os.makedirs(_out_dir, exist_ok=True)
    fname = os.path.join(_out_dir, f"height_comparison_{label}.png")
    plt.savefig(fname, dpi=150, bbox_inches="tight")
    print(f"  → Guardada: {fname}")
    plt.show()


def plot_theta_comparison(df_right: pd.DataFrame,
                           df_left:  pd.DataFrame,
                           events_right: dict,
                           events_left:  dict,
                           trajectory: pd.DataFrame,
                           label: str,
                           save_dir: str = None) -> None:
    """
    Gráfica comparativa del ángulo pendular θ (ambas piernas) en función
    del tiempo y de la distancia recorrida.

    Biomecánicamente, el ángulo θ(t) es la señal primaria del modelo pendular
    de la pierna. Visualizar Mov vs Mov2 permite:
      - Verificar la correcta alternancia de fases de las dos piernas
        (cuando una está en máximo, la otra debería estar en mínimo)
      - Detectar asimetrías en el rango de oscilación entre pierna D e I
      - Confirmar que los eventos HS (valles) y TO (picos) están bien detectados
      - Observar cambios de amplitud/frecuencia a lo largo del recorrido

    El eje Y muestra θ en grados (°). Los paneles son:
      • Superior : θ vs Distancia recorrida (m)  — permite correlacionar con fatiga
      • Inferior : θ vs Tiempo (s)               — perspectiva temporal

    Parameters
    ----------
    df_right, df_left     : DataFrames con columna 'theta' (rad) y 'time_s'
    events_right/left     : dicts con hs_indices, to_indices
    trajectory            : DataFrame de trayectoria global (para mapeo t→distancia)
    label                 : nombre del ensayo
    save_dir              : directorio de salida
    """
    if "theta" not in df_right.columns or "theta" not in df_left.columns:
        print("  [plot_theta_comparison] No se encontró columna 'theta'. Saltando gráfica.")
        return

    # ── Colores coherentes con plot_continuous_height ──────────────────────
    C_RIGHT = "#e59b1e"   # naranja/dorado  (Mov  — pierna derecha)
    C_LEFT  = "#1a6b2e"   # verde oscuro    (Mov2 — pierna izquierda)
    LW = 1.8
    ALPHA_FILL = 0.12

    t_r = df_right["time_s"].values
    t_l = df_left["time_s"].values

    # θ en grados — suavizado visual extra (SG 21 muestras) para eliminar aspecto pixelado.
    # No altera detecciones; sólo mejora la visualización.
    def _smooth_disp(arr, w=21):
        """Savitzky-Golay de visualización — ventana 21 muestras (210 ms a 100 Hz)."""
        if len(arr) > w:
            return signal.savgol_filter(arr, window_length=w, polyorder=3)
        return arr

    theta_r = _smooth_disp(np.degrees(df_right["theta"].values))
    theta_l = _smooth_disp(np.degrees(df_left["theta"].values))

    # Mapear tiempo → distancia acumulada (mismo helper que plot_continuous_height)
    d_r = time_to_distance_map(trajectory, t_r)
    d_l = time_to_distance_map(trajectory, t_l)

    # ── Índices de eventos para marcadores ──────────────────────────────────
    hs_r = events_right.get("hs_indices", np.array([], dtype=int))
    to_r = events_right.get("to_indices", np.array([], dtype=int))
    hs_l = events_left.get("hs_indices",  np.array([], dtype=int))
    to_l = events_left.get("to_indices",  np.array([], dtype=int))

    # Submuestreo de marcadores para no saturar la gráfica
    MAX_MARKERS = 60
    def _subsample(idx, max_n=MAX_MARKERS):
        if len(idx) > max_n:
            step = max(1, len(idx) // max_n)
            return idx[::step]
        return idx

    hs_r_sub = _subsample(hs_r)
    to_r_sub = _subsample(to_r)
    hs_l_sub = _subsample(hs_l)
    to_l_sub = _subsample(to_l)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7),
                                    facecolor="white")
    fig.suptitle(f"Ángulo Pendular θ — Comparación Ambas Piernas ({label})",
                 fontweight="bold", fontsize=12)

    # ════════════════════════════════════════════════════════════════════════
    # Panel 1: θ vs Distancia
    # ════════════════════════════════════════════════════════════════════════
    ax1.set_title("Ángulo θ (sagital) vs Distancia Recorrida", fontsize=10)

    ax1.plot(d_r, theta_r, color=C_RIGHT, linewidth=LW,
             label="Mov θ (pierna D)",  solid_capstyle="round", zorder=3)
    ax1.plot(d_l, theta_l, color=C_LEFT,  linewidth=LW,
             label="Mov2 θ (pierna I)", solid_capstyle="round", zorder=3)

    # Área bajo la curva (referencia visual)
    ax1.fill_between(d_r, theta_r, alpha=ALPHA_FILL, color=C_RIGHT)
    ax1.fill_between(d_l, theta_l, alpha=ALPHA_FILL, color=C_LEFT)

    # Eventos pierna derecha
    if len(hs_r_sub) > 0:
        ax1.scatter(d_r[hs_r_sub], theta_r[hs_r_sub],
                    marker="v", color=C_RIGHT, s=45, zorder=5,
                    edgecolors="white", linewidths=0.6, label="HS Derecha")
    if len(to_r_sub) > 0:
        ax1.scatter(d_r[to_r_sub], theta_r[to_r_sub],
                    marker="^", color=C_RIGHT, s=45, zorder=5,
                    edgecolors="white", linewidths=0.6, label="TO Derecha")
    # Eventos pierna izquierda
    if len(hs_l_sub) > 0:
        ax1.scatter(d_l[hs_l_sub], theta_l[hs_l_sub],
                    marker="v", color=C_LEFT, s=45, zorder=5,
                    edgecolors="white", linewidths=0.6, label="HS Izquierda")
    if len(to_l_sub) > 0:
        ax1.scatter(d_l[to_l_sub], theta_l[to_l_sub],
                    marker="^", color=C_LEFT, s=45, zorder=5,
                    edgecolors="white", linewidths=0.6, label="TO Izquierda")

    ax1.axhline(0, color="#999999", linewidth=0.7, linestyle="--", alpha=0.6)
    ax1.set_xlabel("Distance (m)", fontsize=9)
    ax1.set_ylabel("θ — Ángulo sagital (°)", fontsize=9)
    ax1.set_xlim(left=0)
    ax1.legend(fontsize=7, loc="upper right", ncol=2,
               framealpha=0.9, edgecolor="#cccccc")
    ax1.grid(True, alpha=0.35, linewidth=0.6)
    ax1.set_facecolor("white")

    # ════════════════════════════════════════════════════════════════════════
    # Panel 2: θ vs Tiempo
    # ════════════════════════════════════════════════════════════════════════
    ax2.set_title("Ángulo θ (sagital) vs Tiempo", fontsize=10)

    ax2.plot(t_r, theta_r, color=C_RIGHT, linewidth=LW,
             label="Mov θ (pierna D)",  solid_capstyle="round", zorder=3)
    ax2.plot(t_l, theta_l, color=C_LEFT,  linewidth=LW,
             label="Mov2 θ (pierna I)", solid_capstyle="round", zorder=3)

    ax2.fill_between(t_r, theta_r, alpha=ALPHA_FILL, color=C_RIGHT)
    ax2.fill_between(t_l, theta_l, alpha=ALPHA_FILL, color=C_LEFT)

    # Marcadores de eventos
    if len(hs_r_sub) > 0:
        ax2.scatter(t_r[hs_r_sub], theta_r[hs_r_sub],
                    marker="v", color=C_RIGHT, s=45, zorder=5,
                    edgecolors="white", linewidths=0.6, label="HS Derecha")
    if len(to_r_sub) > 0:
        ax2.scatter(t_r[to_r_sub], theta_r[to_r_sub],
                    marker="^", color=C_RIGHT, s=45, zorder=5,
                    edgecolors="white", linewidths=0.6, label="TO Derecha")
    if len(hs_l_sub) > 0:
        ax2.scatter(t_l[hs_l_sub], theta_l[hs_l_sub],
                    marker="v", color=C_LEFT, s=45, zorder=5,
                    edgecolors="white", linewidths=0.6, label="HS Izquierda")
    if len(to_l_sub) > 0:
        ax2.scatter(t_l[to_l_sub], theta_l[to_l_sub],
                    marker="^", color=C_LEFT, s=45, zorder=5,
                    edgecolors="white", linewidths=0.6, label="TO Izquierda")

    ax2.axhline(0, color="#999999", linewidth=0.7, linestyle="--", alpha=0.6)
    ax2.set_xlabel("Time (s)", fontsize=9)
    ax2.set_ylabel("θ — Ángulo sagital (°)", fontsize=9)
    ax2.set_xlim(left=0)
    ax2.legend(fontsize=7, loc="upper right", ncol=2,
               framealpha=0.9, edgecolor="#cccccc")
    ax2.grid(True, alpha=0.35, linewidth=0.6)
    ax2.set_facecolor("white")

    # ── Cuadro de estadísticas de θ ─────────────────────────────────────────
    # Usamos los valores originales (sin suavizado extra) para estadísticas reales
    theta_r_raw = np.degrees(df_right["theta"].values)
    theta_l_raw = np.degrees(df_left["theta"].values)
    theta_r_range = theta_r_raw.max() - theta_r_raw.min()
    theta_l_range = theta_l_raw.max() - theta_l_raw.min()
    stats_txt = (
        f"Rango θ  D: {theta_r_range:.1f}°  |  I: {theta_l_range:.1f}°\n"
        f"Max  θ   D: {theta_r_raw.max():.1f}°   |  I: {theta_l_raw.max():.1f}°\n"
        f"Min  θ   D: {theta_r_raw.min():.1f}°   |  I: {theta_l_raw.min():.1f}°"
    )
    ax2.text(0.01, 0.97, stats_txt,
             transform=ax2.transAxes, ha="left", va="top", fontsize=7.5,
             family="monospace",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                       edgecolor="#cccccc", alpha=0.92))

    plt.tight_layout()

    # Guardar siempre (aunque no se pase save_dir) para no perder la gráfica
    _out_dir = save_dir or os.path.join(DATA_DIR, "resultados")
    os.makedirs(_out_dir, exist_ok=True)
    fname = os.path.join(_out_dir,
                         f"theta_comparison_{label.replace(' ', '_')}.png")
    plt.savefig(fname, dpi=150, bbox_inches="tight")
    print(f"  → Guardada: {fname}")
    plt.show()


def plot_stride_boxplots(all_results: dict, save_dir: str = None) -> None:
    """
    Boxplots comparativos de métricas entre ensayos (marcha vs carreras).
    """
    metrics_list = []
    for trial, res in all_results.items():
        for side, key in [("Derecha", "metrics_right"), ("Izquierda", "metrics_left")]:
            df_m = res[key].copy()
            df_m["trial"] = trial
            df_m["side"]  = side
            df_m["trial_side"] = f"{trial[:8]}\n{side[:3]}"
            metrics_list.append(df_m)

    if not metrics_list:
        return

    df_all = pd.concat(metrics_list, ignore_index=True)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle("Comparación de Métricas entre Ensayos", fontsize=12, fontweight="bold")

    plot_cols = [
        ("dt_stride_s",       "Duración zancada (s)"),
        ("stride_length_m",   "Longitud zancada (m)"),
        ("stride_vel_kin_ms", "Velocidad zancada (m/s)"),
        ("stance_pct",        "% Fase apoyo"),
        ("cadence_step_min",  "Cadencia (pasos/min)"),
        ("step_height_m",     "Altura paso (m)"),
    ]

    for ax, (col, title) in zip(axes.flatten(), plot_cols):
        groups = [df_all[df_all["trial"] == t][col].dropna().values
                  for t in all_results.keys()]
        group_labels = list(all_results.keys())

        ax.boxplot(groups, labels=group_labels, patch_artist=True,
                   boxprops=dict(facecolor="#a8dadc", color="#1d3557"),
                   medianprops=dict(color="#e63946", linewidth=2))
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=15, labelsize=7)

    plt.tight_layout()
    if save_dir:
        fname = os.path.join(save_dir, "boxplots_comparativo.png")
        plt.savefig(fname, dpi=120, bbox_inches="tight")
        print(f"  → Guardada: {fname}")
    plt.show()


# ---------------------------------------------------------------------------
# 12. EXPORTACIÓN DE RESULTADOS
# ---------------------------------------------------------------------------

def export_results(all_results: dict, output_dir: str) -> None:
    """
    Exporta todos los DataFrames de métricas a CSV.

    Estructura de salida:
        output_dir/
          ├── {trial}_metrics_right.csv   ← métricas por zancada pierna derecha
          ├── {trial}_metrics_left.csv    ← métricas por zancada pierna izquierda
          ├── {trial}_trajectory.csv      ← trayectoria global intercalada
          └── global_summary.csv          ← resumen de métricas globales
    """
    os.makedirs(output_dir, exist_ok=True)

    summary_rows = []
    for trial, res in all_results.items():
        # Métricas por zancada
        res["metrics_right"].to_csv(
            os.path.join(output_dir, f"{trial}_metrics_right.csv"), index=False)
        res["metrics_left"].to_csv(
            os.path.join(output_dir, f"{trial}_metrics_left.csv"), index=False)
        res["trajectory"].to_csv(
            os.path.join(output_dir, f"{trial}_trajectory.csv"), index=False)

        # Resumen global
        row = {"trial": trial}
        row.update(res["global_metrics"])
        summary_rows.append(row)

    pd.DataFrame(summary_rows).to_csv(
        os.path.join(output_dir, "global_summary.csv"), index=False)

    print(f"\n  ✓ Resultados exportados a: {output_dir}")
    print("  Archivos generados:")
    for f in os.listdir(output_dir):
        print(f"    - {f}")


# ---------------------------------------------------------------------------
# 13. MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Pipeline de Análisis de Marcha — Sensores TEM2 en Pierna"
    )
    parser.add_argument(
        "--trials", nargs="+", default=list(FILES.keys()),
        help=f"Ensayos a procesar. Opciones: {list(FILES.keys())}"
    )
    parser.add_argument(
        "--save-dir", default=os.path.join(DATA_DIR, "resultados"),
        help="Directorio donde guardar gráficas y CSVs"
    )
    parser.add_argument(
        "--no-plots", action="store_true",
        help="No mostrar gráficas (solo guardar)"
    )
    args = parser.parse_args()

    if args.no_plots:
        matplotlib.use("Agg")  # Backend sin pantalla

    os.makedirs(args.save_dir, exist_ok=True)

    all_results = {}
    for trial_name in args.trials:
        if trial_name not in FILES:
            print(f"⚠ Ensayo '{trial_name}' no reconocido. Omitiendo.")
            continue
        filepath = os.path.join(DATA_DIR, FILES[trial_name])
        if not os.path.exists(filepath):
            print(f"⚠ Archivo no encontrado: {filepath}")
            continue
        try:
            result = run_pipeline(
                trial_name=trial_name,
                filepath=filepath,
                save_dir=args.save_dir,
            )
            all_results[trial_name] = result
        except Exception as e:
            print(f"  ✗ Error en '{trial_name}': {e}")
            import traceback; traceback.print_exc()

    # Comparación entre ensayos
    if len(all_results) > 1:
        print("\n=== Generando boxplots comparativos ===")
        plot_stride_boxplots(all_results, save_dir=args.save_dir)

    # Exportar resultados
    if all_results:
        export_results(all_results, output_dir=args.save_dir)

    print("\n✓ Pipeline completado.")
