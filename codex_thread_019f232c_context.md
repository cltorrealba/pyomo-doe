# Antecedente de conversación Codex

## Metadatos

- **Thread ID:** `019f232c-2411-7ca1-a59b-4643dad134b0`
- **Transcript:** `C:\Users\ctorrealba\.codex\sessions\2026\07\02\rollout-2026-07-02T10-12-08-019f232c-2411-7ca1-a59b-4643dad134b0.jsonl`
- **Fork de:** `019e9a11-c081-77c3-af6e-9198e048999a`
- **Repositorio:** `cltorrealba/pyomo-doe`
- **Nota:** el transcript fue leído en modo de solo lectura y no fue modificado.

## 1. Objetivo del trabajo

Procesar y auditar los datos experimentales reales de Laboratorio 2026, especialmente los lotes MBDoE 1 y 2, integrando:

- CO2 online y temperatura;
- analítica Y15;
- Oculyze;
- etanol manual;
- eventos operacionales;
- pérdidas de volumen por muestreo.

El trabajo continuó con análisis de calibración e identificabilidad separados por mosto natural/sintético, incorporación del Lote 2 y preparación de material para un agente externo y un kick-off de filtro de Kalman.

## 2. Decisiones tomadas

### Lote 1

- `t=0` se definió usando la primera muestra Oculyze/inoculación real.
- CO2 se cargó desde `CO2_FILT` y `flow_filt_sccm`.
- Y15 se deduplicó y pivotó; se calcularon fructosa y YAN.
- Biomasa se convirtió con `30 pg/célula`.
- PAN/YAN negativos se conservaron como crudos y se censuraron a cero para modelamiento.
- Para F2 se impuso una regla explícita:

  - datos válidos hasta aproximadamente 50 h;
  - YAN para modelo igual a cero entre 50 h y el pulso real en 72.5 h;
  - valores crudos preservados para trazabilidad.

- Volumen extraído:

  - muestra pequeña: `5 mL`;
  - muestra con etanol: `50 mL`.

- CO2 se corrigió por volumen:

  `CO2_corr_2L(t) = CO2_raw(t) * V0 / V(t)`.

- Se distinguieron inputs planificados y ejecutados.
- F3 recibió realmente `0.015 g/L YAN`, no `0.05`; el ratio real/planificado fue `0.3`.

### Comparación por medio y Lote 2

- Se construyeron análisis históricos separados para mosto natural y sintético.
- Los holdouts y perfiles likelihood permanecieron fuera de la calibración correspondiente.
- Para el Lote 2 se generó primero una planilla de etanol y después se integraron Y15, Oculyze, CO2 y pérdidas de volumen.
- CO2 se mantuvo como señal diagnóstica; no se forzó una ecuación flujo–CO2 no validada.

## 3. Resultados alcanzados

### Preview y QC Lote 1

- Oculyze: F1 19 muestras, F2 21, F3 18.
- Y15: F1 19, F2 22, F3 19.
- Volumen final estimado:

  - F1: `1365 mL`;
  - F2: `1260 mL`;
  - F3: `1365 mL`.

- Corrección máxima de CO2:

  - F1/F3: `1.465`;
  - F2: `1.527`.

- Se generaron tablas QC para IDs, conflictos Y15, nitrógeno, etanol, volumen e inputs.

### Estimabilidad histórica por matriz

- Natural:

  - 9 fermentaciones;
  - 11 parámetros libres;
  - WSSE/residuo `19.20`;
  - FIM core `10/11`;
  - FIM secundaria `16/17`;
  - `Kd0` no cierra perfil bilateral; `qN` sí.

- Sintético:

  - 10 fermentaciones;
  - WSSE/residuo `17.16`;
  - FIM core `10/11`;
  - FIM secundaria `13/17`;
  - `Kd0` y `qN` cierran bilateralmente.

- El medio sintético tuvo más información cinética local: `logdet 68.01` frente a `54.79`.
- Los nueve parámetros aromáticos siguieron sin soporte por ausencia de observaciones.

### Lote 2

- 44 muestras Y15.
- 37 mediciones de etanol; 36 utilizadas.
- 49 mediciones de biomasa.
- 52 eventos de muestreo con evidencia analítica.
- Volumen final:

  - F1 `1325 mL`;
  - F2 `1325 mL`;
  - F3 `1425 mL`.

- Combinación de 10 fermentaciones históricas sintéticas con tres del Lote 2:

  - core `11/11`;
  - core + secundarios `17/17`;
  - modelo extendido `17/26`;
  - las nueve direcciones faltantes corresponden exactamente a aromas.

- `Kd0` y `qN` resultaron identificables bilateralmente.
- Commit: `268994c` — `Integrate historical medium and Lot 2 estimability workflows`.
- No se realizó push.

### Handoff para filtro de Kalman

Se recomendaron tres notebooks:

1. `fermentation_lot2_data_preview.executed.ipynb`;
2. `fermentation_estimability_synthetic_prior_plus_lot2.executed.ipynb`;
3. `fermentation_estimability_historical_natural.executed.ipynb`.

Cubren datos de entrada, prior calibrado/incertidumbre y transferencia al medio natural.

## 4. Pendientes

1. Confirmar tiempos reales de los pulsos del Lote 2.
2. Revisar el etanol F2 a 168 h.
3. Formalizar la ecuación de observación del sensor CO2 antes de usar el flujo directamente en calibración.
4. Incorporar datos aromáticos; nueve parámetros permanecen no identificables.
5. Revisar errores estructurados en glucosa F1, YAN residual, biomasa muerta tardía, piruvato, acetaldehído y acetato.
6. Corregir el README de `laboratory_2026`, que aún menciona rutas antiguas en lugar de `MBDoE_2026/`.
7. Los datos auxiliares bajo `raw_data/` no estaban incluidos en el commit.

## 5. Archivos relevantes

- `fermentation_model/laboratory_2026/notebooks/fermentation_lot1_data_preview.executed.ipynb`
- `fermentation_model/laboratory_2026/notebooks/fermentation_lot2_data_preview.executed.ipynb`
- `fermentation_model/laboratory_2026/notebooks/fermentation_estimability_historical_natural.executed.ipynb`
- `fermentation_model/laboratory_2026/notebooks/fermentation_estimability_historical_synthetic.executed.ipynb`
- `fermentation_model/laboratory_2026/notebooks/fermentation_estimability_synthetic_prior_plus_lot2.executed.ipynb`
- `fermentation_model/laboratory_2026/run_estimability_historical_by_medium.py`
- `fermentation_model/laboratory_2026/run_lot2_data_preview.py`
- `fermentation_model/laboratory_2026/results/`
- `fermentation_model/data/Laboratorio 2026/MBDoE_2026/`
- `fermentation_model/data/Laboratorio 2026/Vendimia_2026/`
- `fermentation_model/REPOSITORY_MAP.md`

## Estado de continuidad

El estado más reciente es la calibración sintética histórica + Lote 2 con rango completo para el core y secundarios, pero sin soporte para aromas. El commit `268994c` contiene los flujos y resultados principales y quedó local, sin push.
