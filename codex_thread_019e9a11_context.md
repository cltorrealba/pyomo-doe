# Antecedente de conversación Codex

## Metadatos

- **Thread ID:** `019e9a11-c081-77c3-af6e-9198e048999a`
- **Transcript:** `C:\Users\ctorrealba\.codex\sessions\2026\06\05\rollout-2026-06-05T19-15-21-019e9a11-c081-77c3-af6e-9198e048999a.jsonl`
- **Tipo:** conversación raíz; origen de varios forks posteriores.
- **Repositorio:** `cltorrealba/pyomo-doe`
- **Rama de trabajo observada:** `ctorrealba_fermentation`
- **Nota:** el transcript fue leído en modo de solo lectura y no fue modificado.

## 1. Objetivo del trabajo

La conversación comenzó con la reparación del notebook `fermentation_model_calibration_3_effective_reformulation.ipynb`, que fallaba por conflictos de IPOPT. Después evolucionó hacia un programa completo de:

1. calibración del modelo de fermentación;
2. reducción e identificabilidad estructural/práctica;
3. comparación de multistart, PSO y regularización;
4. diseño MBDoE de una campaña experimental;
5. extensión a glicerol, metabolitos secundarios, O2, CO2 y aromas;
6. incorporación de mosto natural y sintético;
7. planificación operacional de nueve fermentaciones en tres lotes;
8. ingestión y reevaluación de los datos reales obtenidos en el Lote 1.

El objetivo científico era obtener un modelo dinámico suficientemente identificable y predictivo para usarlo posteriormente como restricción en un MPCC y para dirigir nuevas campañas experimentales.

## 2. Decisiones tomadas

### Modelo e identificabilidad

- Se retiraron opciones de IPOPT incompatibles y la UQ problemática basada en `k_aug/PyNumero`.
- Se adoptó inicialmente un modelo identificable de siete parámetros libres:

  - `mu0`, `qN`, `betaG0`, `betaF0`, `qEG`, `qEF`, `iG`.

- Se fijaron parámetros cuando presentaban perfiles planos, límites activos, direcciones FIM casi nulas o sensibilidades numéricamente inestables.
- La selección combinó WSSE, AIC/BIC, FIM, eigenvalores, profile likelihood, Laplace bayesiano, PSO, multistart y log-L2.
- El PSO no justificó liberar parámetros adicionales únicamente para mejorar las curvas.

### Datos y extensiones del modelo

- Se incorporaron datos de mosto natural y sintético.
- Los valores ausentes no debían transformarse en ceros.
- Etanol se convirtió correctamente desde `% v/v`; se eliminaron conversiones que producían valores cercanos a `700 g/L`.
- La biomasa se convirtió usando `30 pg/célula`.
- Se añadieron glicerol, piruvato, acetaldehído, acetato, biomasa viable/muerta y un bloque efectivo de O2.
- Para aromas se consideraron producción, partición gas-líquido, stripping con CO2 y condensado.

### Diseño experimental

- El diseño fue principalmente MBDoE basado en candidatos: simulación de candidatos, cálculo de FIM y selección greedy/híbrida.
- Se distinguió de una optimización continua pura con IPOPT.
- La campaña operacional consideró:

  - nueve fermentaciones en tres lotes de tres;
  - hasta cuatro muestras diarias por fermentación;
  - ventana de laboratorio aproximadamente 10:00–16:00;
  - extracción de `5 mL` en muestra pequeña y cerca de `50 mL` cuando incluía etanol/aromas;
  - sólo una inyección automatizable;
  - registro de tiempos, dosis y condiciones iniciales reales.

- Tras un inicio real del Lote 1 a las 15:30, se decidió reanclar `t=0` a la inoculación real.
- Los pulsos críticos se movieron a las 12:00 para capturar respuesta postpulso. La pérdida FIM frente a las 16:00 era marginal.

### Reevaluación posterior al Lote 1

- Se mantuvo el Lote 2 con F4, F5 y F6.
- El benchmark D-opt puro sólo proponía reemplazar el experimento de stripping CO2/aromas por `synthetic_low_yan_ladder`; se prefirió conservar la campaña híbrida por robustez y valor científico.
- Se corrigió un error de mapeo que impedía que los estados core del Lote 1 entraran realmente al ajuste/FIM.

## 3. Resultados alcanzados

### Modelo base

- El notebook original quedó ejecutable hasta antes de PSO.
- Modelo de siete parámetros:

  - WSSE: `6001.9828`;
  - FIM de rango completo;
  - condición aproximada: `8.38e3`;
  - siete perfiles likelihood bilaterales al 95%;
  - posterior Laplace dominado por datos.

### Campaña y operación

- Se generaron bundles, protocolos, calendarios, tablas de muestreo e inputs.
- La reevaluación numérica condicionada al Lote 1 real confirmó que no convenía reemplazar F4–F6.
- Campaña actual F4–F9 frente a D-opt puro:

  - `logdet`: `178.949` frente a `178.987`;
  - mejor `min_eig`, `trace_inv`, peor reducción de varianza y reducción de `Kd0` para la campaña actual.

### Aporte de información del Lote 1

- Se incorporaron entre 121 y 141 observaciones core por fermentación.
- Residuos totales: `1150 → 1537`.
- Core FIM:

  - `logdet`: `76.70 → 82.70`;
  - rango: `10 → 11`;
  - `trace_inv`: `820.8 → 63.1`.

- FIM extendida: rango `10 → 17`, equivalente a siete nuevas direcciones independientes de información.
- Mejoraron especialmente `kAcStress`, `kAldS_N`, `kPyrS_N`, `kPyrDrain`, `kPyrO2` y `kAcAld`.
- `Kd0` y `qN` cerraron perfiles likelihood bilaterales en los dos conjuntos comparados.
- El WSSE por residual aumentó de aproximadamente `20.75` a `23.59`; la mejora demostrada fue de identificabilidad, no todavía de predicción fuera de muestra.

## 4. Pendientes

1. Completar 22 valores faltantes de etanol del Lote 1:

   - F1: 7;
   - F2: 8;
   - F3: 7.

2. Corregir la ruta prioritaria del script: el libro real está bajo `data/Laboratorio 2026/MBDoE_2026/DOE_Lote_1/`.
3. Reejecutar estimabilidad después de completar etanol.
4. Ejecutar validación predictiva fuera de muestra; la FIM no demuestra por sí sola mejor predicción.
5. Incorporar datos de aromas para identificar sus nueve parámetros todavía no soportados.
6. Actualizar el prior antes de cerrar el Lote 3.
7. Revisar las desviaciones reales de YAN, biomasa y dosis de pulsos antes de interpretar causalmente los ensayos.

## 5. Archivos relevantes

- `fermentation_model/fermentation_model_calibration_3_effective_reformulation.ipynb`
- `fermentation_model/fermentation_lot1_data_preview.executed.ipynb`
- `fermentation_model/fermentation_estimability_old_vs_lot1.ipynb`
- `fermentation_model/fermentation_estimability_old_vs_lot1.executed.ipynb`
- `fermentation_model/run_estimability_old_vs_lot1.py`
- `fermentation_model/results/estimability_old_vs_lot1/estimability_old_vs_lot1_report.md`
- `fermentation_model/results/lot1_campaign_reassessment/`
- `fermentation_model/results/lot1_actual_mbdoe_reassessment/`
- `fermentation_model/results/design_execution_bundle_2026-06-09/`
- `fermentation_model/data/Laboratorio 2026/MBDoE_2026/DOE_Lote_1/DOE_Lote_1_ethanol_manual_entry.xlsx`

## Estado de continuidad

El antecedente más reciente de este hilo es el análisis `historical_only` frente a `historical_plus_lot1`, ya corregido para incorporar los estados core del Lote 1. La siguiente acción concreta era completar el Excel de etanol, corregir la ruta del loader y reejecutar la evaluación antes de cerrar el rediseño del Lote 3.
