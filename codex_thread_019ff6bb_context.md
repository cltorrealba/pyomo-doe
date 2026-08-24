# Antecedente de conversación Codex

## Metadatos

- **Thread ID:** `019ff6bb-6719-7731-9c2e-6fe1494da911`
- **Transcript:** `C:\Users\ctorrealba\.codex\sessions\2026\08\12\rollout-2026-08-12T12-08-24-019ff6bb-6719-7731-9c2e-6fe1494da911.jsonl`
- **Repositorio:** `cltorrealba/pyomo-doe`
- **Ámbito:** calibración CO2 Laboratorio 2026, comparación natural/sintético.
- **Nota:** el transcript fue leído en modo de solo lectura y no fue modificado.

## 1. Objetivo del trabajo

Evaluar el desempeño del modelo `solubility_o2_slow_transition` sobre datos de CO2 de mosto natural y sintético de Laboratorio 2026, calibrar por matriz, realizar validación cruzada con holdouts y diagnosticar el retardo de inicio observado.

El trabajo se amplió para:

- filtrar artefactos de muestreo y sensibilidad en frío;
- alinear CO2, temperatura, química y pulsos;
- incorporar activación biológica gradual;
- estudiar identificabilidad de inicio/duración;
- corregir offsets por sensor;
- diseñar una campaña experimental simple con CO2 manual, O2 disuelto y pH.

## 2. Decisiones tomadas

### Datos y holdouts

- Sintético: DOE-F01–F05 para calibración; DOE-F06 como holdout.
- Natural: LAB004–LAB011 para calibración; LAB012 como holdout.
- Se excluyeron definitivamente `LAB001–LAB003`, `LAB009` y `DOE-F02`.
- El proceso se limitó al intervalo entre primera y última muestra química.
- CO2 se usó sólo en la intersección con la cobertura real del sensor.

### Filtrado y alineación

- Transientes de muestreo se reconstruyeron sólo con recuperación consistente.
- Señales bajas se trataron como censura izquierda, no como cero.
- Filtro final: ventana de muestreo ±1.25 h, mediana robusta de 3 h y Savitzky–Golay de 5 h.
- Temperatura completa se interpoló a timestamps CO2 y entró en cinética y solubilidad.
- Pulsos se aplicaron como saltos con `sample_before_action`.
- Pulsos LAB se obtuvieron mediante densidad 1040 y calendario; dosis `0.14 kg N/m3`.

### Estructura de inicio

- Se diagnosticó que el gate de O2 original retrasaba excesivamente la producción de CO2.
- Se incorporó un intervalo químico de activación y una transición `smoothstep` gradual.
- Parámetros efectivos:

  - `f_s`: fracción del intervalo químico donde comienza la activación;
  - `f_d`: duración relativa de la transición.

- Estos parámetros se interpretaron como activación efectiva, no como tiempo muerto puro del sensor ni como constantes fisiológicas universales.

### Corrección por sensor

- El cero se corrigió por corrida antes del filtrado:

  `y_corr = max(y - Q0.10(primeras 12 h), 0)`.

- Se preservó la asignación física por canal/sensor para los lotes sucesivos.
- Se reconoció que el procedimiento puede remover CO2 biológico temprano y debe reemplazarse por un cero preinoculación medido.

## 3. Resultados alcanzados

### Calibración cruzada inicial

- Todos los escenarios iniciales tuvieron `R2 < 0` y subpredicción sistemática.
- `CO2sat_scale` llegó al límite inferior `0.35` en ambas matrices.
- Holdouts iniciales:

  - DOE-F06/sintético: RMSE `0.670`;
  - DOE-F06/natural: `0.643`;
  - LAB012/natural: `0.417`;
  - LAB012/sintético: `0.489`.

### Filtrado y modelo de activación

- Se corrigieron 199 de 12,694 puntos de 10 min.
- Rugosidad mediana reducida 66%, conservando 99.5–101% de la integral.
- MAE de temperatura usada frente a medida: 0.04–0.26 °C.
- El nuevo modelo redujo:

  - MAE del inicio `33.0 → 14.1 h`;
  - mediana del desfase `+26.2 → +6.4 h`.

- DOE-F06 mejoró de error de inicio `+31.4 → +8.1 h`.
- LAB012 siguió mal, con inicio aproximadamente 32 h demasiado temprano.

### Identificabilidad

Antes de corrección de offsets:

- sintético: `f_s=0.010`, `f_d=0.350`, ambos en límite inferior;
- natural: `f_s=0.670`, `f_d=1.248`.

Después de corregir offsets:

- sintético:

  - inicio/duración permanecen en límites;
  - condición `2.36e11`;
  - leave-one-batch-out siempre regresa a límites.

- natural:

  - `f_s=0.735`, `f_d=1.064`;
  - condición aproximada `130`;
  - correlación local inicio/duración `-0.918`;
  - duración inestable al retirar lotes.

- `O2_qmax` quedó en límite inferior en ambas matrices.
- La transferencia continuó siendo asimétrica; sintético→natural presentó un error de inicio cercano a `-50.7 h`.

### Entregables y commit

- Notebook de 30 celdas, sin errores y con 14 gráficos.
- Commit: `579d6f7c72026550400cca2a7af11de873f11162` — `Add 2026 CO2 cross-matrix calibration analysis`.
- No se realizó push; la rama quedó `ahead 1` en ese cierre.

### Diseño experimental recomendado

- Tres reactores a aproximadamente 15 °C:

  - R1/R2: réplicas basales;
  - R3: O2 inicial 3–4 mg/L mayor, sin superar ~90% de saturación.

- Igualar mosto, volumen, inóculo, nutrientes, headspace, agitación y tuberías.
- Antes de inocular: 2–4 h de baseline para CO2, O2, temperatura y pH.
- Registro continuo ideal cada minuto.
- Química y CO2 manual intensivos alrededor de inoculación, agotamiento de O2 e inicio de salida gaseosa.

## 4. Pendientes

1. Medir el cero de cada sensor antes de inocular; no inferirlo desde horas con actividad biológica posible.
2. Medir O2 disuelto en continuo y CO2 disuelto/manual para separar producción, acumulación y liberación.
3. Obtener química temprana más frecuente.
4. Calibrar la relación pH–CO2 en el propio mosto; pH solo no demuestra saturación.
5. Ejecutar perfiles likelihood, Jacobiano/SVD, bootstrap por fermentación y leave-one-batch-out formal.
6. No interpretar `f_s`, `f_d` u `O2_qmax` como parámetros fisiológicos bien identificados.
7. Validar transferibilidad entre matrices con nueva evidencia; actualmente no está demostrada.

## 5. Archivos relevantes

- `fermentation_model/laboratory_2026/notebooks/co2_solubility_o2_cross_matrix_2026.ipynb`
- `fermentation_model/laboratory_2026/notebooks/co2_solubility_o2_cross_matrix_2026.executed.ipynb`
- `fermentation_model/laboratory_2026/run_co2_matrix_cross_validation_2026.py`
- `fermentation_model/laboratory_2026/results/co2_matrix_cross_validation_2026/`
- `fermentation_model/laboratory_2026/CO2_CROSS_MATRIX_HANDOFF.md`
- `fermentation_model/pilot_2025/run_pilot_2025_co2_solubility_integrated_doe.py`

## Estado de continuidad

El modelo representa mejor el onset, pero inicio, duración, O2 y transferencia continúan parcialmente confundidos. El siguiente paso defendible es la campaña de tres reactores con baseline preinoculación, O2 continuo y CO2 manual/disuelto, no añadir más flexibilidad al modelo actual.
