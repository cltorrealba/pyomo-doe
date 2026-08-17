# Entrega: calibración y transferencia cruzada de CO2, campaña 2026

Actualizado: 2026-08-17

## Ruta de revisión recomendada

Para revisar el análisis sin instalar ni ejecutar nada, abrir primero:

`notebooks/co2_solubility_o2_cross_matrix_2026.executed.ipynb`

El notebook ejecutado contiene las tablas de resultados y 14 figuras embebidas. El
notebook sin ejecutar, el código y todos los resultados tabulares están en el mismo
commit.

El modelo principal de esta iteración es
`solubility_o2_nitrogen_boost_continuous_release`. Se incluyen como comparadores la
liberación con umbral y el modelo histórico `solubility_o2_slow_transition`.

## Figuras prioritarias

Las figuras se encuentran en
`results/co2_matrix_cross_validation_2026/figures/`.

- `sensor_zero_offset_correction.png`: corrección del cero independiente por corrida.
- `sensor_artifact_filter_examples.png`: efecto del filtro de toma de muestras y del
  suavizado.
- `temperature_profiles_model_input.png`: temperatura medida/reconstruida que entra al
  modelo.
- `process_timeline_alignment.png`: alineación entre química, CO2, temperatura y pulsos.
- `initial_gradual_release_comparison.png`: forma de la liberación inicial de CO2.
- `nutrient_pulse_response_comparison.png`: respuesta alrededor del pulso nutricional.
- `heldout_cross_matrix_validation.png`: transferencia mosto sintético a natural y
  natural a sintético.
- `heldout_model_comparison.png`: comparación de modelos en los experimentos de
  validación.
- `activation_identifiability_profiles.png`: perfiles de objetivo para los parámetros
  de inicio y duración.

## Tablas prioritarias

Todas están en `results/co2_matrix_cross_validation_2026/`.

- `fit_parameters.csv`: parámetros calibrados por matriz.
- `validation_summary.csv`: validación interna y transferencia cruzada en LAB012 y
  DOE-F06.
- `model_comparison_validation.csv`: modelo continuo frente al modelo con umbral.
- `sensor_zero_offsets.csv`: offset estimado por sensor y corrida, antes de recorte,
  conversión y filtrado.
- `activation_jacobian_identifiability.csv`: rango y condición local del ajuste.
- `activation_objective_profiles.csv`: perfiles unidimensionales con reoptimización de
  los demás parámetros.
- `activation_leave_one_batch_out_summary.csv`: estabilidad al omitir una fermentación.
- `effective_nutrient_pulses.csv`: pulsos realmente aplicados en el modelo.
- `temperature_alignment_summary.csv`: cobertura y alineación temporal de temperatura.
- `excluded_experiments.csv`: exclusiones solicitadas y su motivo.
- `analysis_manifest.json`: fuentes, supuestos, umbrales y reglas de preprocesamiento.

## Lectura mínima de los resultados

- Se excluyeron LAB001-LAB003, LAB009 y DOE-F02.
- El tiempo de proceso comienza en la primera muestra química y termina en la última.
  CO2 fuera de ese intervalo no participa en el ajuste.
- La temperatura medida/reconstruida entra a la dinámica y a la solubilidad de CO2; el
  setpoint se utiliza como contexto y control de calidad.
- Los pulsos LAB se ubican en el cruce químico de densidad de 1040 g/L cuando existe;
  el calendario queda como auditoría o respaldo. Los pulsos DOE conservan el tiempo de
  proceso registrado.
- En la calibración sintética, los parámetros de activación quedan en sus cotas
  inferiores y el Jacobiano es muy mal condicionado. No deben interpretarse como una
  estimación fisiológica única.
- En la calibración natural, el inicio es más estable que la duración; inicio y duración
  siguen fuertemente correlacionados. La nueva medición de pH temprano está pensada para
  separar activación biológica de latencia/disolución del CO2.

## Reproducción

Desde la raíz del repositorio:

```powershell
python fermentation_model/laboratory_2026/run_co2_matrix_cross_validation_2026.py --timeout 1800
```

Dependencias directas: `numpy`, `pandas`, `scipy`, `matplotlib`, `nbformat` y
`nbclient`. La ejecución vuelve a generar el notebook y los resultados en las rutas
anteriores.

### Insumos externos no incluidos en este commit

El paquete versionado es suficiente para revisar resultados, figuras y trazabilidad.
Para repetir el análisis desde los datos originales también deben entregarse por el
canal de datos aprobado:

1. Los archivos de CO2/temperatura de mosto natural en
   `fermentation_model/data/Laboratorio 2026/raw_data/` (aproximadamente 202 MiB en la
   copia de trabajo actual).
2. El calendario `Fernanda Folch.ics`, actualmente resuelto por el runner desde
   `C:\Users\ctorrealba\Downloads\Fernanda Folch.ics`.

Los datos crudos anteriores no se agregaron automáticamente al commit para evitar
publicar archivos experimentales voluminosos o sensibles sin una decisión explícita.
Los insumos sintéticos procesados y las estimaciones upstream que consume el análisis
ya están versionados en las carpetas `results/lot1_data_preview/`,
`results/lot2_data_preview/`, `results/estimability_historical_natural/` y
`results/estimability_historical_synthetic_plus_lot2/`.

## Asignación de sensores usada

La asignación se realiza por canal físico de adquisición: F1 usa sensor 1, F2 usa
sensor 2 y F3 usa sensor 3. En el lote siguiente, F4 y F5 vuelven a los sensores 1 y
2; F6 conserva la etiqueta de usuario sensor 6, pero el archivo se adquiere por el
canal F3. El offset se estima de forma independiente en cada corrida, aunque se
mantenga el mismo sensor físico.
