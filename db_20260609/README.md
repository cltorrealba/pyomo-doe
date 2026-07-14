# Bundle de diseno experimental para ejecucion

Este bundle contiene la version vigente y operacionalizable de la campana DOE para fermentaciones en reactores de 2 L.

La version recomendada es la carpeta `rec/`, correspondiente al DOE volumen-restringido. La carpeta `ref/` queda solo como referencia historica del diseno anterior sin restriccion de volumen. No ejecutar `ref/` como protocolo experimental.

## Decision vigente

Usar `rec/volume_constrained_doe_report.md` como reporte maestro.

El protocolo anterior tenia demasiadas muestras liquidas si cada muestra de etanol requiere 40 mL y cada muestra de aromas requiere 10 mL. Por eso se rehizo el diseno separando dos tipos de muestra:

- `full_liquid_sample`: 50 mL. Incluye etanol + aromas y panel completo.
- `small_liquid_sample`: 8 mL. No incluye etanol/aromas; se usa para biomasa, G/F, YAN/PAN/NH4, glicerol, piruvato, acetaldehido y acetato.

El criterio de seleccion sigue siendo model-based DOE por FIM extendida, pero el modelo de observacion ahora considera que etanol y aromas solo se observan en muestras `full`.

## Campana recomendada

| Orden | Lote | Inicio | Experimento | Medio | Horizonte | Perfil T |
|---:|---:|---|---|---|---:|---|
| 1 | 1 | 2026-06-15 08:00 | synthetic_lit_SM410_18C_highN_ester | sintetico | 216 h | 18,18,18,18 |
| 2 | 1 | 2026-06-15 08:00 | synthetic_high_biomass_low_N_maintenance | sintetico | 216 h | 16,18,20,18 |
| 3 | 1 | 2026-06-15 08:00 | synthetic_fructose_rich_glucose_pulse | sintetico | 216 h | 17,20,23,20 |
| 4 | 2 | 2026-06-29 08:00 | synthetic_glucose_rich_fructose_pulse | sintetico | 216 h | 17,20,23,20 |
| 5 | 2 | 2026-06-29 08:00 | synthetic_lit_SM410_24C_highN_strip | sintetico | 168 h | 20,24,24,20 |
| 6 | 2 | 2026-06-29 08:00 | synthetic_viable_biomass_step | sintetico | 216 h | 18,22,22,18 |
| 7 | 3 | 2026-07-13 08:00 | synthetic_high_sugar_reference | sintetico | 216 h | 18,20,20,18 |
| 8 | 3 | 2026-07-13 08:00 | natural_glucose_pulse_after_growth | natural | 216 h | 18,20,22,19 |
| 9 | 3 | 2026-07-13 08:00 | synthetic_fast_CO2_aroma_strip_highN | sintetico | 168 h | 20,24,24,20 |

## Volumen

La politica seleccionada fue `full14_small6`.

- Fermentaciones de 216 h: 14 muestras full + 6 small = 748 mL retirados.
- Fermentaciones de 168 h: 12 muestras full + 4 small = 632 mL retirados.
- Volumen inicial asumido: 2000 mL.
- Volumen final estimado minimo: 1252 mL.
- Limite de seguridad usado en el diseno: 1100 mL.

Revisar `rec/volume_audit.csv` antes de ejecutar.

## Stocks

La formula usada para definir volumen de pulso fue:

```text
V_stock = dosis_objetivo * V_reactor_actual / C_stock
```

Concentraciones asumidas:

- N: 20 mg N/mL.
- G: 700 g/L.
- F: 700 g/L.
- E: 789 g/L, no usado en la campana recomendada.
- X: 100 g biomasa seca equivalente/L.

Revisar `rec/stock_concentration_assumptions.csv`. Si se cambia una concentracion de stock, recalcular volumen de pulso y auditoria de volumen.

## Archivos principales

- `rec/volume_constrained_doe_report.md`: reporte maestro.
- `rec/campaign_protocol.csv`: condiciones iniciales y resumen por fermentacion.
- `rec/operational_schedule.csv`: calendario completo de eventos.
- `rec/lot_1_schedule.csv`, `rec/lot_2_schedule.csv`, `rec/lot_3_schedule.csv`: calendarios separados por lote.
- `rec/daily_workload_summary.csv`: carga diaria por lote y tipo de evento.
- `rec/volume_audit.csv`: volumen final estimado por reactor.
- `rec/policy_summary.csv`: comparacion de politicas de muestreo.
- `rec/selected_campaign.csv`: seleccion DOE completa.
- `rec/post_campaign_estimability.csv`: estimabilidad esperada post-campana.
- `rec/plots/`: perfiles de volumen por reactor.
- `rec/fermentation_final_operational_doe_volume_constrained.executed.ipynb`: notebook visual ejecutado de la campana vigente.
- `rec/notebook_plots/`: figuras generadas por el notebook visual, con inputs, muestras, volumen y trayectorias simuladas por fermentacion.

## Material de respaldo

- `diag/secondary_v2_model_report.md`: diagnostico del modelo secundario reducido.
- `diag/secondary_fit_capacity_report.md`: capacidad de ajuste de metabolitos secundarios.
- `ref/final_operational_doe_v2_report.md`: diseno anterior sin restriccion de volumen, solo referencia.
- `ref/pyomo_core_reduced11_checks.csv`: chequeo Pyomo.DoE del bloque core en el diseno previo.
- `scripts/run_final_operational_doe_volume_constrained.py`: script que genero el diseno recomendado.
- `scripts/run_final_operational_doe_v2.py`: script del diseno anterior.

## Restricciones que no deben relajarse sin recalcular DOE

- Etanol y aromas se toman juntos.
- Cada muestra full consume 50 mL.
- Cada muestra small consume 8 mL.
- No programar eventos fuera de ventana laboral.
- Mantener volumen final sobre 1100 mL.
- No convertir small samples en full samples sin auditar volumen.
- Si se elimina un experimento natural, se pierde validacion de transferencia a mosto real.

## Uso sugerido en otro chat

Subir o referenciar el zip completo y partir diciendo:

```text
Usa README.md como contexto. El diseno vigente esta en rec/. Necesito convertir operational_schedule.csv y campaign_protocol.csv en checklist experimental ejecutable por lote, considerando volumen, stocks, muestras full/small, responsables, QC analitico y puntos de decision entre lotes.
```
