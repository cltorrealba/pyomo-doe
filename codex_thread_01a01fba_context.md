# Antecedente de conversación Codex

## Metadatos

- **Thread ID:** `01a01fba-8957-7483-ada2-9fe42e2aeada`
- **Transcript:** `C:\Users\ctorrealba\.codex\sessions\2026\08\20\rollout-2026-08-20T11-11-53-01a01fba-8957-7483-ada2-9fe42e2aeada.jsonl`
- **Fork de:** `019ff6bb-6719-7731-9c2e-6fe1494da911`
- **Repositorio principal:** `cltorrealba/pyomo-doe`
- **Ámbito:** reproducción piloto 2026, CO2, aromas, transferencia de resultados y handoff a clúster UC.
- **Nota:** el transcript fue leído en modo de solo lectura y no fue modificado.

## 1. Objetivo del trabajo

Extender el análisis CO2 de laboratorio a escala piloto 2026 y, sobre esa base, recalibrar el modelo de producción y pérdida aromática. El objetivo incluía:

1. validación cruzada de CO2 entre lotes y escalas;
2. conexión del modelo aromático con el nuevo `rCO2` liberado;
3. auditoría de unidades, diluciones y masas de condensado;
4. evaluación de modelos de transferencia, partición y captura;
5. definición de un gate de validez científica;
6. preparación de una campaña prospectiva;
7. commit/push, transferencia a un artículo y continuidad en el clúster UC mediante Slurm.

## 2. Decisiones tomadas

### CO2 piloto

- Se usaron nueve fermentaciones piloto 2026; seis perfiles CO2 eran modelables.
- Lote 1 se excluyó de calibración CO2 por mala señal.
- Temperatura medida, offsets/ganancias por sensor, pulsos y densidad 1040 se incorporaron explícitamente.
- Se hizo validación interna por lote, entre lotes y entre escalas.
- Se mantuvo como caveat central la mala identificabilidad de parámetros aunque las curvas fueran razonables.

### Modelo de aromas

- Se trabajó inicialmente con acetato de etilo, acetato de isoamilo y octanoato de etilo.
- Holdouts: `26158` y `26211`.
- El nuevo `rCO2` liberado reemplazó al forcing estequiométrico, pero la escala de pérdida compensó el cambio y no mejoró la validación.
- Se declaró que el baseline de tres parámetros era fallido:

  - NRMSE holdout de vino 49–66%;
  - subestimación sistemática;
  - condensados predichos en sólo 1–3.5% de lo observado.

- Se auditó y confirmó el factor de dilución 1:1000 y el cálculo de masa de condensado.
- `alpha` se reinterpretó como escala efectiva de stripping, no eficiencia física del condensador.
- Se corrigió la formulación de transferencia y se evaluaron:

  - NTU histórico;
  - `kLa` líquido corregido;
  - equilibrio UNIFAC;
  - equilibrio Morakul;
  - headspace/reservorio de línea;
  - captura dependiente de etanol;
  - captura completa `eta_A+B=1`.

- Ninguna estructura cumplió el gate de validación; no se autorizó declarar un modelo válido.

### Interpretación física

- La fracción volatilizada estimada resultó consistente con Mouret/Sablayrolles a 18–20 °C:

  - isoamilo ~16.6%;
  - octanoato ~30.1%.

- El conflicto principal quedó en la relación entre masa emitida y masa recuperada en condensado.
- Sin medir aroma después del segundo condensador no puede estimarse la eficiencia total de captura.
- Los factores `eta` calibrados se renombraron conceptualmente como recuperación aparente/efectiva, no eficiencia física.
- Se impuso después `eta_A+B=1` como supuesto de ingeniería y se reservó el condensado como validación independiente.

### Gobernanza y campaña prospectiva

- Se definió un contrato prospectivo con parámetros congelados y sin reajuste.
- El goal se marcó `blocked` al no existir una campaña prospectiva independiente.
- Se decidió no añadir parámetros adicionales a los datos 2026 por riesgo de sobreajuste.

## 3. Resultados alcanzados

### CO2 piloto 2026

- Holdouts internos:

  - `26158`: RMSE `0.098`, onset `+4.3 h`;
  - `26211`: RMSE `0.134`, onset `+12.4 h`.

- Transferencia:

  - lote 3→26158: RMSE `0.179`, onset `-12.2 h`;
  - lote 2→26211: RMSE `0.263`, onset `+29.6 h`.

- Lote 3 presentó condición aproximada `3.9e11` y seis parámetros en límites.
- Notebook ejecutado sin errores y con 12 gráficos.
- Commit `2fb56b7` — `Add pilot-scale CO2 cross-lot validation` fue enviado a `origin/ctorrealba_fermentation`.

### Recalibración aromática con nuevo rCO2

- rCO2 liberado integrado fue aproximadamente 1.82 veces el estequiométrico.
- La escala de pérdida disminuyó y compensó casi completamente el cambio.
- Holdout:

  - acetato de etilo: RMSE `5.547 µg/L`, NRMSE `48.8%`;
  - octanoato de etilo: `48.4 µg/L`, NRMSE `60.7%`;
  - acetato de isoamilo: `723 µg/L`, NRMSE `65.8%`.

### Mejor modelo diagnóstico

- Captura dependiente de etanol + reservorio de línea:

  - octanoato condensado NRMSE `1.188`, mejora `23.4%`;
  - isoamilo condensado NRMSE `0.985`, mejora `25.5%`;
  - vino prácticamente sin deterioro;
  - gate requerido: mejora ≥30% en ambos;
  - resultado: `NO_VALID_MODEL`.

- Parámetros aparentes:

  - octanoato: `eta50=4.00%`, multiplicador etanol `2.34x`, `tau=4.63 h`;
  - isoamilo: `eta50=1.59%`, multiplicador `2.24x`, `tau=5.06 h`.

- Balance mediano:

  - octanoato: 27.38 mg producidos, 8.57 mg emitidos, 1.28 mg recuperados;
  - isoamilo: 328.91 mg producidos, 55.27 mg emitidos, 3.58 mg recuperados.

- Estos valores no demuestran baja eficiencia física del condensador; mezclan pérdida, línea, captura, muestreo y analítica.

### Supuesto de captura completa

- Con `eta_A+B=1`, el balance cerró, pero el modelo siguió en `NO_VALID_MODEL`:

  - octanoato: vino NRMSE `0.723`, condensado `1.230`;
  - isoamilo: vino `1.009`, condensado `0.922`;
  - reservorio llevado al límite `tau=72 h`;
  - pérdida severa de identificabilidad.

- El handoff final describe como baseline operativo un modelo calibrado sólo contra vino, captura completa como supuesto y condensado como prueba externa:

  - NRMSE LORO vino: octanoato `0.419`, isoamilo `0.475`;
  - cierre de condensado sobrepredicho: `4.288x` y `6.570x`.

### Handoff y clúster

- Clúster:

  - host `cluster-uc`;
  - repo `/home/cltorrealba/pyomo-doe`;
  - rama `ctorrealba_fermentation`;
  - Python `/home/cltorrealba/pyomo-doe/.venv/bin/python`;
  - 23 tests Linux aprobados.

- Handoff completo:

  - `/home/cltorrealba/pyomo-doe/fermentation_model/pilot_2026/HANDOFF_AROMA_MODEL_2026.md`;
  - commit remoto `a126672`;
  - árbol limpio.

## 4. Pendientes

### Científicos

1. Resolver la sobrepredicción del cierre de condensado `4.288x/6.570x` sin recalibrar la producción biológica contra condensado.
2. Ejecutar modelos anidados de transferencia mediante Slurm.
3. Obtener campaña prospectiva sellada con:

   - `process.csv`;
   - `events.csv`;
   - `liquid_aroma.csv`;
   - `outlet_gas.csv`;
   - `condensate.csv`;
   - `trap_standards.csv`.

4. Medir simultáneamente vino, gas de salida y condensado.
5. Medir flujo de gas, volumen y %EtOH de condensado y calibraciones del tren de captura.
6. Mantener parámetros congelados para una validación prospectiva real.
7. No declarar un modelo válido mientras falle el gate independiente.

### Operacionales

- Ejecutar los cálculos pesados con Slurm, no en el nodo de acceso.
- Continuar desde `HANDOFF_AROMA_MODEL_2026.md` en la sesión remota.
- Para el repositorio del artículo se indicó la ruta `/home/cltorrealba/DC_dFVB_2026_fermentation_si_2026`, rama `paper/fermentation-si-2026`, con datos restringidos fuera de Git.

## 5. Archivos relevantes

- `fermentation_model/pilot_2026/notebooks/pilot_2026_co2_solubility_cross_lot_validation.executed.ipynb`
- `fermentation_model/pilot_2026/run_co2_solubility_cross_lot_validation_2026.py`
- `fermentation_model/pilot_2026/notebooks/pilot_2026_aroma_co2_release_recalibration.executed.ipynb`
- `fermentation_model/pilot_2026/notebooks/pilot_2026_aroma_final_nested_validation.executed.ipynb`
- `fermentation_model/pilot_2026/notebooks/pilot_2026_aroma_assumed_complete_capture_validation.executed.ipynb`
- `fermentation_model/pilot_2026/notebooks/pilot_2026_aroma_wine_calibrated_complete_capture.executed.ipynb`
- `fermentation_model/pilot_2026/results/aroma_final_nested_validation_2026/`
- `fermentation_model/pilot_2026/results/aroma_wine_calibrated_complete_capture_2026/`
- `fermentation_model/pilot_2026/adaptive_design/aroma_prospective_validation_contract.json`
- `fermentation_model/pilot_2026/protocols/aroma_prospective_validation_2026.md`
- `fermentation_model/pilot_2026/HANDOFF_AROMA_MODEL_2026.md`
- `publication_transfer/DC_dFVB_2026_fermentation_si_2026/`

## Estado de continuidad

La continuidad oficial está en el clúster UC, commit remoto `a126672`, mediante `HANDOFF_AROMA_MODEL_2026.md`. El modelo de producción se calibra sólo contra vino; `eta_A+B=1` es un supuesto, no una eficiencia medida; el condensado queda como validación externa. El gate científico permanece `NO_VALID_MODEL` hasta resolver el cierre o disponer de una campaña prospectiva.
