# Antecedente de conversación Codex

## Metadatos

- **Thread ID:** `019f53ba-6fa4-7ab1-9341-7aeb311a40cd`
- **Transcript consultado:** `C:\Users\ctorrealba\.codex\sessions\2026\07\11\rollout-2026-07-11T20-29-22-019f53ba-6fa4-7ab1-9341-7aeb311a40cd.jsonl`
- **Repositorio:** `https://github.com/cltorrealba/pyomo-doe.git`
- **Rama observada:** `ctorrealba_fermentation`
- **Directorio de trabajo:** `C:\Users\ctorrealba\OneDrive - Viña Concha y Toro S.A\Documentos\Doctorado\Artículos\Artículo - Estimación_dFBA\pyomo-doe`
- **Nota:** el transcript fue leído en modo de solo lectura y no fue modificado.

## 1. Objetivo del trabajo

El objetivo evolucionó en varias etapas:

1. Resolver incompatibilidades de IPOPT y ejecutar `fermentation_model_calibration_3_effective_reformulation.ipynb` hasta antes de PSO.
2. Obtener un modelo de fermentación con buen ajuste e identificabilidad estructural y práctica, usando FIM, eigenvalores, profile likelihood y una aproximación bayesiana de Laplace.
3. Diseñar experimentos MBDoE para mejorar parámetros débiles mediante perfiles de temperatura, nutrientes y, inicialmente, composición de mosto sintético.
4. Extender el modelo con glicerol, piruvato, acetaldehído, acetato, CO2, aromas y pérdidas al condensador.
5. Transferir el framework a datos de escala piloto 2025.
6. Corregir la liberación demasiado temprana de CO2 incorporando un estado macroscópico de oxígeno, transición hacia anaerobiosis, efecto Crabtree y CO2 respiratorio dentro del notebook piloto existente.

## 2. Decisiones tomadas

### Calibración e identificabilidad inicial

- Se eliminaron opciones de IPOPT incompatibles.
- La UQ basada en `k_aug/PyNumero` fue reemplazada por diferencias finitas.
- La base inicial aceptada fue un modelo identificable de siete parámetros:

  - Libres: `mu0`, `qN`, `betaG0`, `betaF0`, `qEG`, `qEF`, `iG`.
  - Fijos: `sN`, `qXG`, `qXF`, `sG`, `sF`, `iE`, `Kd0`, `m0`.

- La selección de parámetros y modelos combinó:

  - ajuste WSSE;
  - límites activos;
  - robustez numérica;
  - FIM y eigenvalores;
  - profile likelihood;
  - aproximación de Laplace;
  - AIC/BIC;
  - multistart, PSO y regularización log-L2.

- El PSO piloto no encontró una cuenca claramente mejor. No se justificó liberar parámetros únicamente para mejorar visualmente las curvas.

### Aromas y partición gas-líquido

- Se aplicó selección cinética por estado. Para ethyl acetate quedó seleccionada la estructura `ea_ethanol_nlimited`.
- El equilibrio gas-líquido se representó mediante Antoine + UNIFAC.
- Se incorporó una mezcla agua-etanol-azúcar, usando `G+F` como glucosa equivalente.
- Las columnas aromáticas se interpretaron de la siguiente forma:

  - `*_total = retenido en vino + condensado`;
  - `*_condensado` representa pérdida acumulada capturada.

### Datos piloto y CO2

- El MBDoE piloto se restringió a palancas realistas para mosto natural: temperatura y nutrientes.
- Curación de datos CO2:

  - `25150` y `25151`: excluidos.
  - `25170`: conservado desde su inicio operacional.
  - `25171`: redefinido finalmente con `t=0` en `ING25-SB012-Pre reinóculo (25171)`, alrededor de 114 h del reloj original.

- El modelo de liberación de CO2 evolucionó desde liberación instantánea hacia:

  1. `threshold`;
  2. solubilidad explícita;
  3. `solubility_o2_slow_transition`.

### Bloque macroscópico de oxígeno

La última estructura seleccionada:

- inicializa O2 cerca de saturación en mosto fresco;
- usa O2 inicial bajo en `25171`;
- modula la fracción fermentativa mediante una transición dependiente de O2 con piso Crabtree;
- incluye una contribución respiratoria pequeña de CO2;
- mantiene fijado el parámetro de transición/consumo de O2, porque al estimarlo alcanzaba el límite inferior.

### Organización del proyecto

`fermentation_model/pilot_2025` fue reorganizado de modo que:

- el flujo vigente permanece en la raíz de `pilot_2025`;
- los módulos auxiliares están en `support/`;
- notebooks, scripts y resultados antiguos están en `legacy/`.

## 3. Resultados alcanzados

### Modelo base

- El notebook original quedó ejecutable hasta antes de PSO.
- Para el modelo de siete parámetros:

  - WSSE: `6001.9828`;
  - FIM de rango completo;
  - número de condición: `8.38e3`;
  - profile likelihood satisfactorio para los siete parámetros;
  - posterior de Laplace dominado por los datos.

### Framework extendido

- Se generaron campañas MBDoE para mosto sintético y natural.
- Se incorporaron restricciones reales de personal, horarios, volumen de muestreo y pulsos.
- Se construyó el framework piloto con ocho fermentaciones:

  - `25026`;
  - `25027`;
  - `25085`;
  - `25086`;
  - `25150`;
  - `25151`;
  - `25170`;
  - `25171`.

- El contexto piloto cargaba:

  - modelo aromático `ea_ethanol_nlimited`;
  - modelo secundario `secondary_full_chem_o2fixed`;
  - 51 parámetros en el contexto extendido.

### Resultado final del bloque O2/CO2

- Modelo seleccionado: `solubility_o2_slow_transition`.
- Correlación CO2 para `25170`: `0.843`, frente a aproximadamente `0.43` con solubilidad simple.
- Correlación CO2 para `25171`: `0.931`.
- FIM focalizada CO2/O2:

  - rango numérico: `2/2`;
  - eigenvalor relativo mínimo: `0.307`;
  - número de condición: `3.26`.

Las direcciones `kCO2_release_h` y `CO2sat_scale` quedaron informadas por los datos online de CO2.

## 4. Pendientes

1. Recalcular la FIM global completa con todos los parámetros aromáticos después de incorporar O2. Cada evaluación de `aroma_residual` tomaba aproximadamente 63 segundos, por lo que en el último trabajo se utilizó una FIM focalizada CO2/O2.
2. Validar el bloque macroscópico de O2. No representa explícitamente toda la respiración, los rendimientos aerobios, la regulación redox, esteroles o lípidos.
3. Validar parámetros fijados desde literatura usando mediciones de DO o nuevos experimentos.
4. Resolver problemas de ajuste o estimabilidad todavía presentes, especialmente en:

   - `ethyl_acetate_condensate`;
   - `alpha_EA_loss`;
   - `k_EA_XE` y otros parámetros asociados a EA;
   - `kAldRed`;
   - `kAcAssim` y `kAcStress`;
   - `kPyrS_stat`;
   - algunas direcciones asociadas a `gammaF0`.

5. Recalibrar conjuntamente el modelo ampliado y comprobar que la mejora de CO2 no degrade el ajuste de etanol, azúcares, biomasa, secundarios o aromas.
6. Actualizar la campaña MBDoE con nuevos datos experimentales, ya que el diseño actual es local al punto nominal empleado.

## 5. Archivos y repositorios relevantes

### Flujo piloto vigente

- `fermentation_model/pilot_2025/pilot_2025_co2_solubility_integrated_doe.executed.ipynb`
- `fermentation_model/pilot_2025/pilot_2025_co2_solubility_integrated_doe.ipynb`
- `fermentation_model/pilot_2025/run_pilot_2025_co2_solubility_integrated_doe.py`
- `fermentation_model/pilot_2025/results/co2_solubility_integrated_doe/`
- `fermentation_model/pilot_2025/README.md`
- `fermentation_model/pilot_2025/support/`
- `fermentation_model/pilot_2025/legacy/`

### Antecedentes metodológicos

- `fermentation_model/fermentation_model_calibration_3_effective_reformulation.ipynb`
- `fermentation_model/run_fit_strategy_analysis.py`
- `fermentation_model/run_doe_experiment_design.py`
- `fermentation_model/fermentation_extended_campaign_doe.ipynb`
- `fermentation_model/fermentation_aroma_campaign_doe.ipynb`
- `fermentation_model/aroma_partition_unifac.py`

### Datos

- `fermentation_model/data/Piloto 2025/Calibration_data_vl3.xlsx`
- `fermentation_model/data/Piloto 2025/Sensores CO2/`

## Estado de continuidad

El antecedente vigente para continuar el trabajo es el flujo `pilot_2025_co2_solubility_integrated_doe`, utilizando `solubility_o2_slow_transition` como última estructura seleccionada. La mejora del ajuste de CO2 está demostrada localmente, pero la recalibración y evaluación global del modelo extendido, particularmente del bloque aromático, permanecen abiertas.
