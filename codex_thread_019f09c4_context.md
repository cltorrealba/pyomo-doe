# Antecedente de conversación Codex

## Metadatos

- **Thread ID:** `019f09c4-4819-75e2-9d65-86ee984d4ee8`
- **Transcript:** `C:\Users\ctorrealba\.codex\sessions\2026\06\27\rollout-2026-06-27T11-48-06-019f09c4-4819-75e2-9d65-86ee984d4ee8.jsonl`
- **Fork de:** `019e9a11-c081-77c3-af6e-9198e048999a`
- **Repositorio:** `cltorrealba/pyomo-doe`
- **Nota:** el transcript fue leído en modo de solo lectura y no fue modificado. Este resumen reconoce el contexto heredado y enfatiza el trabajo propio del fork.

## 1. Objetivo del trabajo

Reutilizar el framework de calibración, identificabilidad y MBDoE desarrollado a escala laboratorio para construir un flujo homólogo con datos piloto 2025. El objetivo específico fue:

1. integrar proceso, CO2 y aromas de escala piloto;
2. seleccionar estructuras cinéticas capaces de explicar cada estado;
3. incorporar partición agua–etanol–azúcar y stripping con CO2;
4. corregir el inicio excesivamente temprano de liberación de CO2;
5. agregar un estado macroscópico de O2;
6. dejar un workflow piloto reproducible y ordenar el repositorio.

## 2. Decisiones tomadas

### Datos piloto

- Se trabajó con ocho fermentaciones piloto 2025.
- CO2 de `25150` y `25151` fue excluido.
- `25171` se reindexó finalmente desde `ING25-SB012-Pre reinóculo (25171)`, aproximadamente a 114 h del reloj original.
- `25170` se mantuvo, aunque su ajuste inicialmente era débil.
- `*_total` se interpretó como aroma retenido en vino más condensado; `*_condensado` como masa acumulada capturada.

### Modelo de aromas y CO2

- Para ethyl acetate se seleccionó `ea_ethanol_nlimited`.
- Se usó Antoine + UNIFAC con agua, etanol y `G+F` como glucosa equivalente.
- El modelo de CO2 pasó por las variantes `instant`, `threshold`, `lag_threshold` y `solubility_scaled`.
- La liberación simulada temprana de `25170` se atribuyó a que el modelo primario producía etanol/CO2 desde `t=0` sin estar regulado por O2.
- No se adoptó un modelo aeróbico completo por falta de observaciones y riesgo de confusión paramétrica.
- Se implementó un bloque efectivo O2→anaerobiosis con:

  - O2 inicial cercano a saturación para mosto fresco;
  - O2 inicial bajo para `25171`;
  - piso Crabtree;
  - CO2 respiratorio pequeño;
  - parámetros de O2 mayoritariamente fijados desde literatura.

- La estructura final seleccionada fue `solubility_o2_slow_transition`.

### Organización

- Primero se ordenó `pilot_2025`, separando flujo vigente, `support/` y `legacy/`.
- Después se realizó una limpieza conservadora de `fermentation_model`:

  - 31 notebooks/runners sustituidos fueron archivados;
  - no se eliminaron resultados aún consumidos;
  - se generaron README, mapa de dependencias y manifiesto de archivo.

- Se reconoció que esta limpieza no implementaba todavía la arquitectura final `src/ + campaigns/ + configs/ + run_id`.

## 3. Resultados alcanzados

### Modelo piloto integrado

- Contexto cargado: `ea_ethanol_nlimited`, `secondary_full_chem_o2fixed`, 51 parámetros.
- Antes del bloque O2, `solubility_scaled` produjo:

  - `kCO2_release_h = 0.232 h^-1`;
  - `CO2sat_scale = 0.491`;
  - correlación `25170 = 0.427`;
  - correlación `25171 = 0.926`.

- Después de incorporar O2:

  - modelo seleccionado: `solubility_o2_slow_transition`;
  - correlación `25170 = 0.843`;
  - correlación `25171 = 0.931`;
  - FIM focalizada CO2/O2 de rango `2/2`;
  - eigenvalor relativo mínimo `0.307`;
  - condición `3.26`.

- El modelo con `qO2` libre llegó al límite inferior; se conservó la variante slow-transition con ese parámetro fijado.
- Persistió un mal ajuste importante de `ethyl_acetate_condensate`.

### Organización del repositorio

- `fermentation_model/pilot_2025` quedó con un único workflow vigente.
- Se crearon:

  - `fermentation_model/README.md`;
  - `fermentation_model/REPOSITORY_MAP.md`;
  - `fermentation_model/results/README.md`;
  - `fermentation_model/legacy/development_2026/ARCHIVE_MANIFEST.md`.

- Validaciones:

  - 29 archivos Python activos con sintaxis válida;
  - 10 notebooks vigentes como JSON válido y sin outputs de error;
  - sin imports activos hacia módulos archivados.

- Commit creado: `816bf51` — `Organize fermentation model workflow artifacts`.
- No se realizó push.

## 4. Pendientes

1. Ejecutar la FIM global completa después de O2. La FIM focalizada no evalúa la confusión con todos los parámetros aromáticos.
2. Resolver el ajuste de ethyl acetate y otros parámetros débiles de secundarios/aromas.
3. Validar el bloque O2 con DO medido; es un modelo efectivo, no una descripción aeróbica completa.
4. Implementar la arquitectura final propuesta:

   - `src/fermentation/`;
   - `campaigns/laboratory_2025`, `pilot_2025`, `laboratory_2026`, `pilot_2026`;
   - `data/raw/` y `data/processed/`;
   - configuraciones versionadas;
   - `results/<campaign>/<run_id>/`;
   - manifiestos, checksums y tests de paridad.

5. Migrar primero `laboratory_2026` de extremo a extremo y luego las otras campañas.
6. El directorio vacío `fermentation_model/fermentation_model/` permaneció por restricciones de OneDrive.
7. Publicar el commit sólo cuando corresponda; el hilo terminó dejándolo para otro agente.

## 5. Archivos relevantes

- `fermentation_model/pilot_2025/pilot_2025_co2_solubility_integrated_doe.ipynb`
- `fermentation_model/pilot_2025/pilot_2025_co2_solubility_integrated_doe.executed.ipynb`
- `fermentation_model/pilot_2025/run_pilot_2025_co2_solubility_integrated_doe.py`
- `fermentation_model/pilot_2025/results/co2_solubility_integrated_doe/`
- `fermentation_model/pilot_2025/README.md`
- `fermentation_model/pilot_2025/support/`
- `fermentation_model/pilot_2025/legacy/`
- `fermentation_model/REPOSITORY_MAP.md`
- `fermentation_model/legacy/development_2026/ARCHIVE_MANIFEST.md`

## Estado de continuidad

El workflow científico vigente es `pilot_2025_co2_solubility_integrated_doe` con `solubility_o2_slow_transition`. El ordenamiento realizado fue una limpieza conservadora preliminar, no la migración arquitectónica definitiva. El commit local de organización quedó listo para un agente posterior.
