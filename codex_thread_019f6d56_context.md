# Antecedente de conversación Codex

## Metadatos

- **Thread ID:** `019f6d56-7a1f-7c52-8d68-62a32d8bd07e`
- **Transcript:** `C:\Users\ctorrealba\.codex\sessions\2026\07\16\rollout-2026-07-16T19-50-16-019f6d56-7a1f-7c52-8d68-62a32d8bd07e.jsonl`
- **Repositorio:** `cltorrealba/pyomo-doe`
- **Ámbito:** Piloto 2026, diseño adaptativo Wave 1.
- **Nota:** el transcript fue leído en modo de solo lectura y no fue modificado.

## 1. Objetivo del trabajo

Construir, recalificar y operacionalizar un workflow híbrido MBDoE para Piloto 2026 con gobernanza fail-closed. El proceso debía integrar datos, validar el adaptador/FIM, calibrar el modelo, ejecutar búsqueda Sobol–PSO, refinar con IPOPT, optimizar muestreo y producir un plan físicamente ejecutable sólo si todos los gates eran aprobados.

El objetivo evolucionó hacia:

- cobertura operacional entre perfiles frío, ancla y cálido;
- selección robusta de nutrición N76/N80;
- ejecución en tres tanques;
- reforecast adaptativo después de descubrir que el inóculo real era sólo 5–10% del supuesto;
- rediseño condicionado a las primeras 24 h ya ejecutadas.

## 2. Decisiones tomadas

### Gobernanza y gates

- Se separó estrictamente mejora computacional de autorización física.
- Si la calibración o la calificación operacional fallaban, no se emitía schedule ejecutable.
- Datos raw, manifests y auditorías no podían rebaselinarse silenciosamente.
- La muestra coincidente con una acción usa `sample_before_action`.

### Diseño numérico

- FIM con escalas nominales congeladas y validación de paso FD/malla.
- Pulsos representados como saltos causales.
- Ensemble robusto de 64 miembros.
- Exploración Sobol multisemilla + PSO con checkpoints y reanudación.
- Refinamiento local con IPOPT y aceptación contra el objetivo real full-ensemble.
- Comparación conjunta de candidatos originales, continuados y refinados.
- Selección lexicográfica con seguridad operacional dentro de 0.5% del mejor score.

### Cobertura operacional

- Se compararon cuatro estrategias de cobertura.
- Se seleccionó la Estrategia II:

  - ancla;
  - frío/temprano;
  - cálido/tardío.

- Se adoptó N76 en los tres tanques por margen de inventario y pérdida de información menor a 0.5% frente a N80.
- La primera campaña oficial quedó como `READY FOR OWNER PHYSICAL RELEASE REVIEW`, pero no como autorización automática.

### Cambio de inóculo y reforecast

- El diseño original asumía `X0 = 0.02244 g/L` y no incluía incertidumbre en `X0`.
- Reducir el inóculo a 5–10% invalidó FIM, predicciones aromáticas y calendario del run previo para la ejecución real.
- Se congelaron las primeras 24 h ya ejecutadas y se recalculó con:

  - YAN `210 mg/L`;
  - `1.7 millones células/mL`;
  - `17 °Brix = 156.29 g/L` de azúcares;
  - inicio real `21-jul-2026 16:00`.

## 3. Resultados alcanzados

### Primera integración y calibración

- Dataset validado:

  - 1,728 observaciones primarias;
  - 11,575 inputs de temperatura;
  - 6,119 observaciones CO2;
  - 45 eventos;
  - 234 observaciones de aromas en vino;
  - 198 observaciones de condensado;
  - ESS CO2 `206.57` en la corrida inicial.

- La primera calibración falló de forma controlada por auditoría/gobernanza:

  - 53 errores y 1 warning;
  - revisión QC pendiente;
  - definiciones MassView y YAN pendientes;
  - conversión Oculyze y escalas de error pendientes.

### Recalificación y cobertura

- Se corrigieron adaptador/FIM, causalidad de pulsos, representación de dosis, rutas largas, checkpoints y selección de candidatos.
- La corrida oficial de cobertura final quedó en:

  - run `wave1_operational_coverage_final/20260719T180451Z_9b8093`;
  - veredicto `PASS — READY FOR OWNER PHYSICAL RELEASE REVIEW`;
  - 64 miembros;
  - convergencia final ≤0.1%;
  - margen mínimo de secado `32.7205 h`;
  - temperatura robusta `15–27 °C`;
  - perfiles ancla, frío y cálido diferenciados.

- Rama `wave1_operational_coverage`:

  - remoto base `9ecc162`;
  - HEAD local informado `3d35310`;
  - tres commits por delante del remoto en ese cierre.

### Impacto del bajo inóculo

- Contraprueba con perfiles previos:

  - retraso mediano de secado entre aproximadamente 11 y 28 h;
  - un miembro del ensemble frío no secó antes de 504 h;
  - el diseño no podía considerarse vigente sin recalificación.

### Reforecast adaptativo final

Diseño recomendado desde `t=24 h`:

| Perfil | Setpoint | Nutrición |
|---|---:|---|
| Ancla | 18.0 °C | 23-jul 09:00 |
| Frío | 16.4 °C | ya ejecutada 21-jul 18:00 |
| Cálido | 25.6 °C | 23-jul 09:00 |

Cada N76 para 230 L:

- `87.4 g` SpringFerm Xtrem;
- `43.7 g` DAP;
- `76 mg YAN/L`.

Resultado:

- candidato `ef9288e679a442b9`;
- 63/64 miembros completan, probabilidad `98.44%`;
- margen mínimo antes de secado `88.7 h`;
- gates térmico, secado, muestreo y FIM: `PASS`;
- 79/79 pruebas específicas aprobadas;
- estado `PASS_CONDITIONAL — NOT READY FOR EXECUTION`.

## 4. Pendientes

1. Confirmar que `1.7 millones células/mL` corresponde a células viables, no totales.
2. Confirmar temperaturas reales de las primeras 24 h.
3. Confirmar qué muestras, además de la basal, ya fueron tomadas.
4. Registrar por tanque masa/volumen de inóculo, viabilidad, lote y hora real.
5. Obtener release físico explícito del owner antes de ejecutar el plan recalculado.
6. Verificar estado Git/publicación: varias etapas quedaron sólo locales o por delante del remoto.
7. Resolver auditorías históricas raw/metadata sin rebaselinarlas de forma implícita.

## 5. Archivos relevantes

- `fermentation_model/pilot_2026/adaptive_design/README.md`
- `fermentation_model/pilot_2026/adaptive_design/pilot_mbdoe_adapter.py`
- `fermentation_model/pilot_2026/adaptive_design/hybrid_optimizer.py`
- `fermentation_model/pilot_2026/adaptive_design/local_refinement.py`
- `fermentation_model/pilot_2026/adaptive_design/run_wave1_hybrid_search.py`
- `fermentation_model/pilot_2026/adaptive_design/optimize_wave1_sampling_and_plots.py`
- `fermentation_model/pilot_2026/adaptive_design/design_constraints.json`
- `fermentation_model/pilot_2026/adaptive_design/wave1_mbdoe_config.json`
- `fermentation_model/pilot_2026/adaptive_design/campaign_state.json`
- `fermentation_model/pilot_2026/adaptive_design/coverage_design_config.json`
- `fermentation_model/pilot_2026/results/adaptive_design_2026/wave1_operational_coverage_final/20260719T180451Z_9b8093/`
- `fermentation_model/pilot_2026/results/adaptive_design_2026/wave1_adaptive_reforecast/20260722T202428Z_2eecec/`

## Estado de continuidad

El estado final es un reforecast adaptativo condicionado a las primeras 24 h y a las nuevas condiciones iniciales. El diseño supera los gates computacionales, pero permanece condicionado a confirmar viabilidad del inóculo, temperaturas ejecutadas y muestras ya tomadas antes de cualquier liberación física.
