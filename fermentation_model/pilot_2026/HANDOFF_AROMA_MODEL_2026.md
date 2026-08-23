# Traspaso técnico — modelo de producción y pérdida aromática, piloto 2026

## 1. Propósito de este documento

Este archivo es la fuente de contexto obligatoria para continuar el trabajo en
el cluster UC sin reinterpretar decisiones científicas ya discutidas. Antes de
editar o ejecutar el modelo, el agente debe leer este documento y verificar los
artefactos citados.

El objetivo activo es desarrollar un modelo defendible de producción y pérdida
aromática para las fermentaciones piloto 2026, usando el modelo de CO2 validado
como fuente de `rCO2`, separando producción biológica de transferencia/captura,
y evaluando generalización mediante leave-one-run-out (LORO). No existe todavía
un modelo completo aceptado.

## 2. Estado reproducible

- Repositorio en cluster: `/home/cltorrealba/pyomo-doe`
- Rama: `ctorrealba_fermentation`
- Commit base del traspaso: `bbadc16`
- Entorno Python: `/home/cltorrealba/pyomo-doe/.venv/bin/python`
- Prueba en Linux: 23 tests aprobados.
- El checkout del cluster estaba limpio al realizar este traspaso.
- Las calibraciones pesadas deben enviarse por Slurm; no se deben ejecutar en el
  nodo de acceso.

Comprobación mínima:

```bash
cd /home/cltorrealba/pyomo-doe
source .venv/bin/activate
git branch --show-current
git rev-parse --short HEAD
git status --short
python -m pytest \
  fermentation_model/tests/test_pilot_2026_aroma_calibration.py \
  fermentation_model/tests/test_pilot_2026_aroma_nested_models.py \
  fermentation_model/tests/test_pilot_2026_aroma_prospective_validation.py -q
```

## 3. Pregunta científica y separación de dominios

Se deben distinguir tres bloques:

1. Producción biológica del aroma en el líquido.
2. Transferencia gas-líquido y volatilización conducida por CO2.
3. Recuperación/observación del aroma volatilizado en los condensadores A+B.

La analítica química manda. En la estrategia actual, los parámetros de
producción se estiman **sólo con química de vino**. El condensado combinado A+B
no participa en la función objetivo y se usa como prueba independiente de
cierre. Esta separación evita que una recuperación desconocida altere los
rendimientos biológicos.

## 4. Supuesto de captura solicitado por el usuario

Se impone

```text
eta_AB = 1
```

es decir, toda la masa que el modelo predice que abandona el líquido se asigna
al condensado combinado A+B. Éste es un **supuesto de diseño de ingeniería**, no
una eficiencia medida ni identificada.

No existe medición de aroma en la corriente gaseosa que sale después del segundo
condensador. Por lo tanto, con los datos actuales no es posible determinar de
manera independiente cuánto aroma:

- permanece en A+B;
- queda en la línea;
- sale al ambiente;
- o fue sobrepredicho por el submodelo de transferencia.

Regla no negociable: no presentar un parámetro de captura ajustado a partir del
MIX A+B como eficiencia física del sistema. Sin medición de salida o balance
independiente, ese parámetro sólo sería un factor efectivo de observación.

## 5. Datos piloto usados

Reactores/corridas incluidos en el análisis aromático actual:

```text
26157, 26158, 26159, 26210, 26211, 26212
```

Los datos armonizados están versionados bajo:

```text
fermentation_model/pilot_2026/results/adaptive_design_2026/model_dataset/
20260717T125004Z_b03413/
```

La calibración piloto upstream utilizada está en:

```text
fermentation_model/pilot_2026/results/adaptive_design_2026/
baseline_calibration/20260717T153612Z_b5b67b/
```

Entradas relevantes:

- química en vino con tiempo de proceso;
- temperatura medida y setpoint;
- `rCO2` reconstruido por el modelo de CO2 piloto validado;
- consumo de azúcares y estados biológicos del modelo upstream;
- tiempos de pulsos nutricionales;
- muestras de condensado MIX A+B por intervalo.

La señal de proceso se interpola a la malla del modelo. Las observaciones
químicas mantienen sus timestamps y operadores de observación. Los pulsos se
usan con sus tiempos registrados y se muestran en los gráficos. No desplazar
tiempos químicos para mejorar visualmente el ajuste.

## 6. Unidades y reconstrucción del condensado

- Aroma en vino: `ug/L`.
- Concentración analítica del MIX A+B: `ug/L`.
- Volúmenes de condensadores: `volume_a_ml` y `volume_b_ml`.
- Volumen total: `total_condensate_ml = volume_a_ml + volume_b_ml`.
- Masa observada por intervalo:

```text
captured_mass_ug = mix_concentration_ug_l * total_condensate_ml / 1000
```

Para valores bajo LOQ se usa el operador censurado y la misma conversión con el
LOQ corregido. La concentración del MIX es una sola medición del combinado; no
hay concentración específica de cada condensador.

La integración está implementada en:

```text
fermentation_model/pilot_2026/run_data_integration.py
```

## 7. Estructura del modelo actual

Para cada compuesto se calculan rendimientos de formación asociados a la fase
de crecimiento y a la fase estacionaria. De forma esquemática:

```text
P(t) = [Y_growth * phi(t) + Y_stationary * (1 - phi(t))] * q_sugar(t)
```

El candidato con pulso agrega una respuesta transitoria posterior al pulso
nutricional. La dinámica en el líquido tiene la forma:

```text
dC_liq/dt = P(t) - k_loss(t) * C_liq(t)
```

La volatilización utiliza:

- partición dependiente de temperatura y etanol basada en Morakul et al.;
- caudal gaseoso obtenido de `rCO2` y densidad de CO2;
- temperatura y etanol variables en el tiempo.

En el candidato de equilibrio actualmente mostrado:

```text
q_gas(t) = rCO2(t) / rho_CO2(T)
k_loss(t) = K_gas/liquid(T, ethanol) * q_gas(t)
dM_captured/dt = V_liquid * k_loss(t) * C_liq(t)
```

Como `eta_AB=1`, la masa volatilizada y la masa capturada predicha son iguales.
El balance numérico cierra, pero eso no demuestra que el modelo físico describa
correctamente el tren real.

Código principal:

```text
fermentation_model/pilot_2026/adaptive_design/pilot_aroma_calibration.py
fermentation_model/pilot_2026/adaptive_design/pilot_aroma_nested_models.py
```

## 8. Iteraciones realizadas y lecciones

### 8.1 Escala empírica de pérdida

Un parámetro `alpha` absorbía simultáneamente errores de partición,
transferencia, producción y recuperación. No era una propiedad física
identificable.

### 8.2 Dependencia de etanol y transferencia dinámica

Se incorporó explícitamente el cambio temporal de volatilidad por concentración
de etanol, temperatura, `rCO2` y gas turnover. Esto corrigió una omisión
conceptual importante, pero no resolvió por sí solo el condensado.

### 8.3 Modelos con captura ajustada y reservorio de línea

El mejor modelo efectivo anterior fue
`ethanol_capture_plus_line_reservoir`. Sus ajustes LORO aproximados fueron:

```text
Octanoato de etilo: vino 0.421; condensado 1.188
Acetato de isoamilo: vino 0.475; condensado 0.985
```

Estimó factores de captura cercanos a 4.0 % y 1.6 %, con multiplicadores por
etanol alrededor de 2.3 y taus de línea alrededor de 4.6–5.1 h. Esos factores
no pueden interpretarse como eficiencias físicas porque no existe medición de la
corriente de salida.

### 8.4 Calibración conjunta imponiendo captura completa

Al ajustar simultáneamente vino y condensado con `eta_AB=1`, el reservorio llevó
`tau` al límite superior de 72 h y degradó fuertemente el ajuste al vino. Esto
demuestra incompatibilidad entre supuestos/datos; no demuestra baja eficiencia
del condensador.

### 8.5 Estrategia actual: vino solamente y cierre externo

La producción se calibra sólo con vino. La masa de condensado se predice sin
participar en el ajuste. Ésta es la comparación científicamente más limpia con
los datos disponibles.

## 9. Resultado actual

Veredicto formal:

```text
NO_VALID_MODEL
```

El mejor candidato diagnóstico es:

```text
assumed_complete_capture_equilibrium
```

Parámetros all-data del candidato basal:

| Compuesto | Formación crecimiento (ug/g azúcar) | Formación estacionaria (ug/g azúcar) |
|---|---:|---:|
| Octanoato de etilo | 1.09527 | 0.74149 |
| Acetato de isoamilo | 8.05407 | 9.67236 |

Métricas LORO primarias:

| Compuesto | NRMSE vino | Umbral vino | NRMSE condensado independiente | Umbral condensado |
|---|---:|---:|---:|---:|
| Octanoato de etilo | 0.41933 | 0.60 | 4.74485 | 1.05509 |
| Acetato de isoamilo | 0.47511 | 0.60 | 6.82270 | 0.96933 |

El bloque de producción pasa el gate interno en vino. La cadena completa falla
el gate independiente de condensado.

Cierre agregado del candidato basal:

| Compuesto | Observado (ug) | Predicho (ug) | Predicho/observado |
|---|---:|---:|---:|
| Octanoato de etilo | 7,519.18 | 32,243.18 | 4.288 |
| Acetato de isoamilo | 22,070.90 | 145,006.99 | 6.570 |

La sobrepredicción ocurre en prácticamente todos los intervalos. Existen
corridas sin muestra de condensado para algunos compuestos; los gráficos ahora
las rotulan explícitamente y no muestran paneles blancos ambiguos.

## 10. Pulsos nutricionales

Se comparó producción basal contra una respuesta retardada posterior al pulso.
El modelo retardado no mejoró la validación LORO:

```text
Octanoato de etilo: NRMSE vino 0.42482, 1.31 % peor que basal
Acetato de isoamilo: NRMSE vino 0.47881, 0.78 % peor que basal
```

No afirmar que el pulso no tiene efecto biológico. La conclusión admisible es
que los datos actuales no respaldan parámetros adicionales de pulso con
capacidad predictiva entre corridas. El pulso y cambios de temperatura del
protocolo siguen parcialmente confundidos.

## 11. Comparación con Mouret/Sablayrolles

Base bibliográfica ya utilizada:

- Mouret et al. 2014, DOI `10.1016/j.foodres.2014.02.044`.
- Morakul et al. 2011, DOI `10.1016/j.procbio.2011.01.034`.
- Mouret et al. 2012, DOI `10.1016/j.lwt.2012.04.031`.

El grupo de Sablayrolles mide concentración en gas mediante GC en línea y
calcula la pérdida como la integral de `C_gas * Q_gas`. La trampa Tenax sirve
como dispositivo analítico, no como supuesto de recuperación completa de un
tren de condensación.

Fracciones de pérdida publicadas aproximadamente a 18 °C:

```text
Acetato de isoamilo: 14–19 %
Octanoato de etilo: 26–34 %
```

La escala de pérdida gaseosa de nuestras iteraciones anteriores era congruente
con esos rangos. La discrepancia actual está en el cierre entre volatilización
predicha y MIX A+B observado, no necesariamente en la escala total de pérdida
desde el vino.

## 12. Artefactos que se deben inspeccionar

Notebook principal ejecutado:

```text
fermentation_model/pilot_2026/notebooks/
pilot_2026_aroma_wine_calibrated_complete_capture.executed.ipynb
```

Tiene 9 celdas, 5 de código ejecutadas, 0 errores y 7 gráficos embebidos.

Configuración y runner:

```text
fermentation_model/pilot_2026/adaptive_design/
aroma_wine_calibrated_complete_capture_config.json

fermentation_model/pilot_2026/
run_aroma_wine_calibrated_complete_capture_2026.py
```

Resultados principales:

```text
fermentation_model/pilot_2026/results/
aroma_wine_calibrated_complete_capture_2026/
```

Archivos prioritarios:

```text
analysis_summary.json
validation_gate.json
condensate_closure_summary.csv
loro_metrics.csv
all_data_parameter_estimates.csv
capture_contract_audit.csv
figures/04_liquid_*.png
figures/05_condensate_*.png
```

El archivo
`adaptive_design/aroma_prospective_validation_contract.json` pertenece a una
etapa anterior y bloquea el mejor modelo efectivo con captura ajustada. Debe
preservarse como trazabilidad, pero **no** debe citarse como el modelo actual ni
como validación de `eta_AB=1`.

## 13. Gates de aceptación

Un candidato completo sólo puede llamarse válido si cumple simultáneamente:

- NRMSE LORO en vino <= 0.60 para ambos compuestos;
- NRMSE independiente de condensado <= 1.05509 para octanoato de etilo;
- NRMSE independiente de condensado <= 0.96933 para acetato de isoamilo;
- convergencia de todos los fits;
- fracción de parámetros en límites <= 0.25;
- error relativo de balance de masa <= `1e-8`;
- sin usar condensado para calibrar los parámetros de producción;
- confirmación prospectiva posterior antes de una afirmación externa fuerte.

No suavizar ni cambiar estos gates después de ver el resultado sin declarar una
nueva versión preregistrada.

## 14. Siguiente plan recomendado

1. Auditar nuevamente unidades y escala de `K(T, ethanol)`, `q_gas` y
   `k_loss`, con pruebas numéricas de dimensiones y casos manuales.
2. Verificar aplicabilidad a 230 L de la formulación tomada de Morakul/Mouret,
   separando equilibrio de partición y limitación de transferencia.
3. Auditar por intervalo y reactor la reconstrucción del MIX A+B, límites de
   los intervalos, volúmenes, muestras faltantes y valores censurados.
4. Ejecutar sensibilidad con incertidumbre bibliográfica explícita para la
   partición/transferencia. No introducir un factor libre que simplemente
   reproduzca la antigua eficiencia efectiva.
5. Mantener calibración de producción sólo en vino y condensado como validación
   externa durante las comparaciones de transferencia.
6. Comparar modelos anidados por LORO y parsimonia. No agregar flexibilidad si no
   mejora ambos dominios y la estabilidad de parámetros.
7. Si ninguna formulación defendible cierra el balance, documentar que el diseño
   actual no identifica el mecanismo y planificar medición prospectiva de gas de
   salida, tercer capturador o balance manual/GC independiente.

## 15. Instrucción de arranque para el nuevo agente

Usar este texto al abrir el chat remoto:

```text
Lee completamente fermentation_model/pilot_2026/HANDOFF_AROMA_MODEL_2026.md.
Luego verifica rama, HEAD, working tree, los artefactos de resultados y los 23
tests indicados. Continúa el objetivo científico descrito allí sin reinterpretar
eta_AB=1 como eficiencia medida, sin usar condensado para recalibrar producción
y sin declarar un modelo válido mientras falle el gate independiente de cierre.
Antes de modificar código, presenta un diagnóstico concreto del término de
transferencia que podría explicar la sobrepredicción 4.288x/6.570x y un plan de
comparación anidada. Ejecuta cargas pesadas mediante Slurm.
```

