# Decision tecnica - A02

## Decision

`validacion_tecnica_condicionada`

## Alcance

Modelo aceptable para planificacion y diseno secuencial con base independiente vendimia 2025.

## Condiciones de uso

Ejecutar primer bloque mixto, incluir una prueba natural, recalibrar y recomputar perfiles/FIM antes de cerrar bloques siguientes.

## Criterios evaluados

| criterio                                | umbral_o_regla                                                                            | valor_observado                            | estado      | interpretacion                                                                                                    |
|:----------------------------------------|:------------------------------------------------------------------------------------------|:-------------------------------------------|:------------|:------------------------------------------------------------------------------------------------------------------|
| Base independiente 2025 cargada         | Debe contener lotes naturales y sinteticos 2025 homologados.                              | natural=9; synthetic=10                    | cumple      | La base independiente contiene ambos dominios de validacion.                                                      |
| Robustez theta base vs mejor multistart | Mejor multistart no debe deteriorar la mayoria de combinaciones medio/estado.             | 13/14 combinaciones mejoran o empatan      | cumple      | La alternativa multistart mejora casi todos los estados/medios; el deterioro principal reportado es Xd sintetico. |
| Factibilidad de predicciones DOE        | No activar banderas duras en etanol, biomasa, biomasa muerta, glicerol o azucar residual. | issue_count_total=0                        | cumple      | Las simulaciones seleccionadas quedan en rangos plausibles.                                                       |
| Prueba de transferencia a mosto natural | El set de validacion/diseno debe incluir al menos un candidato natural.                   | n_candidatos_naturales=2                   | cumple      | Existe prueba natural para transferencia de matriz, pero debe mantenerse en el primer bloque.                     |
| Sensibilidad a theta                    | Registrar desacuerdo base-vs-best para priorizar experimentos informativos.               | max_theta_disagreement_weighted_rmse=2.168 | advertencia | Alto desacuerdo implica informacion util, pero no cierre automatico de parametros debiles.                        |
| Redundancia de trayectorias             | Pares con distancia <=0.25 y correlacion >=0.90 requieren seleccion manual.               | pares_redundantes=2                        | advertencia | No ejecutar pares redundantes en el primer bloque; elegir el que ataque mejor la brecha.                          |
| Identificabilidad practica              | Registrar parametros no identificables y tratarlos como brecha tecnica.                   | 4/17 parametros identificables al 95%      | brecha      | La validacion permite uso secuencial, no cierre definitivo de todos los parametros.                               |

## Brechas y acciones

| brecha                                                | evidencia                                                                                                             | impacto                                             | accion_recomendada                                                                             | estado                                  |
|:------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------|:-----------------------------------------------------------------------------------------------|:----------------------------------------|
| Identificabilidad practica incompleta                 | Perfiles de verosimilitud con direcciones one-sided o no identificables.                                              | No declarar cierre parametrico definitivo.          | Usar campana secuencial: ejecutar bloque 1, recalibrar, recomputar FIM/perfiles y rerankear.   | controlada_para_validacion_condicionada |
| Transferencia de medio requiere prueba natural        | El reporte de transferencia indica que la matriz natural debe mantenerse para distinguir transferencia vs estructura. | No usar solo mosto sintetico como validacion final. | Mantener `natural_cold_hot_switch` o `natural_glucose_pulse_after_growth` en el primer bloque. | controlada_con_bloque_mixto             |
| Redundancia entre algunos candidatos                  | Pares con baja distancia RMS y alta correlacion de trayectoria.                                                       | Baja eficiencia si ambos se ejecutan temprano.      | No poner ambos en el primer bloque; escoger segun parametro debil objetivo.                    | controlada_por_protocolo                |
| Deterioro menor en Xd sintetico bajo mejor multistart | Curve validation report: synthetic Xd es la unica combinacion con deterioro leve.                                     | Xd permanece observable sensible/ruidoso.           | Revisar calidad de conteos de biomasa muerta y priorizar medidas consistentes en bloque 1.     | advertencia                             |
| Par candidato redundante detectado                    | synthetic_late_ethanol_death_probe vs synthetic_ethanol_inhibition_challenge; distancia=0.203; correlacion=0.966      | Puede aportar informacion duplicada.                | Ejecutar solo uno en bloque inicial.                                                           | advertencia                             |
| Par candidato redundante detectado                    | natural_glucose_pulse_after_growth vs natural_cold_hot_switch; distancia=0.243; correlacion=0.943                     | Puede aportar informacion duplicada.                | Ejecutar solo uno en bloque inicial.                                                           | advertencia                             |

## Texto sugerido para informe

Se valido el modelo contra una base independiente de vendimia 2025 compuesta por fermentaciones en mosto natural y sintetico. La comparacion predicho-observado fue evaluada mediante metricas de error por lote, estado y medio, y se contrastaron las predicciones bajo theta base y una solucion best-multistart. Las simulaciones seleccionadas no activaron banderas duras de factibilidad y el analisis de curvas mostro mejora o estabilidad en la mayoria de estados/medios. Se identifican brechas de identificabilidad practica y redundancia entre algunos candidatos, por lo que la decision tecnica es validar el modelo para uso secuencial condicionado: ejecutar un primer bloque mixto con al menos una prueba natural, recalibrar y recomputar FIM/perfiles antes de cerrar la campana completa.
