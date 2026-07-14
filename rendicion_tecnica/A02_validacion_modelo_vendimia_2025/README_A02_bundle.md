# Bundle A02 - Validacion del modelo con datos de vendimia 2025

Actividad: OE2-2.14 - Validacion del modelo con datos de vendimia 2025 (Actividad 2.14)

Fecha de cierre documental: 2026-06-23

## Proposito

Este bundle organiza la evidencia para el anexo solicitado en la imagen: base independiente 2025, protocolo de validacion, predicho vs observado, metricas, brechas y decision tecnica.

## Tabla del anexo

| ID anexo   | OE   |   Actividad unica | Titulo propuesto                                                  | Nombre de archivo sugerido                                   | Contenido minimo                                                                                                                      | Evidencia fuente                                                                                                                 |
|:-----------|:-----|------------------:|:------------------------------------------------------------------|:-------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------|
| A02        | OE2  |              2.14 | Validacion del modelo con datos de vendimia 2025 (Actividad 2.14) | PI-4497_IA5_A02_Act-2.14_Validacion-modelo-vendimia-2025.pdf | Documentar la base independiente 2025, protocolo de validacion, comparacion predicho-observado, metricas, brechas y decision tecnica. | Base independiente vendimia 2025, predicho vs. observado, metricas de error, criterio de validacion, brechas y decision tecnica. |

## Archivos principales

- `ANEXO_A02_borrador.md`: texto base del anexo.
- `01_tabla_anexo_A02.csv`: tabla con las columnas de la imagen.
- `02_matriz_cumplimiento_A02.csv`: mapeo contenido minimo -> evidencia.
- `03_protocolo_validacion_A02.md`: protocolo y criterio de validacion.
- `04_decision_tecnica_A02.md`: decision, brechas y condiciones.
- `08_version_codigo.md`: version de codigo y hashes.
- `00_manifest.csv`: inventario completo con SHA256.

## Tablas generadas

- `tables/03_base_independiente_2025_resumen_lotes.csv`
- `tables/03_base_independiente_2025_observaciones.csv`
- `tables/03_base_independiente_2025_diccionario.csv`
- `tables/04_predicho_observado_metricas.csv`
- `tables/04_metricas_error_por_medio_estado.csv`
- `tables/04_residuales_agregados_base_vs_best.csv`
- `tables/05_criterios_validacion.csv`
- `tables/06_brechas_y_acciones.csv`
- `tables/07_decision_tecnica.csv`

## Resumen cuantitativo

- Lotes independientes 2025 resumidos: 19
- Filas de metricas por medio/estado: 14
- Criterios en estado `cumple`: 4/7

## Uso recomendado

1. Usar `ANEXO_A02_borrador.md` como cuerpo narrativo.
2. Insertar la tabla `01_tabla_anexo_A02.csv` en la matriz general de anexos.
3. Usar `02_matriz_cumplimiento_A02.csv` para demostrar que el contenido minimo esta cubierto.
4. Adjuntar `source_reports/`, `source_tables/` y `figures/` como respaldo.
5. Convertir el anexo final con el nombre sugerido: `PI-4497_IA5_A02_Act-2.14_Validacion-modelo-vendimia-2025.pdf`.

Para regenerar:

```powershell
python scripts/build_A02_validation_bundle.py
```
