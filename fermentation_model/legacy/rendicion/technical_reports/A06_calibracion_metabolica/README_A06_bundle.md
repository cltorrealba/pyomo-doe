# Bundle A06 - Calibracion del submodelo metabolico

Actividad: OE3-3.13 - Calibracion del submodelo metabolico del Gemelo Digital (Actividad 3.13)

Fecha de cierre documental: 2026-06-26

## Proposito

Este bundle organiza evidencia para el anexo solicitado: dataset, estrategia de calibracion, parametros metabolicos, funcion objetivo, metricas, incertidumbre/identificabilidad, resultados e interpretacion.

## Tabla del anexo

| ID anexo   | OE   | Actividad unica   | Titulo propuesto                                                         | Nombre de archivo sugerido                                    | Contenido minimo                                                                                                 | Evidencia fuente                                                                                                          |
|:-----------|:-----|:------------------|:-------------------------------------------------------------------------|:--------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------|
| A06        | OE3  | 3.13              | Calibracion del submodelo metabolico del Gemelo Digital (Actividad 3.13) | PI-4497_IA5_A06_Act-3.13_Calibracion-submodelo-metabolico.pdf | Documentar dataset, estrategia de calibracion, parametros metabolicos, metricas, incertidumbre e interpretacion. | Dataset de calibracion, parametros metabolicos, funcion objetivo, metricas, incertidumbre/identificabilidad y resultados. |

## Archivos principales

- `ANEXO_A06_borrador.md`: texto base del anexo.
- `01_tabla_anexo_A06.csv`: tabla con las columnas solicitadas.
- `02_matriz_cumplimiento_A06.csv`: mapeo contenido minimo -> evidencia.
- `03_metodologia_calibracion_A06.md`: dataset, estrategia y funcion objetivo.
- `04_resultados_interpretacion_A06.md`: parametros, metricas, incertidumbre e interpretacion.
- `08_version_codigo.md`: version de codigo y estado git.
- `00_manifest.csv`: inventario completo con SHA256.

## Resumen

- Archivos del bundle: 134
- Filas del dataset integrado: 520
- Cumplimiento documental: 7/7 items cubiertos
- WSSE secondary_v2_reduced_o2fixed: 5756.427666344939
- Parametros debiles/no identificables reportados: 14

## Uso recomendado

1. Usar `ANEXO_A06_borrador.md` como cuerpo narrativo.
2. Insertar `01_tabla_anexo_A06.csv` en la matriz general de anexos.
3. Usar `02_matriz_cumplimiento_A06.csv` para justificar cobertura.
4. Adjuntar `data_sources/`, `source_reports/`, `source_tables/`, `figures/`, `code_sources/` y `notebooks/`.
5. Convertir el anexo final con el nombre sugerido: `PI-4497_IA5_A06_Act-3.13_Calibracion-submodelo-metabolico.pdf`.

Para regenerar:

```powershell
python scripts/build_A06_metabolic_calibration_bundle.py
```
