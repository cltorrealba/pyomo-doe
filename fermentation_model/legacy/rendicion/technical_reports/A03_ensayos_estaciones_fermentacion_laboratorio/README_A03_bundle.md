# Bundle A03 - Ensayos en estaciones de fermentacion a escala laboratorio

Actividad: OE3-3.10 - Ensayos en estaciones de fermentacion a escala laboratorio (Actividad 3.10)

Fecha de cierre documental: 2026-06-23

## Proposito

Este bundle organiza la evidencia para el anexo solicitado: diseno experimental, tratamientos, equipos, ejecucion, registros de proceso, incidencias, resultados y cumplimiento, asociado al archivo `mosto_sintetico_vl3.xlsx`.

## Tabla del anexo

| ID anexo   | OE   | Actividad unica   | Titulo propuesto                                                            | Nombre de archivo sugerido                                   | Contenido minimo                                                                                                                | Evidencia fuente                                                                                                                                      |
|:-----------|:-----|:------------------|:----------------------------------------------------------------------------|:-------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------|
| A03        | OE3  | 3.10              | Ensayos en estaciones de fermentacion a escala laboratorio (Actividad 3.10) | PI-4497_IA5_A03_Act-3.10_Ensayos-estaciones-fermentacion.pdf | Documentar diseno experimental, tratamientos, equipos, ejecucion, registros de proceso, incidencias, resultados y cumplimiento. | Diseno experimental, matriz de tratamientos, protocolos, IDs de fermentacion, registros de temperatura/nutrientes/densidad/CO2, fotos y datos crudos. |

## Archivos principales

- `ANEXO_A03_borrador.md`: texto base del anexo.
- `01_tabla_anexo_A03.csv`: tabla con las columnas de la imagen.
- `02_matriz_cumplimiento_A03.csv`: mapeo contenido minimo -> evidencia.
- `03_protocolo_y_contexto_A03.md`: protocolo, fuentes y nota CO2.
- `04_resultados_y_cumplimiento_A03.md`: resultados, incidencias y cumplimiento.
- `08_version_codigo.md`: version de codigo y hashes.
- `00_manifest.csv`: inventario completo con SHA256.

## Resumen

- Cubadas documentadas: 10
- Filas homologadas: 86
- Estado registro CO2: `no_disponible_justificado`
- Cumplimiento: 8/8 items cubiertos

## Uso recomendado

1. Usar `ANEXO_A03_borrador.md` como cuerpo narrativo.
2. Insertar `01_tabla_anexo_A03.csv` en la matriz general de anexos.
3. Usar `02_matriz_cumplimiento_A03.csv` para justificar cobertura.
4. Adjuntar `data_sources/`, `reference_documents/`, `source_tables/`, `tables/` y `figures/`.
5. Convertir el anexo final con el nombre sugerido: `PI-4497_IA5_A03_Act-3.10_Ensayos-estaciones-fermentacion.pdf`.

Para regenerar:

```powershell
python scripts/build_A03_fermentation_station_bundle.py
```
