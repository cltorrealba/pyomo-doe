# Bundle A02 - Balance de masa y perdidas aromaticas

Actividad: OE2-2.13 - Desarrollar y calibrar modelos de balance de masa para cuantificar perdidas aromaticas

Fecha de cierre documental: 2026-06-22

## Proposito

Este bundle organiza la evidencia tecnica para cubrir el 100% de las categorias solicitadas en la fila de rendicion: ecuaciones/supuestos, datos de calibracion, parametros, cierre de masa, metricas de ajuste, version de codigo y analisis de sensibilidad.

La cobertura es documental y tecnica. Si una revision externa exige validacion experimental directa de gas o condensado, usar este paquete como base y anexar esas mediciones como una capa adicional.

## Archivos principales

- `ANEXO_A02_borrador.md`: texto base listo para transformar en anexo.
- `01_matriz_cumplimiento_A02.csv`: mapeo requisito -> evidencia.
- `02_fila_propuesta_planilla_A02.csv`: fila sugerida para la planilla de seguimiento.
- `03_ecuaciones_y_supuestos_A02.md`: ecuaciones, supuestos y referencias a codigo.
- `08_version_codigo.md`: commit, estado git y hashes.
- `00_manifest.csv`: inventario completo del bundle con rutas y SHA256.

## Tablas generadas

- `tables/02_resumen_datos_fuente.csv`: resumen de libros Excel, hojas, columnas y variables.
- `tables/02_diccionario_variables_clave.csv`: diccionario de estados/variables.
- `tables/04_metricas_calibracion.csv`: metricas de ajuste y capacidad de ajuste.
- `tables/05_parametros_cineticos.csv`: parametros cineticos calibrados.
- `tables/05_parametros_particion_aromas.csv`: parametros de particion aromaticos UNIFAC/surrogado.
- `tables/06_cierre_masa_aromatico.csv`: cierre de masa por candidato y especie.
- `tables/07_sensibilidad_estimabilidad.csv`: FIM, estimabilidad, reduccion de varianza y eigen-resumen.

## Resumen cuantitativo

- Filas de resumen de datos fuente: 44
- Filas de metricas/calibracion: 62
- Filas de cierre de masa aromatico: 54
- Error relativo maximo de cierre de masa aromatico por integral numerica: 0.8547%
- Error relativo maximo de balance de perdida aromatico por integral numerica: 1.9533%
- Filas de sensibilidad/estimabilidad: 97

## Como usarlo para construir el anexo

1. Abrir `ANEXO_A02_borrador.md` y usarlo como cuerpo narrativo del anexo.
2. Insertar o convertir a tablas formales los CSV de `tables/`.
3. Adjuntar como respaldo los reportes de `source_reports/`, las figuras de `figures/` y los datos en `data_sources/`.
4. Usar `01_matriz_cumplimiento_A02.csv` para demostrar que cada requisito minimo tiene una evidencia concreta.
5. Usar `02_fila_propuesta_planilla_A02.csv` para completar la planilla de control de rendicion.
6. Si se actualizan datos o resultados, regenerar con:

```powershell
python scripts/build_A02_evidence_bundle.py
```
