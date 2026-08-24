# Transferencia de datos Piloto 2025 para el artículo DFVB

Este directorio contiene una derivada apta para publicación de los datos de
fermentación piloto 2025 usados en este repositorio. La salida principal tiene
siete procesos homogéneos de desarrollo/calibración y documenta, por separado,
un octavo proceso excluido por inóculo inicial no viable y reinoculación.

No se exportan IDs industriales, IDs de muestra ni fechas calendario. Los
procesos y campañas se identifican mediante `process_token` y `campaign_token`.
La relación privada se debe mantener fuera del repositorio público usando
`private/process_id_map.template.csv`.

## Contenido

- `data/process_observations.csv`: química, temperatura y biomasa de los siete procesos elegibles.
- `data/aroma_observations.csv`: aromas total-equivalente, condensado-equivalente y retenido derivado.
- `data/nitrogen_events.csv`: registros de `pulso_nut`, conservando explícitamente que su unidad no está resuelta.
- `data/process_registry.csv`: proceso → campaña, rango térmico, ventana analítica, uso previo y exclusiones.
- `data/measurement_dictionary.csv`: unidades y operadores recomendados.
- `data/quality_checks.csv`: pruebas y advertencias de calidad.
- `docs/DECISIONS_AND_ANSWERS.md`: respuestas a las dudas del traspaso.
- `docs/MEASUREMENT_OPERATORS.md`: correspondencia dato–estado del modelo.
- `docs/VALIDATION_SCOPE.md`: lenguaje de validación permitido y protocolo confirmatorio.

## Reproducibilidad

El exportador de fuente sólo es reproducible dentro del repositorio que contiene
el workbook restringido:

```powershell
python publication_transfer/DC_dFVB_2026_fermentation_si_2026/scripts/build_transfer_package.py
python publication_transfer/DC_dFVB_2026_fermentation_si_2026/scripts/validate_transfer_package.py
```

El primer comando requiere `pandas`, un motor de lectura de Excel y acceso al
workbook original; se conserva sólo como trazabilidad y no debe ejecutarse en el
repositorio del artículo si el origen restringido no está montado. El validador
no necesita el libro fuente y sí debe ejecutarse después de cada copia.

## Alcance

Los siete procesos elegibles ya fueron utilizados por los pipelines de ajuste o
selección de modelo del repositorio fuente. Por tanto, esta entrega permite
replay, calibración, sensibilidad y validación cruzada interna, pero no una
afirmación de validación independiente. Esa afirmación requiere una campaña
nueva, sellada y adquirida después del freeze del modelo.
