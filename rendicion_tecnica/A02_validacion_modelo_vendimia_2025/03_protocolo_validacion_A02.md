# Protocolo de validacion - A02

## Objetivo

Documentar la validacion del modelo con base independiente de vendimia 2025, cubriendo datos fuente, protocolo, comparacion predicho-observado, metricas de error, brechas y decision tecnica.

## Base independiente

La base independiente 2025 se compone de:

- `fermentation_model/data/mosto_natural_xthiol.xlsx`
- `fermentation_model/data/mosto_sintetico_vl3.xlsx`

La base historica `Calibration_data_vl3.xlsx` se incluye como referencia de calibracion previa y auditoria, pero la evidencia principal de vendimia 2025 esta en los libros natural/sintetico homologados y en `results/new_must_data_loading`.

## Estados validados

- Biomasa viable (`X`)
- Biomasa muerta (`Xd`)
- Nitrogeno asimilable (`N`)
- Glucosa (`G`)
- Fructosa (`F`)
- Etanol (`E`)
- Glicerol (`Gly`)

## Procedimiento

1. Cargar y normalizar los datos 2025 con `new_must_data_loader.py`.
2. Resumir la base por lote, medio, horizonte, condiciones iniciales y numero de observaciones.
3. Simular curvas con theta base y theta best-multistart.
4. Comparar predicho vs observado por lote/estado mediante RMSE y weighted RMSE.
5. Evaluar robustez base-vs-best con `best_minus_base_weighted_rmse`.
6. Revisar factibilidad de predicciones DOE: etanol, biomasa, biomasa muerta, glicerol y azucar residual.
7. Detectar sensibilidad a theta y redundancias de trayectorias.
8. Declarar brechas y decision tecnica.

## Criterio de validacion usado en este bundle

- La base debe contener lotes naturales y sinteticos de vendimia 2025.
- La comparacion predicho-observado debe estar disponible por medio/estado.
- Las simulaciones seleccionadas no deben activar banderas duras de factibilidad.
- La alternativa best-multistart debe mejorar o empatar la mayoria de combinaciones medio/estado.
- El set recomendado debe incluir una prueba natural para transferencia a mosto real.
- Las brechas de identificabilidad se documentan como condicion de uso secuencial, no como bloqueo del anexo.

## Decision operacional derivada

Validacion tecnica condicionada: el modelo es utilizable para diseno y planificacion secuencial con datos de vendimia 2025. El primer bloque recomendado debe ser mixto y debe incluir transferencia a mosto natural; despues del bloque 1 se debe recalibrar y rerankear.
