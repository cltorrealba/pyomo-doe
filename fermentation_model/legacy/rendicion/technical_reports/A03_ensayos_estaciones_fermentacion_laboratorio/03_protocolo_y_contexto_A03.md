# Protocolo y contexto experimental - A03

## Objetivo

Documentar los ensayos en estaciones de fermentacion a escala laboratorio asociados a la generacion del archivo `mosto_sintetico_vl3.xlsx`.

## Fuentes principales

- `fermentation_model/data/mosto_sintetico_vl3.xlsx`
- PDF de referencia sobre sistema/estacion de vinificacion de laboratorio.
- PDF de referencia sobre experimentos a escala laboratorio, vendimia 2024.

## Diseno experimental

El archivo contiene un diseno CCD de dos factores resumido en la hoja `Diseno_CCD`:

- Cubadas: MS007 a MS016.
- Factor termico: temperatura operacional de la columna `Temperatura`.
- Factor de adicion: columna `Adicion`, tambien mapeada a `pulso_nut` para compatibilidad con el modelo.
- Registros de proceso: fecha/hora, horas, Brix, densidad, temperatura medida y operacional, nutrientes, azucares, biomasa/viabilidad, etanol, glicerol, metabolitos y aromas totales.

## Protocolo de documentacion del bundle

1. Copiar el Excel fuente y exportar hojas clave a CSV.
2. Resumir matriz de tratamientos por cubada.
3. Resumir ejecucion temporal por cubada.
4. Auditar cobertura de registros de proceso.
5. Registrar incidencias/flags de calidad.
6. Consolidar resultados finales por cubada.
7. Declarar cumplimiento y brechas.

## CO2

No se incluyen registros de CO2 para estos ensayos. La razon operacional indicada para este anexo es que los sensores de CO2 aun no estaban instalados durante la ejecucion de estos ensayos. Por tanto, CO2 se documenta como registro no disponible con justificacion, no como omision de procesamiento.

## Documentos PDF de referencia localizados

| document_key                 | filename                                                                |   pages | status   |
|:-----------------------------|:------------------------------------------------------------------------|--------:|:---------|
| anexo_3_sistema_vinificacion | Anexo 3 - Nuevo sistema de vinificación laboratorio.pdf                 |      10 | copied   |
| anexo_8_experimentos_2024    | Anexo 8 - Informe Experimentos a Escala Laboratorio – Vendimia 2024.pdf |      12 | copied   |
