# Resultados y cumplimiento - A03

## Resultado tecnico

Se consolido la evidencia documental de los ensayos en estaciones de fermentacion a escala laboratorio asociados al archivo `mosto_sintetico_vl3.xlsx`. El bundle incluye diseno CCD, matriz de tratamientos, datos crudos y homologados, IDs de fermentacion, registros de proceso, flags de calidad, resultados finales por cubada y documentos de referencia de equipos/estacion.

## Cumplimiento

| criterio                   | estado                    | valor                         | comentario                                                                                       |
|:---------------------------|:--------------------------|:------------------------------|:-------------------------------------------------------------------------------------------------|
| Cobertura contenido minimo | cumple                    | 8/8 items cubiertos           | Todos los items se cubren documentalmente; CO2 se cubre como ausencia justificada.               |
| Registros CO2              | no_disponible_justificado | columnas_presentes=           | No hay registros CO2 en estos ensayos porque los sensores aun no estaban instalados.             |
| Datos crudos preservados   | cumple                    | Excel fuente y CSV exportados | La hoja 00_Resumen indica que no se modificaron valores crudos; flags solo recomiendan revision. |

## Resultados por cubada

| Cubada   |   n_registros |   t_final_h |   temperatura_setpoint |   densidad_inicial |   densidad_final |   GLUCOSE_inicial |   GLUCOSE_final |   FRUCTOSE_inicial |   FRUCTOSE_final |   YAN_inicial |   YAN_final |   ETANOL_final |   GLYCEROL_final |   Viability_final |   Peso_Seco_final |   Ethyl_Acetate_total_final |   isoamil_acetate_total_final |   octanoate_de_etilo_total_final |
|:---------|--------------:|------------:|-----------------------:|-------------------:|-----------------:|------------------:|----------------:|-------------------:|-----------------:|--------------:|------------:|---------------:|-----------------:|------------------:|------------------:|----------------------------:|------------------------------:|---------------------------------:|
| MS007    |             9 |       213   |                     15 |             1098.6 |           1020.9 |            125.34 |           14.64 |             120.52 |            49.08 |           224 |         119 |          10.04 |             6.76 |             63.57 |               nan |                         nan |                           nan |                              nan |
| MS008    |             9 |       213   |                     15 |             1097.1 |           1014.7 |            124.28 |            9.99 |             119.63 |            42.07 |           218 |          60 |          10.64 |             6.91 |             69.76 |               nan |                         nan |                           nan |                              nan |
| MS009    |             9 |       213   |                     23 |             1095.4 |            989.8 |            120.89 |            1.13 |             114.19 |             0.89 |           236 |          34 |          14.02 |             7.24 |             48.48 |               nan |                         nan |                           nan |                              nan |
| MS010    |             9 |       213.5 |                     19 |             1097.9 |            992   |            124.65 |            0.96 |             113.21 |             3.92 |           221 |         nan |          14.04 |             7.99 |             71.54 |               nan |                         nan |                           nan |                              nan |
| MS011    |             9 |       213.5 |                     19 |             1097.2 |            996.4 |            123.57 |            2.48 |             119.77 |            11.8  |           218 |         nan |          13.35 |             7.93 |             70.96 |               nan |                         nan |                           nan |                              nan |
| MS012    |             9 |       213.5 |                     23 |             1095.3 |            988.8 |            122.56 |            2.97 |             114.23 |            -1.44 |           223 |         nan |          14.2  |             7.68 |             56.97 |               nan |                         nan |                           nan |                              nan |
| MS013    |             8 |       213   |                     19 |             1098.4 |            993.9 |            125.07 |            1.09 |             116.2  |             8.43 |           218 |          90 |          13.9  |             8.08 |             67.21 |               nan |                         nan |                           nan |                              nan |
| MS014    |             8 |       213   |                     19 |             1101.2 |            995.6 |            127.15 |            1.54 |             126.5  |            12.28 |           226 |          98 |          13.9  |             7.92 |             67.19 |               nan |                         nan |                           nan |                              nan |
| MS015    |             8 |       213   |                     15 |             1098.3 |           1016.3 |            121.91 |           11.15 |             118.49 |            43    |           220 |         111 |          10.8  |             7.39 |             71.76 |               nan |                         nan |                           nan |                              nan |
| MS016    |             8 |       208   |                     23 |             1098.5 |            989   |            125.22 |            0.6  |             116.88 |             0.12 |           214 |         129 |          14.5  |             7.86 |             54.08 |               nan |                         nan |                           nan |                              nan |

## Incidencias registradas

| Codigo   | Cubada   |   Horas | Columna   |   Valor | Criterio   | Accion                                                                        |
|:---------|:---------|--------:|:----------|--------:|:-----------|:------------------------------------------------------------------------------|
| MS009-2  | MS009    |    22.5 | Brix      |   45797 | Brix > 40  | Valor crudo preservado; revisar si corresponde a error de carga/fecha serial. |

## Nota CO2

No hay registros de CO2 en estos ensayos porque los sensores aun no estaban instalados. Esta salvedad queda trazada en `tables/06_cobertura_registros_proceso.csv` y en la matriz de cumplimiento.
