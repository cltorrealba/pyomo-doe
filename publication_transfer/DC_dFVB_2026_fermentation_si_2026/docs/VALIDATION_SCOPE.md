# Alcance de validación

## Ruta admisible con los datos actuales

Los siete procesos elegibles se usaron previamente en los pipelines de ajuste,
selección estructural o evaluación del repositorio fuente. Deben etiquetarse como
**desarrollo/calibración**. Con ellos se puede reportar:

- replay o bondad de ajuste retrospectiva;
- validación cruzada interna, idealmente dejando fuera campañas pareadas completas y no filas aleatorias;
- bootstrap por proceso/campaña;
- sensibilidad a masa por célula, partición de N y unidad/composición del pulso;
- desempeño por estado y no sólo una función de pérdida agregada.

No debe usarse “validación independiente”, “validación externa” ni
“confirmación prospectiva” para estos siete procesos. Reservar ahora uno de
ellos sería un holdout post-hoc porque el desarrollo ya tuvo acceso a sus datos.
El proceso con reinoculación tampoco es un candidato de validación: viola el
protocolo homogéneo y se mantiene sólo como análisis de robustez/incidente.

## Ruta confirmatoria

Para una afirmación prospectiva se debe:

1. congelar código, parámetros, operadores de medición, criterios de exclusión y métricas;
2. registrar el hash/commit del freeze;
3. ejecutar una campaña nueva después del freeze;
4. ocultar sus resultados al equipo de calibración hasta cerrar el análisis;
5. aplicar el pipeline sin reajuste;
6. reportar todos los procesos, fallos y exclusiones predefinidas.

La campaña sellada debe medir al menos temperatura, glucosa/fructosa, PAN,
amoníaco, peso seco o recuento viable, etanol y tiempos/composición real de los
pulsos. Si se evalúan aromas, debe incluir líquido y condensado con LOD/LOQ. Si
el artículo formula conclusiones sobre transferencia de CO2, se requiere además
un sensor CO2 válido y un criterio de curation fijado antes de abrir los datos.

## Set reservado actual

No hay procesos retrospectivos reservados para evaluación independiente. El
campo `independent_validation` de `process_registry.csv` es falso para todos los
procesos por diseño.
