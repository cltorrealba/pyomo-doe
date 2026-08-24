# Informe de validación del paquete

## Evaluación global: compartible con caveats

El paquete es reproducible y estructuralmente consistente para el uso de los
siete procesos como desarrollo/calibración. No está listo para convertir
`pulso_nut` en una entrada cuantitativa de N ni para sostener validación
independiente.

Fuente auditada: `Calibration_data_vl3.xlsx`, campaña piloto 2025. Revisión y
empaquetado: 2026-08-20.

## Verificaciones realizadas

- 8 procesos fuente; 7 incluidos y 1 excluido por incidente de reinoculación.
- 249 filas proceso-tiempo, sin duplicados y sin claves nulas.
- El conteo por proceso coincide exactamente con el registro de campaña.
- 344 filas proceso-tiempo-especie de aromas, con 8 especies y sin duplicados.
- 14 registros nutricionales: una carga inicial y un pulso de proceso por cada proceso elegible.
- 110 observaciones completas de PAN/amoníaco/YAN en los siete procesos cierran exactamente con `YAN = PAN + 0.82·AMMONIA`.
- PAN representa 60,2–74,2 % del YAN inicial; promedio simple 68,4 %. Esto no respalda 50/50 como supuesto central.
- No se detectaron IDs industriales, IDs de muestra ni fechas calendario en los CSV públicos.
- No se detectaron valores de aroma retenido menores que total menos condensado fuera de tolerancia numérica.

## Issues y caveats

1. **Alto — unidad/composición de `pulso_nut` sin resolver.** El libro no indica
   unidad, producto ni equivalencia de N. Los siete registros durante proceso
   tienen valor 100, ocurren a 24, 36 o 42 h y a densidades 1036–1059 kg/m³.
   No deben transformarse en un salto de N hasta recuperar la bitácora.
2. **Alto — alcance de validación.** Los siete procesos ya fueron utilizados en
   ajuste o selección de modelo; sólo permiten evaluación interna retrospectiva.
3. **Medio — operador PAN para DFVB.** PAN está en mg N/L. El modelo DFVB separa
   aminoácidos en base masa/molar; dividir PAN por 1000 no produce masa de la
   mezcla de aminoácidos.
4. **Medio — biomasa.** La conversión de células viables a g/L depende de 30
   versus 35,307 pg/célula entre repositorios. Debe estimarse o analizarse por
   sensibilidad.
5. **Medio — LOD/LOQ ausente.** Hay 44 ceros de amoníaco y valores cero de
   azúcares sin flags de censura. Además, existen dos pesos secos negativos
   (mínimo −0,0145 g/L). Se preservaron; no deben recortarse silenciosamente.
6. **Bajo — condensado incompleto.** Hay 56 filas de aroma total sin medición de
   condensado. En ellas se compara `total` contra retenido + capturado, no se
   inventa el pool retenido.

## Estado del handoff

El checkout del artículo está accesible por SSH en el clúster y el paquete puede
ser incorporado mediante commit local. La privacidad y URL final del remoto
siguen pendientes; por eso no corresponde hacer push. La plantilla MDPI y los
datos originales permanecen fuera de Git.

## Veredicto

Los datos y operadores documentados son aptos para compartir como paquete de
desarrollo/calibración. Antes de ejecutar el modelo DFVB se deben resolver la
unidad/composición de `pulso_nut` y el operador PAN-N → mezcla de aminoácidos.
Para una afirmación confirmatoria se requiere una nueva campaña sellada.
