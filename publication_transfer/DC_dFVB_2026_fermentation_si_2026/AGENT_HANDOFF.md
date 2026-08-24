# Handoff para el agente de `paper/fermentation-si-2026`

## Alcance autorizado

- Trabajar sólo en la rama `paper/fermentation-si-2026` del checkout
  `/home/cltorrealba/DC_dFVB_2026_fermentation_si_2026`.
- Se permiten commits locales; no hacer push.
- No agregar el directorio local `repro/article_fermentation_MPDI/`.
- No publicar workbooks originales, IDs industriales, tabla privada completada,
  claves HMAC ni rutas privadas.
- La privacidad/URL final del remoto está pendiente.
- Este handoff no autoriza jobs Slurm. Antes de crear un goal o ejecutar jobs se
  requiere confirmación explícita de la política exacta para `N=0`.

## Qué contiene el paquete

- Siete procesos piloto 2025 elegibles como desarrollo/calibración.
- Un octavo proceso aparece sólo en `process_registry.csv` como exclusión por
  inóculo no viable y reinoculación; sus observaciones no están en los datasets
  principales.
- 249 filas de proceso/química, 344 filas de aromas y 14 eventos nutricionales.
- Tokens públicos `P25_XX` y `C25_XX`; no hay IDs industriales ni fechas absolutas.

## Decisiones científicas que el modelo debe respetar

1. `PAN` está en mg N/L.
2. `AMMONIA` está en mg NH3/L.
3. `YAN = PAN + 0.82 * AMMONIA` cierra exactamente en el origen.
4. No sumar PAN y AMMONIA crudos ni tratarlos como la misma unidad.
5. La distribución inicial media es aproximadamente 68/32 PAN-N/amoniacal-N;
   50/50 sólo puede ser sensibilidad sobre base N-equivalente.
6. `pulso_nut` no tiene unidad/composición demostrada. No convertirlo a salto de
   N hasta recuperar producto, dosis y secuencia pre/post-muestra.
7. `Viability` se interpreta como millones de células viables/mL. Comparar con
   biomasa mediante un operador de masa/célula, no directamente.
8. `w_total` de aroma representa retenido + condensado acumulado; no es una tasa
   instantánea de síntesis.
9. Densidad es un proxy empírico de G+F y no un estado mecanístico directo.
10. Ceros y negativos se preservan porque el origen no incluye LOD/LOQ.

## Implicancia específica para DFVB

El modelo DFVB separa amoníaco y PAN, pero PAN-N/1000 produce g N/L, no g/L de
una mezcla de aminoácidos. Si los estados DFVB representan masa o moles de la
mezcla, implementar un operador PAN-N → aminoácidos basado en contenido de N y
mantener esa composición como hipótesis/sensibilidad del modelo.

Los siete procesos ya fueron usados por pipelines previos. Reportar replay,
bootstrap o validación cruzada por campaña; no llamarlos validación independiente.

## Primeros pasos recomendados

1. Ejecutar `python scripts/validate_transfer_package.py` dentro del paquete.
2. Leer `docs/DECISIONS_AND_ANSWERS.md` y `docs/MEASUREMENT_OPERATORS.md`.
3. Integrar un loader que conserve por separado `pan_n_mg_l`,
   `ammonia_compound_mg_l` y `ammonia_n_mg_l`.
4. Mantener `pulso_nut` fuera del balance hasta resolver sus metadatos.
5. Antes de cualquier cálculo con `N=0`, pedir/registrar la política aprobada y
   agregar pruebas de frontera; no inferirla desde este paquete.
