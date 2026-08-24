# Decisiones y respuestas para el traspaso

## 1. Base de unidades de nitrógeno

La base se resolvió por cierre interno del libro fuente. En las 134 filas donde
coexisten los tres ensayos se cumple exactamente:

\[
YAN_{mg\,N/L}=PAN_{mg\,N/L}+0.82\,AMMONIA_{mg\,NH_3/L}.
\]

Por tanto:

- `PAN` está expresado como **mg N/L**.
- `AMMONIA` está expresado como **mg de NH3/L**, no como mg N/L.
- `YAN` está expresado como **mg N/L**.
- La columna derivada `ammonia_n_mg_l = 0.82 × AMMONIA` es la que puede sumarse a PAN.

El factor 0,82 es el factor redondeado N/NH3 usado por el propio libro. No se
debe sumar `PAN + AMMONIA` sin conversión. El cargador vigente del repositorio
fuente etiqueta ambos componentes como si compartieran base y, además, trata
`pulso_nut` como mg/L de N. Esa interpretación no está sustentada por los
metadatos del libro y no se propagó a esta transferencia.

### `pulso_nut`

El libro no contiene unidad, comentario de celda, producto, composición ni base
N-equivalente para `pulso_nut`. La serie registra una carga inicial y un registro
durante el proceso, pero no permite decidir si el número corresponde a mg/L,
g/hL, dosis de producto o mg N/L equivalente. En consecuencia:

- el valor se exporta como `recorded_amount`;
- la unidad queda como `unresolved_in_source_workbook`;
- no debe convertirse ni inyectarse al estado de N hasta recuperar la bitácora o ficha del producto;
- los registros en `t=0` se consideran carga inicial y no un pulso dinámico si el estado inicial ya es post-adición.

### ¿50 % amonio / 50 % PAN?

No es una aproximación central congruente con los datos iniciales. En los siete
procesos elegibles, PAN representa 60,2–74,2 % del YAN inicial y el componente
amoniacal, ya convertido a N, representa 25,8–39,8 %. Si falta el desglose, una
prior operacional cercana a 70/30 PAN-N/amoniacal-N representa mejor estas
condiciones iniciales. Para un pulso de producto de composición desconocida, ni
70/30 ni 50/50 deben fijarse sin ficha técnica: 50/50 puede usarse únicamente
como escenario de sensibilidad y siempre sobre base **N-equivalente**, nunca
mezclando masas de PAN y NH3.

## 2. Mapping experimental

`data/process_registry.csv` entrega:

- proceso público → campaña pública;
- mínimo, media y máximo de temperatura;
- ventana entre primera y última observación analítica;
- número de cargas iniciales y pulsos durante el proceso;
- uso previo en ajuste/selección;
- rol de submission y motivo de exclusión.

`data/nitrogen_events.csv` entrega tiempo, densidad observada y valor registrado
para cada evento. Los siete procesos principales son `development_calibration`.
El proceso con reinoculación queda como `excluded_process_incident`; no es un
set reservado de evaluación. Actualmente no existe un proceso retrospectivo que
pueda declararse independiente.

## 3. Compatibilidad con el repositorio DFVB

La copia local auditada de DFVB usa estados separados de amoníaco y PAN en g/L,
y reinicia el estado amoniacal en el segundo tramo de simulación. Para esta nueva
base se requieren dos precauciones:

1. `AMMONIA / 1000` produce g NH3/L y es compatible con un estado de masa de NH3.
2. `PAN / 1000` produce g N/L, no g/L de una mezcla de aminoácidos. Si el modelo
   representa masa de aminoácidos, se necesita un operador que convierta PAN-N
   a la mezcla asumida mediante su contenido de N; no basta dividir por 1000.

La composición PAN de seis aminoácidos usada por el modelo DFVB es una hipótesis
del modelo, no una medición de composición del ensayo. Debe declararse como tal
y evaluarse por sensibilidad.

## 4. Identificadores

Los CSV públicos sólo contienen tokens. La tabla privada puede completarse a
partir de `private/process_id_map.template.csv` y, si se desea, protegerse con
HMAC usando una clave mantenida fuera del repositorio. Nunca se debe commitear la
tabla completa ni la clave.
