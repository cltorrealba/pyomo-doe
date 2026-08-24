# Operadores de medición recomendados

## Viabilidad y biomasa

En el libro Piloto 2025, `Viability` se comporta como concentración de células
viables en millones de células/mL, no como fracción o porcentaje. No debe
compararse directamente con el estado de biomasa seca `X`.

El repositorio fuente usa 30 pg/célula:

\[
X_{proxy}[g/L]=0.030\,C_{viable}[10^6\,cell/mL].
\]

La copia local de DFVB usa 35,307 pg/célula. Esa diferencia de 17,7 % demuestra
que el factor es un operador de medición, no una constante universal. Para el
artículo se recomienda estimarlo o propagar un rango 30–35,3 pg/célula, y usar
el peso seco como observable separado cuando exista. No se debe mezclar
porcentaje de viabilidad con concentración viable.

## PAN, amoníaco y estado agregado de N

Si el modelo conserva un único estado `N` en base N:

\[
h_N(x)=PAN_N+0.82\,NH3_{compound}=YAN.
\]

PAN y amoníaco no son dos réplicas del mismo estado; son componentes con
cinéticas potencialmente diferentes. Si el modelo los separa, cada estado debe
usar su unidad física y su operador correspondiente. El cierre reportado en
`yan_closure_residual_mg_l` sirve para verificar cada observación.

## Pulsos de nutriente

El valor `pulso_nut` no se usa como salto de N hasta confirmar dosis de producto,
unidad, contenido de N y partición PAN/amoniacal. Una vez confirmados esos datos,
el operador correcto es un evento discreto:

\[
N_i(t_p^+)=N_i(t_p^-)+\Delta N_i,
\]

aplicado después de la muestra pre-pulso si ésa fue la secuencia operacional.
La tabla actual no resuelve por sí sola si la muestra de la fila fue pre o
post-adición; esa convención debe recuperarse de la bitácora.

## Aromas

Los campos `*_total` se interpretan como concentración total equivalente en
vino más condensado. Los campos `*_condensado` son la concentración equivalente
acumulada en el condensador. Por tanto:

\[
w_{total}=w_{retained}+w_{condensate},\qquad
w_{retained}=w_{total}-w_{condensate}.
\]

`w_total` no debe compararse con la tasa instantánea de síntesis. Debe compararse
con la suma de los estados retenido y capturado/acumulado. La tasa de síntesis
alimenta esos estados mediante un balance que incluya pérdida/stripping. Sin
condensado o gas, síntesis y pérdida no son separables de forma identificable.

## Densidad

La densidad no es un estado de azúcar. En el repositorio fuente se ajustó, para
mosto natural, el operador empírico:

\[
G+F\,[g/L]\approx -2119+2.128\,\rho\,[kg/m^3],
\]

con R² 0,9963 y RMSE 3,54 g/L. Puede usarse como proxy con esa incertidumbre,
pero las mediciones directas de glucosa y fructosa prevalecen.

## Temperatura

La temperatura se entrega al modelo como input exógeno, interpolado linealmente
entre los tiempos de observación. No se sincronizan por igualdad exacta de reloj:
el estado se integra sobre una función `T(t)` construida con el perfil y se
evalúa en los tiempos químicos. La derivada pública sólo incluye tiempo relativo;
no exporta fechas calendario.

## Ausentes, ceros y LOD/LOQ

Los vacíos permanecen ausentes. El libro no contiene flags ni límites LOD/LOQ,
por lo que un cero no se recodifica automáticamente como censurado. Los valores
negativos tampoco se recortan silenciosamente: deben revisarse contra el método
analítico. Sólo cuando exista un límite conocido corresponde usar una
verosimilitud censurada o un intervalo `[0, LOD]`.
