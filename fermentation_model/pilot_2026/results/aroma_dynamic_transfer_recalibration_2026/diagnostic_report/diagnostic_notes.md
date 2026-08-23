# Notas de diagnóstico: transferencia y pérdida aromática

Fecha de corte: 2026-08-21.

## Fuentes locales

- `aroma_prediction_metrics.csv`: comparación de RMSE/NRMSE de calibración y holdout.
- `aroma_condensate_predictions.csv`: masa observada o censurada y masa predicha por intervalo MIX.
- `aroma_transfer_diagnostics.csv`: `rCO2`, recambio gaseoso, `K`, eficiencia de transferencia, `kLa`, temperatura, etanol y azúcar.
- `aroma_wine_predictions.csv`: observaciones y predicciones en vino.
- `aroma_parameter_estimates.csv`: parámetros del modelo dinámico y del modelo empírico previo.

## Cálculos auditados

1. La implementación dinámica usa `epsilon = 1-exp(-kLa/qgas)` y `lambda = K*qgas*epsilon`.
2. Para un `kLa` definido sobre la fase líquida y un gas entrante sin aroma, el NTU coherente con el balance `kLa*(CL-CG/K)` es `kLa/(K*qgas)`. La alternativa preferida es resolver explícitamente los balances líquido-headspace y no depender de una fórmula NTU cerrada.
3. Para 26157-MIX-03, octanoato de etilo, la implementación ajustada entrega, en promedio entre 91 y 115 h: `K=0.00519`, `qgas=0.287 h-1`, `kLa=0.00310 h-1` y `epsilon=1.08 %`.
4. En ese MIX se observaron 614.56 ug y se predijeron 5.05 ug. Escalar sólo por la penalización de transferencia da aproximadamente 467 ug; es un diagnóstico de orden de magnitud, no una recalibración.
5. El techo de stripping se calculó como `Mmax = eta_trap*VL*integral(K*qgas*CL_obs dt)`, interpolando linealmente la química líquida. En 26157-MIX-03 da 611.2 ug; en 26211-MIX-04 da 654 ug frente a 2146 ug observados.
6. Ocho de 19 intervalos observados de octanoato superan ese techo. El resultado puede cambiar si hubo máximos líquidos entre muestras químicas; por eso el techo se interpreta como prueba de inconsistencia del operador con los datos disponibles, no como violación termodinámica demostrada.
7. El perfil residual es tardío: entre MIX-02 y MIX-04 el `rCO2` mediano cae de 0.74 a 0.27 g/L/h, mientras la pérdida observada de octanoato sube de 1.62 a 27.37 ug/h.
8. La relación observado/predicho mediana es 34.6x en calibración y 44.0x en holdout para octanoato; 26.9x y 39.3x para acetato de isoamilo.
9. El NRMSE de condensado en holdout prácticamente no cambia al agregar `kLa(E)`: octanoato 1.6193 a 1.6197; isoamilo 1.3018 a 1.2994.

## Interpretación y límites

- `PASS_CONDITIONAL` indica convergencia y propagación de direcciones débiles; no certifica adecuación predictiva ni validez física.
- El acetato de etilo está censurado en los 33 intervalos de condensado, por lo que su pérdida no es identificable con esta campaña.
- Temperatura, programa y reactor/línea están parcialmente confundidos; las asociaciones residuales no permiten causalidad.
- La eficiencia fija de trampa no ha sido validada para cada caudal, composición de gas y temperatura de los dos condensadores.
- No se debe interpretar el `kLa` ni su multiplicador de etanol como parámetros físicos hasta corregir la base del NTU y separar producción, transferencia y captura.

## Criterios de aceptación de la siguiente iteración

- Comparar por validación dejando fuera un reactor/lote completo.
- Reducir NRMSE de condensado de octanoato e isoamilo al menos 30 % sin deteriorar RMSE de vino más de 10 %.
- Eliminar tendencia sistemática de observado/predicho con MIX, etanol, temperatura y volumen de condensado.
- Mantener parámetros dentro de rangos físicos y estables ante cambios razonables de pesos y ventanas de muestreo.
- No introducir escalas libres por MIX; cada nuevo estado o parámetro debe corresponder a una medición o hipótesis contrastable.
