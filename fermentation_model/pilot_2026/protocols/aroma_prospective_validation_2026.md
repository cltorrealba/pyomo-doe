# Protocolo prospectivo para validar producción y pérdida aromática

## Objetivo

Validar, sin recalibración, el modelo `ethanol_capture_plus_line_reservoir` para octanoato de etilo y acetato de isoamilo. La prueba separa tres cantidades que hoy están confundidas: producción en el líquido, emisión por el gas y recuperación del tren de condensado.

El contrato legible por máquina está en `adaptive_design/aroma_prospective_validation_contract.json`. Sus nombres de archivo, columnas y unidades son obligatorios.

## Diseño mínimo

- Tres reactores piloto independientes.
- Dos reactores repiten exactamente temperatura, inoculación y pulso nutricional.
- El tercero puede usar la secuencia térmica alternativa, pero sus predicciones se calculan con parámetros congelados.
- No se permite ajustar parámetros después de ver los aromas.
- Los IDs y el script de evaluación se congelan antes de recibir resultados analíticos.

## Calibración del tren de captura

Antes, a mitad y después de la campaña, inyectar estándares gaseosos de ambos compuestos en tres niveles. Para cada nivel realizar tres réplicas y registrar masa alimentada, masa recuperada, flujo, temperatura y fracción de etanol del portador.

La recuperación debe recalcularse como `recovered_mass_ug / input_mass_ug`. El CV por nivel debe ser ≤15 %. Un estándar fuera de criterio invalida los condensados comprendidos entre las dos sesiones de calibración que lo rodean.

## Muestreo durante fermentación

Tomar líquido y gas de salida pareados a −2, 0, +1, +3, +6, +12 y +24 h respecto del pulso. Agregar puntos al 10, 25, 50, 75 y 90 % del avance por CO₂ acumulado si no coinciden con los anteriores.

Recolectar condensado en intervalos de máximo 12 h entre −6 y +24 h del pulso, y máximo 24 h fuera de esa ventana. En cada MIX registrar volumen real y fracción volumétrica de etanol, no sólo concentración aromática.

Cada determinación de aroma debe tener al menos dos réplicas analíticas. Conservar cromatogramas, integración, dilución, LOD, LOQ y motivo de exclusión.

## Mediciones obligatorias

- `process.csv`: temperatura, setpoint, rCO₂, etanol, azúcar y volumen líquido.
- `events.csv`: pulso nutricional y cambios de temperatura.
- `liquid_aroma.csv`: concentración en vino.
- `outlet_gas.csv`: concentración en gas y flujo durante cada intervalo.
- `condensate.csv`: concentración, volumen, etanol y masa capturada.
- `trap_standards.csv`: recuperación independiente del tren.

La masa de condensado debe satisfacer, dentro de 1 %, `concentration_ug_l × condensate_volume_l × dilution_factor`.

## Decisión predefinida

El modelo se ejecuta una vez con los parámetros congelados. Deben pasar simultáneamente:

- NRMSE de vino ≤0,55 para cada compuesto;
- NRMSE de condensado ≤1,15 para octanoato y ≤1,05 para isoamilo;
- NRMSE de masa gaseosa integrada ≤0,60;
- sesgo absoluto/media observada ≤0,25;
- balance producción–líquido–gas dentro de 20 %;
- al menos dos de tres reactores mejoran frente a captura constante;
- todos los gates de calidad del contrato.

Si un criterio falla, el resultado es `FAIL`. Los parámetros no se reajustan con esta campaña; cualquier modelo nuevo debe declararse como una fase de desarrollo distinta.
