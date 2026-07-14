# Ecuaciones y supuestos tecnicos - A02

Actividad: OE2-2.13 - Desarrollar y calibrar modelos de balance de masa para cuantificar perdidas aromaticas

## Modelo base de fermentacion

El estado dinamico base utilizado para calibracion y diseno experimental es:

`X, Xd, N, G, F, E, Gly`

con:

- `X`: biomasa viable.
- `Xd`: biomasa muerta.
- `N`: nitrogeno asimilable.
- `G`: glucosa.
- `F`: fructosa.
- `E`: etanol.
- `Gly`: glicerol.

Las ecuaciones estan implementadas en `fermentation_model/run_new_must_glycerol_estimability_doe.py`:

- `X_balance`: linea 1266.
- `N_balance`: linea 1276.
- `G_balance`: linea 1281.
- `F_balance`: linea 1288.
- `E_balance`: linea 1295.
- `Gly_balance`: linea 1300.

Forma resumida:

```text
dX/dt   = (mu - kd) X + X_input
dXd/dt  = kd X
dN/dt   = -qN phi_N X + N_input
dG/dt   = -(qXG phi_N + qEG phi_G + m_G) X + G_input
dF/dt   = -(qXF phi_N + qEF phi_F + m_F) X + F_input
dE/dt   = (betaG phi_G + betaF phi_F) X + E_input
dGly/dt = (gammaG0 phi_G + gammaF0 phi_F) X
```

Los factores `phi_G` y `phi_F` incorporan temperatura, saturacion de sustrato, interaccion glucosa/fructosa e inhibicion por etanol.

## Capa de perdidas aromaticas

La capa aromatica esta implementada en `fermentation_model/run_aroma_campaign_doe.py`, desde `add_aroma_layer` en linea 385.

Especies incluidas:

- `ethyl_acetate`
- `isoamyl_acetate`
- `ethyl_octanoate`

Forma resumida:

```text
Q_CO2(t)        = conversion_CO2 * CO2_rate(t)
K_lg(i,t)       = K20_i exp(aT_i (T-20) + aE_i (E-50) + aS_i (G+F-100))
r_syn(i,t)      = (k_growth_i phi_growth + k_stationary_i (1-phi_growth)) sugar_uptake(t)
r_loss(i,t)     = alpha_i K_lg(i,t) Q_CO2(t) A_liq(i,t)
dA_liq(i,t)/dt  = r_syn(i,t) - r_loss(i,t)
dA_loss(i,t)/dt = r_loss(i,t)
A_cond(i)       = eta_i A_loss(i,t_final)
```

Referencias de codigo:

- `aroma_loss_rate`: linea 468.
- `aroma_liquid_balance`: linea 475.
- `aroma_loss_balance`: linea 479.
- `aroma_condensate_measurement`: linea 489.

## Particion gas-liquido

El coeficiente `K_lg` proviene de una regresion log-lineal ajustada sobre calculos UNIFAC de dilucion infinita. La funcion base `unifac_partition_K` esta en `fermentation_model/aroma_partition_unifac.py`, linea 223.

El modo operacional usa:

```text
log(K_lg) = log(K20) + temp_slope (T - 20) + ethanol_slope (E - 50) + sugar_slope (G + F - 100)
```

## Supuestos tecnicos principales

- La fructosa se aproxima como glucosa-equivalente en la composicion liquida UNIFAC cuando no existe asignacion directa disponible.
- Las perdidas aromaticas se calculan por stripping proporcional a `Q_CO2`, `K_lg` y concentracion liquida `A_liq`.
- Las eficiencias de trampa terminal son fijas por especie en esta version del modelo.
- Los parametros de particion se fijan por defecto para evitar confundir equilibrio gas-liquido con eficiencia de trampa.
- El condensado terminal aporta cierre de masa acumulado, no dinamica temporal de gas.
- La calibracion base de fermentacion usa datos VL3/natural/sintetico; la capa aromatica queda documentada como extension mecanistica con DOE/FIM y cierre de masa simulado.
