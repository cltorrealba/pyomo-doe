# Secondary-metabolite modeling notes

## Data-driven interpretation

The review plots should be read as a model-structure diagnostic. A state is a good candidate for direct calibration only when it has repeated time trajectories, coverage over the phase where the state changes, and enough perturbation across media or inputs. A state should remain an auxiliary or regularized output when it is sparsely measured, mostly terminal, or confounded with an unmeasured loss process.

Current evidence supports adding pyruvate as a transient state. Acetaldehyde and acetic acid are also informative, but their dynamics should be coupled to redox and oxygen exposure instead of being fitted as independent Monod products. Liquid aromas should be treated as net liquid outputs unless gas or condensate measurements are available for the same experiment.

PAN and ammonium should not be collapsed blindly into a single nitrogen pool. A minimal upgrade is to keep two external nitrogen states and let the growth/uptake terms use an effective nitrogen availability:

$$N_{eff}=w_{NH4}N_{NH4}+w_{PAN}N_{PAN}$$

with $w_{NH4} \geq w_{PAN}$ or literature-informed priors. This preserves the fast ammonium depletion and slower PAN behavior visible in the data.

## Acetaldehyde, acetate, and pyruvate

The observed trajectories should not be forced into simple Monod-like accumulation. These compounds are intracellular/extracellular exchange nodes and often show transient peaks, lag-phase behavior, and redox-stress dependence.

Recommended first dynamic layer:

$$r_S=r_G+r_F$$

$$\frac{dPyr}{dt}=k_{Pyr,S}r_S+k_{Pyr,O}g_{O_2}(O_2)X-k_{Pyr,drain}PyrX$$

$$\frac{dAcAld}{dt}=k_{Ald,Pyr}PyrX+k_{Ald,S}r_S-k_{Ald,red}AcAldXh_{ana}(O_2)-k_{Ac,Ald}AcAldXg_{O_2}(O_2)$$

$$\frac{dAcetate}{dt}=k_{Ac,Ald}AcAldXg_{O_2}(O_2)+k_{Ac,stress}Xh_E(E)-k_{Ac,assim}AcetateXh_N(N_{eff})$$

For the ODE calibration phase, use parsimonious empirical rates rather than a full biochemical network. Candidate drivers are sugar uptake rate, biomass, ethanol inhibition, nitrogen limitation, and an oxygen exposure factor.

The practical rule is to estimate only the terms supported by the data. For example, if oxygen is not measured, fix or tightly regularize the oxygen-gated terms and estimate the sugar-driven pyruvate production plus one drain term first.

## Net liquid aroma outputs

With liquid-only aroma measurements, the identifiable quantity is the net liquid concentration:

$$\frac{dA_i^L}{dt}=r_{i,prod}(z,t)-k_{i,loss}(T,E,CO_2)A_i^L-k_{i,deg}A_i^L$$

$$y_{i,liq}(t)=A_i^L(t)+\epsilon_i(t)$$

If a final condensate measurement exists, it can constrain the time-integrated loss:

$$A_i^C(t_f)=\int_0^{t_f} k_{i,loss}(T,E,CO_2)A_i^L(t)dt$$

Here `r_prod` includes biological production and `k_loss` lumps volatilization, stripping, adsorption, and gas-liquid transfer. Without condensate/gas measurements, production and loss are not separately identifiable. Use literature or UNIFAC/Henry estimates to fix/regularize `k_loss`, and estimate only a small number of production parameters.

For MPCC use, the robust output is therefore the predicted liquid concentration band, not a uniquely decomposed biological synthesis flux. The decomposition becomes reliable only when liquid, gas/condensate, and CO2-rate information are combined.

## Oxygen

The initial oxygen pulse should be represented, even if not fully controlled. A minimal state is dissolved oxygen plus an oxygen-exposure integral:

$$\frac{dO_2}{dt}=k_La(O_2^*(T,E)-O_2)-q_{O_2}Xf_S(S)f_N(N_{eff})$$

$$\frac{d\Omega}{dt}=O_2$$

$$\mu_{eff}=\mu_{anaer}(S,N,T,E)\,f_{lag}(\Omega)$$

If `k_La` is not controlled, treat it as experiment-specific or fixed from a short oxygen measurement protocol. The oxygen state should mainly gate lag-phase growth, sterol/lipid limitation, acetaldehyde/acetate behavior, and aroma-side pathways. Do not allow unconstrained oxygen parameters to explain late fermentation unless oxygen is measured.

A practical measurement protocol is sparse but early: dissolved oxygen at inoculation, 2 h, 6 h, 12 h, and 24 h for a subset of fermentations. That is enough to estimate or validate an initial oxygen exposure lump, but not enough to identify a detailed aeration/transfer model.