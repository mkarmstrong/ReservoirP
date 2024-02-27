# Reservoir Pressure in R
Reservoir excess pressure analysis in R

This R script includes functions for calculating reservoir pressure variables from pressure waveforms measured continuously over one cardiac cycle.

Diastolic pressure (dicrotic notch to end) is fit with an exponential with a fixed asymptote of 25 mmHg (per Kottenburg-Assenmacher 2009 & Schipke 2003). The systolic component of the reservoir pressure is estimated by finding continuity with fitted diastolic pressure at the dicrotic notch. Note that these calculations differ slightly from the original calculations from Professors K. Parker & A. Hughes.

List of values returned:

| Variable  | Description (units)                  |
|-----------|--------------------------------------|
| rp_amp    | reservoir pressure amplitude (mmHg)  |
| ep_amp    | excess pressure amplitude (mmHg)     |
| rp_int    | reservoir pressure integral (mmHg*s) |
| ep_int    | excess pressure integral (mmHg*s)    |
| sys_k     | systolic rate constant (s)           |
| dia_k     | diastolic rate constant (s)          |
| crit_p    | critical closing pressure (mmHg)     |
|kratio     | tau ratio (unitless)                 |



![Example 1](rp_1.png)

![Example 2](rp_2.png)
