**Date:** 2025-05-12  Update: NEE diagnostics and conservation check to 3-layer SOM model

This version of model introduces a new layer for explicit layer-specific modeling for SOM decomposition and NEE diagnostics and conservation checks

**Key point**
- Implemented 3-layer structure: SOM and litter now dimensioned by soil layers (`nlayers = 3`).
- Mass conservation checks at each timestep
- Consistent with variable timestep lengths

---

**Date:** 2025-05-09  Update: NEE diagnostics and conservation check to 2-layer SOM model

This version of model introduces explicit layer-specific modeling for SOM decomposition and NEE diagnostics and conservation checks

**Key point**
- Net Ecosystem Exchange (NEE) diagnostic for comparison with eddy flux data
- Mass conservation checks at each timestep
- Consistent with variable timestep lengths

## Outputs

- `Diagnostics.csv`: Yearly output including SOM per layer, input, respiration, and NEE
- Console summary of total SOM, cumulative NEE, and conservation check (`SOM + NEE ≈ 0`)

 Total SOM (kg C/m2)     :    150.00000
 Cumulative NEE (kg C/m2):   -150.00000
 Total (SOM + NEE)       :      0.00000

---

**Date:** 2025-05-09  Update: Two-layer SOM structure with layer speific calculations and mass-conservation

This version of model introduces explicit layer-specific modeling for SSOM decomposition.

**Key point**
- Implemented 2-layer structure: SOM and litter now dimensioned by soil layers (`nlayers = 2`).
- Introduced a layer aware loop to compute litter input, respiration, and SOM update independently for each layer.
- Replaced scalar SOM tracking with array-based structure and used sum(SOM(:)) for total system mass.
- Added `resp_total` to aggregate respiration across layers and included it in the mass conservation calculation.
- Verified timestep-wise mass conservation using `mass_start = sum(SOM(:)) + dt * input_rate` and `mass_end = sum(SOM(:)) + dt * resp_total`.

Results consistency
- SOM stock and respiration increase gradually over time in each layer and asymptotically reach equilibrium.
- Output remains yearly, with results saved to Diagnostics.csv.
- Behavior converges similarly to the original 1-layer model but now with layer-resolved dynamics and mass conservation.

---

Update: Numerical stability of resp for multiple dt and diagnostic output

- Ran toymodel_v0 with multiple timestep values (dt = 1.0, 0.5, 0.25, 0.1, 0.01)
- Fixed output to ensure respiration and input are reported as annual fluxes (not per timestep)
- Saved results as diagnostics

| dt     | Description       |
|--------|-------------------|
| `1.0`  | Annual            |
| `0.5`  | Semi-annual       |
| `0.25` | Quarterly         |
| `0.1`  | Monthly           |
| `0.01` | Near-daily        |

## Results

The result is shown in the graph below. The notable points from this implementation are:
- The **updated model** with parameters (`dt = 1.0/ 0.5/ 0.25/ 0.1/ 0.01`) mirrors that SOM dynamics and respiration are preserved for all reasonable dt, allowing both mass conservation and numerical stability across.

The greatest differences across time steps (dt) were observed in:
- **SOM** at Year 143, with a maximum difference of **0.19177** across dt values (<0.2%).
- **Respired Carbon** at Year 137, with a maximum difference of **0.00135** across dt values (<0.21%).

The attached figure shows for multiple dt where 

- **Left panel**: SOM pool (kg C m⁻²) over time.
- **Right panel**: Annual respiration flux (kg C m⁻² yr⁻¹).

**Multiple dt plot**
![SOM and respiration across various dt](Plots/Plot_v6a.jpg)

**Zoomed at the greatest diff**
![SOM and respiration in their maximum diff across various dt](Plots/Plot_v6b.jpg)
 
---

**Date:** 2025-05-08  Update: Assessing mass conservation and numerical stability 

| dt     | Description       |
|--------|-------------------|
| `1.0`  | Annual            |
| `0.5`  | Semi-annual       |
| `0.25` | Quarterly         |
| `0.1`  | Monthly           |
| `0.01` | Near-daily        |

## Results

The result is shown in the graph below. The notable points from this implementation are:
- The **updated model** with parameters (`dt = 1.0/ 0.5/ 0.25/ 0.1/ 0.01`) mirrors that SOM dynamics are preserved for all reasonable dt, allowing both mass conservation and numerical stability across.
- 
The greatest differences across time steps (dt) were observed in:
- **SOM** at Year 143, with a maximum difference of **0.19177** units across dt values.
- **Respired Carbon** at Year 1745, with a maximum difference of **1.0395** units across dt values.

The attached figure shows for multiple dt where 

- **Left panel with blue line**: SOM pool (kg C m⁻²) over time.
- **Right panel with red line**: Annual respiration flux (kg C m⁻² yr⁻¹).

**Multiple dt plot**
![SOM and respiration across various dt](Plots/Plot_v5.jpg)

**600 years**
![SOM and respiration over first 600 years across various dt](Plots/Plot_v4a.jpg)

**Full period**
![SOM and respiration over full time scale across various dt](Plots/Plot_v4.jpg)

---

**Date:** 2025-05-07  Update: Mass-Conserving Variable Timestep Version

toymodel_v0 is the first version of a simple SOM decomposition model designed to simulate long-term carbon accumulation in peat soils. This version is deliberately minimal, developed to start the model framework with just the core components: **one layer, one SOM pool, one carbon input, and one decay (respiration) output**.

---

## Model Description

This is a **single-pool, first-order kinetic model** for soil organic matter (SOM), using the standard carbon mass balance:

dS/dt = I - k * S

Where:
- `S` = SOM stock (kg C m⁻²)
- `I` = input flux (kg C m⁻² yr⁻¹)
- `k` = decay constant (/yr)

Integration is performed with a 1-year timestep.

---

## Parameters and Units

| Parameter     | Value         | Units             | Description                                               |
|---------------|---------------|-------------------|-----------------------------------------------------------|
| `input_rate`  | 1.05          | kg C/m²/year      | Annual NPP-derived litter input (Packer et al., 2017)     |
| `k_decay`     | 0.007         | /year             | First-order decay constant (derived from peat data)       |
| `dt`          | 1.0           | years             | Timestep length                                           |
| `nsteps`      | 6000          | unitless          | Number of timesteps (= 6000 years / dt)                   |
| `eps`         | 1.0e-8        | kg C/m²           | Tolerance for mass conservation per timestep              |

---

## Mass Conservation

The model now uses the update formulation
SOM = SOM + dt * dSOM
where `dSOM` is calculated as the net SOM rate of change per year (kg C m⁻² yr⁻¹), and `dt` is the timestep fraction.

The mass conservation check is performed after the SOM update, using:
mass_start = SOM_before + litter
mass_end   = SOM_after + resp

If `abs(mass_end - mass_start) > eps`, the model stops with diagnostic output. This ensures no loss/gain of carbon outside defined fluxes.

These values are placed at the top and bottom of the timestep loop, respectively, to clearly bracket the state update and ensure transparency in mass accounting.

Units used in the model:
`dSOM`: kg C m⁻² yr⁻¹ (annual rate of SOM change)
`litter`: kg C m⁻² timestep⁻¹
`resp`: kg C m⁻² timestep⁻¹

---

## Output

The model prints one line of output per year and `dt` can be changed to improve numerical accuracy without changing the output structure.

---

## Version Notes

This version satisfies the following:
- Timestep-wise carbon mass conservation
- Variable timestep control
- Fully removes accumulated total tracking
- Explicit SOM pool updates
- Suitable for testing long-term peat C dynamics over millennia

---

Date: **2025-05-06**

## Parameters and Justification

| Parameter      | Value            | Description & Source                                       |
|----------------|------------------|------------------------------------------------------------|
| `input_rate`   | 1.05 kg C/m²/yr   | Based on observed total NPP for Fenland reedbeds (Packer et al., 2017) |
| `k_decay`      | 0.007 /yr     | Derived numerically to match a carbon stock of 150 kg C/m² over 6000 years |
| `SOM_init`     | 0.0 kg C/m²      | Assumes peat accumulation starts from bare mineral surface |
| `nyears`       | 6000              | Matches estimated Fenland peat development timescale (Waller, 1994) |

---

## Output

The model prints annual values of:
- SOM pool size (kg C m⁻²)
- Input (kg C m⁻² yr⁻¹)
- Respired carbon (kg C m⁻² yr⁻¹)

This output allows for quick verification that the model behaves as expected and conserves mass.

---

## Results

The result is shown in the graph below. The notable points from this implementation are:

- The **original model** (`k = 0.05 yr⁻¹`) converges to an equilibrium of by ~90 years.
- The **updated model** with parameters (`k = 0.007 yr⁻¹`, `input_rate = 1.05 kg C/m²/yr`, `SOM = 0.0 kg C/m²`, `nyears = 6000`) allows equilibrium to be approached within the 6000-year window, without assuming the system is already at steady state.

The attached figure shows 1) Original model and 2) Updated model:

- **Left panel**: SOM pool (kg C m⁻²) over time.
- **Right panel**: Annual respiration flux (kg C m⁻² yr⁻¹).

**Original model**
![SOM over time and respiration](Plots/Plot_v0.jpg)

**Updated model**
![SOM over time and respiration with updated parameters](Plots/Plot_3.jpg)

---

## Decay Constant Derivation (Updated)

The decay constant was derived to satisfy the analytical solution:

S(t) = (I / k) * (1 - exp(-k * t))

Where:
- I = 1.05 kg C m⁻² yr⁻¹  (input rate based on NPP)
- S(t = 6000) = 150 kg C m⁻²  (observed peat C stock)
- S₀ = 0  (initial condition))

Solving numerically gives:
- **k = 0.007 yr⁻¹**

This decay rate provides a realistic simplification for single-pool dynamics, consistent with gradual but ongoing accumulation observed in modern Fenland peatlands.

## Mass Conservation Check (Timestep-wise, Double Precision)

To ensure numerical accuracy and physical realism, the model now includes a **timestep-wise mass conservation test**. 

### Implementation
At each timestep, the model checks that:
Beginning mass = SOM + Input
Ending mass = SOM + Respiration

If the absolute difference between beginning and ending mass exceeds a small tolerance `eps = 1.0e-8`, the model stops and prints diagnostics.

### Floating round Precision
- The check uses **double precision arithmetic** from **float** throughout the model to minimize floating-point rounding errors.

---

Date: **2025-04-30**

## Parameters and Justification

| Parameter      | Value            | Description & Source                                       |
|----------------|------------------|------------------------------------------------------------|
| `input_rate`   | 0.2 kg C/m²/yr   | Chosen for testing. Roughly reflects a modest litter input rate for a semi-productive organic soil, e.g., in a partially rewet or degraded fen. This value is conservative; natural Fenland inputs are closer to 500–1000 g C/m²/yr (Packer et al., 2017; Stout, 1971). |
| `k_decay`      | 0.000167 /yr     | **Updated decay constant** based on field data (see below). Much slower decay. |
| `SOM_init`     | 1.0 kg C/m²      | Arbitrary initial condition for testing. As expected, the model converges to a steady-state stock of I/k = 4.0 kg C/m². |
| `nyears`       | 100              | Run length chosen to approach equilibrium (~98% of steady state is reached in 88 years with this k). |

---

## Output

The model prints annual values of:
- SOM pool size (kg C m⁻²)
- Input (kg C m⁻² yr⁻¹)
- Respired carbon (kg C m⁻² yr⁻¹)

This output allows for quick verification that the model behaves as expected and conserves mass.

---

## Results

The result is shown in the graph below. The notable points from this implementation are:

- The **original model** (`k = 0.05 yr⁻¹`) converges to an equilibrium of by ~90 years.
- The **updated model** (`k = 0.000167 yr⁻¹`) accumulates slowly and remains far from equilibrium after 100 years.

The attached figure shows 1) Original model and 2) Updated `k_decay` model:

- **Left panel**: SOM pool (kg C m⁻²) over time.
- **Right panel**: Annual respiration flux (kg C m⁻² yr⁻¹).

**Original model**
![SOM over time and respiration](Plots/Plot_v0.jpg)

**Updated model**
![SOM over time and respiration with updated k](Plots/Plot_v0_with_updated_k.jpg)

---

## Decay Constant Derivation (Updated)

From field estimates:

- **Peat C stock** ≈ 150 kg C/m²  
- **Long-term C accumulation** ≈ 25 g C/m²/yr
- Using: k = I / S = 0.025 / 150 ≈ 1.67 × 10⁻⁴ yr⁻¹

This value matches literature estimates for long-term decay in catotelm peat (Yu, 2011; Clymo, 1984).

---


Date: **2025-04-30**

## Parameters and Justification

| Parameter      | Value       | Description & Source                                      |
|----------------|-------------|------------------------------------------------------------|
| `input_rate`   | 0.2 kg C/m²/yr | Chosen for testing. Roughly reflects a modest litter input rate for a semi-productive organic soil, e.g., in a partially rewet or degraded fen. This value is conservative; natural Fenland inputs are closer to 500–1000 g C/m²/yr (Packer et al., 2017; Stout, 1971). |
| `k_decay`      | 0.05 /yr     | Represents a moderate SOM turnover rate, corresponding to a half-life of ~14 years. This is faster than passive peat but slower than active microbial pools. Chosen as a simplified representation of combined aerobic + anaerobic decay. |
| `SOM_init`     | 1.0 kg C/m²  | Arbitrary initial condition for testing. As expected, the model converges to a steady-state stock of I/k = 4.0 kg C/m². |
| `nyears`       | 100          | Run length chosen to approach equilibrium (~98% of steady state is reached in 88 years with this k). |

---

## Results

The result is shown in the graph below. The notable points from this implementation are:

- **SOM (Soil Organic Matter)** increases over time and asymptotically approaches a steady-state value of **4.0 kg C m⁻²**.
- **Respiration flux** also increases over time and stabilizes at the same rate as the input (0.2 kg C m⁻² yr⁻¹), confirming the system reaches equilibrium.
- The model reaches approximately **98% of equilibrium** by **year 88**, making 100 years a useful reference for steady-state behavior under these parameters.

The attached figure shows:

- **Left panel**: SOM pool (kg C m⁻²) over time.
- **Right panel**: Annual respiration flux (kg C m⁻² yr⁻¹).

![SOM over time and respiration](Plots/Plot_v0.jpg)

---

## Purpose

This version is a **toy model**, used to:
- Check that the code framework is working as expected.
- Explore basic dynamics of SOM accumulation and loss.
- Serve as a baseline before adding complexity (e.g., multiple pools, depth layers, or environmental controls).

It is not yet parameterized for real-world Fenland peatlands. That will require updated values for input, decay rate, and time horizon based on observational data (see Yu, 2011; Waller, 1994; Peacock et al., 2019).

---

## References

- **Packer, J.G., et al. (2017)**. Biological Flora of the British Isles: *Phragmites australis*. *Journal of Ecology*, 105(4), 1123–1145.
- **Stout, J.D. (1971)**. Aspects of the microbiology and oxidation of Wicken Fen soil. *Soil Biology & Biochemistry*, 3(1), 9–25.
- **Yu, Z. (2011)**. Holocene carbon flux histories of the world’s peatlands. *The Holocene*, 21(5), 761–774.
- **Waller, M.P. (1994)**. *Flandrian Environmental Change in Fenland*. East Anglian Archaeology Monograph 70.
- **Peacock, M., et al. (2019)**. The full carbon balance of a rewetted cropland fen and a conservation-managed fen. *Agric. Ecosyst. Environ.*, 269, 1–12.
- **Clymo, R.S. (1984)**. The limits to peat bog growth. *Phil. Trans. R. Soc. B*, 303(1117), 605–654.
