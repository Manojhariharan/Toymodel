2025-04-30

peat_v0 is the first version of a simple SOM decomposition model designed to simulate long-term carbon accumulation in peat soils. This version is deliberately minimal, developed to start the model framework with just the core components: **one layer, one SOM pool, one carbon input, and one decay (respiration) output**.

---

## Model Description

This is a **single-pool, first-order kinetic model** for soil organic matter (SOM), using the standard carbon mass balance:


Where:
- `S` = SOM stock (kg C m⁻²)
- `I` = input flux (kg C m⁻² yr⁻¹)
- `k` = decay constant (/yr)

Integration is performed with a 1-year timestep.

---

## Parameters and Justification

| Parameter      | Value       | Description & Source                                      |
|----------------|-------------|------------------------------------------------------------|
| `input_rate`   | 0.2 kg C/m²/yr | Chosen for testing. Roughly reflects a modest litter input rate for a semi-productive organic soil, e.g., in a partially rewet or degraded fen. This value is conservative; natural Fenland inputs are closer to 500–1000 g C/m²/yr (Packer et al., 2017; Stout, 1971). |
| `k_decay`      | 0.05 /yr     | Represents a moderate SOM turnover rate, corresponding to a half-life of ~14 years. This is faster than passive peat but slower than active microbial pools. Chosen as a simplified representation of combined aerobic + anaerobic decay. |
| `SOM_init`     | 1.0 kg C/m²  | Arbitrary initial condition for testing. As expected, the model converges to a steady-state stock of I/k = 4.0 kg C/m². |
| `nyears`       | 100          | Run length chosen to approach equilibrium (~98% of steady state is reached in 88 years with this k). |

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
