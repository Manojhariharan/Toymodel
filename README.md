# Toy Model

## Overview

This model simulates the vertical redistribution of Soil Organic Matter (SOM) and mineral soil mass across a fixed set of soil layers. It is designed to:
- Accept SOM inputs at the top layer.
- Adjust SOM-to-mineral ratios while conserving total mass.
- Maintain constant interface depths for all layers **except the bottom one**, which grows as peat accumulates.

## Model Features

- Tracks SOM and mineral mass over `nlayers` soil layers.
- Redistributes excess mass to maintain consistent interface depths.
- Bottom interface depth increases with cumulative SOM input.
- Validates **mass conservation** of SOM and mineral soil at each timestep.
- Outputs initial and final interface depths, total masses, and conservation status.

## File Structure

### `toymodel.f90`
Main Fortran source code implementing the peat accumulation model.

### `debug.txt`
Output from the model.

### `README.md`
This documentation file.

## How It Works

1. **Initialization**:
   - Each layer starts with 1.0 g/cm2 SOM and 1.0 g/cm2 mineral soil.
   - SOM density = 1.0 g/cm3, Mineral density = 1.5 g/cm3.
   - 1.0 g/cm2 of SOM is added per timestep to the top layer.

2. **Redistribution Logic**:
   - Each timestep, new SOM is added to layer 1.
   - If the layer exceeds its thickness, SOM and mineral mass are redistributed downward.
   - This ensures upper interface depths are preserved.

3. **Conservation Checks**:
   - Total SOM and mineral mass before and after each timestep are compared.
   - Warnings are issued if mass is not conserved within a tolerance (1e-5).
   - A final summary is printed after the simulation.

## Compilation and Running

### Using `gfortran`:

```bash
gfortran -o peat peat.f90
./peat
