# Slosh Tank Simulation (OpenFOAM)

## Overview

This case models **free-surface sloshing in a partially filled rectangular tank**, a canonical problem in fluid dynamics with direct applications to spacecraft propellant dynamics, stability, and control.

The simulation uses a **two-phase Volume of Fluid (VOF) formulation** to capture the interface between liquid and gas. A **time-dependent lateral acceleration input** is applied to excite the fluid, producing a controlled slosh response.

This case is intended as:

* A **numerical experiment** for slosh dynamics
* A **validation environment** for reduced-order modeling
* A **data-generation tool** for estimation and control, such as EKF, MHE, and LQR development

---

## Governing Equations

The flow is governed by the **incompressible Navier--Stokes equations**.

### Continuity Equation

The incompressibility constraint is given by

$$
\nabla \cdot \mathbf{U} = 0
$$

where $\mathbf{U}$ is the velocity field.

### Momentum Equation

The momentum equation is written as

$$
\frac{\partial \left(\rho \mathbf{U}\right)}{\partial t}
+
\nabla \cdot \left(\rho \mathbf{U}\mathbf{U}\right)
=
-\nabla p
+
\nabla \cdot \left(\mu \nabla \mathbf{U}\right)
+
\rho \mathbf{g}
+
\rho \mathbf{a}(t)
$$

where:

* $\mathbf{U}$ is the velocity field.
* $p$ is the pressure.
* $\rho$ is the density.
* $\mu$ is the dynamic viscosity.
* $\mathbf{g}$ is the gravitational acceleration.
* $\mathbf{a}(t)$ is the imposed time-dependent acceleration forcing.

---

## Multiphase Model: Volume of Fluid Method

The liquid-gas interface is tracked using the water volume-fraction field,

$$
\alpha_{\text{water}} \in [0,1]
$$

where:

* $\alpha_{\text{water}} = 1$ corresponds to liquid.
* $\alpha_{\text{water}} = 0$ corresponds to gas.
* $0 < \alpha_{\text{water}} < 1$ corresponds to the liquid-gas interface.

The phase-fraction transport equation is written as

$$
\frac{\partial \alpha}{\partial t}
+
\nabla \cdot \left(\mathbf{U}\alpha\right)
+
\nabla \cdot \left[\mathbf{U}_c \alpha \left(1-\alpha\right)\right]
=
0
$$

where the last term is an **interface-compression term** used to maintain a sharp free surface.

---

## Physical Assumptions

* Incompressible flow
* Newtonian fluids
* Laminar regime
* Immiscible phases, air and water
* No phase change
* Rigid tank walls
* Closed domain with no inflow or outflow

---

## Problem Configuration

* Geometry: rectangular tank
* Fill fraction: approximately $60\%$
* Gravity: lunar-like acceleration, $g = -1.62~\text{m/s}^2$
* Excitation: lateral acceleration pulse

---

## Forcing Model

The slosh motion is induced using a time-dependent body acceleration,

$$
\mathbf{a}(t) = a_x(t)\,\hat{\mathbf{i}}
$$

where $a_x(t)$ is defined using a smooth temporal forcing profile.

This produces:

* Initial perturbation of the free surface
* Excitation of the dominant slosh modes
* Subsequent free-decay response

---

## Numerical Method

The simulation is performed using:

```text
interFoam
```

### Key Features

* Finite volume discretization
* VOF interface tracking
* PIMPLE algorithm, which combines features of PISO and SIMPLE
* Second-order spatial schemes

---

## Case Structure

```text
SIM_Slosh_60/
├── 0/           # Initial and boundary conditions
├── constant/    # Mesh, fluid properties, forcing
├── system/      # Numerical setup
├── cleanCase.sh # Case reset script
```

---

## Key Files

### 0/

* `alpha.water`: phase-fraction field
* `U`: velocity field
* `p_rgh`: modified pressure field

---

### constant/

#### polyMesh/

Structured mesh with **non-uniform vertical refinement**:

* Fine near the free surface
* Coarser near the bottom of the tank

#### phaseProperties

Defines the multiphase system.

#### physicalProperties.*

Defines density and viscosity.

#### momentumTransport

Defines the laminar transport model.

#### fvModels

Defines the **time-dependent acceleration forcing**.

---

### system/

#### controlDict

Defines time integration, run controls, and output settings.

#### fvSchemes

Defines the discretization schemes.

#### fvSolution

Defines the linear solvers and pressure-velocity coupling controls.

#### setFieldsDict

Defines the initial liquid region.

---

## Important Numerical Concepts

### Pressure Formulation

OpenFOAM solves for the modified pressure variable,

$$
p_{rgh} = p - \rho g h
$$

This improves numerical stability in gravity-dominated flows.

---

### Interface Compression

The interface compression strength is controlled using:

```text
cAlpha
```

This parameter balances:

* Interface sharpness
* Numerical stability

---

### Courant Number Control

The Courant number is defined as

$$
Co = \frac{U \Delta t}{\Delta x}
$$

and is maintained below a selected threshold to improve numerical stability.

---

## Execution Procedure

### 1. Load OpenFOAM

```bash
source ../../scripts/load_openfoam.sh
```

---

### 2. Navigate to the Case Directory

```bash
cd sims/SIM_Slosh_60
```

---

### 3. Reset the Case

```bash
./cleanCase.sh

```
---
### 4.Build liquidCOM Logging Library 
```bash
./../../scripts/setup_liquidCOM.sh 
```
---

### 5. Generate the Mesh

```bash
blockMesh
```

---

### 6. Initialize the Fluid Region

```bash
setFields
```

---

### 7. Decompose the Domain

```bash
decomposePar -force -time 0
```

---

### 8. Run the Simulation

Run the simulation in parallel using:

```bash
mpirun -np 16 interFoam -parallel
```
Change the number 16 to match the number of cores you want to use, based on your machine.

Run the simulation in parallel while displaying the solver output in the terminal and saving the full log file:

```bash
mpirun -np 16 interFoam -parallel 2>&1 | tee log.interFoam
```

---

### 9. Reconstruct the Solution

```bash
reconstructPar
```

To reconstruct only the `alpha.water` field and move the time directories into a `data` folder, use:

```bash
rm -rf data && reconstructPar -fields '(alpha.water)' && mkdir -p data && find . -maxdepth 1 -type d \( -regex './[0-9.]+' -o -name 'processor*' \) -exec mv {} data/ \;
```

To merge the liquid center-of-mass output into one file, use:

```bash
cd postProcessing/liquidCOM

awk 'FNR==1 && NR!=1 {next} {print}' */liquidCOM.dat > liquidCOM_all.dat
```

To write cell centers and cell volumes, use:

```bash
postProcess -func writeCellCentres -parallel -time 0
postProcess -func writeCellVolumes -parallel -time 0
```

---

## Visualization

Open the case in ParaView using:

```bash
paraFoam -builtin
```

### In ParaView

* Select `alpha.water`.
* Apply a **Contour** filter.
* Use the contour value:

```text
alpha.water = 0.5
```

This approximates the liquid-gas free surface.

---

## Animation

To create an animation from exported image frames, use:

```bash
ffmpeg -framerate 30 -i SLOSH.%04d.png \
-vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
-c:v libx264 -pix_fmt yuv460p slosh.mp4
```

---

## Expected Results

The simulation should exhibit:

* Smooth free-surface oscillation
* Dominant first-mode slosh behavior
* Gradual decay after the forcing input is removed

---

## Relevance to the Project

This simulation provides:

* Time-resolved free-surface dynamics
* Wall-force data
* Liquid center-of-mass motion
* System response to a known excitation input

These data are used to develop:

* Reduced-order slosh models
* Parameter estimation algorithms
* Control-oriented force representations

---

## Summary

```text
Mesh → Initialize → Solve → Reconstruct → Visualize → Analyze
```

This case forms the **computational foundation for slosh modeling and control development**.