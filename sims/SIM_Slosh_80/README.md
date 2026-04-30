# Slosh Tank Simulation (OpenFOAM)

## Overview

This case models **free-surface sloshing in a partially filled rectangular tank**, a canonical problem in fluid dynamics with direct applications to spacecraft propellant dynamics, stability, and control.

The simulation uses a **two-phase Volume of Fluid (VOF) formulation** to capture the interface between liquid and gas. A **time-dependent lateral acceleration input** is applied to excite the fluid, producing a controlled slosh response.

This case is intended as:

* A **numerical experiment** for slosh dynamics
* A **validation environment** for reduced-order modeling
* A **data-generation tool** for estimation and control (e.g., EKF, LQR)

---

## Governing Equations

The flow is governed by the **incompressible Navier–Stokes equations**:

### Continuity (mass conservation)

[
\nabla \cdot \mathbf{U} = 0
]

### Momentum equation

[
\frac{\partial (\rho \mathbf{U})}{\partial t}

* \nabla \cdot (\rho \mathbf{U} \mathbf{U})
  = -\nabla p + \nabla \cdot (\mu \nabla \mathbf{U}) + \rho \mathbf{g} + \rho \mathbf{a}(t)
  ]

where:

* (\mathbf{U}) = velocity
* (p) = pressure
* (\rho) = density
* (\mu) = dynamic viscosity
* (\mathbf{g}) = gravitational acceleration
* (\mathbf{a}(t)) = imposed time-dependent acceleration (forcing)

---

## Multiphase Model (VOF)

The interface is tracked using a scalar field:

[
\alpha_{\text{water}} \in [0,1]
]

* (\alpha = 1): liquid
* (\alpha = 0): gas
* (0 < \alpha < 1): interface

The transport equation for the phase fraction is:

[
\frac{\partial \alpha}{\partial t}

* \nabla \cdot (\mathbf{U} \alpha)
* \nabla \cdot \left( \mathbf{U}_c \alpha (1-\alpha) \right) = 0
  ]

where the last term is an **interface compression term** used to maintain a sharp free surface.

---

## Physical Assumptions

* Incompressible flow
* Newtonian fluids
* Laminar regime
* Immiscible phases (air and water)
* No phase change
* Rigid tank walls
* Closed domain (no inflow/outflow)

---

## Problem Configuration

* Geometry: rectangular tank
* Fill fraction: ~80%
* Gravity: lunar-like ((g = -1.62, \text{m/s}^2))
* Excitation: lateral acceleration pulse

---

## Forcing Model

The slosh is induced via a **time-dependent body acceleration**:

[
\mathbf{a}(t) = a_x(t),\hat{i}
]

where (a_x(t)) is defined using a smooth temporal profile.

This produces:

* initial perturbation of the free surface
* excitation of dominant slosh modes
* subsequent free decay

---

## Numerical Method

The simulation is performed using:

```text
interFoam
```

### Key Features

* Finite volume discretization
* VOF interface tracking
* PIMPLE algorithm (PISO + SIMPLE hybrid)
* Second-order spatial schemes

---

## Case Structure

```
SIM_Slosh_80/
├── 0/           # Initial and boundary conditions
├── constant/    # Mesh, fluid properties, forcing
├── system/      # Numerical setup
├── cleanCase.sh # Case reset script
```

---

## Key Files

### 0/

* `alpha.water`: phase fraction
* `U`: velocity field
* `p_rgh`: modified pressure field

---

### constant/

#### polyMesh/

Structured mesh with **non-uniform vertical refinement**:

* fine near free surface
* coarse at bottom

#### phaseProperties

Defines multiphase system.

#### physicalProperties.*

Defines density and viscosity.

#### momentumTransport

Laminar model.

#### fvModels

Defines **time-dependent acceleration forcing**.

---

### system/

#### controlDict

Time integration and output.

#### fvSchemes

Discretization schemes.

#### fvSolution

Linear solvers and coupling.

#### setFieldsDict

Initial condition for liquid region.

---

## Important Numerical Concepts

### Pressure formulation

OpenFOAM solves for:

[
p_{rgh} = p - \rho g h
]

This improves numerical stability in gravity-dominated flows.

---

### Interface compression

Controlled via:

```
cAlpha
```

Balances:

* sharp interface
* numerical stability

---

### Courant number control

[
Co = \frac{U \Delta t}{\Delta x}
]

Maintained below a threshold to ensure stability.

---

## Execution Procedure

### 1. Load OpenFOAM

```
source ../../scripts/load_openfoam.sh
```

---

### 2. Navigate to case

```
cd sims/SIM_Slosh_80
```

---

### 3. Reset case

```
./cleanCase.sh

rm -rf dynamicCode/liquidCOM
```

---

### 4. Generate mesh

```
blockMesh
```

---

### 5. Initialize fluid

```
setFields
```

---

### 6. Decompose domain

```
decomposePar -force -time 0
```

---

### 7. Run simulation

```
mpirun -np 16 interFoam -parallel 

For only progress :
mpirun -np 16 interFoam -parallel 2>&1 \
| tee log.interFoam \
| grep --line-buffered "Progress:" \
| sed -E 's/.*Progress: ([0-9.]+%).*/Progress: \1/'

```

---

### 8. Reconstruct solution

```
reconstructPar

rm -rf data && reconstructPar  -fields '(alpha.water)' && mkdir -p data && find . -maxdepth 1 -type d \( -regex './[0-9.]+' -o -name 'processor*' \) -exec mv {} data/ \;

cd postProcessing/liquidCOM

# merge all files into one
awk 'FNR==1 && NR!=1 {next} {print}' */liquidCOM.dat > liquidCOM_all.dat

postProcess -func writeCellCentres -parallel -time 0
postProcess -func writeCellVolumes -parallel -time 0
```

---

## Visualization

```
paraFoam -builtin
```

### In ParaView

* Select `alpha.water`
* Apply **Contour** at:

```
alpha.water = 0.5
```

This approximates the free surface.

---

## Animation

```
ffmpeg -framerate 30 -i SLOSH.%04d.png \
-vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
-c:v libx264 -pix_fmt yuv420p slosh.mp4
```

---

## Expected Results

The simulation should exhibit:

* smooth free-surface oscillation
* dominant first-mode slosh
* gradual decay after excitation

---

## Relevance to Project

This simulation provides:

* time-resolved free-surface dynamics
* force and motion data
* system response to known excitation

These are used to develop:

* reduced-order slosh models
* parameter estimation algorithms
* control-oriented representations

---

## Summary

```
Mesh → Initialize → Solve → Reconstruct → Visualize → Analyze
```

This case forms the **computational foundation for slosh modeling and control development**.
