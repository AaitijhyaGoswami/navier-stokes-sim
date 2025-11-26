# 2D Navier–Stokes Fluid Simulation  
*A Solo Computational Physics Project (Java)*

## Overview
This project implements a simplified 2D Navier–Stokes solver in Java to model incompressible fluid flow.  
The simulation numerically integrates the governing equations of fluid dynamics using a grid-based finite-difference approach, enabling interactive visualization of velocity fields, vorticity, and advection–diffusion behavior.

The goal is to build an educational, modular, and extensible fluid simulation engine from scratch, inspired by classic real-time solvers used in graphics and physics research.

---

## Governing Equations
The simulation numerically solves the incompressible Navier–Stokes equations:

\[
\frac{\partial \mathbf{u}}{\partial t}
+ (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{F}
\]

\[
\nabla \cdot \mathbf{u} = 0
\]

Where  
- **u** — velocity field  
- **p** — pressure  
- **ν** — kinematic viscosity  
- **F** — external forces  

---

## Features (Planned / In Progress)
- **Stable fluid solver** based on Stam’s semi-Lagrangian method  
- **Pressure projection** using iterative solvers  
- **Advection & diffusion** using finite-difference approximations  
- **Velocity + pressure grids** on a uniform 2D lattice  
- **Interactive UI** for injecting dye and momentum  
- **Visualizations**:
  - velocity vector field  
  - scalar dye field  
  - vorticity heatmap  

---

## Implementation Details

### Numerical Techniques
- Semi-Lagrangian advection  
- Gauss–Seidel relaxation for diffusion & pressure Poisson equation  
- Projection step to enforce incompressibility  
- Explicit time integration  

### Architecture
