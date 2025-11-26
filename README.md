# 2D Navier–Stokes Fluid Simulation  
*A Solo Computational Physics Project (Java)*

## Overview
This project implements a simplified solver for the **2D incompressible Navier–Stokes equations** in Java.  
The goal is to build a minimal but accurate educational simulation of fluid flow, using stable numerical techniques that balance physical correctness with real-time visualization.

---

## Governing Equations

The simulation models an **incompressible Newtonian fluid** using the Navier–Stokes PDE system:

### **Momentum Equation**
\[
\frac{\partial \mathbf{u}}{\partial t}
+ (\mathbf{u} \cdot \nabla)\mathbf{u}
= -\frac{1}{\rho}\nabla p
+ \nu \nabla^2 \mathbf{u}
+ \mathbf{F}
\]

where  
- \(\mathbf{u}(x, y, t)\) = velocity field  
- \(p(x, y, t)\) = pressure field  
- \(\rho\) = density (assumed constant)  
- \(\nu\) = kinematic viscosity  
- \(\mathbf{F}\) = external force field  

### **Incompressibility Constraint**
\[
\nabla \cdot \mathbf{u} = 0
\]

This constraint is enforced via a **projection step**, solving the Poisson equation for pressure:

\[
\nabla^2 p = \rho \, \nabla \cdot \mathbf{u}^{*}
\]

where \(\mathbf{u}^{*}\) is the intermediate velocity before enforcing incompressibility.

---

## Numerical Method
The solver follows the classic **Stable Fluids** approach:

1. **Add Forces**  
   \[
   \mathbf{u} \leftarrow \mathbf{u} + \Delta t \, \mathbf{F}
   \]

2. **Diffuse**  
   Solve  
   \[
   \frac{\partial \mathbf{u}}{\partial t} = \nu \nabla^2 \mathbf{u}
   \]  
   using Gauss–Seidel relaxation.

3. **Advect**  
   Semi-Lagrangian backward particle tracing:  
   \[
   \mathbf{u}(x, y, t + \Delta t)
   = \mathbf{u}(x - \mathbf{u} \Delta t, \; y - \mathbf{u} \Delta t, \; t)
   \]

4. **Project**  
   Enforce  
   \[
   \nabla \cdot \mathbf{u} = 0
   \]  
   by solving the Poisson equation for pressure and subtracting the gradient from velocity.

This produces a stable, visually smooth simulation that remains robust even under large timesteps.

---

## Features (Planned / In Progress)
- 2D incompressible Navier–Stokes solver (finite-difference)
- Semi-Lagrangian advection
- Diffusion & viscosity modeling
- Pressure projection using iterative solvers
- Dye field for visualizing flow
- Velocity arrows & vorticity map visualizations
- Real-time interactive rendering (Java Swing / JavaFX)

---

## Tools & Requirements
- **Java 17+**
- Optional: Gradle or Maven for build automation
- No external simulation libraries  
  *(all numerical methods implemented from scratch for educational value)*

---

## Learning Objectives
- Understand the mathematical structure of Navier–Stokes equations  
- Implement numerical PDE solvers in Java  
- Explore advection, diffusion, and incompressibility constraints  
- Visualize vector fields and fluid behavior  
- Develop modular scientific codebases  

---

## References
- Jos Stam — *“Stable Fluids,”* SIGGRAPH 1999  
- Robert Bridson — *Fluid Simulation for Computer Graphics*  
- Chorin — *A Numerical Method for Solving Incompressible Flow Problems*  

---

## License
MIT License
