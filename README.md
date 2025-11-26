# 2D Navier–Stokes Fluid Simulation  
*A Solo Computational Physics Project (Java)*

## Overview
This project is a simplified educational implementation of a 2D incompressible fluid simulator in Java.  
It numerically models how velocity, pressure, viscosity, and external forces interact on a 2D grid.

---

## Physical Parameters Used in the Simulation

### Velocity Field
- **u(x, y, t):** Horizontal and vertical velocity components at each grid cell  
- **Units:** arbitrary (simulation units)

### Pressure Field
- **p(x, y, t):** Scalar pressure value at each grid cell  
- Used to enforce incompressibility

### Density
- **rho:** Fluid density  
- Assumed constant (typically set to 1.0 for simulation)

### Viscosity
- **nu:** Kinematic viscosity  
- Controls how quickly the fluid diffuses and smooths out  
- Typical simulation range: 0.0001 – 0.1

### External Forces
- **F:** External forces applied to the fluid  
- Examples: user dragging mouse, injected momentum, gravity-like effects

### Time Step
- **dt:** Size of each simulation step  
- Smaller dt gives more accurate results  
- Typical values: 0.01 – 0.1

### Grid Resolution
- **N:** Number of cells along each axis (NxN grid)  
- Higher N = more accuracy, slower simulation  
- Typical beginner value: 64 or 128

---

## Numerical Components Used
- **Advection:** Moving velocity/dye along the flow  
- **Diffusion:** Smoothing due to viscosity  
- **Projection:** Adjusting the velocity field to make divergence = 0  
- **Boundary Conditions:** No-slip or reflective walls  
- **Gauss–Seidel iteration:** For solving diffusion and pressure steps  

---

## Requirements
- Java 17+  
- Optional: JavaFX or Swing for animation and visualization  

---

## References
- Jos Stam — “Stable Fluids” (1999)  
- Bridson — *Fluid Simulation for Computer Graphics*  

---

## License
MIT License

