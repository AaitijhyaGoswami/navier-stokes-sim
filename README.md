# Navier-Stokes Fluid Simulation

A Java-based 2D fluid simulation implementing the Navier-Stokes equations using Jos Stam's stable fluids method.

## Overview

This project simulates incompressible fluid flow in a 2D grid. It models:
- **Velocity field** (u, v components) - how the fluid moves
- **Density/scalar field** - visual representation of fluid particles
- **Viscosity** - fluid's resistance to flow
- **Diffusion** - spreading of density through the fluid

The simulation uses a stable numerical method that ensures the simulation remains stable even with large time steps.

## Physics

The Navier-Stokes equations describe the motion of viscous fluid substances:

```
∂u/∂t + (u·∇)u = -∇p + ν∇²u + f
∇·u = 0
```

Where:
- `u` is the velocity field
- `p` is pressure
- `ν` is kinematic viscosity
- `f` represents external forces
- The second equation enforces incompressibility (mass conservation)

## Features

- **Diffusion step**: Spreads velocity and density based on viscosity and diffusion coefficients
- **Advection step**: Moves quantities through the velocity field
- **Projection step**: Enforces incompressibility by removing divergence from velocity field
- **Boundary conditions**: Proper handling of fluid at domain boundaries
- **Visualization**: ASCII-based visualization of density field
- **File output**: Saves simulation results to file for analysis

## Requirements

- Java 11 or higher
- Maven 3.6 or higher

## Building

Build the project using Maven:

```bash
mvn clean package
```

This will:
1. Compile all source files
2. Run unit tests
3. Create an executable JAR in `target/fluid-sim-1.0-SNAPSHOT.jar`

## Running the Simulation

Run the main simulation:

```bash
java -cp target/fluid-sim-1.0-SNAPSHOT.jar com.navierstokes.FluidSimulation
```

Or using Maven:

```bash
mvn exec:java -Dexec.mainClass="com.navierstokes.FluidSimulation"
```

## Output

The simulation will:
1. Display progress to the console with ASCII visualization of the density field
2. Save detailed results to `fluid_output.txt`
3. Print final statistics including:
   - Total and average density
   - Maximum density values
   - Average and maximum velocity

### Visualization Legend

The ASCII visualization uses these characters to represent density levels:
- ` ` (space) - No density (< 0.1)
- `.` - Very low density (0.1 - 1.0)
- `:` - Low density (1.0 - 5.0)
- `+` - Medium-low density (5.0 - 10.0)
- `*` - Medium density (10.0 - 20.0)
- `#` - High density (20.0 - 40.0)
- `@` - Very high density (> 40.0)

## Testing

Run the unit tests:

```bash
mvn test
```

The test suite includes:
- Grid initialization tests
- Density and velocity manipulation tests
- Simulation step validation
- Conservation law verification
- Multi-step simulation stability tests

## Code Structure

```
src/main/java/com/navierstokes/
├── FluidGrid.java         - Core simulation engine implementing Navier-Stokes solver
└── FluidSimulation.java   - Main application with demo and visualization

src/test/java/com/navierstokes/
└── FluidGridTest.java     - Unit tests for the fluid grid
```

## Configuration

Key simulation parameters in `FluidSimulation.java`:

- `GRID_WIDTH` / `GRID_HEIGHT`: Simulation grid dimensions (default: 64x64)
- `VISCOSITY`: Fluid viscosity coefficient (default: 0.0001)
- `DIFFUSION`: Density diffusion rate (default: 0.0001)
- `DT`: Time step size (default: 0.1)
- `ITERATIONS`: Number of simulation steps (default: 100)

## Algorithm Details

The simulation uses Jos Stam's stable fluids method with these key steps each frame:

1. **Add forces**: Apply external forces and sources
2. **Diffuse velocity**: Spread momentum due to viscosity
3. **Project velocity**: Make velocity field divergence-free
4. **Advect velocity**: Move velocity through itself
5. **Project velocity**: Enforce incompressibility again
6. **Diffuse density**: Spread density through diffusion
7. **Advect density**: Move density through velocity field

The projection step uses a Poisson solver to compute pressure and remove divergence, ensuring the fluid behaves as incompressible.

## References

- Stam, Jos. "Stable Fluids." SIGGRAPH 99 (1999)
- Bridson, Robert. "Fluid Simulation for Computer Graphics." (2015)

## License

See LICENSE file for details.