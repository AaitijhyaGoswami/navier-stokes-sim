package com.navierstokes;

/**
 * Represents a 2D grid for fluid simulation using the Navier-Stokes equations.
 * This implementation uses Jos Stam's stable fluids method.
 */
public class FluidGrid {
    private final int width;
    private final int height;
    private final int size;
    
    // Velocity fields
    private double[] u;  // x-component of velocity
    private double[] v;  // y-component of velocity
    private double[] uPrev;
    private double[] vPrev;
    
    // Density/scalar field
    private double[] density;
    private double[] densityPrev;
    
    // Simulation parameters
    private double viscosity;
    private double diffusion;
    private double dt;  // time step
    
    /**
     * Creates a new fluid grid.
     * @param width Grid width
     * @param height Grid height
     * @param viscosity Fluid viscosity (affects momentum diffusion)
     * @param diffusion Density diffusion rate
     * @param dt Time step for simulation
     */
    public FluidGrid(int width, int height, double viscosity, double diffusion, double dt) {
        this.width = width;
        this.height = height;
        this.size = width * height;
        this.viscosity = viscosity;
        this.diffusion = diffusion;
        this.dt = dt;
        
        // Initialize arrays
        this.u = new double[size];
        this.v = new double[size];
        this.uPrev = new double[size];
        this.vPrev = new double[size];
        this.density = new double[size];
        this.densityPrev = new double[size];
    }
    
    /**
     * Advances the simulation by one time step.
     */
    public void step() {
        // Diffuse velocity
        diffuse(1, uPrev, u, viscosity);
        diffuse(2, vPrev, v, viscosity);
        
        // Project to maintain incompressibility
        project(uPrev, vPrev, u, v);
        
        // Advect velocity
        advect(1, u, uPrev, uPrev, vPrev);
        advect(2, v, vPrev, uPrev, vPrev);
        
        // Project again
        project(u, v, uPrev, vPrev);
        
        // Diffuse and advect density
        diffuse(0, densityPrev, density, diffusion);
        advect(0, density, densityPrev, u, v);
    }
    
    /**
     * Adds density at a specific location.
     */
    public void addDensity(int x, int y, double amount) {
        int index = IX(x, y);
        if (index >= 0 && index < size) {
            density[index] += amount;
        }
    }
    
    /**
     * Adds velocity at a specific location.
     */
    public void addVelocity(int x, int y, double amountX, double amountY) {
        int index = IX(x, y);
        if (index >= 0 && index < size) {
            u[index] += amountX;
            v[index] += amountY;
        }
    }
    
    /**
     * Diffusion step using Gauss-Seidel relaxation.
     */
    private void diffuse(int b, double[] x, double[] x0, double diff) {
        double a = dt * diff * (width - 2) * (height - 2);
        linearSolve(b, x, x0, a, 1 + 6 * a);
    }
    
    /**
     * Projection step to enforce incompressibility (mass conservation).
     */
    private void project(double[] velocX, double[] velocY, double[] p, double[] div) {
        // Calculate divergence
        for (int j = 1; j < height - 1; j++) {
            for (int i = 1; i < width - 1; i++) {
                div[IX(i, j)] = -0.5 * (
                    velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] +
                    velocY[IX(i, j + 1)] - velocY[IX(i, j - 1)]
                ) / width;
                p[IX(i, j)] = 0;
            }
        }
        
        setBoundary(0, div);
        setBoundary(0, p);
        linearSolve(0, p, div, 1, 6);
        
        // Subtract gradient of pressure
        for (int j = 1; j < height - 1; j++) {
            for (int i = 1; i < width - 1; i++) {
                velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * width;
                velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * height;
            }
        }
        
        setBoundary(1, velocX);
        setBoundary(2, velocY);
    }
    
    /**
     * Advection step - moves quantities through velocity field.
     */
    private void advect(int b, double[] d, double[] d0, double[] velocX, double[] velocY) {
        double dtx = dt * (width - 2);
        double dty = dt * (height - 2);
        
        for (int j = 1; j < height - 1; j++) {
            for (int i = 1; i < width - 1; i++) {
                double x = i - dtx * velocX[IX(i, j)];
                double y = j - dty * velocY[IX(i, j)];
                
                if (x < 0.5) x = 0.5;
                if (x > width - 1.5) x = width - 1.5;
                int i0 = (int) x;
                int i1 = i0 + 1;
                
                if (y < 0.5) y = 0.5;
                if (y > height - 1.5) y = height - 1.5;
                int j0 = (int) y;
                int j1 = j0 + 1;
                
                double s1 = x - i0;
                double s0 = 1 - s1;
                double t1 = y - j0;
                double t0 = 1 - t1;
                
                d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                              s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
            }
        }
        
        setBoundary(b, d);
    }
    
    /**
     * Solves linear system using Gauss-Seidel method.
     */
    private void linearSolve(int b, double[] x, double[] x0, double a, double c) {
        double cRecip = 1.0 / c;
        for (int k = 0; k < 20; k++) {
            for (int j = 1; j < height - 1; j++) {
                for (int i = 1; i < width - 1; i++) {
                    x[IX(i, j)] = (x0[IX(i, j)] +
                        a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] +
                             x[IX(i, j + 1)] + x[IX(i, j - 1)])) * cRecip;
                }
            }
            setBoundary(b, x);
        }
    }
    
    /**
     * Sets boundary conditions.
     */
    private void setBoundary(int b, double[] x) {
        for (int i = 1; i < width - 1; i++) {
            x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, height - 1)] = b == 2 ? -x[IX(i, height - 2)] : x[IX(i, height - 2)];
        }
        
        for (int j = 1; j < height - 1; j++) {
            x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
            x[IX(width - 1, j)] = b == 1 ? -x[IX(width - 2, j)] : x[IX(width - 2, j)];
        }
        
        // Corners
        x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, height - 1)] = 0.5 * (x[IX(1, height - 1)] + x[IX(0, height - 2)]);
        x[IX(width - 1, 0)] = 0.5 * (x[IX(width - 2, 0)] + x[IX(width - 1, 1)]);
        x[IX(width - 1, height - 1)] = 0.5 * (x[IX(width - 2, height - 1)] + x[IX(width - 1, height - 2)]);
    }
    
    /**
     * Converts 2D coordinates to 1D array index.
     */
    private int IX(int x, int y) {
        return x + y * width;
    }
    
    // Getters
    public int getWidth() {
        return width;
    }
    
    public int getHeight() {
        return height;
    }
    
    public double getDensity(int x, int y) {
        return density[IX(x, y)];
    }
    
    public double getVelocityX(int x, int y) {
        return u[IX(x, y)];
    }
    
    public double getVelocityY(int x, int y) {
        return v[IX(x, y)];
    }
    
    /**
     * Fades density over time.
     */
    public void fadeDensity(double fadeAmount) {
        for (int i = 0; i < size; i++) {
            density[i] *= (1 - fadeAmount);
        }
    }
}
