package com.navierstokes;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Main class for running the Navier-Stokes fluid simulation.
 * This demonstrates a 2D fluid simulation with density and velocity fields.
 */
public class FluidSimulation {
    
    private static final int GRID_WIDTH = 64;
    private static final int GRID_HEIGHT = 64;
    private static final double VISCOSITY = 0.0001;
    private static final double DIFFUSION = 0.0001;
    private static final double DT = 0.1;
    private static final int ITERATIONS = 100;
    
    public static void main(String[] args) {
        System.out.println("Navier-Stokes Fluid Simulation");
        System.out.println("==============================");
        System.out.println("Grid size: " + GRID_WIDTH + "x" + GRID_HEIGHT);
        System.out.println("Viscosity: " + VISCOSITY);
        System.out.println("Diffusion: " + DIFFUSION);
        System.out.println("Time step: " + DT);
        System.out.println("Iterations: " + ITERATIONS);
        System.out.println();
        
        // Create fluid grid
        FluidGrid fluid = new FluidGrid(GRID_WIDTH, GRID_HEIGHT, VISCOSITY, DIFFUSION, DT);
        
        // Initialize with some density sources and velocity
        initializeFluid(fluid);
        
        // Run simulation
        System.out.println("Running simulation...");
        for (int i = 0; i < ITERATIONS; i++) {
            // Add continuous input
            if (i % 10 == 0) {
                addSources(fluid);
            }
            
            // Step simulation
            fluid.step();
            
            // Fade density to prevent accumulation
            fluid.fadeDensity(0.01);
            
            // Print progress
            if (i % 20 == 0) {
                System.out.println("Iteration " + i + "/" + ITERATIONS);
                printDensitySnapshot(fluid);
            }
        }
        
        System.out.println("\nSimulation complete!");
        
        // Save final state to file
        try {
            saveToFile(fluid, "fluid_output.txt");
            System.out.println("Results saved to fluid_output.txt");
        } catch (IOException e) {
            System.err.println("Error saving results: " + e.getMessage());
        }
        
        // Print final statistics
        printStatistics(fluid);
    }
    
    /**
     * Initializes the fluid with some initial conditions.
     */
    private static void initializeFluid(FluidGrid fluid) {
        int centerX = GRID_WIDTH / 2;
        int centerY = GRID_HEIGHT / 2;
        
        // Add density source at center
        for (int i = -3; i <= 3; i++) {
            for (int j = -3; j <= 3; j++) {
                fluid.addDensity(centerX + i, centerY + j, 100.0);
            }
        }
        
        // Add some initial velocity to create swirling motion
        for (int i = 0; i < GRID_WIDTH; i++) {
            for (int j = 0; j < GRID_HEIGHT; j++) {
                double dx = i - centerX;
                double dy = j - centerY;
                double dist = Math.sqrt(dx * dx + dy * dy);
                
                if (dist > 0 && dist < 15) {
                    // Create circular velocity field
                    fluid.addVelocity(i, j, -dy * 0.1, dx * 0.1);
                }
            }
        }
    }
    
    /**
     * Adds continuous sources to the simulation.
     */
    private static void addSources(FluidGrid fluid) {
        // Add density at multiple points
        fluid.addDensity(GRID_WIDTH / 4, GRID_HEIGHT / 2, 50.0);
        fluid.addDensity(3 * GRID_WIDTH / 4, GRID_HEIGHT / 2, 50.0);
        
        // Add some turbulence
        fluid.addVelocity(GRID_WIDTH / 4, GRID_HEIGHT / 2, 2.0, 0.0);
        fluid.addVelocity(3 * GRID_WIDTH / 4, GRID_HEIGHT / 2, -2.0, 0.0);
    }
    
    /**
     * Prints a visual snapshot of the density field.
     */
    private static void printDensitySnapshot(FluidGrid fluid) {
        System.out.println("\nDensity field:");
        
        // Sample every 4th cell for compact display
        for (int j = 0; j < GRID_HEIGHT; j += 4) {
            for (int i = 0; i < GRID_WIDTH; i += 4) {
                double density = fluid.getDensity(i, j);
                char c = getDensityChar(density);
                System.out.print(c);
            }
            System.out.println();
        }
        System.out.println();
    }
    
    /**
     * Converts density value to a character for visualization.
     */
    private static char getDensityChar(double density) {
        if (density < 0.1) return ' ';
        if (density < 1.0) return '.';
        if (density < 5.0) return ':';
        if (density < 10.0) return '+';
        if (density < 20.0) return '*';
        if (density < 40.0) return '#';
        return '@';
    }
    
    /**
     * Saves the fluid state to a file.
     */
    private static void saveToFile(FluidGrid fluid, String filename) throws IOException {
        try (PrintWriter writer = new PrintWriter(new FileWriter(filename))) {
            writer.println("# Navier-Stokes Fluid Simulation Output");
            writer.println("# Grid size: " + fluid.getWidth() + "x" + fluid.getHeight());
            writer.println("# Format: x y density velocityX velocityY");
            writer.println();
            
            for (int j = 0; j < fluid.getHeight(); j++) {
                for (int i = 0; i < fluid.getWidth(); i++) {
                    double density = fluid.getDensity(i, j);
                    double vx = fluid.getVelocityX(i, j);
                    double vy = fluid.getVelocityY(i, j);
                    
                    // Only write non-zero values to keep file smaller
                    if (density > 0.01 || Math.abs(vx) > 0.001 || Math.abs(vy) > 0.001) {
                        writer.printf("%d %d %.6f %.6f %.6f%n", i, j, density, vx, vy);
                    }
                }
            }
        }
    }
    
    /**
     * Prints statistics about the fluid state.
     */
    private static void printStatistics(FluidGrid fluid) {
        double totalDensity = 0;
        double maxDensity = 0;
        double totalVelocity = 0;
        double maxVelocity = 0;
        
        for (int j = 0; j < fluid.getHeight(); j++) {
            for (int i = 0; i < fluid.getWidth(); i++) {
                double density = fluid.getDensity(i, j);
                double vx = fluid.getVelocityX(i, j);
                double vy = fluid.getVelocityY(i, j);
                double vel = Math.sqrt(vx * vx + vy * vy);
                
                totalDensity += density;
                maxDensity = Math.max(maxDensity, density);
                totalVelocity += vel;
                maxVelocity = Math.max(maxVelocity, vel);
            }
        }
        
        int gridSize = fluid.getWidth() * fluid.getHeight();
        
        System.out.println("\nFinal Statistics:");
        System.out.println("-----------------");
        System.out.printf("Total density: %.2f%n", totalDensity);
        System.out.printf("Average density: %.4f%n", totalDensity / gridSize);
        System.out.printf("Maximum density: %.2f%n", maxDensity);
        System.out.printf("Average velocity: %.4f%n", totalVelocity / gridSize);
        System.out.printf("Maximum velocity: %.4f%n", maxVelocity);
    }
}
