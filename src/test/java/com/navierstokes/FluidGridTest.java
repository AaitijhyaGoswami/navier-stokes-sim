package com.navierstokes;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Unit tests for FluidGrid class.
 */
public class FluidGridTest {
    
    private static final double EPSILON = 1e-6;
    
    @Test
    public void testGridInitialization() {
        FluidGrid grid = new FluidGrid(10, 10, 0.0001, 0.0001, 0.1);
        
        assertEquals(10, grid.getWidth());
        assertEquals(10, grid.getHeight());
        
        // Check all cells are initially zero
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                assertEquals(0.0, grid.getDensity(i, j), EPSILON);
                assertEquals(0.0, grid.getVelocityX(i, j), EPSILON);
                assertEquals(0.0, grid.getVelocityY(i, j), EPSILON);
            }
        }
    }
    
    @Test
    public void testAddDensity() {
        FluidGrid grid = new FluidGrid(10, 10, 0.0001, 0.0001, 0.1);
        
        grid.addDensity(5, 5, 100.0);
        
        assertTrue(grid.getDensity(5, 5) > 0);
        assertEquals(100.0, grid.getDensity(5, 5), EPSILON);
    }
    
    @Test
    public void testAddVelocity() {
        FluidGrid grid = new FluidGrid(10, 10, 0.0001, 0.0001, 0.1);
        
        grid.addVelocity(5, 5, 2.0, 3.0);
        
        assertEquals(2.0, grid.getVelocityX(5, 5), EPSILON);
        assertEquals(3.0, grid.getVelocityY(5, 5), EPSILON);
    }
    
    @Test
    public void testFadeDensity() {
        FluidGrid grid = new FluidGrid(10, 10, 0.0001, 0.0001, 0.1);
        
        grid.addDensity(5, 5, 100.0);
        double initialDensity = grid.getDensity(5, 5);
        
        grid.fadeDensity(0.1);
        
        double fadedDensity = grid.getDensity(5, 5);
        assertTrue(fadedDensity < initialDensity);
        assertEquals(90.0, fadedDensity, EPSILON);
    }
    
    @Test
    public void testSimulationStep() {
        FluidGrid grid = new FluidGrid(20, 20, 0.0001, 0.0001, 0.1);
        
        // Add some initial density and velocity
        grid.addDensity(10, 10, 100.0);
        grid.addVelocity(10, 10, 1.0, 0.0);
        
        // Run one simulation step
        grid.step();
        
        // After simulation, density should have spread
        // We can't predict exact values, but we can verify the simulation runs
        double totalDensity = 0;
        for (int i = 0; i < 20; i++) {
            for (int j = 0; j < 20; j++) {
                totalDensity += grid.getDensity(i, j);
            }
        }
        
        // Total density should be conserved (approximately)
        assertTrue(totalDensity > 0);
    }
    
    @Test
    public void testDensityConservation() {
        FluidGrid grid = new FluidGrid(30, 30, 0.0, 0.0, 0.1);
        
        // Add density
        grid.addDensity(15, 15, 100.0);
        
        // Calculate initial total
        double initialTotal = 0;
        for (int i = 0; i < 30; i++) {
            for (int j = 0; j < 30; j++) {
                initialTotal += grid.getDensity(i, j);
            }
        }
        
        // Run simulation without diffusion
        grid.step();
        
        // Calculate final total
        double finalTotal = 0;
        for (int i = 0; i < 30; i++) {
            for (int j = 0; j < 30; j++) {
                finalTotal += grid.getDensity(i, j);
            }
        }
        
        // With no diffusion, density should be approximately conserved
        assertEquals(initialTotal, finalTotal, 1.0);
    }
    
    @Test
    public void testMultipleSteps() {
        FluidGrid grid = new FluidGrid(20, 20, 0.0001, 0.0001, 0.1);
        
        grid.addDensity(10, 10, 50.0);
        grid.addVelocity(10, 10, 0.5, 0.5);
        
        // Run multiple steps - should not crash or produce NaN values
        for (int i = 0; i < 10; i++) {
            grid.step();
        }
        
        // Verify no NaN values
        for (int i = 0; i < 20; i++) {
            for (int j = 0; j < 20; j++) {
                assertFalse(Double.isNaN(grid.getDensity(i, j)));
                assertFalse(Double.isNaN(grid.getVelocityX(i, j)));
                assertFalse(Double.isNaN(grid.getVelocityY(i, j)));
            }
        }
    }
    
    @Test
    public void testGetterBoundsChecking() {
        FluidGrid grid = new FluidGrid(10, 10, 0.0001, 0.0001, 0.1);
        
        grid.addDensity(5, 5, 100.0);
        grid.addVelocity(5, 5, 2.0, 3.0);
        
        // Test valid coordinates
        assertEquals(100.0, grid.getDensity(5, 5), EPSILON);
        assertEquals(2.0, grid.getVelocityX(5, 5), EPSILON);
        assertEquals(3.0, grid.getVelocityY(5, 5), EPSILON);
        
        // Test out of bounds coordinates - should return 0.0
        assertEquals(0.0, grid.getDensity(-1, 5), EPSILON);
        assertEquals(0.0, grid.getDensity(5, -1), EPSILON);
        assertEquals(0.0, grid.getDensity(10, 5), EPSILON);
        assertEquals(0.0, grid.getDensity(5, 10), EPSILON);
        
        assertEquals(0.0, grid.getVelocityX(-1, 5), EPSILON);
        assertEquals(0.0, grid.getVelocityX(10, 5), EPSILON);
        assertEquals(0.0, grid.getVelocityY(5, -1), EPSILON);
        assertEquals(0.0, grid.getVelocityY(5, 10), EPSILON);
    }
}
