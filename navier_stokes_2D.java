import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.*;
import java.io.File;

public class NavierStokes2D {

    // Simulation resolution
    static final int N = 128;            // Grid size
    static final int ITER = 10;          // Solver iterations
    static final int FRAMES = 300;       // Number of frames to save
    static final double DIFF = 0.0001;   // Diffusion
    static final double VISC = 0.0001;   // Viscosity
    static final double DT = 0.1;        // Timestep

    // Fluid fields
    static double[][] u = new double[N][N];     // velocity x
    static double[][] v = new double[N][N];     // velocities y
    static double[][] uPrev = new double[N][N];
    static double[][] vPrev = new double[N][N];

    static double[][] dens = new double[N][N];
    static double[][] densPrev = new double[N][N];

    public static void main(String[] args) throws Exception {

        // Make sure frames folder exists
        new File("frames").mkdirs();

        for (int frame = 1; frame <= FRAMES; frame++) {

            // --- Add some density + velocity in the center ---
            int cx = N / 2;
            int cy = N / 2;

            for (int i = -5; i <= 5; i++) {
                for (int j = -5; j <= 5; j++) {
                    densPrev[cx + i][cy + j] = 50.0;
                    uPrev[cx + i][cy + j] = 20 * Math.sin(frame * 0.1);
                    vPrev[cx + i][cy + j] = 20 * Math.cos(frame * 0.1);
                }
            }

            step();

            saveFrame(frame);

            System.out.println("Saved frame " + frame);
        }

        System.out.println("Done.");
    }

    // ----------------------------
    // Main Navier–Stokes solver step
    // ----------------------------
    static void step() {
        diffuse(1, uPrev, u, VISC, DT);
        diffuse(2, vPrev, v, VISC, DT);

        project(uPrev, vPrev, u, v);

        advect(1, u, uPrev, uPrev, vPrev, DT);
        advect(2, v, vPrev, uPrev, vPrev, DT);

        project(u, v, uPrev, vPrev);

        diffuse(0, densPrev, dens, DIFF, DT);
        advect(0, dens, densPrev, u, v, DT);

        // Clear prev fields for next add forces
        clear(uPrev);
        clear(vPrev);
        clear(densPrev);
    }

    // ----------------------------
    // Util functions
    // ----------------------------
    static void clear(double[][] x) {
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                x[i][j] = 0;
    }

    static void setBnd(int b, double[][] x) {
        for (int i = 1; i < N - 1; i++) {
            x[i][0] = b == 2 ? -x[i][1] : x[i][1];
            x[i][N - 1] = b == 2 ? -x[i][N - 2] : x[i][N - 2];
            x[0][i] = b == 1 ? -x[1][i] : x[1][i];
            x[N - 1][i] = b == 1 ? -x[N - 2][i] : x[N - 2][i];
        }
        x[0][0] = 0.33 * (x[1][0] + x[0][1]);
        x[0][N - 1] = 0.33 * (x[1][N - 1] + x[0][N - 2]);
        x[N - 1][0] = 0.33 * (x[N - 2][0] + x[N - 1][1]);
        x[N - 1][N - 1] = 0.33 * (x[N - 2][N - 1] + x[N - 1][N - 2]);
    }

    static void linSolve(int b, double[][] x, double[][] x0, double a, double c) {
        for (int k = 0; k < ITER; k++) {
            for (int i = 1; i < N - 1; i++) {
                for (int j = 1; j < N - 1; j++) {
                    x[i][j] = (x0[i][j] + a * (
                            x[i - 1][j] + x[i + 1][j] +
                            x[i][j - 1] + x[i][j + 1]
                    )) / c;
                }
            }
            setBnd(b, x);
        }
    }

    static void diffuse(int b, double[][] x, double[][] x0, double diff, double dt) {
        double a = dt * diff * (N - 2) * (N - 2);
        linSolve(b, x, x0, a, 1 + 4 * a);
    }

    static void advect(int b, double[][] d, double[][] d0, double[][] u, double[][] v, double dt) {
        int i0, j0, i1, j1;
        double x, y, s0, t0, s1, t1;

        double dt0 = dt * (N - 2);

        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {

                x = i - dt0 * u[i][j];
                y = j - dt0 * v[i][j];

                if (x < 0.5) x = 0.5;
                if (x > N - 1.5) x = N - 1.5;

                if (y < 0.5) y = 0.5;
                if (y > N - 1.5) y = N - 1.5;

                i0 = (int) x;
                i1 = i0 + 1;
                j0 = (int) y;
                j1 = j0 + 1;

                s1 = x - i0; s0 = 1 - s1;
                t1 = y - j0; t0 = 1 - t1;

                d[i][j] =
                        s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1]) +
                        s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);
            }
        }

        setBnd(b, d);
    }

    static void project(double[][] u, double[][] v, double[][] p, double[][] div) {
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                div[i][j] = -0.5 * (
                        u[i + 1][j] - u[i - 1][j] +
                        v[i][j + 1] - v[i][j - 1]
                ) / N;
                p[i][j] = 0;
            }
        }
        setBnd(0, div);
        setBnd(0, p);

        linSolve(0, p, div, 1, 4);

        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                u[i][j] -= 0.5 * N * (p[i + 1][j] - p[i - 1][j]);
                v[i][j] -= 0.5 * N * (p[i][j + 1] - p[i][j - 1]);
            }
        }

        setBnd(1, u);
        setBnd(2, v);
    }

    // ----------------------------
    // Save grayscale density PNG
    // ----------------------------
    static void saveFrame(int frameNumber) throws Exception {
        BufferedImage img = new BufferedImage(N, N, BufferedImage.TYPE_INT_RGB);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {

                int c = (int) Math.min(255, dens[i][j] * 255.0);
                if (c < 0) c = 0;

                int rgb = (c << 16) | (c << 8) | c;
                img.setRGB(i, j, rgb);
            }
        }

        File out = new File("frames/frame_" + String.format("%04d", frameNumber) + ".png");
        ImageIO.write(img, "png", out);
    }
}
