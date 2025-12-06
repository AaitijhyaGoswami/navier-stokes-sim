import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.util.Arrays;

/**
 * Simple 2D Navier-Stokes (Stable Fluids) simulation in plain Java + Swing.
 * Single-file demo suitable for pushing to GitHub.
 *
 * Controls:
 *  - Drag mouse with left button to add velocity
 *  - Drag mouse with right button to add density (smoke)
 *
 * Compile:
 *  javac NavierStokes2D.java
 * Run:
 *  java NavierStokes2D
 */
public class NavierStokes2D extends JPanel implements ActionListener, MouseListener, MouseMotionListener {
    // Simulation parameters
    static final int N = 128; // grid size (NxN)
    static final int SCALE = 4; // pixel scale for display
    static final float DIFF = 0.00001f; // diffusion rate
    static final float VISC = 0.00005f; // viscosity
    static final float DT = 0.5f; // timestep

    // Flattened arrays, index i + (N+2)*j
    final int size = (N + 2) * (N + 2);
    float[] s = new float[size], density = new float[size];
    float[] Vx = new float[size], Vy = new float[size];
    float[] Vx0 = new float[size], Vy0 = new float[size];

    // Rendering
    BufferedImage img;
    Timer timer;

    // Mouse
    int mx = -1, my = -1;
    boolean leftDown = false, rightDown = false;

    public NavierStokes2D() {
        setPreferredSize(new Dimension((N + 2) * SCALE, (N + 2) * SCALE));
        img = new BufferedImage((N + 2) * SCALE, (N + 2) * SCALE, BufferedImage.TYPE_INT_ARGB);
        addMouseListener(this);
        addMouseMotionListener(this);
        timer = new Timer(16, this);
        timer.start();
    }

    // Helper: index in flattened arrays
    private int IX(int i, int j) {
        return i + (N + 2) * j;
    }

    // Add sources
    void addDensity(int x, int y, float amount) {
        int i = Math.max(1, Math.min(N, x));
        int j = Math.max(1, Math.min(N, y));
        density[IX(i, j)] += amount;
    }

    void addVelocity(int x, int y, float amountX, float amountY) {
        int i = Math.max(1, Math.min(N, x));
        int j = Math.max(1, Math.min(N, y));
        Vx[IX(i, j)] += amountX;
        Vy[IX(i, j)] += amountY;
    }

    // Time step
    void step() {
        int iter = 20;
        // velocity step
        addSource(Vx, Vx0);
        addSource(Vy, Vy0);
        Arrays.fill(Vx0, 0f);
        Arrays.fill(Vy0, 0f);

        diffuse(1, Vx0, Vx, VISC, DT, iter);
        diffuse(2, Vy0, Vy, VISC, DT, iter);

        project(Vx0, Vy0, Vx, Vy, iter);

        advect(1, Vx, Vx0, Vx0, Vy0, DT);
        advect(2, Vy, Vy0, Vx0, Vy0, DT);

        project(Vx, Vy, Vx0, Vy0, iter);

        // density step
        addSource(density, s);
        Arrays.fill(s, 0f);

        diffuse(0, s, density, DIFF, DT, iter);
        advect(0, density, s, Vx, Vy, DT);
    }

    void addSource(float[] x, float[] s) {
        for (int i = 0; i < x.length; i++) x[i] += DT * s[i];
    }

    void diffuse(int b, float[] x, float[] x0, float diff, float dt, int iter) {
        float a = dt * diff * N * N;
        linearSolve(b, x, x0, a, 1 + 4 * a, iter);
    }

    void linearSolve(int b, float[] x, float[] x0, float a, float c, int iter) {
        int Np2 = N + 2;
        for (int k = 0; k < iter; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
                }
            }
            setBnd(b, x);
        }
    }

    void project(float[] velocX, float[] velocY, float[] p, float[] div, int iter) {
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                div[IX(i, j)] = -0.5f * (velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] + velocY[IX(i, j + 1)] - velocY[IX(i, j - 1)]) / N;
                p[IX(i, j)] = 0;
            }
        }
        setBnd(0, div);
        setBnd(0, p);

        linearSolve(0, p, div, 1, 4, iter);

        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
                velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
            }
        }
        setBnd(1, velocX);
        setBnd(2, velocY);
    }

    void advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float dt) {
        float dt0 = dt * N;
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                float x = i - dt0 * velocX[IX(i, j)];
                float y = j - dt0 * velocY[IX(i, j)];
                if (x < 0.5f) x = 0.5f;
                if (x > N + 0.5f) x = N + 0.5f;
                int i0 = (int) Math.floor(x);
                int i1 = i0 + 1;
                if (y < 0.5f) y = 0.5f;
                if (y > N + 0.5f) y = N + 0.5f;
                int j0 = (int) Math.floor(y);
                int j1 = j0 + 1;
                float s1 = x - i0;
                float s0 = 1 - s1;
                float t1 = y - j0;
                float t0 = 1 - t1;
                d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                              s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
            }
        }
        setBnd(b, d);
    }

    void setBnd(int b, float[] x) {
        for (int i = 1; i <= N; i++) {
            x[IX(0, i)] = (b == 1) ? -x[IX(1, i)] : x[IX(1, i)];
            x[IX(N + 1, i)] = (b == 1) ? -x[IX(N, i)] : x[IX(N, i)];
            x[IX(i, 0)] = (b == 2) ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N + 1)] = (b == 2) ? -x[IX(i, N)] : x[IX(i, N)];
        }
        x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
        x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
        x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
    }

    // Rendering to BufferedImage
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        // draw density field as grayscale
        for (int j = 0; j < N + 2; j++) {
            for (int i = 0; i < N + 2; i++) {
                float d = density[IX(i, j)];
                int c = (int) (Math.min(1f, d) * 255);
                int rgb = (255 << 24) | (c << 16) | (c << 8) | c;
                // fill SCALExSCALE block
                for (int y = 0; y < SCALE; y++) {
                    for (int x = 0; x < SCALE; x++) {
                        img.setRGB(i * SCALE + x, j * SCALE + y, rgb);
                    }
                }
            }
        }
        g.drawImage(img, 0, 0, null);

        // Optionally draw velocity vectors
        g.setColor(Color.RED);
        int stride = 8;
        for (int j = 1; j <= N; j += stride) {
            for (int i = 1; i <= N; i += stride) {
                int x = i * SCALE;
                int y = j * SCALE;
                float vx = Vx[IX(i, j)];
                float vy = Vy[IX(i, j)];
                int ex = (int) (x + vx * 10);
                int ey = (int) (y + vy * 10);
                g.drawLine(x, y, ex, ey);
            }
        }
    }

    // Main loop
    @Override
    public void actionPerformed(ActionEvent e) {
        // if mouse dragging, add sources
        if (mx >= 0 && my >= 0) {
            int gx = mx / SCALE;
            int gy = my / SCALE;
            if (rightDown) addDensity(gx, gy, 100f * DT);
            if (leftDown) {
                // use mouse movement to create velocity
                // naive small random / directional velocity
                addVelocity(gx, gy, (float) (Math.random() - 0.5) * 50f * DT, (float) (Math.random() - 0.5) * 50f * DT);
            }
        }
        step();
        repaint();
    }

    // Mouse events
    @Override public void mousePressed(MouseEvent e) {
        mx = e.getX(); my = e.getY();
        if (SwingUtilities.isLeftMouseButton(e)) leftDown = true;
        if (SwingUtilities.isRightMouseButton(e)) rightDown = true;
    }
    @Override public void mouseReleased(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e)) leftDown = false;
        if (SwingUtilities.isRightMouseButton(e)) rightDown = false;
    }
    @Override public void mouseMoved(MouseEvent e) { mx = e.getX(); my = e.getY(); }
    @Override public void mouseDragged(MouseEvent e) { mx = e.getX(); my = e.getY(); }
    @Override public void mouseClicked(MouseEvent e) {}
    @Override public void mouseEntered(MouseEvent e) {}
    @Override public void mouseExited(MouseEvent e) {}

    public static void main(String[] args) {
        JFrame frame = new JFrame("2D Navier-Stokes (Stable Fluids) - Java");
        NavierStokes2D panel = new NavierStokes2D();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setContentPane(panel);
        frame.pack();
        frame.setResizable(false);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);

        // Optional: add some initial density/velocity
        for (int j = N/4; j < N/2; j++) {
            for (int i = N/4; i < N/2; i++) {
                panel.density[panel.IX(i, j)] = 50f;
                panel.Vx[panel.IX(i, j)] = 5f;
            }
        }
    }
}
