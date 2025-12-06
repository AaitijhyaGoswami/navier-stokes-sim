import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.util.Arrays;

public class NavierStokes2DSmooth extends JPanel implements ActionListener, MouseListener, MouseMotionListener {
    static final int N = 256;       // grid resolution
    static final int SCALE = 3;     // visual scaling
    float DIFF = 0.000005f;
    float VISC = 0.0001f;
    float DT = 0.5f;
    float DENSITY_AMOUNT = 500f;
    float VELOCITY_AMOUNT = 50f;

    final int size = (N+2)*(N+2);
    float[] s = new float[size], density = new float[size];
    float[] Vx = new float[size], Vy = new float[size];
    float[] Vx0 = new float[size], Vy0 = new float[size];

    BufferedImage img;
    Timer timer;

    int mx=-1,my=-1;
    boolean leftDown=false,rightDown=false;

    public NavierStokes2DSmooth() {
        setPreferredSize(new Dimension((N+2)*SCALE + 200, (N+2)*SCALE));
        img = new BufferedImage(N+2, N+2, BufferedImage.TYPE_INT_ARGB);

        addMouseListener(this);
        addMouseMotionListener(this);

        timer = new Timer(16, this);
        timer.start();
    }

    private int IX(int i, int j) { return i + (N+2)*j; }

    void addDensity(int x,int y,float amount){
        int i = Math.max(1, Math.min(N, x));
        int j = Math.max(1, Math.min(N, y));
        density[IX(i,j)] += amount;
    }

    void addVelocity(int x,int y,float amountX,float amountY){
        int i = Math.max(1, Math.min(N, x));
        int j = Math.max(1, Math.min(N, y));
        Vx[IX(i,j)] += amountX;
        Vy[IX(i,j)] += amountY;
    }

    void step(){
        int iter = 20;

        addSource(Vx,Vx0);
        addSource(Vy,Vy0);
        Arrays.fill(Vx0,0f); Arrays.fill(Vy0,0f);

        diffuse(1,Vx0,Vx,VISC,DT,iter);
        diffuse(2,Vy0,Vy,VISC,DT,iter);
        project(Vx0,Vy0,Vx,Vy,iter);
        advect(1,Vx,Vx0,Vx0,Vy0,DT);
        advect(2,Vy,Vy0,Vx0,Vy0,DT);
        project(Vx,Vy,Vx0,Vy0,iter);

        addSource(density,s);
        Arrays.fill(s,0f);
        diffuse(0,s,density,DIFF,DT,iter);
        advect(0,density,s,Vx,Vy,DT);
    }

    void addSource(float[] x,float[] s){for(int i=0;i<x.length;i++) x[i]+=DT*s[i];}
    void diffuse(int b,float[] x,float[] x0,float diff,float dt,int iter){
        float a=dt*diff*N*N;
        linearSolve(b,x,x0,a,1+4*a,iter);
    }
    void linearSolve(int b,float[] x,float[] x0,float a,float c,int iter){
        for(int k=0;k<iter;k++){
            for(int j=1;j<=N;j++)
                for(int i=1;i<=N;i++)
                    x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)] + x[IX(i+1,j)] + x[IX(i,j-1)] + x[IX(i,j+1)]))/c;
            setBnd(b,x);
        }
    }

    void project(float[] velocX,float[] velocY,float[] p,float[] div,int iter){
        for(int j=1;j<=N;j++)
            for(int i=1;i<=N;i++){
                div[IX(i,j)] = -0.5f*(velocX[IX(i+1,j)]-velocX[IX(i-1,j)] + velocY[IX(i,j+1)]-velocY[IX(i,j-1)])/N;
                p[IX(i,j)] = 0;
            }
        setBnd(0,div); setBnd(0,p);
        linearSolve(0,p,div,1,4,iter);
        for(int j=1;j<=N;j++)
            for(int i=1;i<=N;i++){
                velocX[IX(i,j)] -= 0.5f*(p[IX(i+1,j)]-p[IX(i-1,j)])*N;
                velocY[IX(i,j)] -= 0.5f*(p[IX(i,j+1)]-p[IX(i,j-1)])*N;
            }
        setBnd(1,velocX); setBnd(2,velocY);
    }

    void advect(int b,float[] d,float[] d0,float[] velocX,float[] velocY,float dt){
        float dt0 = dt*N;
        for(int j=1;j<=N;j++){
            for(int i=1;i<=N;i++){
                float x=i - dt0*velocX[IX(i,j)];
                float y=j - dt0*velocY[IX(i,j)];
                x=Math.max(0.5f,Math.min(N+0.5f,x));
                y=Math.max(0.5f,Math.min(N+0.5f,y));
                int i0=(int)Math.floor(x), i1=i0+1;
                int j0=(int)Math.floor(y), j1=j0+1;
                float s1=x-i0,s0=1-s1,t1=y-j0,t0=1-t1;
                d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)] + t1*d0[IX(i0,j1)]) + s1*(t0*d0[IX(i1,j0)] + t1*d0[IX(i1,j1)]);
            }
        }
        setBnd(b,d);
    }

    void setBnd(int b,float[] x){
        for(int i=1;i<=N;i++){
            x[IX(0,i)] = (b==1)? -x[IX(1,i)]:x[IX(1,i)];
            x[IX(N+1,i)] = (b==1)? -x[IX(N,i)]:x[IX(N,i)];
            x[IX(i,0)] = (b==2)? -x[IX(i,1)]:x[IX(i,1)];
            x[IX(i,N+1)] = (b==2)? -x[IX(i,N)]:x[IX(i,N)];
        }
        x[IX(0,0)] = 0.5f*(x[IX(1,0)]+x[IX(0,1)]);
        x[IX(0,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0,N)]);
        x[IX(N+1,0)] = 0.5f*(x[IX(N,0)]+x[IX(N+1,1)]);
        x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N)]);
    }

    protected void paintComponent(Graphics g){
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;

        // draw density into 1:1 pixel image
        for(int j=0;j<N+2;j++){
            for(int i=0;i<N+2;i++){
                float d = density[IX(i,j)];
                int c = Math.min(255,(int)d);
                int rgb = (255<<24)|(c<<16)|(c<<8)|c;
                img.setRGB(i,j,rgb);
            }
        }

        // smooth upscale
        g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
        g2.drawImage(img, 0, 0, (N+2)*SCALE, (N+2)*SCALE, null);

        // draw velocity vectors
        g2.setColor(Color.RED);
        int stride = 8;
        for(int j=1;j<=N;j+=stride)
            for(int i=1;i<=N;i+=stride){
                int x=i*SCALE,y=j*SCALE;
                float vx=Vx[IX(i,j)], vy=Vy[IX(i,j)];
                int ex=(int)(x+vx*10), ey=(int)(y+vy*10);
                g2.drawLine(x,y,ex,ey);
            }

        // draw cursor indicator
        if(mx >= 0 && my >= 0){
            if(leftDown){
                // velocity: blue arrow
                float vx = (float)(Math.random()-0.5)*VELOCITY_AMOUNT*DT;
                float vy = (float)(Math.random()-0.5)*VELOCITY_AMOUNT*DT;
                g2.setColor(Color.BLUE);
                g2.drawLine(mx, my, mx + (int)(vx*10), my + (int)(vy*10));
            }
            if(rightDown){
                // density: white semi-transparent circle
                g2.setColor(new Color(255,255,255,100));
                int radius = (int)(DENSITY_AMOUNT*0.01);
                g2.fillOval(mx - radius, my - radius, radius*2, radius*2);
            }
        }
    }

    @Override
    public void actionPerformed(ActionEvent e){
        if(mx>=0 && my>=0){
            int gx=mx/SCALE, gy=my/SCALE;
            if(rightDown) addDensity(gx,gy,DENSITY_AMOUNT*DT);
            if(leftDown) addVelocity(gx,gy,(float)(Math.random()-0.5)*VELOCITY_AMOUNT*DT,(float)(Math.random()-0.5)*VELOCITY_AMOUNT*DT);
        }
        step();
        repaint();
    }

    @Override public void mousePressed(MouseEvent e){ mx=e.getX(); my=e.getY(); if(SwingUtilities.isLeftMouseButton(e)) leftDown=true; if(SwingUtilities.isRightMouseButton(e)) rightDown=true;}
    @Override public void mouseReleased(MouseEvent e){ if(SwingUtilities.isLeftMouseButton(e)) leftDown=false; if(SwingUtilities.isRightMouseButton(e)) rightDown=false;}
    @Override public void mouseMoved(MouseEvent e){ mx=e.getX(); my=e.getY();}
    @Override public void mouseDragged(MouseEvent e){ mx=e.getX(); my=e.getY();}
    @Override public void mouseClicked(MouseEvent e){}
    @Override public void mouseEntered(MouseEvent e){}
    @Override public void mouseExited(MouseEvent e){}

    public static void main(String[] args){
        JFrame frame = new JFrame("2D Navier-Stokes Smooth Demo");
        NavierStokes2DSmooth panel = new NavierStokes2DSmooth();

        JPanel controlPanel = new JPanel();
        controlPanel.setLayout(new GridLayout(0,1));
        frame.setLayout(new BorderLayout());
        frame.add(panel, BorderLayout.CENTER);
        frame.add(controlPanel, BorderLayout.EAST);

        JSlider diffSlider = new JSlider(1,100,5);
        diffSlider.setBorder(BorderFactory.createTitledBorder("Diffusion"));
        diffSlider.addChangeListener(e -> panel.DIFF = diffSlider.getValue()/1e6f);
        controlPanel.add(diffSlider);

        JSlider viscSlider = new JSlider(1,200,100);
        viscSlider.setBorder(BorderFactory.createTitledBorder("Viscosity"));
        viscSlider.addChangeListener(e -> panel.VISC = viscSlider.getValue()/1e6f);
        controlPanel.add(viscSlider);

        JSlider dtSlider = new JSlider(1,100,50);
        dtSlider.setBorder(BorderFactory.createTitledBorder("Time Step"));
        dtSlider.addChangeListener(e -> panel.DT = dtSlider.getValue()/100f);
        controlPanel.add(dtSlider);

        JSlider densSlider = new JSlider(50,1000,500);
        densSlider.setBorder(BorderFactory.createTitledBorder("Density Amount"));
        densSlider.addChangeListener(e -> panel.DENSITY_AMOUNT = densSlider.getValue());
        controlPanel.add(densSlider);

        JSlider velSlider = new JSlider(10,200,50);
        velSlider.setBorder(BorderFactory.createTitledBorder("Velocity Strength"));
        velSlider.addChangeListener(e -> panel.VELOCITY_AMOUNT = velSlider.getValue());
        controlPanel.add(velSlider);

        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setResizable(false);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
}
