import jplots.*;

JPlot plt;

void setup() {
    size(600,600,P2D);
    
    int n = 256;
    float[] xx = new float[n+1];
    float[] yy = new float[n+1];
    float[][] zz = new float[n+1][n+1];
    for(int j=0; j<=n; j++) {
        yy[j] = j - 0.5f*n;
        for(int i=0; i<=n; i++) {
            if(j==0)
                xx[i] = i - 0.5f*n;
            zz[j][i] = xx[i]*yy[j];
        }
    }
    
    plt = new JPlot(this);
    plt.debug(true);
    plt.figure(2d,2d);
    plt.contour(xx,yy,zz,8,"linewidth",3f);
}

void draw() {
    background(0xffffffff);
    int s = min(width,height);
    image(plt.show(),0,0,s,s);
}
