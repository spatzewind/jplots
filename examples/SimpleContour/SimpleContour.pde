import jplots.*;
import jplots.colour.*;
import jplots.layer.*;

int n = 257, m = 257;
JPlot plt;

void setup() {
    size(600,600,P2D);
    
    plt = new JPlot(this);
    plt.figure(4d,4d);
    newPlot();
}

void draw() {
    background(0xffffffff);
    int s = min(width,height);
    image(plt.show(),0,0,s,s);
    plt.debug(false);
}

void keyReleased() {
    if(key=='n') newPlot();
    if(key=='r') {
        plt.redraw(true);
        plt.debug(true);
    }
}

void newPlot() {
    noiseSeed(System.currentTimeMillis());
    float[] xx = new float[n+1];
    float[] yy = new float[m+1];
    float[][] zz = new float[m+1][n+1];
    for(int j=0; j<=m; j++) {
        yy[j] = j - 0.5*m;
        for(int i=0; i<=n; i++) {
            if(j==0)
                xx[i] = i + 0.5; //338.84736; // - 0.5*n;
            zz[j][i] = noise(1.2345*xx[i]/n+19.632, 3.6224*yy[j]/m-8635.3257);
        }
    }
    
    if(noise(0.87598f,42.49374f,.98758f)>0.5f) {
        int klen = (int) random(1,sqrt(n)+5);
        for(int k=0; k<klen; k++) {
            int i = (int) (xx.length*10*noise(k*0.1653+k*k*0.038265, 89376.625354/k)) % xx.length;
            int j = (int) (yy.length*10*noise(k/(0.0235*k), k*k*0.07352+8274.375264, 842646526.077345/(k+k*k*0.376254+92.8264))) % yy.length;
            if(noise(0.453686*k+78.2376, 0.263543*k*k-76267.026365)>0.5f) {
                for(int j2=max(0,j-2); j2<min(j+3,yy.length); j2++)
                    for(int i2=max(0,i-2); i2<min(i+3,xx.length); i2++)
                        zz[j2][i2] = Float.NaN;
            } else {
                for(int j2=-4; j2<5; j2++)
                    for(int i2=-4; i2<5; i2++)
                        if(i2*i2+j2*j2<20)
                            zz[max(0,min(yy.length-1,j+j2))][max(0,min(xx.length-1,i+i2))] = Float.NaN;
            }
        }
    }
    plt.clear();
    plt.contourp(xx,yy,zz,10,"lines",true,"linewidth",3f,"linecolor",0xff7f0000);
    //plt.setAsTimeAxis('x', "dayssince1.1.2000 00:00:00", "proleptic", "dd.mm.'yy");
    plt.colourbar();
}
