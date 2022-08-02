import jplots.*;
import jplots.colour.*;
import jplots.layer.*;
import jplots.maths.*;

JPlot plt;

void setup() {
    size(600,600,P2D);
    
    //prepare data:
    float[] x = new float[10];
    float[] y = new float[10];
    float[][] z = new float[10][10];
    for(int j=0; j<10; j++) {
    	for(int i=0; i<10; i++) {
    		z[j][i] = noise(0.7345*i+129.632, 0.6224*j-8635.3257);
    	}
    	x[j] = j;
    	y[j] = j;
    }
    int num_levels = 10;
    //create plot
    plt = new JPlot(this);
    plt.figure(4d,4d);
    //filled contours
    plt.contourp(x, y, z, num_levels);
    //draw contour lines
    plt.contour(x, y, z, num_levels);
    //add a colourbar
    plt.colourbar();
}

void draw() {
    background(0xffffffff);
    int s = min(width,height);
    image(plt.show(),0,0,s,s);
}
