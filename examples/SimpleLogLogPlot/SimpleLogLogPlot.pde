import jplots.*;

JPlot plt;
PImage img;

void setup() {
  size(800,400, P2D);
  smooth();
  
  //create data
  int count = 220;
  float[] x = new float[count];
  float[] y = new float[count];
  float[] z = new float[count];
  float[] r = new float[count];
  for(int i=0; i<count; i++) {
  	x[i] = 0.5+i/10.0;
  	y[i] = exp(2.1+2.0*sin(15.0*log(x[i])/log(10)));
  	z[i] = 0.7+100.0*(1.0+cos(3.0*log(x[i])/log(10)));
  	r[i] = pow(10,random(0.5,2.5));
  }
  PFont font = createFont("",60);
  
  //create plot
  plt = new JPlot(this);
  plt.debug(true);
  plt.subplots(6.0,3.0, 1,2);
  JAxis[] ax = plt.ga();
  //draw curves on normal plot
  ax[0].plot(x,y);
  ax[0].plot(x,z,color(255,0,0),3.0,",",null);
  ax[0].scatter(x,r,color(0,127,255),3.0,"[]",null);
  ax[0].setFont(font);
  ax[0].setGrid();
  
  //draw curves on log-log plot
  ax[1].plot(x,y);
  ax[1].plot(x,z,color(255,0,0),3.0,",",null);
  ax[1].scatter(x,r,color(0,127,255),3.0,"[]",null);
  ax[1].setLogarithmicAxis('b'); //set both (b) axes to logarithmic
  ax[1].setFont(font);
  ax[1].setGrid();
  
  //create image
  img = plt.show();
}

void draw() {
  background(255);
  image(img,0,0,width,height);
}
