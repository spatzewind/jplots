import pplots.*;

PPlot plt;
PImage img;

void setup() {
  size(512,512, P2D);
  smooth();
  
  //create data
  float[] x = new float[100];
  float[] y = new float[100];
  for(int i=0; i<100; i++) {
  	x[i] = 0.1*i;
  	y[i] = sin(0.1*i);
  }
  PFont font = createFont("",40);
  
  //create plot
  plt = new PPlot(this);
  plt.debug(true);
  plt.figure();
  plt.plot(x,y);
  plt.setFont(font);
  img = plt.show();
}

void draw() {
  background(255);
  image(img,0,0,width,height);
}