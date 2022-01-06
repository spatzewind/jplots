import jplots.*;

JPlot plt;
PImage img;

void setup() {
  size(512,512, P2D);
  smooth();
  
  //create data
  float[] x = new float[101];
  float[] y = new float[101];
  float[] z = new float[101];
  float[] r = new float[101];
  for(int i=0; i<101; i++) {
  	x[i] = i/10.0;
  	y[i] = sin(x[i]);
  	z[i] = 1+cos(2*x[i]);
  	r[i] = random(-0.5,1.5);
  }
  PFont font = createFont("",60);
  
  //create plot
  plt = new JPlot(this);
  plt.debug(true);
  plt.figure(3.0,3.0);
  plt.plot(x,y,color(0,0,0),3.0,"-","label","black line");
  plt.plot(x,z,color(255,0,0),3.0,",","label","red dashes");
  plt.scatter(x,r,color(0,127,255),3.0,"[]","label","blue boxes");
  plt.setFont(font);
  plt.setYRange(-1.0,3.5);
  plt.legend();
  img = plt.show();
}

void draw() {
  background(255);
  image(img,0,0,width,height);
}
