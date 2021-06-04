import jplots.*;
import jplots.transform.*;

JPlot plt;
PImage img;

void setup() {
  size(1440,720, P2D);
  smooth();
  
  //create data
  float[] longitude = new float[361];
  float[] curve1 = new float[361];
  float[] curve2 = new float[361];
  float[] curve3 = new float[361];
  for(int i=0; i<361; i++) {
  	longitude[i] = i-180;
  	curve1[i] = 60+20*sin(radians(longitude[i]));
  	curve2[i] = 15*cos(radians(2*longitude[i]));
  	curve3[i] = random(-70,-20);
  }
  
  //create plot
  plt = new JPlot(this);
  plt.debug(true);
  plt.subplots(4d,2d, 1,2);     //create figure twice big than high with 2 subplots
  JAxis[] ax = plt.ga();        //access JAxis objects
  //set first subplot to ortographic projection with center at 0°N, 0°E
  ax[0].setGeoProjection(new OrthographicJProjection(0,0,0,true));
  ax[1].setGrid();
  //set stereographic projection with center at 0°N, 0°E and true scale at 45°N
  ax[1].setGeoProjection(new StereographicJProjection(0,0,0,true));
  ax[1].setGrid();
  //create parameters:
  Object[] params = {"transform", new PlateCarreeJProjection(true)};
  //plot curves on both subplots:
  for(int i=0; i<2; i++) {
    ax[i].plot(longitude,curve1,color(0,0,0),3.0,"-",params);
    ax[i].plot(longitude,curve2,color(255,0,0),3.0,",",params);
    ax[i].scatter(longitude,curve3,color(0,127,255),3.0,"<>",params);
  }
  img = plt.show();
}

void draw() {
  background(255);
  image(img,0,0,width,height);
}