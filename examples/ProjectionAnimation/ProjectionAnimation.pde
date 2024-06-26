import jplots.*;
import jplots.transform.*;

JPlot plt;
PImage img;

void setup() {
  size(1440,720, P2D);
  smooth();
  //frameRate(1);
  
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
  plt.subplots(4d,2d, 1,2);     //create figure twice big than high with 2 subplots
  JAxis[] ax = plt.ga();        //access JAxis objects
  //set first subplot to ortographic projection with center at 0°N, 0°E
  ax[0].setGeoProjection(new OrthographicJProjection(0,0,0,true));
  ax[0].predefImgShow("etopo1");
  ax[0].setGrid();
  ax[0].coastLines();
  //set stereographic projection with center at 0°N, 0°E and true scale at 45°N
  ax[1].setGeoProjection(new StereographicJProjection(0,0,0,true));
  ax[1].predefImgShow("etopo1");
  ax[1].setGrid();
  ax[1].coastLines();
  //create parameters:
  Object[] params = {"transform", new PlateCarreeJProjection(true)};
  //plot curves on both subplots:
  for(int i=0; i<2; i++) {
    ax[i].plot(longitude,curve1,color(0,0,0),3.0,"-",params);
    ax[i].plot(longitude,curve2,color(255,0,0),3.0,";",params);
    ax[i].scatter(longitude,curve3,color(0,0,127),3.0,"<>",params);
  }
}

void draw() {
  background(255);
  int a = (1*frameCount) % 360;
  plt.ga(0).setGeoProjection(new OrthographicJProjection(a,2*a,30*sin(radians(a)),true));
  plt.ga(0).setGrid();
  plt.ga(1).setGeoProjection(new StereographicJProjection(a,2*a,30*sin(radians(a)),true));
  plt.ga(1).setGrid();
  image(plt.show(),0,0,width,height);
}
