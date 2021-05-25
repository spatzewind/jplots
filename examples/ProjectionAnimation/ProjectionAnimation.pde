import com.metzner.enrico.JKriging.helper.FormatHelper;
import pplots.*;
import pplots.transform.*;

PPlot plt;
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
  plt = new PPlot(this);
  //plt.debug(true);
  plt.subplots(4d,2d, 1,2);     //create figure twice big than high with 2 subplots
  PAxis[] ax = plt.ga();        //access PAxis objects
  //set first subplot to ortographic projection with center at 0°N, 0°E
  ax[0].setGeoProjection(new POrthographicProjection(0,0,0,true));
  ax[0].predefImgShow("earth2");
  ax[0].setGrid();
  //set stereographic projection with center at 0°N, 0°E and true scale at 45°N
  ax[1].setGeoProjection(new PStereographicProjection(0,0,0,true));
  ax[1].predefImgShow("jupiter");
  ax[1].setGrid();
  //create parameters:
  Object[] params = {"transform", new PPlateCarreeProjection(true)};
  //plot curves on both subplots:
  for(int i=0; i<2; i++) {
    ax[i].plot(longitude,curve1,color(0,0,0),3.0,"-",params);
    ax[i].plot(longitude,curve2,color(255,0,0),3.0,",",params);
    ax[i].scatter(longitude,curve3,color(0,127,255),3.0,"<>",params);
  }
}

void draw() {
  background(255);
  int a = (1*frameCount) % 360;
  PAxis[] ax = plt.ga();
  ax[0].setGeoProjection(new POrthographicProjection(a,2*a,30*sin(radians(a)),true));
  ax[1].setGeoProjection(new PStereographicProjection(a,2*a,30*sin(radians(a)),true));
  //println("\nAngle: "+a);
  //FormatHelper.printMat(System.out, ((POrthographicProjection) ax[0].getGeoProjection()).getRotmatP2X());
  image(plt.show(),0,0,width,height);
  //println(frameCount);
}
