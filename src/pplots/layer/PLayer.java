package pplots.layer;

import pplots.PAxis;
import pplots.PPlot;
import pplots.PPlotMath;
import pplots.shapes.PGroupShape;
import pplots.transform.PIdentityProjection;
import pplots.transform.PProjection;
import processing.core.PGraphics;
import processing.core.PImage;

public abstract class PLayer {

	//package private variables
	double lw;
	double minX,maxX,minY,maxY;
	double angleMode;
	PImage img;
	PProjection inputProj;
	
	public PLayer() {
		lw = 2d;
		minX = -1d;
		maxX =  1d;
		minY = -1d;
		maxY =  1d;
		angleMode = 1d;
		img = null;
		inputProj = new PIdentityProjection();
	}
	
	public void setRange(double xmin, double xmax, double ymin, double ymax) {
		minX = xmin; maxX = xmax;
		minY = ymin; maxY = ymax;
	}
	public double[] getRange() {
		return new double[] {minX,maxX,minY,maxY};
	}
	public void setSourceProjection(PProjection in_proj) {
		inputProj = in_proj;
	}
	public void angleMode(String angle_type) {
		if(angle_type==null) {
			System.err.println("[ERROR] Cannot read anglemode from null!");
			return;
		}
		if("degrees".equals(angle_type.toLowerCase())) { angleMode = PPlotMath.DEG_TO_RAD; return; };
		if("gon".equals(angle_type.toLowerCase())) { angleMode = Math.PI/200d; return; };
		if("radians".equals(angle_type.toLowerCase())) { angleMode = 1d; return; };
		System.err.println("[ERROR] Cannot interprete anglemode \""+angle_type+"\". Reset to RADIANS");
		angleMode = 1d;
	}

	public abstract void createRasterImg(PPlot plot, PGraphics g);
	public abstract void createVectorImg(PAxis ax, int layernum, PGroupShape s);
}
