package jplots.layer;

import jplots.JAxis;
import jplots.JPlot;
import jplots.color.JColourtable;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.transform.IdentityJProjection;
import jplots.transform.JProjection;
import processing.core.PGraphics;
import processing.core.PImage;

public abstract class JPlotsLayer {

	//package private variables
	boolean invertAxisX, invertAxisY, logAxisX, logAxisY;
	boolean drawLines, singleLineColour;
	int lc;
	int[] lcs;
	double lw;
	double minX,maxX,minY,maxY;
	double angleMode;
	String ls;
	PImage img;
	JProjection inputProj;
	JColourtable colourtable;
	
	public JPlotsLayer() {
		lw = 2d;
		lc = 0xff000000;
		singleLineColour = true;
		lcs = new int[] { 0xff000000 };
		minX = -1d;
		maxX =  1d;
		minY = -1d;
		maxY =  1d;
		angleMode = 1d;
		img = null;
		inputProj = new IdentityJProjection();
		colourtable = JColourtable.pctables.get("default");
		invertAxisX = false;
		invertAxisY = false;
		logAxisX = false;
		logAxisY = false;
	}
	
	public void setRange(double xmin, double xmax, double ymin, double ymax) {
		minX = xmin; maxX = xmax;
		minY = ymin; maxY = ymax;
	}
	public double[] getRange() {
		return new double[] {minX,maxX,minY,maxY};
	}
	public void setSourceProjection(JProjection in_proj) {
		inputProj = in_proj;
	}
	public void setColourtable(JColourtable ct) {
		colourtable = ct;
	}
	public void angleMode(String angle_type) {
		if(angle_type==null) {
			System.err.println("[ERROR] Cannot read anglemode from null!");
			return;
		}
		if("degrees".equals(angle_type.toLowerCase())) { angleMode = JPlotMath.DEG_TO_RAD; return; };
		if("gon".equals(angle_type.toLowerCase())) { angleMode = Math.PI/200d; return; };
		if("radians".equals(angle_type.toLowerCase())) { angleMode = 1d; return; };
		System.err.println("[ERROR] Cannot interprete anglemode \""+angle_type+"\". Reset to RADIANS");
		angleMode = 1d;
	}
	public void lines(boolean l) {
		drawLines = l;
	}
	public void setLineColour(int _lc) {
		lc = _lc; singleLineColour=true; }
	public void setLineColour(int[] _lcs) {
		lcs = _lcs; singleLineColour=false; }
	public void setLineWidth(double lwd) {
		lw = lwd; }
	public void setLineStyle(String lst) {
		ls = lst; }
	public void invert(String which, boolean _invert) {
		boolean ix = false, iy = false;
		if("both".equalsIgnoreCase(which)) { ix=true; iy=true; }
		if("x".equalsIgnoreCase(which)) ix=true;
		if("y".equalsIgnoreCase(which)) iy=true;
		if(ix) invertAxisX = _invert;
		if(iy) invertAxisY = _invert;
	}
	
	public abstract void createRasterImg(JPlot plot, PGraphics g);
	public abstract void createVectorImg(JAxis ax, int layernum, JGroupShape s);
}
