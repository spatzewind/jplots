package jplots.layer;

import jplots.JAxis;
import jplots.JPlot;
import jplots.colour.JColourtable;
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
	String label;
	PImage img;
	JProjection inputProj;
	JColourtable colourtable;
	Object parallelArray;
	
	public JPlotsLayer() {
		lw = 2d;
		lc = 0xff000000;
		ls = "";
		singleLineColour = true;
		lcs = new int[] { 0xff000000 };
		minX = -1d;
		maxX =  1d;
		minY = -1d;
		maxY =  1d;
		angleMode = 1d;
		img = null;
		label = "";
		inputProj = new IdentityJProjection();
		colourtable = JColourtable.pctables.get("default");
		invertAxisX = false;
		invertAxisY = false;
		logAxisX = false;
		logAxisY = false;
		parallelArray = null;
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
	public int getColour() {
		return singleLineColour ? lc : 0; }
	public void setLineWidth(double lwd) {
		lw = lwd; }
	public double getLineWidth() {
		return lw; }
	public void setStyle(String lst) {
		ls = lst; }
	public String getStyle() {
		return ls; }
	public void invert(String which, boolean _invert) {
		boolean ix = false, iy = false;
		if("both".equalsIgnoreCase(which)) { ix=true; iy=true; }
		if("x".equalsIgnoreCase(which)) ix=true;
		if("y".equalsIgnoreCase(which)) iy=true;
		if(ix) invertAxisX = _invert;
		if(iy) invertAxisY = _invert;
	}
	public void addParallelArray(Object parallel) {
		parallelArray = parallel; }
	public void setLabel(String _l) {
		label = _l; }
	public String getLabel() {
		return label; }
	
	public JColourtable getColourtable() { return colourtable; }
	
	public abstract void createRasterImg(JPlot plot, PGraphics g);
	public abstract void createVectorImg(JAxis ax, int layernum, JGroupShape s);
}
