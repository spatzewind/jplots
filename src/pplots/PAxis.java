package pplots;

import java.util.ArrayList;
import java.util.List;

import pplots.layer.PLayer;
import pplots.layer.PXYLayer;
import processing.core.PApplet;
import processing.core.PFont;
import processing.core.PShape;

public class PAxis {

	private PPlot pplot;
	private boolean rangeFixed;
	private int px, py, pw, ph;
	private double minX,maxX,minY,maxY;
	private List<PLayer> layers;
	private PFont pfont;
	
	PAxis(PPlot plot, int pos_x, int pos_y, int width, int height) {
		pplot = plot;
		px = pos_x;
		py = pos_y;
		pw = width;
		ph = height;
		layers = new ArrayList<PLayer>();
		rangeFixed = false;
		if(plot.isDebug())
			System.out.println("[DEBUG] created PAxis-object: x/y="+pos_x+"/"+pos_y+" w/h="+width+"/"+height);
	}
	
	public PAxis setPositionAndSize(int pos_x, int pos_y, int width, int height) {
		px = pos_x;
		py = pos_y;
		pw = width;
		ph = height;
		if(pplot.isDebug())
			System.out.println("[DEBUG] resize PAxis-object: x/y="+pos_x+"/"+pos_y+" w/h="+width+"/"+height);
		return this;
	}


	//************************************
	//**** PUBLIC ************************
	//************************************

	//....
	public void contour(double[] x, double[] y, double[][] z) {
		this.contour(x, y, z, 10, (Object[])null); }
	public void contour(double[] x, double[] y, double[][] z, int levels, Object... params) {
		
	}
	
	public void plot(float[] x, float[] y) {
		this.plot(x, y, 0xff000000, 3f, "-", (Object[])null); }
	public void plot(float[] x, float[] y, int colour, float linewidth, String linestyle, Object... params) {
		layers.add(new PXYLayer(x, y, colour, linewidth, linestyle)); updateRange(layers.get(layers.size()-1)); }
	public void plot(double[] x, double[] y) {
		this.plot(x, y, 0xff000000, 3d, "-", (Object[])null); }
	public void plot(double[] x, double[] y, int colour, double linewidth, String linestyle, Object... params) {
		layers.add(new PXYLayer(x, y, colour, linewidth, linestyle)); updateRange(layers.get(layers.size()-1)); }
	
	public void scatter(double[] x, double[] y) {
		this.scatter(x, y, 0xff000000, 1d, "c", (Object[])null); }
	public void scatter(double[] x, double[] y, int colour, double iconsize, String symbol, Object... params) {
		
	}
	
	//....
	public void setRange(float xmin, float xmax, float ymin, float ymax) {
		minX = xmin; maxX = xmax; minY = ymin; maxY = ymax; rangeFixed = true; }
	public void setRange(double xmin, double xmax, double ymin, double ymax) {
		minX = xmin; maxX = xmax; minY = ymin; maxY = ymax; rangeFixed = true; }
	
	//....
	public void setFont(PFont font) {
		pfont = font; }
	


	//************************************
	//**** GETTER ************************
	//************************************
	
	public PPlot getPlot() { return pplot; }
	public int[] getSize() { return new int[] {px,py,pw,ph}; }


	//************************************
	//**** PACKAGE PRIVATE ***************
	//************************************
	
	PShape createPlot(PApplet applet, int w, int h) {
		PShape graph = applet.createShape(PShape.GROUP);
		for(PLayer layer: layers) {
			layer.setRange(minX,maxX,minY,maxY);
			layer.createVectorImg(pplot, graph, px,py, pw,ph);
		}
		return graph;
	}


	//************************************
	//**** PRIVATE ***********************
	//************************************
	
	private void updateRange(PLayer layer) {
		if(rangeFixed) return;
		double[] r = layer.getRange();
		if(r[0]<minX) minX = r[0];
		if(r[1]>maxX) maxX = r[1];
		if(r[2]<minY) minY = r[2];
		if(r[3]>maxY) maxY = r[3];
	}

}
