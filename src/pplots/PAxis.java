package pplots;

import java.util.ArrayList;
import java.util.List;

import pplots.gfx.PGroupShape;
import pplots.gfx.PLineShape;
import pplots.gfx.PPlotShape;
import pplots.gfx.PTextShape;
import pplots.layer.PLayer;
import pplots.layer.PScatterLayer;
import pplots.layer.PXYLayer;
import processing.core.PApplet;
import processing.core.PFont;

public class PAxis {

	private PPlot pplot;
	private boolean rangeFixed, xAxOn, yAxOn, xGrdOn, yGrdOn;
	private int px, py, pw, ph;
	private double minX,maxX,minY,maxY;
	private double txtsize;
	private List<PLayer> layers;
	private PFont pfont;
	
	PAxis(PPlot plot, int pos_x, int pos_y, int width, int height) {
		pplot = plot;
		px = pos_x;
		py = pos_y;
		pw = width;
		ph = height;
		layers = new ArrayList<PLayer>();
		if(plot.isDebug())
			System.out.println("[DEBUG] created PAxis-object: x/y="+pos_x+"/"+pos_y+" w/h="+width+"/"+height);
		defaults();
	}
	
	private void defaults() {
		rangeFixed = false;
		xAxOn = true;
		yAxOn = true;
		xGrdOn = false;
		yGrdOn = false;
		minX = -1d;
		maxX =  1d;
		minY = -1d;
		maxY =  1d;
		txtsize = 10d;
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
		this.plot(x, y, 0xff000000, 3f, "-", (Object)null); }
	public void plot(float[] x, float[] y, int colour, float linewidth, String linestyle, Object... params) {
		layers.add(new PXYLayer(x, y, colour, linewidth, linestyle)); updateRange(layers.get(layers.size()-1)); }
	public void plot(double[] x, double[] y) {
		this.plot(x, y, 0xff000000, 3d, "-", (Object)null); }
	public void plot(double[] x, double[] y, int colour, double linewidth, String linestyle, Object... params) {
		layers.add(new PXYLayer(x, y, colour, linewidth, linestyle)); updateRange(layers.get(layers.size()-1)); }

	public void scatter(float[] x, float[] y) {
		this.scatter(x, y, 0xff000000, 1f, "c", (Object)null); }
	public void scatter(float[] x, float[] y, int colour, float iconsize, String symbol, Object... params) {
		layers.add(new PScatterLayer(x, y, colour, iconsize, symbol)); updateRange(layers.get(layers.size()-1)); }
	public void scatter(double[] x, double[] y) {
		this.scatter(x, y, 0xff000000, 1d, "c", (Object)null); }
	public void scatter(double[] x, double[] y, int colour, double iconsize, String symbol, Object... params) {
		layers.add(new PScatterLayer(x, y, colour, iconsize, symbol)); updateRange(layers.get(layers.size()-1)); }
	
	//....
	public void setRange(float xmin, float xmax, float ymin, float ymax) {
		minX = xmin; maxX = xmax; minY = ymin; maxY = ymax; rangeFixed = true; }
	public void setRange(double xmin, double xmax, double ymin, double ymax) {
		minX = xmin; maxX = xmax; minY = ymin; maxY = ymax; rangeFixed = true; }
	public void setAxis(String axis, String which, boolean onoff) {
		boolean setX=false, setY=false;
		boolean setAx=false, setGrd=false;
		if("both".equals(axis.toLowerCase())) { setX=true; setY=true; }
		if("x".equals(axis.toLowerCase())) setX=true;
		if("y".equals(axis.toLowerCase())) setY=true;
		if("both".equals(which.toLowerCase())) { setAx=true; setGrd=true; }
		if("a".equals(which.toLowerCase()) || "axis".equals(which.toLowerCase())) setAx=true;
		if("g".equals(which.toLowerCase()) || "grid".equals(which.toLowerCase())) setGrd=true;
		if(setAx) { if(setX) xAxOn=onoff; if(setY) yAxOn=onoff; }
		if(setGrd) { if(setX) xGrdOn=onoff; if(setY) yGrdOn=onoff; }
	}
	
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
	
	PGroupShape createPlot(PApplet applet, int w, int h) {
		if(pplot.isDebug())
			System.out.println("[DEBUG] PAxis-object: min/max={x:"+minX+"/"+maxX+", y:"+minY+"/"+maxY+
				"} with "+layers.size()+" layers");
		PGroupShape graph = new PGroupShape();
		if(xAxOn || xGrdOn) graph.addChild(createXAxis());
		if(yAxOn || yGrdOn) graph.addChild(createYAxis());
		if(xAxOn || yAxOn) {
			PPlotShape.stroke(0xff000000); PPlotShape.strokeWeight(3f);
			graph.addChild(new PLineShape(px, py, px+pw, py));
			graph.addChild(new PLineShape(px,py+ph,px+pw,py+ph));
			graph.addChild(new PLineShape(px,py,px,py+ph));
			graph.addChild(new PLineShape(px+pw,py,px+pw,py+ph));
		}
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
		if(layers.size()==1) {
			minX = r[0];
			maxX = r[1];
			minY = r[2];
			maxY = r[3];
		} else {
			if(r[0]<minX) minX = r[0];
			if(r[1]>maxX) maxX = r[1];
			if(r[2]<minY) minY = r[2];
			if(r[3]>maxY) maxY = r[3];
		}
	}
	private PGroupShape createXAxis() {
		PGroupShape axisgrid = new PGroupShape();
		double[] ticks = PMath.optimalLinearTicks(minX, maxX);
		double[] tcpos = PMath.dlerp(ticks,minX,maxX,px,px+pw);
		if(pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for(int t=2; t<ticks.length; t++) {
				tickStr += ", "+PApplet.nf((float)ticks[t],0,2);
				posStr  += ", "+PApplet.nf((float)tcpos[t],0,2);
			}
			System.out.println("[DEBUG] PAxis-object: Xtickfactors={p10: "+ticks[0]+", f: "+ticks[1]+"}");
			System.out.println("[DEBUG] PAxis-object: Xtickval={"+tickStr.substring(2)+"}");
			System.out.println("[DEBUG] PAxis-object: Xtickpos={"+posStr.substring(2)+"}");
		}
		if(xGrdOn) {
			PPlotShape.stroke(0xff999999); PPlotShape.strokeWeight(2f);
			for(int t=2; t<ticks.length; t++)
				if(ticks[t]>=Math.min(minX, maxX) && ticks[t]<=Math.max(minX, maxX))
					axisgrid.addChild(new PLineShape((float)tcpos[t],py,(float)tcpos[t],py+ph));
		}
		if(xAxOn) {
			PPlotShape.stroke(0xff000000); PPlotShape.strokeWeight(2f);
			double vf = 1d/(ticks[0]*ticks[1]);
			for(int t=2; t<ticks.length; t++)
				if(ticks[t]>=Math.min(minX, maxX) && ticks[t]<=Math.max(minX, maxX)) {
					axisgrid.addChild(new PLineShape((float)tcpos[t],py+ph,(float)tcpos[t],py+1.1f*ph));
					//axisgrid.addChild(ap.createShape(PShape.TEXT, "H", (float)tcpos[t],py-0.1f*ph,(float)tcpos[t],py));
					axisgrid.addChild(new PTextShape(""+(ticks[t]*vf), (float)tcpos[t], py+1.11f*ph, (float)txtsize, PApplet.CENTER, PApplet.TOP, 0xff000000));
				}
		}
		return axisgrid;
	}
	private PGroupShape createYAxis() {
		PGroupShape axisgrid = new PGroupShape();
		double[] ticks = PMath.optimalLinearTicks(minY, maxY);
		double[] tcpos = PMath.dlerp(ticks,minY,maxY,py+ph,py);
		if(pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for(int t=2; t<ticks.length; t++) {
				tickStr += ", "+PApplet.nf((float)ticks[t],0,2);
				posStr  += ", "+PApplet.nf((float)tcpos[t],0,2);
			}
			System.out.println("[DEBUG] PAxis-object: Ytickfactors={p10: "+ticks[0]+", f: "+ticks[1]+"}");
			System.out.println("[DEBUG] PAxis-object: Ytickval={"+tickStr.substring(2)+"}");
			System.out.println("[DEBUG] PAxis-object: Ytickpos={"+posStr.substring(2)+"}");
		}
		if(yGrdOn) {
			PPlotShape.stroke(0xff999999); PPlotShape.strokeWeight(2f);
			for(int t=2; t<ticks.length; t++)
				if(ticks[t]>=Math.min(minY, maxY) && ticks[t]<=Math.max(minY, maxY))
					axisgrid.addChild(new PLineShape(px,(float)tcpos[t],px+pw,(float)tcpos[t]));
		}
		if(yAxOn) {
			PPlotShape.stroke(0xff000000); PPlotShape.strokeWeight(2f);
			double vf = 1d/(ticks[0]*ticks[1]);
			for(int t=2; t<ticks.length; t++)
				if(ticks[t]>=Math.min(minY, maxY) && ticks[t]<=Math.max(minY, maxY)) {
					axisgrid.addChild(new PLineShape(px-0.1f*pw,(float)tcpos[t],px,(float)tcpos[t]));
					//axisgrid.addChild(ap.createShape(PShape.TEXT, "H", (float)tcpos[t],py-0.1f*ph,(float)tcpos[t],py));
					axisgrid.addChild(new PTextShape(""+(ticks[t]*vf), px-0.11f*pw, (float)tcpos[t], (float)txtsize, PApplet.RIGHT, PApplet.CENTER, 0xff000000));
				}
		}
		return axisgrid;
	}

}
