package jplots;

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import jplots.color.ColourSequenceJColourtable;
import jplots.color.JColourtable;
import jplots.color.LinearSegmentedJColourtable;
import jplots.helper.FileLoader;
import jplots.layer.JContourLayer;
import jplots.layer.JImageLayer;
import jplots.layer.JLineLayer;
import jplots.layer.JPlotsLayer;
import jplots.layer.JScatterLayer;
import jplots.layer.JXYLayer;
import jplots.layer.JShapesLayer;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JTextShape;
import jplots.transform.IdentityJProjection;
import jplots.transform.JProjection;
import processing.core.PApplet;
import processing.core.PFont;
import processing.core.PImage;

public class JAxis {
	
	private static Map<String, PImage> loadedPreDefImgs = new HashMap<>();

	protected JPlot pplot;
	private boolean xRangeFix,yRangeFix, isGeoAxis;
	private boolean xAxOn, yAxOn, xGrdOn, yGrdOn, xTkOn, yTkOn, xAxInv, yAxInv;
	protected int px, py, pw, ph;
	private double minX,maxX,minY,maxY;
	protected double txtsize;
	private String titleX, titleY, titleP;
	private List<JPlotsLayer> layers;
	protected PFont pfont;
	private JProjection projection;
	
	public JAxis(JPlot plot, int pos_x, int pos_y, int width, int height) {
		pplot = plot;
		px = pos_x;
		py = pos_y;
		pw = width;
		ph = height;
		layers = new ArrayList<JPlotsLayer>();
		if(plot.isDebug())
			System.out.println("[DEBUG] created PAxis-object: x/y="+pos_x+"/"+pos_y+" w/h="+width+"/"+height);
		defaults();
	}
	
	private void defaults() {
		xRangeFix = false;
		xAxOn = true;
		yAxOn = true;
		xGrdOn = false;
		yGrdOn = false;
		xTkOn = true;
		yTkOn = true;
		xAxInv = false;
		yAxInv = false;
		minX = -1d;
		maxX =  1d;
		minY = -1d;
		maxY =  1d;
		txtsize = 300d*10d/72d;
		titleX = "";
		titleY = "";
		titleP = "";
		layers.clear();
		projection = new IdentityJProjection();
	}


	//************************************
	//**** PUBLIC ************************
	//************************************

	//....
	public void contour(float[] x, float[] y, float[][] z) {
		this.contour(x, y, z, 10, (Object[])null); }
	public void contour(float[] x, float[] y, float[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, JColourtable.pctables.get("default"), 2.0f, true, false, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contour(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0f, true, false, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contour(float[][] x, float[][] y, float[][] z) {
		this.contour(x, y, z, 10, (Object[])null); }
	public void contour(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, JColourtable.pctables.get("default"), 2.0f, true, false, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contour(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0f, true, false, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contour(double[] x, double[] y, double[][] z) {
		this.contour(x, y, z, 10, (Object[])null); }
	public void contour(double[] x, double[] y, double[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, JColourtable.pctables.get("default"), 2.0d, true, false, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contour(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0d, true, false, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contour(double[][] x, double[][] y, double[][] z) {
		this.contour(x, y, z, 10, (Object[])null); }
	public void contour(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, JColourtable.pctables.get("default"), 2.0d, true, false, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contour(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0d, true, false, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}

	public void contourf(float[] x, float[] y, float[][] z) {
		this.contourf(x, y, z, 10, (Object[])null); }
	public void contourf(float[] x, float[] y, float[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, JColourtable.pctables.get("default"), 2.0f, false, true, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourf(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0f, false, true, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourf(float[][] x, float[][] y, float[][] z) {
		this.contourf(x, y, z, 10, (Object[])null); }
	public void contourf(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, JColourtable.pctables.get("default"), 2.0f, false, true, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourf(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0f, false, true, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourf(double[] x, double[] y, double[][] z) {
		this.contourf(x, y, z, 10, (Object[])null); }
	public void contourf(double[] x, double[] y, double[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, JColourtable.pctables.get("default"), 2.0d, false, true, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourf(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0d, false, true, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourf(double[][] x, double[][] y, double[][] z) {
		this.contourf(x, y, z, 10, (Object[])null); }
	public void contourf(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, JColourtable.pctables.get("default"), 2.0d, false, true, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourf(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0d, false, true, false);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}

	public void contourp(float[] x, float[] y, float[][] z) {
		this.contourf(x, y, z, 10, (Object[])null); }
	public void contourp(float[] x, float[] y, float[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, JColourtable.pctables.get("default"), 2.0f, false, true, true);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourp(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0f, false, true, true);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourp(float[][] x, float[][] y, float[][] z) {
		this.contourf(x, y, z, 10, (Object[])null); }
	public void contourp(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, JColourtable.pctables.get("default"), 2.0f, false, true, true);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourp(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0f, false, true, true);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourp(double[] x, double[] y, double[][] z) {
		this.contourf(x, y, z, 10, (Object[])null); }
	public void contourp(double[] x, double[] y, double[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, JColourtable.pctables.get("default"), 2.0d, false, true, true);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourp(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0d, false, true, true);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourp(double[][] x, double[][] y, double[][] z) {
		this.contourf(x, y, z, 10, (Object[])null); }
	public void contourp(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, JColourtable.pctables.get("default"), 2.0d, false, true, true);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}
	public void contourp(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0d, false, true, true);
		layers.add(cnl); readParams(cnl, params); updateRange(cnl);
	}

	public void plot(float[] x, float[] y) {
		this.plot(x, y, 0xff000000, 3f, "-", (Object)null); }
	public void plot(float[] x, float[] y, int colour, float linewidth, String linestyle, Object... params) {
		JPlotsLayer xyl = new JXYLayer(x, y, colour, linewidth, linestyle); layers.add(xyl);
		readParams(xyl, params); updateRange(xyl); }
	public void plot(double[] x, double[] y) {
		this.plot(x, y, 0xff000000, 3d, "-", (Object)null); }
	public void plot(double[] x, double[] y, int colour, double linewidth, String linestyle, Object... params) {
		JPlotsLayer xyl = new JXYLayer(x, y, colour, linewidth, linestyle); layers.add(xyl);
		readParams(xyl, params); updateRange(xyl); }

	public void scatter(float[] x, float[] y) {
		this.scatter(x, y, 0xff000000, 1f, "c", (Object)null); }
	public void scatter(float[] x, float[] y, int colour, float iconsize, String symbol, Object... params) {
		JPlotsLayer scl = new JScatterLayer(x, y, colour, iconsize, symbol); layers.add(scl);
		readParams(scl, params); updateRange(scl); }
	public void scatter(double[] x, double[] y) {
		this.scatter(x, y, 0xff000000, 1d, "c", (Object)null); }
	public void scatter(double[] x, double[] y, int colour, double iconsize, String symbol, Object... params) {
		JPlotsLayer scl = new JScatterLayer(x, y, colour, iconsize, symbol); layers.add(scl);
		readParams(scl, params); updateRange(scl); }

	public void axhline(float y) {
		axhline(y, 0xff000000, 3f, "-"); }
	public void axhline(float y, int colour, float linewidth, String linestyle) {
		JPlotsLayer xyl = new JLineLayer(y, 'h', colour, linewidth, linestyle); layers.add(xyl);
		updateRange(xyl, "y"); }
	public void axhline(double y) {
		axhline(y, 0xff000000, 3f, "-"); }
	public void axhline(double y, int colour, double linewidth, String linestyle) {
		JPlotsLayer xyl = new JLineLayer(y, 'h', colour, linewidth, linestyle); layers.add(xyl);
		updateRange(xyl, "y"); }
	public void axvline(float x) {
		axvline(x, 0xff000000, 3f, "-"); }
	public void axvline(float x, int colour, float linewidth, String linestyle) {
		JPlotsLayer xyl = new JLineLayer(x, 'v', colour, linewidth, linestyle); layers.add(xyl);
		updateRange(xyl, "x"); }
	public void axvline(double x) {
		axvline(x, 0xff000000, 3f, "-"); }
	public void axvline(double x, int colour, double linewidth, String linestyle) {
		JPlotsLayer xyl = new JLineLayer(x, 'v', colour, linewidth, linestyle); layers.add(xyl);
		updateRange(xyl, "x"); }
	
	public void colourbar() {
		pplot.colourbar(this); }
	public void colourbar(String name) {
		pplot.colourbar(this, name); }
	
	public void coastLines() {
		coastLines(110); }
	public void coastLines(int resolution) {
		JPlotsLayer shl = new JShapesLayer(FileLoader.loadResourceShapeFile("/data/ne_"+resolution+"m_coastline"), "line");
		layers.add(shl); }
	
	/**
	 * predefined images are used as background images in plot, especially with geographical projections
	 * <p>
	 *     earth&nbsp; -- NASA image of earths surface<br>
	 *     earth2 -- ETOPO-like image of earths surface<br>
	 *     etopo1 -- ETOPO image from <a href="https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/image/">NOAA</a><br>
	 *     mercury,venus,mars,jupiter,saturn,uranus,neptune -- NASA images of surfaces of planets in our solar system<br>
	 *     hsbpalette -- rainbow palette of all possible hsb-colors
	 * </p>
	 * 
	 * @param predefined_image one name of above list of predefined images
	 */
	public void predefImgShow(String predefined_image) {
		JPlotsLayer iml;
		if(loadedPreDefImgs.containsKey(predefined_image)) {
			iml = new JImageLayer(loadedPreDefImgs.get(predefined_image));
		} else {
			iml = new JImageLayer(loadPreDefImg(pplot.getApplet(), predefined_image));
		}
		layers.add(0, iml);
	}
	public void imgShow(PImage img) {
		JPlotsLayer iml = new JImageLayer(img); layers.add(iml);
		//updateRange(iml);
	}

	/**
	 * removes all plotting infos
	 * also all configuration will be reseted
	 */
	public void clear() {
		defaults();
	}
	
	//....
	public JAxis setPositionAndSize(int pos_x, int pos_y, int width, int height) {
		px = pos_x;
		py = pos_y;
		pw = width;
		ph = height;
		if(pplot.isDebug())
			System.out.println("[DEBUG] resize PAxis-object: x/y="+pos_x+"/"+pos_y+" w/h="+width+"/"+height);
		return this;
	}
	public JAxis setXRange(float xmin, float xmax) { minX = xmin; maxX = xmax; xRangeFix = true; return this; }
	public JAxis setXRange(double xmin, double xmax) { minX = xmin; maxX = xmax; xRangeFix = true; return this; }
	public JAxis setYRange(float ymin, float ymax) { minY = ymin; maxY = ymax; yRangeFix = true; return this; }
	public JAxis setYRange(double ymin, double ymax) { minY = ymin; maxY = ymax; yRangeFix = true; return this; }
	public JAxis setRange(float xmin, float xmax, float ymin, float ymax) {
		minX = xmin; maxX = xmax; minY = ymin; maxY = ymax; xRangeFix = true; yRangeFix = true; return this; }
	public JAxis setRange(double xmin, double xmax, double ymin, double ymax) {
		minX = xmin; maxX = xmax; minY = ymin; maxY = ymax; xRangeFix = true; yRangeFix = true; return this; }
	public void setAxis(String axis, String which, boolean onoff) {
		boolean setX=false, setY=false;
		boolean setAx=false, setGrd=false, setTck=false;
		if("both".equals(axis.toLowerCase())) { setX=true; setY=true; }
		if("x".equals(axis.toLowerCase())) setX=true;
		if("y".equals(axis.toLowerCase())) setY=true;
		String w = which.toLowerCase();
		if("all".equals(w)) { setAx=true; setGrd=true; }
		if("a".equals(w) || "axis".equals(w)) setAx=true;
		if("g".equals(w) || "grid".equals(w)) setGrd=true;
		if("t".equals(w) || "tick".equals(w) || "ticks".equals(w)) setTck=true;
		if(setAx) { if(setX) xAxOn=onoff; if(setY) yAxOn=onoff; }
		if(setGrd) { if(setX) xGrdOn=onoff; if(setY) yGrdOn=onoff; }
		if(setTck) { if(setX) xTkOn=onoff; if(setY) yTkOn=onoff; }
	}
	public void setGrid() {
		setGrid("both", true); }
	public void setGrid(String axis, boolean onoff) {
		if(isGeoAxis) {
			
		} else {
			boolean setX=false, setY=false;
			if("both".equals(axis.toLowerCase())) { setX=true; setY=true; }
			if("x".equals(axis.toLowerCase())) setX=true;
			if("y".equals(axis.toLowerCase())) setY=true;
			if(setX) xGrdOn=onoff;
			if(setY) yGrdOn=onoff;
		}
	}
	public JAxis setGeoProjection(JProjection proj) {
		projection = proj;
		isGeoAxis = true;
		double[] rr = projection.defaultMapExtend();
		setRange(rr[0], rr[1], rr[2], rr[3]);
		xAxOn = false;  yAxOn = false;
		xGrdOn = false; yGrdOn = false;
		pplot.redraw(true);
		return this;
	}

	//....
	public void setFont(PFont font) {
		pfont = font; }
	public void setTextSize(double ts) {
		txtsize = JPlot.dpi*ts/72d; }
	public void setXTitle(String xtitle) {
		titleX = xtitle; }
	public void setYTitle(String ytitle) {
		titleY = ytitle; }
	public void setTitle(String _title) {
		titleP = _title; }


	//************************************
	//**** GETTER ************************
	//************************************
	
	public JPlot getPlot() { return pplot; }
	public int[] getSize() { return new int[] {px,py,pw,ph}; }
	public double[] getRange() { return new double[] {minX,maxX,minY,maxY}; }
	public boolean isGeoAxis() { return isGeoAxis; }
	public JProjection getGeoProjection() { return projection; }
	public JPlotsLayer getLayer(int layernum) { return layers.get(layernum); }
	public List<JPlotsLayer> getLayers() { return layers; }
	public void addJPlotsLayer(JPlotsLayer layer) { layers.add(layer); }
	public PFont getFont() { return pfont; }


	//************************************
	//**** PACKAGE PRIVATE ***************
	//************************************
	
	public JGroupShape createPlot(PApplet applet, int w, int h) {
		if(isGeoAxis) {
			double r = 0.5d * Math.max((maxX-minX)/pw, (maxY-minY)/ph);
			double xm = 0.5d*(minX+maxX);
			double ym = 0.5d*(minY+maxY);
			minX = xm - r*pw; maxX = xm + r*pw;
			minY = ym - r*ph; maxY = ym + r*ph;
		}
		if(pplot.isDebug())
			System.out.println("[DEBUG] JAxis-object: min/max={x:"+minX+"/"+maxX+", y:"+minY+"/"+maxY+
				"} with "+layers.size()+" layer"+(layers.size()>1?"s":""));
		JGroupShape graph = new JGroupShape();
		if(isGeoAxis) {
			projection.drawBorder(this, graph);
		} else {
			if(xAxOn || xGrdOn) graph.addChild(createXAxis());
			if(yAxOn || yGrdOn) graph.addChild(createYAxis());
			if(xAxOn || yAxOn) {
				if(xAxOn) {
					graph.addChild(new JLineShape(px, py,    px+pw, py   ));
					graph.addChild(new JLineShape(px, py+ph, px+pw, py+ph));
				}
				if(yAxOn) {
					graph.addChild(new JLineShape(px,    py, px,    py+ph));
					graph.addChild(new JLineShape(px+pw, py, px+pw, py+ph));
				}
			}
		}
		for(int l=0; l<layers.size(); l++) {
			JPlotsLayer layer = layers.get(l);
			layer.setRange(minX,maxX,minY,maxY);
			if(xAxInv || yAxInv)
				layer.invert(xAxInv ? (yAxInv ? "both" : "x") : "y", true);
			layer.createVectorImg(this, l, graph);
		}
		if(titleP.length()>0) {
			if(pplot.isDebug())
				System.out.println("[DEBUG] JAxis: add title \""+titleP+"\" to graphic.");
			graph.addChild(createTitle());
		}
		return graph;
	}
	public JProjection getProjection() {
		return projection; }

	//************************************
	//**** PRIVATE ***********************
	//************************************
	
	private void readParams(JPlotsLayer layer, Object... params) {
		if(params==null)
			return;
		int o=0;
		while(o<params.length) {
			if(params[o] instanceof String) {
				String p = ((String) params[o]).toLowerCase();
				boolean isunread = true;
				if(isunread && ("tf".equals(p) || "transform".equals(p)) && o+1<params.length) {
					layer.setSourceProjection((JProjection) params[o+1]);
					o++; isunread=false;
				}
				if(isunread && ("am".equals(p) || "anglemode".equals(p)) && o+1<params.length) {
					layer.angleMode((String) params[o+1]);
					o++; isunread=false;
				}
				if(isunread && ("l".equals(p) || "lines".equals(p)) && o+1<params.length) {
					layer.lines((boolean)params[o+1]);
					o++; isunread=false;
				}
				if(isunread && ("ls".equals(p) || "linestyle".equals(p)) && o+1<params.length) {
					layer.setLineStyle((String)params[o+1]);
					o++; isunread=false;
				}
				if(isunread && ("lw".equals(p) || "linewidth".equals(p)) && o+1<params.length) {
					layer.setLineWidth(params[o+1] instanceof Float ? (float)params[o+1] : (double)params[o+1]);
					o++; isunread=false;
				}
				if(isunread && ("lc".equals(p) || "linecolor".equals(p) || "linecolour".equals(p)) && o+1<params.length) {
					if(params[o+1] instanceof Integer) {
						layer.setLineColour((int)params[o+1]);
						o++; isunread=false;
					} else if(params[o+1] instanceof int[]) {
						layer.setLineColour((int[])params[o+1]);
						o++; isunread=false;
					}
				}
				if(isunread && ("lb".equals(p) || "label".equals(p))) {
					if(params[o+1] instanceof String) {
						layer.setLabel((String) params[o+1]);
						o++; isunread=false;
					}
				}
				if(isunread && ("ct".equals(p) || "colortable".equals(p) || "colourtable".equals(p)) && o+1<params.length) {
					if(params[o+1] instanceof ColourSequenceJColourtable) {
						layer.setColourtable((ColourSequenceJColourtable) params[o+1]);
						o++; isunread=false;
					} else
					if(params[o+1] instanceof LinearSegmentedJColourtable) {
						layer.setColourtable((LinearSegmentedJColourtable) params[o+1]);
						o++; isunread=false;
					} else
					if(params[o+1] instanceof JColourtable){
						layer.setColourtable((JColourtable) params[o+1]);
						o++; isunread=false;
					}
				}
				if(isunread && ("z".equals(p)) && o+1<params.length) {
					layer.addParallelArray(params[o+1]); o++; isunread=false; }
				if(isunread && "invertxaxis".equals(p)) {
					layer.invert("x", true); xAxInv=true; isunread=false; }
				if(isunread && "invertyaxis".equals(p)) {
					layer.invert("y", true); yAxInv=true; isunread=false; }
			} else {
				System.err.println("[ERROR] Cannot interprete param "+o+": "+params[o]);
			}
			o++;
		}
	}
	
	private void updateRange(JPlotsLayer layer) {
		double[] r = layer.getRange();
		double xmin=r[0],xmax=r[1], ymin=r[2],ymax=r[3];
		if(r[1]-r[0]<1.e-20) {
			double xm=0.5d*(r[0]+r[1]); double xr = Math.max(1.0e-10d, Math.abs(xm)*1.0e-10d);
			xmin = xm-xr; xmax = xm+xr;
		}
		if(r[3]-r[2]<1.e-20) {
			double ym=0.5d*(r[2]+r[3]); double yr = Math.max(1.0e-10d, Math.abs(ym)*1.e-10d);
			ymin = ym-yr; ymax = ym+yr;
		}
		if(layers.size()==1) {
			if(!xRangeFix) {
				minX = xmin;
				maxX = xmax;
			}
			if(!yRangeFix) {
				minY = ymin;
				maxY = ymax;
			}
		} else {
			if(!xRangeFix) {
				if(xmin<minX) minX = xmin;
				if(xmax>maxX) maxX = xmax;
			}
			if(!yRangeFix) {
				if(ymin<minY) minY = ymin;
				if(ymax>maxY) maxY = ymax;
			}
		}
	}
	private void updateRange(JPlotsLayer layer, String axis) {
		double[] r = layer.getRange();
		double xmin=r[0],xmax=r[1], ymin=r[2],ymax=r[3];
		if(r[1]-r[0]<1.e-20) {
			double xm=0.5d*(r[0]+r[1]); double xr = Math.max(1.0e-10d, Math.abs(xm)*1.0e-10d);
			xmin = xm-xr; xmax = xm+xr;
		}
		if(r[3]-r[2]<1.e-20) {
			double ym=0.5d*(r[2]+r[3]); double yr = Math.max(1.0e-10d, Math.abs(ym)*1.e-10d);
			ymin = ym-yr; ymax = ym+yr;
		}
		if(layers.size()==1) {
			if(!xRangeFix && axis.equals("x")) {
				minX = xmin;
				maxX = xmax;
			}
			if(!yRangeFix && axis.equals("y")) {
				minY = ymin;
				maxY = ymax;
			}
		} else {
			if(!xRangeFix && axis.equals("x")) {
				if(xmin<minX) minX = xmin;
				if(xmax>maxX) maxX = xmax;
			}
			if(!yRangeFix && axis.equals("y")) {
				if(ymin<minY) minY = ymin;
				if(ymax>maxY) maxY = ymax;
			}
		}
	}
	private JGroupShape createXAxis() {
		JGroupShape axisgrid = new JGroupShape();
		//first estimate of ticks
		double[] oticks = JPlotMath.optimalLinearTicks(minX, maxX);
		double vf = 1d/(oticks[0]);
		int decimal = (int) (1000d*oticks[1]+0.5d);
		decimal = decimal%100==0 ? 1 : decimal%10==0 ? 2 : 3;
		double tmlen = 0d, tmlc = 0d;
		pplot.getGraphic().textSize(200);
		pplot.getGraphic().textAlign(PApplet.LEFT,PApplet.TOP);
		//create tickmark strings and calc mean tickmark text width
		for(int t=2; t<oticks.length; t++) {
			String tm = PApplet.nf((float)(oticks[t]*vf),0,decimal).replace(",",".");
			tmlen += pplot.getGraphic().textWidth(tm) / 200f;
			tmlc += 1d;
		}
		tmlen *= this.txtsize / (oticks.length-2);
		//with upper bound of number of ticks
		int tickcount = Math.max(2, (int) (pw/(1.2d*tmlen)+0.99999999d));
		if(pplot.isDebug())
			System.out.println("[DEBUG] JAxis-object: tmlen="+tmlen+" -> tickcount approx. "+tickcount);
		//create new ticks
		double[] ticks = JPlotMath.optimalLinearTicks(minX, maxX, tickcount);
		double[] tcpos = JPlotMath.dlerp(ticks,minX,maxX,px,px+pw);
		String[] tickmark = new String[ticks.length];
		if(xAxInv)
			for(int t=0; t<tcpos.length; t++)
				tcpos[t] = 2*px+pw-tcpos[t];
		vf = 1d/(ticks[0]); decimal = (int) (1000d*ticks[1]+0.5d);
		decimal = decimal%1000==0 ? 0 : decimal%100==0 ? 1 : decimal%10==0 ? 2 : 3;
		for(int t=0; t<ticks.length; t++)
			tickmark[t] = PApplet.nf((float)(ticks[t]*vf),0,decimal).replace(",",".");
		if(pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for(int t=2; t<ticks.length; t++) {
				tickStr += ", "+tickmark[t];
				posStr  += ", "+PApplet.nf((float)tcpos[t],0,2).replace(",",".");
			}
			System.out.println("[DEBUG] JAxis-object: Xtickfactors={p10: "+ticks[0]+", f: "+ticks[1]+"}");
			System.out.println("[DEBUG] JAxis-object: Xtickval={"+tickStr.substring(2)+"}");
			System.out.println("[DEBUG] JAxis-object: Xtickpos={"+posStr.substring(2)+"}");
		}
		if(xGrdOn) {
			JPlotShape.stroke(0xff999999); JPlotShape.strokeWeight(2f);
			for(int t=2; t<ticks.length; t++)
				if(ticks[t]>=Math.min(minX, maxX) && ticks[t]<=Math.max(minX, maxX))
					axisgrid.addChild(new JLineShape((float)tcpos[t],py,(float)tcpos[t],py+ph));
		}
		if(xAxOn) {
			if(xTkOn) {
				JPlotShape.stroke(0xff000000); JPlotShape.strokeWeight(2f);
				for(int t=2; t<ticks.length; t++)
					if(ticks[t]>=Math.min(minX, maxX) && ticks[t]<=Math.max(minX, maxX)) {
						axisgrid.addChild(new JLineShape((float)tcpos[t],py+ph,(float)tcpos[t],py+1.02f*ph));
						//axisgrid.addChild(ap.createShape(PShape.TEXT, "H", (float)tcpos[t],py-0.1f*ph,(float)tcpos[t],py));
						axisgrid.addChild(new JTextShape(tickmark[t], (float)tcpos[t], py+1.03f*ph, (float)txtsize, PApplet.CENTER, PApplet.TOP, 0xff000000, 0));
					}
			}
			if(titleX.length()>0) {
				if(pplot.isDebug())
					System.out.println("[DEBUG] JAxi-object: add x-axis title \""+titleX+"\"");
				axisgrid.addChild(new JTextShape(titleX, px+0.5f*pw, py+1.04f*ph+(float)txtsize, (float)(1.1d*txtsize), PApplet.CENTER, PApplet.TOP, 0xff000000, 0));
			}
		}
		return axisgrid;
	}
	private JGroupShape createYAxis() {
		JGroupShape axisgrid = new JGroupShape();
		double[] ticks = JPlotMath.optimalLinearTicks(minY, maxY);
		double[] tcpos = JPlotMath.dlerp(ticks,minY,maxY,py+ph,py);
		if(yAxInv)
			for(int t=0; t<tcpos.length; t++)
				tcpos[t] = 2*py+ph-tcpos[t];
		if(pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for(int t=2; t<ticks.length; t++) {
				tickStr += ", "+PApplet.nf((float)ticks[t],0,2);
				posStr  += ", "+PApplet.nf((float)tcpos[t],0,2);
			}
			System.out.println("[DEBUG] JAxis-object: Ytickfactors={p10: "+ticks[0]+", f: "+ticks[1]+"}");
			System.out.println("[DEBUG] JAxis-object: Ytickval={"+tickStr.substring(2)+"}");
			System.out.println("[DEBUG] JAxis-object: Ytickpos={"+posStr.substring(2)+"}");
		}
		if(yGrdOn) {
			JPlotShape.stroke(0xff999999); JPlotShape.strokeWeight(2f);
			for(int t=2; t<ticks.length; t++)
				if(ticks[t]>=Math.min(minY, maxY) && ticks[t]<=Math.max(minY, maxY))
					axisgrid.addChild(new JLineShape(px,(float)tcpos[t],px+pw,(float)tcpos[t]));
		}
		if(yAxOn) {
			JPlotShape.stroke(0xff000000); JPlotShape.strokeWeight(2f);
			float tw = 0f;
			if(yTkOn) {
				double vf = 1d/(ticks[0]);
				int decimal = (int) (1000d*ticks[1]+0.5d);
				decimal = decimal%1000==0 ? 0 : decimal%100==0 ? 1 : decimal%10==0 ? 2 : 3;
				for(int t=2; t<ticks.length; t++)
					if(ticks[t]>=Math.min(minY, maxY) && ticks[t]<=Math.max(minY, maxY)) {
						axisgrid.addChild(new JLineShape(px-0.02f*pw,(float)tcpos[t],px,(float)tcpos[t]));
						String tmstr = PApplet.nf((float)(ticks[t]*vf),0,decimal).replace(",",".");
						tw = Math.max(tw, (float)txtsize * pplot.getGraphic().textWidth(tmstr)/pplot.getGraphic().textSize);
						axisgrid.addChild(new JTextShape(tmstr, px-0.03f*pw, (float)tcpos[t], (float)txtsize,
								PApplet.RIGHT, PApplet.CENTER, 0xff000000, 0));
					}
			}
			if(titleY.length()>0)
				axisgrid.addChild(new JTextShape(titleY, px-0.03f*pw-tw, py+0.5f*ph, (float)(1.1d*txtsize), PApplet.CENTER, PApplet.BOTTOM, 0xff000000, JPlotShape.ROTATE_COUNTERCLOCKWISE));
		}
		return axisgrid;
	}
	private JTextShape createTitle() {
		return new JTextShape(titleP, px+0.5f*pw, py-0.04f*ph, (float)(1.3d*txtsize), PApplet.CENTER, PApplet.BOTTOM, 0xff000000, 0f);
	}
	
	private static PImage loadPreDefImg(PApplet applet, String name) {
		BufferedImage bimg;
		try {
			bimg = ImageIO.read(JPlot.class.getResourceAsStream("/data/"+name+".png"));
			PImage limg = applet.createImage(bimg.getWidth(), bimg.getHeight(), PApplet.ARGB);
			bimg.getRGB(0, 0, bimg.getWidth(), bimg.getHeight(), limg.pixels, 0, bimg.getWidth());
			loadedPreDefImgs.put(name, limg);
			return limg;
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.err.println("The image \""+name+".png\" is missing in the jar or inaccessible. Contact the author of ##library.name## ##library.prettyVersion## if the image should be there.");
		return null;
	}
}
