package jplots.axes;

import java.util.List;

import jplots.JPlot;
import jplots.JPlotConstants;
import jplots.layer.JPlotsLayer;
import jplots.shapes.JGroupShape;
import jplots.shapes.JPlotShape;
import jplots.transform.IdentityJProjection;
import jplots.transform.JProjection;
import processing.core.PApplet;
import processing.core.PFont;

public abstract class JAxis {
	
	public final static int NONE   = JPlotConstants.NONE;
	public final static int TOP    = JPlotConstants.TOP;
	public final static int LEFT   = JPlotConstants.LEFT;
	public final static int RIGHT  = JPlotConstants.RIGHT;
	public final static int BOTTOM = JPlotConstants.BOTTOM;
	public final static int BOTH   = JPlotConstants.BOTH;
	public final static int CENTER = JPlotConstants.CENTER;
	
	public final static float ROTATE_CLOCKWISE        = JPlotConstants.ROTATE_CLOCKWISE;
	public final static float ROTATE_COUNTERCLOCKWISE = JPlotConstants.ROTATE_COUNTERCLOCKWISE;
	
	//reference
	protected JPlot pplot;
	//style
	protected int px, py, pw, ph;
	protected double txtsize;
	protected PFont pfont;
	//content
	protected boolean xAxOn, yAxOn, xGrdOn, yGrdOn, xTkOn, yTkOn;
	protected int xDrawSide, yDrawSide;
	protected String titleP;
	protected JProjection projection;
	
	public JAxis(JPlot plot, int pos_x, int pos_y, int width, int height) {
		pplot = plot;
		px = pos_x;
		py = pos_y;
		pw = width;
		ph = height;
	}
	public JAxis(JAxis src_axis) {
		pplot = src_axis.getPlot();
		px    = src_axis.px;
		py    = src_axis.py;
		pw    = src_axis.pw;
		ph    = src_axis.ph;
	}
	
	
	
	
	// ******************
	// **** ABSTRACT ****
	// ******************
	
	public abstract JAxis copy();
	public abstract JGroupShape createPlot(PApplet applet, int w, int h);
	public abstract JPlotShape createPlotOnlyAxes(PApplet applet, int w, int h);
	public abstract int[] getSize();
	public abstract double[] getRange();
	public abstract boolean isXlogAxis();
	public abstract boolean isYlogAxis();
	public abstract boolean isGeoAxis();
	public abstract void setGeoProjection(JProjection proj);
	public abstract List<JPlotsLayer> getLayers();
	
	
	
	
	// ****************
	// **** PUBLIC ****
	// ****************
	
	/**
	 * removes all plotting infos also all configuration will be reseted
	 */
	public void clear() {
		defaults();
	}
	
	public JPlot getPlot() {
		return pplot;
	}
	public double getTextSize() {
		return txtsize;
	}
	public JProjection getGeoProjection() {
		return projection;
	}
	public boolean isXGridVisible() {
		return xGrdOn;
	}
	public boolean isYGridVisible() {
		return yGrdOn;
	}
	public PFont getFont() {
		return pfont;
	}
	public void printInfo() {
		System.out.print(
			 this.getClass().getSimpleName()+":\n"
			+"    pos/size: "+px+"|"+py+" / "+pw+"|"+ph+"\n"
			+"    font/txtsize: "+pfont+"/"+txtsize+"\n"
			+"    x-Axis: log="+isXlogAxis()+" draw="+xAxOn+" ticks="+xTkOn+" grid="+xGrdOn+"\n"
			+"    y-Axis: log="+isYlogAxis()+" draw="+yAxOn+" ticks="+yTkOn+" grid="+yGrdOn+"\n"
			+"    proj:   geo="+isGeoAxis()+" which="+projection+"\n"
			+"    title:  "+(titleP!=null?"\""+titleP+"\"":null)+"\n"
		);
	}
	
	public JColourbar colourbar() {
		return pplot.colourbar(this, "", "neither");
	}
	public JColourbar colourbar(String name) {
		return pplot.colourbar(this, name, "neither");
	}
	public JColourbar colourbar(String name, String extent) {
		return pplot.colourbar(this, name, extent);
	}
	
	public void setFont(PFont font) {
		pfont = font;
	}
	public void setTextSize(double ts) {
		txtsize = JPlot.dpi * ts / 72d;
	}
	public void showXAxisSide(String side) {
		if(side.equalsIgnoreCase("none"))   { yDrawSide = NONE;   yAxOn = false; return; }
		if(side.equalsIgnoreCase("top"))    { yDrawSide = TOP;    yAxOn = true;  return; }
		if(side.equalsIgnoreCase("bottom")) { yDrawSide = BOTTOM; yAxOn = true;  return; }
		if(side.equalsIgnoreCase("both"))   { yDrawSide = BOTH;   yAxOn = true;  return; }
		System.err.println("Unknown value "+side+" for side.");
	}
	public void showXAxisSide(int side) {
		if(side==NONE)   { yDrawSide = NONE;   yAxOn = false; return; }
		if(side==TOP)    { yDrawSide = TOP;    yAxOn = true;  return; }
		if(side==BOTTOM) { yDrawSide = BOTTOM; yAxOn = true;  return; }
		if(side==BOTH)   { yDrawSide = BOTH;   yAxOn = true;  return; }
		System.err.println("Unknown value "+side+" for side.");
	}
	public void showYAxisSide(String side) {
		if(side.equalsIgnoreCase("none"))  { xDrawSide = NONE;  xAxOn = false; return; }
		if(side.equalsIgnoreCase("left"))  { xDrawSide = LEFT;  xAxOn = true;  return; }
		if(side.equalsIgnoreCase("right")) { xDrawSide = RIGHT; xAxOn = true;  return; }
		if(side.equalsIgnoreCase("both"))  { xDrawSide = BOTH;  xAxOn = true;  return; }
		System.err.println("Unknown value "+side+" for side.");
	}
	public void showYAxisSide(int side) {
		if(side==NONE)  { xDrawSide = NONE;  xAxOn = false; return; }
		if(side==LEFT)  { xDrawSide = LEFT;  xAxOn = true;  return; }
		if(side==RIGHT) { xDrawSide = RIGHT; xAxOn = true;  return; }
		if(side==BOTH)  { xDrawSide = BOTH;  xAxOn = true;  return; }
		System.err.println("Unknown value "+side+" for side.");
	}
	public void setTitle(String _title) {
		titleP = _title;
	}
	
	
	
	
	// *******************
	// **** PROTECTED ****
	// *******************
	
	protected void defaults() {
		xAxOn = true;
		yAxOn = true;
		xGrdOn = false;
		yGrdOn = false;
		xTkOn = true;
		yTkOn = true;
		titleP = "";
		txtsize = 300d * 10d / 72d;
		xDrawSide = BOTTOM;
		yDrawSide = LEFT;
		projection = new IdentityJProjection();
	}
	
	public JAxis setPositionAndSize(int pos_x, int pos_y, int width, int height) {
		px = pos_x;
		py = pos_y;
		pw = width;
		ph = height;
		if (pplot.isDebug())
			System.out.println("[DEBUG] resize PAxis-object: x/y=" + px + "/" + py + " w/h=" + pw + "/" + ph);
		return this;
	}
	

	public void setAxis(String axis, String which, boolean onoff) {
		boolean setX = false, setY = false;
		boolean setAx = false, setGrd = false, setTck = false;
		if ("both".equals(axis.toLowerCase())) {
			setX = true;
			setY = true;
		}
		if ("x".equals(axis.toLowerCase()))
			setX = true;
		if ("y".equals(axis.toLowerCase()))
			setY = true;
		String w = which.toLowerCase();
		if ("all".equals(w)) {
			setAx = true;
			setGrd = true;
		}
		if ("a".equals(w) || "axis".equals(w))
			setAx = true;
		if ("g".equals(w) || "grid".equals(w))
			setGrd = true;
		if ("t".equals(w) || "tick".equals(w) || "ticks".equals(w))
			setTck = true;
		if (setAx) {
			if (setX)
				xAxOn = onoff;
			if (setY)
				yAxOn = onoff;
		}
		if (setGrd) {
			if (setX)
				xGrdOn = onoff;
			if (setY)
				yGrdOn = onoff;
		}
		if (setTck) {
			if (setX)
				xTkOn = onoff;
			if (setY)
				yTkOn = onoff;
		}
	}
	public void setGrid() {
		setGrid("both", true);
	}
	public void setGrid(String axis, boolean onoff) {
		boolean setX = false, setY = false;
		if ("both".equals(axis.toLowerCase())) {
			setX = true;
			setY = true;
		}
		if ("x".equals(axis.toLowerCase()))
			setX = true;
		if ("y".equals(axis.toLowerCase()))
			setY = true;
		if (setX)
			xGrdOn = onoff;
		if (setY)
			yGrdOn = onoff;
	}
	
}
