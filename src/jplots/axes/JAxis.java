package jplots.axes;

import java.util.List;

import jplots.JPlot;
import jplots.JPlotConstants;
import jplots.layer.JPlotsLayer;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLatexShape;
import jplots.shapes.JTextShape;
import jplots.transform.IdentityJProjection;
import jplots.transform.JProjection;
import processing.core.PApplet;
import processing.core.PFont;

public abstract class JAxis implements JPlotConstants {
	
	//reference
	protected JPlot pplot;
	//style
	protected int px, py, pw, ph;
	protected double txtsize, tickscale;
	protected PFont pfont;
	//content
	protected boolean xAxOn, yAxOn, xGrdOn, yGrdOn, xTkOn, yTkOn, xTkLbOn, yTkLbOn;
	protected int xDrawSide, yDrawSide;
	protected String titleP;
	protected JProjection projection;
	
	public JAxis(JPlot plot, int pos_x, int pos_y, int width, int height) {
		pplot = plot;
		px = pos_x;
		py = pos_y;
		pw = width;
		ph = height;
		defaults();
	}
	public JAxis(JAxis src_axis) {
		pplot = src_axis.getPlot();
		px    = src_axis.px;
		py    = src_axis.py;
		pw    = src_axis.pw;
		ph    = src_axis.ph;
		defaults();
	}
	/**
	 * set default values for
	 * <ul>
	 * <li>xAxOn/yAxOn: \t show x/y-axis, (true)</li>
	 * <li>xGrdOn/yGrdOn: \t show grid lines (false)</li>
	 * <li>xTkOn/yTkOn: \t show tick marks (true)</li>
	 * <li>xTkLbOn/yTkLbOn: \t show tick labes (true)</li>
	 * <li>xDrawSide/yDrawSide: \t side to put ticks/labels/axis-title (BOTTOM/LEFT)</li>
	 * <li>titleP: \t title of the subplot ("")</li>
	 * <li>txtsize: \t reference size for text (12pt by 300dpi)</li>
	 * <li>projection: \t geographic projection (identity = no projection)</li>
	 * </ul>
	 */
	private void defaults() {
		xAxOn = true;
		yAxOn = true;
		xGrdOn = false;
		yGrdOn = false;
		xTkOn = true;
		yTkOn = true;
		xTkLbOn = true;
		yTkLbOn = true;
		titleP = "";
		txtsize = JPlot.dpi * 10d / 72d;
		tickscale = JPlot.dpi * 10d / 450d;
		xDrawSide = BOTTOM;
		yDrawSide = LEFT;
		projection = new IdentityJProjection();
	}
	
	
	
	
	// ******************
	// **** ABSTRACT ****
	// ******************
	
	public abstract JAxis copy();
	public abstract JGroupShape createPlot(PApplet applet, int w, int h);
	public abstract JGroupShape createPlotOnlyAxes(PApplet applet, int w, int h);
	public abstract int[] getSize();
	public abstract double[] getRange();
	public abstract AxisScale getScaleX();
	public abstract AxisScale getScaleY();
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
	public AxisScale getScaleX(int subaxis) {
		return getScaleX();
	}
	public AxisScale getScaleY(int subaxis) {
		return getScaleY();
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
			+"    x-Axis: draw="+xAxOn+" ticks="+xTkOn+" tickLabels="+xTkLbOn+" grid="+xGrdOn+"\n"
			+"    y-Axis: draw="+yAxOn+" ticks="+yTkOn+" tickLabels="+yTkLbOn+" grid="+yGrdOn+"\n"
			+"    proj:   geo="+isGeoAxis()+" which="+projection+"\n"
			+"    title:  "+(titleP!=null?"\""+titleP+"\"":null)+"\n"
		);
	}
	
	public JColourbar colourbar() {
		return pplot.colourbar(this, "", "neither", 1);
	}
	public JColourbar colourbar(String name) {
		return pplot.colourbar(this, name, "neither", 1);
	}
	public JColourbar colourbar(String name, String extent) {
		return pplot.colourbar(this, name, extent, 1);
	}
	public JColourbar colourbar(String name, String extent, int orientation) {
		return pplot.colourbar(this, name, extent, orientation);
	}
	
	public void setFont(PFont font) {
		pfont = font;
	}
	public void setTextSize(double ts) {
		txtsize = JPlot.dpi * ts / 72d;
	}
	public void setTickScale(double ts) {
		tickscale = JPlot.dpi * ts / 450d;
	}
	public void showXAxisSide(String side) {
		if(side.equalsIgnoreCase("none"))   { xDrawSide = NONE;   xAxOn = false; return; }
		if(side.equalsIgnoreCase("top"))    { xDrawSide = TOP;    xAxOn = true;  return; }
		if(side.equalsIgnoreCase("bottom")) { xDrawSide = BOTTOM; xAxOn = true;  return; }
		if(side.equalsIgnoreCase("both"))   { xDrawSide = BOTH;   xAxOn = true;  return; }
		System.err.println("Unknown value "+side+" for side.");
	}
	public void showXAxisSide(int side) {
		if(side==NONE)   { xDrawSide = NONE;   xAxOn = false; return; }
		if(side==TOP)    { xDrawSide = TOP;    xAxOn = true;  return; }
		if(side==BOTTOM) { xDrawSide = BOTTOM; xAxOn = true;  return; }
		if(side==BOTH)   { xDrawSide = BOTH;   xAxOn = true;  return; }
		System.err.println("Unknown value "+side+" for side.");
	}
	public void showYAxisSide(String side) {
		if(side.equalsIgnoreCase("none"))  { yDrawSide = NONE;  yAxOn = false; return; }
		if(side.equalsIgnoreCase("left"))  { yDrawSide = LEFT;  yAxOn = true;  return; }
		if(side.equalsIgnoreCase("right")) { yDrawSide = RIGHT; yAxOn = true;  return; }
		if(side.equalsIgnoreCase("both"))  { yDrawSide = BOTH;  yAxOn = true;  return; }
		System.err.println("Unknown value "+side+" for side.");
	}
	public void showYAxisSide(int side) {
		if(side==NONE)  { yDrawSide = NONE;  yAxOn = false; return; }
		if(side==LEFT)  { yDrawSide = LEFT;  yAxOn = true;  return; }
		if(side==RIGHT) { yDrawSide = RIGHT; yAxOn = true;  return; }
		if(side==BOTH)  { yDrawSide = BOTH;  yAxOn = true;  return; }
		System.err.println("Unknown value "+side+" for side.");
	}
	public void setTitle(String _title) {
		titleP = _title;
	}
	
	
	
	
	// *******************
	// **** PROTECTED ****
	// *******************
	
	
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
		boolean setAx = false, setGrd = false, setTck = false, setTkLab = false;
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
		if ("l".equals(w) || "label".equals(w) || "ticklabel".equals(w) || "labels".equals(w) || "ticklabels".equals(w))
			setTkLab = true;
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
			if (setX) {
				xTkOn = onoff;
				if(!onoff) xTkLbOn = onoff;
			}
			if (setY) {
				yTkOn = onoff;
				if(!onoff) yTkLbOn = onoff;
			}
		}
		if (setTkLab) {
			if (setX)
				xTkLbOn = onoff;
			if (setY)
				yTkLbOn = onoff;
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
	
	
	
	protected void addAxisText(JGroupShape gs, char ax, String txt, double posX, double posY, double ts, int alignX, int alignY, int col, double rotation, String style) {
		float posx = (float)posX, posy = (float)posY, rot = (float)rotation, tsf = (float) ts;
		if(ax=='x') {
			if(xDrawSide==BOTTOM || xDrawSide==BOTH) {
				if(JPlot.supportLatex) gs.addChild(new JLatexShape(txt, posx, posy, tsf, alignX, alignY, col, rot, style));
				else                   gs.addChild(new JTextShape(txt, posx, posy, tsf, alignX, alignY, col, rot, style));
			}
			if(xDrawSide==TOP || xDrawSide==BOTH) {
				posy = 2*py+ph - posy;
//				int alignx = alignX==LEFT?RIGHT : alignX==RIGHT?LEFT : alignX;
				int aligny = alignY==TOP?BOTTOM : alignY==BOTTOM?TOP : alignY;
				if(JPlot.supportLatex) gs.addChild(new JLatexShape(txt, posx, posy, tsf, alignX, aligny, col, rot, style));
				else                   gs.addChild(new JTextShape(txt, posx, posy, tsf, alignX, aligny, col, rot, style));
			}
		}
		if(ax=='y') {
			if(yDrawSide==LEFT || yDrawSide==BOTH) {
				if(JPlot.supportLatex) gs.addChild(new JLatexShape(txt, posx, posy, tsf, alignX, alignY, col, rot, style));
				else                   gs.addChild(new JTextShape(txt, posx, posy, tsf, alignX, alignY, col, rot, style));
			}
			if(yDrawSide==RIGHT || yDrawSide==BOTH) {
				posx = 2*px+pw - posx;
				int alignx = alignX==LEFT?RIGHT : alignX==RIGHT?LEFT : alignX;
				int aligny = alignY==TOP?BOTTOM : alignY==BOTTOM?TOP : alignY;
				if(JPlot.supportLatex) gs.addChild(new JLatexShape(txt, posx, posy, tsf, alignx, aligny, col, rot, style));
				else                   gs.addChild(new JTextShape(txt, posx, posy, tsf, alignx, aligny, col, rot, style));
			}
		}
	}
}
