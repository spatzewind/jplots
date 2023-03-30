package jplots.layer;

import jplots.JPlot;
import jplots.axes.JAxis;
import jplots.colour.JColourtable;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.transform.IdentityJProjection;
import jplots.transform.JProjection;
import processing.core.PGraphics;
import processing.core.PImage;

public abstract class JPlotsLayer {
	
	protected boolean invertAxisX, invertAxisY, logAxisX, logAxisY;
	protected boolean drawLines, singleLineColour, singleFillColour;
	protected int lc, pc;
	protected int[] lcs, pcs;
	protected double lw;
	protected double minX, maxX, minY, maxY, minZ, maxZ;
	protected double angleMode;
	protected String ls, label;
	protected PImage img;
	protected JProjection inputProj;
	protected JColourtable colourtable;
	protected Object parallelArray;
	
	
	public JPlotsLayer() {
		lw = 2d;
		lc = 0xff000000;
		pc = 0xff7f7f7f;
		ls = "-";
		singleLineColour = true;
		lcs = new int[] { 0xff000000 };
		pcs = new int[] { 0xff7f7f7f };
		minX = -1d;
		maxX = 1d;
		minY = -1d;
		maxY = 1d;
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
	
	// ******************
	// **** ABSTRACT ****
	// ******************
	
	public abstract void createRasterImg(JPlot plot, PGraphics g);
	public abstract void createVectorImg(JAxis ax, int layernum, JGroupShape s);
	
	// ****************
	// **** PUBLIC ****
	// ****************
	
	public void setRange(double xmin, double xmax, double ymin, double ymax) {
		minX = xmin;
		maxX = xmax;
		minY = ymin;
		maxY = ymax;
	}
	public double[] getRange() {
		return new double[] { minX, maxX, minY, maxY };
	}
	public void setZRange(double zmin, double zmax) {
		minZ = zmin;
		maxZ = zmax;
	}

	public void setSourceProjection(JProjection in_proj) {
		inputProj = in_proj;
	}
	
	public void setColourtable(JColourtable ct) {
		colourtable = ct;
	}
	public JColourtable getColourtable() {
		return colourtable;
	}
	
	public void angleMode(String angle_type) {
		if (angle_type == null) {
			System.err.println("[ERROR] Cannot read anglemode from null!");
			return;
		}
		if ("degrees".equals(angle_type.toLowerCase())) {
			angleMode = JPlotMath.DEG_TO_RAD;
			return;
		}
		if ("gon".equals(angle_type.toLowerCase())) {
			angleMode = Math.PI / 200d;
			return;
		}
		if ("radians".equals(angle_type.toLowerCase())) {
			angleMode = 1d;
			return;
		}
		System.err.println("[ERROR] Cannot interprete anglemode \"" + angle_type + "\". Reset to RADIANS");
		angleMode = 1d;
	}
	
	public void lines(boolean l) {
		drawLines = l;
	}
	public void setLineColour(int _lc) {
		lc = _lc;
		singleLineColour = true;
	}
	public void setLineColours(int[] _lcs) {
		lcs = _lcs;
		singleLineColour = false;
	}
	public int getLineColour() {
		return singleLineColour ? lc : 0;
	}
	
	public void setLineWidth(double lwd) {
		lw = lwd;
	}
	public double getLineWidth() {
		return lw;
	}
	
	public void setFillColour(int _pc) {
		pc = _pc;
		singleFillColour = true;
	}
	public void setFillColours(int[] _pcs) {
		pcs = _pcs;
		singleFillColour = false;
	}
	public int getFillColour() {
		return singleFillColour ? pc : 0;
	}

	public void setStyle(String lst) {
		ls = lst;
	}
	public String getStyle() {
		return ls;
	}

	public void invert(String which, boolean _invert) {
		boolean ix = false, iy = false;
		if ("both".equalsIgnoreCase(which)) {
			ix = true;
			iy = true;
		}
		if ("x".equalsIgnoreCase(which))
			ix = true;
		if ("y".equalsIgnoreCase(which))
			iy = true;
		if (ix)
			invertAxisX = _invert;
		if (iy)
			invertAxisY = _invert;
	}

	public void addParallelArray(Object parallel) {
		parallelArray = parallel;
	}

	public void setLabel(String _l) {
		label = _l;
	}
	public String getLabel() {
		return label;
	}
	
}
