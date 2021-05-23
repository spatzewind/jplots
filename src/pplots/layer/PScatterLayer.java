package pplots.layer;

import pplots.PMath;
import pplots.PPlot;
import pplots.gfx.PEllipseShape;
import pplots.gfx.PGroupShape;
import pplots.gfx.PLineShape;
import pplots.gfx.PPlotShape;
import pplots.gfx.PRectShape;
import processing.core.PApplet;
import processing.core.PGraphics;

public class PScatterLayer extends PLayer {
	
	private double[] xarrayx, yarrayy;
	private int col;
	private double lw;
	private String ls;
	private double minX,maxX,minY,maxY;
	
	public PScatterLayer(float[] x, float[] y, int colour, float linewidth, String linestyle) {
		xarrayx = new double[x.length];
		for(int i=0; i<x.length; i++)
			xarrayx[i] = x[i];
		yarrayy = new double[y.length];
		for(int i=0; i<y.length; i++)
			yarrayy[i] = y[i];
		minX = PMath.dmin(xarrayx);
		maxX = PMath.dmax(xarrayx);
		minY = PMath.dmin(yarrayy);
		maxY = PMath.dmax(yarrayy);
		col = colour;
		lw = linewidth;
		ls = linestyle;
	}
	public PScatterLayer(double[] x, double[] y, int colour, double linewidth, String linestyle) {
		xarrayx = x;
		yarrayy = y;
		minX = PMath.dmin(xarrayx);
		maxX = PMath.dmax(xarrayx);
		minY = PMath.dmin(yarrayy);
		maxY = PMath.dmax(yarrayy);
		col = colour;
		lw = linewidth;
		ls = linestyle;
	}
	
	
	@Override
	public void setRange(double xmin, double xmax, double ymin, double ymax) {
		minX = xmin; maxX = xmax;
		minY = ymin; maxY = ymax;
	}
	
	@Override
	public double[] getRange() {
		return new double[] {minX,maxX,minY,maxY};
	}

	@Override
	public void createRasterImg(PPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(PPlot plot, PGroupShape s, int x, int y, int w, int h) {
		PApplet ap = plot.getApplet();
		if("o".equals(ls) || "@".equals(ls)) {
			PPlotShape.fill(col); } else { PPlotShape.noFill(); }
		PPlotShape.stroke(col); ap.strokeWeight((float)lw);
		PGroupShape xyShape = new PGroupShape();
		double xs = w/(maxX-minX), ys = h/(maxY-minY);
		for(int i=0; i<xarrayx.length; i++)
			if(Double.isFinite(xarrayx[i]) && Double.isFinite(yarrayy[i])) {
				double x1 = x+xs*(xarrayx[i]-minX);
				double y1 = y+ys*(maxY-yarrayy[i]);
				if(x1<x || x1>x+w) continue;
				if(y1<y || y1>y+h) continue;
				if("c".equals(ls) || "o".equals(ls)) {
					xyShape.addChild(new PEllipseShape((float)x1, (float)y1, (float)(6*lw), (float)(6*lw))); }
				if("[]".equals(ls) || "@".equals(ls)) {
					xyShape.addChild(new PRectShape((float)(x1-2*lw), (float)(y1-2*lw), (float)(x1+2*lw), (float)(y1+2*lw))); }
				if("x".equals(ls)) {
					xyShape.addChild(new PLineShape((float)(x1-2*lw), (float)(y1-2*lw), (float)(x1+2*lw), (float)(y1+2*lw)));
					xyShape.addChild(new PLineShape((float)(x1-2*lw), (float)(y1+2*lw), (float)(x1+2*lw), (float)(y1-2*lw))); }
				if("+".equals(ls)) {
					xyShape.addChild(new PLineShape((float)(x1-3*lw), (float)y1, (float)(x1+3*lw), (float)y1));
					xyShape.addChild(new PLineShape((float)x1, (float)(y1-3*lw), (float)x1, (float)(y1+3*lw))); }
				if("<>".equals(ls)) {
					xyShape.addChild(new PLineShape((float)x1, (float)(y1-3*lw), (float)(x1+3*lw), (float)y1));
					xyShape.addChild(new PLineShape((float)(x1+3*lw), (float)y1, (float)x1, (float)(y1+3*lw)));
					xyShape.addChild(new PLineShape((float)x1, (float)(y1+3*lw), (float)(x1-3*lw), (float)y1));
					xyShape.addChild(new PLineShape((float)(x1-3*lw), (float)y1, (float)x1, (float)(y1-3*lw))); }
			}
		s.addChild(xyShape);
	}

}
