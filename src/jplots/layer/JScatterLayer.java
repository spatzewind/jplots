package jplots.layer;

import jplots.JAxis;
import jplots.JPlot;
import jplots.maths.JPlotMath;
import jplots.shapes.JEllipseShape;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JRectShape;
import processing.core.PGraphics;

public class JScatterLayer extends JPlotsLayer {
	
	private double[] xarrayx, yarrayy;
	private int col;
	private double lw;
	private String ls;
	
	public JScatterLayer(float[] x, float[] y, int colour, float linewidth, String linestyle) {
		xarrayx = new double[x.length];
		for(int i=0; i<x.length; i++)
			xarrayx[i] = x[i];
		yarrayy = new double[y.length];
		for(int i=0; i<y.length; i++)
			yarrayy[i] = y[i];
		minX = JPlotMath.dmin(xarrayx);
		maxX = JPlotMath.dmax(xarrayx);
		minY = JPlotMath.dmin(yarrayy);
		maxY = JPlotMath.dmax(yarrayy);
		col = colour;
		lw = linewidth;
		ls = linestyle;
	}
	public JScatterLayer(double[] x, double[] y, int colour, double linewidth, String linestyle) {
		xarrayx = x;
		yarrayy = y;
		minX = JPlotMath.dmin(xarrayx);
		maxX = JPlotMath.dmax(xarrayx);
		minY = JPlotMath.dmin(yarrayy);
		maxY = JPlotMath.dmax(yarrayy);
		col = colour;
		lw = linewidth;
		ls = linestyle;
	}
	

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		if("o".equals(ls) || "@".equals(ls)) {
			JPlotShape.fill(col); } else { JPlotShape.noFill(); }
		JPlotShape.stroke(col); JPlotShape.strokeWeight((float)lw);
		JGroupShape xyShape = new JGroupShape();
		double xs = p[2]/(maxX-minX), ys = p[3]/(maxY-minY);
		for(int i=0; i<xarrayx.length; i++)
			if(Double.isFinite(xarrayx[i]) && Double.isFinite(yarrayy[i])) {
				double[] xy = inputProj.fromPROJtoLATLON(xarrayx[i], yarrayy[i], false);
				if(ax.isGeoAxis())
					xy = ax.getGeoProjection().fromLATLONtoPROJ(xy[0], xy[1], false);
				double x1 = p[0]+xs*(xy[0]-minX);
				double y1 = p[1]+ys*(maxY-xy[1]);
				if(x1<p[0] || x1>p[0]+p[2]) continue;
				if(y1<p[1] || y1>p[1]+p[3]) continue;
				if("()".equals(ls) || "o".equals(ls)) {
					xyShape.addChild(new JEllipseShape((float)x1, (float)y1, (float)(6*lw), (float)(6*lw))); }
				if("[]".equals(ls) || "@".equals(ls)) {
					xyShape.addChild(new JRectShape((float)(x1-2*lw), (float)(y1-2*lw), (float)(x1+2*lw), (float)(y1+2*lw))); }
				if("x".equals(ls)) {
					xyShape.addChild(new JLineShape((float)(x1-2*lw), (float)(y1-2*lw), (float)(x1+2*lw), (float)(y1+2*lw)));
					xyShape.addChild(new JLineShape((float)(x1-2*lw), (float)(y1+2*lw), (float)(x1+2*lw), (float)(y1-2*lw))); }
				if("+".equals(ls)) {
					xyShape.addChild(new JLineShape((float)(x1-3*lw), (float)y1, (float)(x1+3*lw), (float)y1));
					xyShape.addChild(new JLineShape((float)x1, (float)(y1-3*lw), (float)x1, (float)(y1+3*lw))); }
				if("<>".equals(ls)) {
					xyShape.addChild(new JLineShape((float)x1, (float)(y1-3*lw), (float)(x1+3*lw), (float)y1));
					xyShape.addChild(new JLineShape((float)(x1+3*lw), (float)y1, (float)x1, (float)(y1+3*lw)));
					xyShape.addChild(new JLineShape((float)x1, (float)(y1+3*lw), (float)(x1-3*lw), (float)y1));
					xyShape.addChild(new JLineShape((float)(x1-3*lw), (float)y1, (float)x1, (float)(y1-3*lw))); }
			}
		s.addChild(xyShape);
	}

}
