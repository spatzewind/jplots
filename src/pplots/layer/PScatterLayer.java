package pplots.layer;

import pplots.PAxis;
import pplots.PPlotMath;
import pplots.shapes.PEllipseShape;
import pplots.shapes.PGroupShape;
import pplots.shapes.PLineShape;
import pplots.shapes.PPlotShape;
import pplots.shapes.PRectShape;
import pplots.PPlot;
import processing.core.PGraphics;

public class PScatterLayer extends PLayer {
	
	private double[] xarrayx, yarrayy;
	private int col;
	private double lw;
	private String ls;
	
	public PScatterLayer(float[] x, float[] y, int colour, float linewidth, String linestyle) {
		xarrayx = new double[x.length];
		for(int i=0; i<x.length; i++)
			xarrayx[i] = x[i];
		yarrayy = new double[y.length];
		for(int i=0; i<y.length; i++)
			yarrayy[i] = y[i];
		minX = PPlotMath.dmin(xarrayx);
		maxX = PPlotMath.dmax(xarrayx);
		minY = PPlotMath.dmin(yarrayy);
		maxY = PPlotMath.dmax(yarrayy);
		col = colour;
		lw = linewidth;
		ls = linestyle;
	}
	public PScatterLayer(double[] x, double[] y, int colour, double linewidth, String linestyle) {
		xarrayx = x;
		yarrayy = y;
		minX = PPlotMath.dmin(xarrayx);
		maxX = PPlotMath.dmax(xarrayx);
		minY = PPlotMath.dmin(yarrayy);
		maxY = PPlotMath.dmax(yarrayy);
		col = colour;
		lw = linewidth;
		ls = linestyle;
	}
	

	@Override
	public void createRasterImg(PPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(PAxis ax, int layernum, PGroupShape s) {
		int[] p = ax.getSize();
		if("o".equals(ls) || "@".equals(ls)) {
			PPlotShape.fill(col); } else { PPlotShape.noFill(); }
		PPlotShape.stroke(col); PPlotShape.strokeWeight((float)lw);
		PGroupShape xyShape = new PGroupShape();
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
