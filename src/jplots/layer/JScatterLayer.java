package jplots.layer;

import jplots.JPlot;
import jplots.axes.JAxis;
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
		for (int i = 0; i < x.length; i++)
			xarrayx[i] = x[i];
		yarrayy = new double[y.length];
		for (int i = 0; i < y.length; i++)
			yarrayy[i] = y[i];
		minX = JPlotMath.dmin(xarrayx);
		maxX = JPlotMath.dmax(xarrayx);
		minY = JPlotMath.dmin(yarrayy);
		maxY = JPlotMath.dmax(yarrayy);
		setLineColour(colour);
		col = colour;
		lw = linewidth;
		setStyle(linestyle);
		ls = linestyle;
	}

	public JScatterLayer(double[] x, double[] y, int colour, double linewidth, String linestyle) {
		xarrayx = x;
		yarrayy = y;
		minX = JPlotMath.dmin(xarrayx);
		maxX = JPlotMath.dmax(xarrayx);
		minY = JPlotMath.dmin(yarrayy);
		maxY = JPlotMath.dmax(yarrayy);
		setLineColour(colour);
		col = colour;
		lw = linewidth;
		setStyle(linestyle);
		ls = linestyle;
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		double Xin = ax.isXlogAxis() ? Math.log10(minX) : minX, Xax = ax.isXlogAxis() ? Math.log10(maxX) : maxX;
		double Yin = ax.isYlogAxis() ? Math.log10(minY) : minY, Yax = ax.isYlogAxis() ? Math.log10(maxY) : maxY;
		if ("o".equals(ls) || "@".equals(ls)) {
			JPlotShape.fill(col);
		} else {
			JPlotShape.noFill();
		}
		JPlotShape.stroke(col);
		JPlotShape.strokeWeight((float) lw);
		JGroupShape xyShape = new JGroupShape();
		double xs = p[2] / (Xax - Xin), ys = p[3] / (Yax - Yin);
		for (int i = 0; i < xarrayx.length; i++)
			if (Double.isFinite(xarrayx[i]) && Double.isFinite(yarrayy[i])) {
				double[] xy = inputProj.fromPROJtoLATLON(ax.isXlogAxis() ? Math.log10(xarrayx[i]) : xarrayx[i],
						ax.isYlogAxis() ? Math.log10(yarrayy[i]) : yarrayy[i], false, true);
				if (ax.isGeoAxis())
					xy = ax.getGeoProjection().fromLATLONtoPROJ(xy[0], xy[1], false, true);
				double x1 = p[0] + xs * (invertAxisX ? Xax - xy[0] : xy[0] - Xin);
				double y1 = p[1] + ys * (invertAxisY ? xy[1] - Yin : Yax - xy[1]);
				if (x1 < p[0] || x1 > p[0] + p[2] || y1 < p[1] || y1 > p[1] + p[3])
					continue;
				if ("()".equals(ls) || "o".equals(ls)) {
					xyShape.addChild(new JEllipseShape((float) x1, (float) y1, (float) (6 * lw), (float) (6 * lw)));
				}
				if ("[]".equals(ls) || "@".equals(ls)) {
					xyShape.addChild(new JRectShape((float) (x1 - 2 * lw), (float) (y1 - 2 * lw), (float) (x1 + 2 * lw),
							(float) (y1 + 2 * lw)));
				}
				if ("x".equals(ls)) {
					xyShape.addChild(new JLineShape((float)lw, col, (float) (x1 - 2 * lw), (float) (y1 - 2 * lw), (float) (x1 + 2 * lw),
							(float) (y1 + 2 * lw)));
					xyShape.addChild(new JLineShape((float)lw, col, (float) (x1 - 2 * lw), (float) (y1 + 2 * lw), (float) (x1 + 2 * lw),
							(float) (y1 - 2 * lw)));
				}
				if ("+".equals(ls)) {
					xyShape.addChild(
							new JLineShape((float)lw, col, (float) (x1 - 3 * lw), (float) y1, (float) (x1 + 3 * lw), (float) y1));
					xyShape.addChild(
							new JLineShape((float)lw, col, (float) x1, (float) (y1 - 3 * lw), (float) x1, (float) (y1 + 3 * lw)));
				}
				if ("<>".equals(ls)) {
					xyShape.addChild(
							new JLineShape((float)lw, col, (float) x1, (float) (y1 - 3 * lw), (float) (x1 + 3 * lw), (float) y1));
					xyShape.addChild(
							new JLineShape((float)lw, col, (float) (x1 + 3 * lw), (float) y1, (float) x1, (float) (y1 + 3 * lw)));
					xyShape.addChild(
							new JLineShape((float)lw, col, (float) x1, (float) (y1 + 3 * lw), (float) (x1 - 3 * lw), (float) y1));
					xyShape.addChild(
							new JLineShape((float)lw, col, (float) (x1 - 3 * lw), (float) y1, (float) x1, (float) (y1 - 3 * lw)));
				}
			}
		s.addChild(xyShape);
	}

}
