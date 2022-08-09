package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import jplots.JAxis;
import jplots.JPlot;
import jplots.maths.AffineBuilder;
import jplots.maths.JDLine;
import jplots.maths.JDPoint;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import processing.core.PGraphics;

public class JXYLayer extends JPlotsLayer {

	private double[] xarrayx, yarrayy;
	private int col;

	public JXYLayer(float[] x, float[] y, int colour, float linewidth, String linestyle) {
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
		maxZ = Double.NaN;
		minZ = Double.NaN;
		setLineColour(colour);
		col = colour;
		lw = linewidth;
		setStyle(linestyle);
		ls = linestyle;
	}
	
	public JXYLayer(double[] x, double[] y, int colour, double linewidth, String linestyle) {
		xarrayx = x;
		yarrayy = y;
		minX = JPlotMath.dmin(xarrayx);
		maxX = JPlotMath.dmax(xarrayx);
		minY = JPlotMath.dmin(yarrayy);
		maxY = JPlotMath.dmax(yarrayy);
		maxZ = Double.NaN;
		minZ = Double.NaN;
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
		double xs = p[2] / (Xax - Xin), ys = p[3] / (Yax - Yin);
		AffineBuilder affine = new AffineBuilder().scale(invertAxisX ? -1d : 1d, invertAxisY ? 1d : -1d)
				.translate(invertAxisX ? Xax : -Xin, invertAxisY ? -Yin : Yax).scale(xs, ys).translate(p[0], p[1]);
		double[][] affmat = affine.getMatrix();
		
		JDPoint[] points = new JDPoint[xarrayx.length];
		for(int i=0; i<xarrayx.length; i++) {
			double x = ax.isXlogAxis() ? Math.log10(xarrayx[i]) : xarrayx[i];
			double y = ax.isYlogAxis() ? Math.log10(yarrayy[i]) : yarrayy[i];
			double[] xy = {x,y};
			if (ax.isGeoAxis()) {
				xy = inputProj.fromPROJtoLATLON(xy[0], xy[1], false, false);
				xy = ax.getGeoProjection().fromLATLONtoPROJ(xy[0], xy[1], false, false);
			}
			points[i] = new JDPoint(xy[0],xy[1]);
		}
		JDLine line = new JDLine(points);
		JGroupShape xyShape = new JGroupShape();;
		double[] pa = null;
		if (parallelArray != null) {
			pa = JPlotMath.toDoubleArray1D(parallelArray);
		}
		if(pa!=null) {
			int[] cols = new int[pa.length];
			if(pa.length==1) {
				cols[0] = colourtable.getColour(0.5d);
			} else {
				double pmin = JPlotMath.dmin(pa), pmax = JPlotMath.dmax(pa);
				if(!Double.isNaN(minZ)) pmin = minZ;
				if(!Double.isNaN(maxZ)) pmax = maxZ;
				for(int i = 0; i < cols.length; i++)
					cols[i] = colourtable.getColour(pa[i], pmin, pmax);
			}
			if(cols.length>=points.length) {
				line.addVertexValues(cols);
			} else {
				line.addLineValues(cols);
			}
		}
		JDLine[] lineArr = new JDLine[0];
		if(ax.isGeoAxis()) {
			lineArr = ax.getGeoProjection().splitByMapBorder(line).toArray(new JDLine[0]);
		} else {
			lineArr = new JDLine[] {line};
		}
		List<JDLine> lines = new ArrayList<>();
		for(JDLine l: lineArr) {
			l.affine(affmat);
			lines.addAll(l.intersectsAABB(p[0],p[0]+p[2], p[1],p[1]+p[3]));
		}
		
		double lln = 1d, llf = 0d, lpn = 0d, lpf = 0d;
		if ("-".equals(ls)) { lln = 1000 * lw; llf = 0; lpn = 0; lpf = 0; }
		if (".".equals(ls)) { lln = 0; llf = 0; lpn = 1 * lw; lpf = 3 * lw; }
		if (",".equals(ls)) { lln = 8 * lw; llf = 7 * lw; lpn = 0; lpf = 0; }
		if (";".equals(ls)) { lln = 8 * lw; llf = 3 * lw; lpn = 1 * lw; lpf = 3 * lw; }
		
		JPlotShape.stroke(lc);
		JPlotShape.strokeWeight((float)lw);
		for(JDLine l: lines)
			drawSingleLine(xyShape, l, lln, llf, lpn, lpf, pa!=null);
		
		s.addChild(xyShape);
	}
	
	public double[] getZRange() {
		return new double[] {minZ, maxZ};
	}
	
	private void drawSingleLine(JGroupShape linesh, JDLine line, double lln, double llf, double lpn, double lpf, boolean useColorGradient) {
		if(line.getPoints().length<2) return;
		int lc = JPlotShape.strokeColour;
		float lw = JPlotShape.strokeWeight;
		if(lln>999d*lw) { //continuous line
			if(useColorGradient) {
				float[] ff = line.getCoordsAsFloats();
				int[] cc = line.getLineColors();
				for(int i=0; i<cc.length; i++)
					linesh.addChild(new JLineShape(lw, cc[i], ff[2*i], ff[2*i+1], ff[2*i+2], ff[2*i+3]));
			} else {
				linesh.addChild(new JLineShape(lw, lc, line.getCoordsAsFloats()));
			}
			return;
		}
		JDPoint[] pnts = line.getPoints();
		int[] cols = line.allColors();
		int li = 0;
		double lpos = 0d, ldif = 0d, loff = 0d;
		for (int i=0,j=1; j<pnts.length; i=j++) {
			JDPoint lvs = pnts[i]; int cs = cols[2*i+1];
			JDPoint lve = pnts[j]; int ce = cols[2*j];
			double x1 = lvs.x, x2 = lve.x;
			double y1 = lvs.y, y2 = lve.y;
			double dx = x2 - x1, dy = y2 - y1;
			double l = Math.sqrt(dx * dx + dy * dy);
			dx /= l; dy /= l;
			lpos = 0d;
			while (lpos < l) {
				ldif = 0d;
				switch (li) {
					case 0:
						if (lln == 0d) break;
						ldif = Math.min(l - lpos, lln - loff);
						break;
					case 1:
						if (llf == 0d)
							break;
						ldif = Math.min(l - lpos, llf - loff);
						break;
					case 2:
						if (lpn == 0d)
							break;
						ldif = Math.min(l - lpos, lpn - loff);
						break;
					case 3:
						if (lpf == 0d)
							break;
						ldif = Math.min(l - lpos, lpf - loff);
						break;
				}
				float	xf1 = (float) (x1 + lpos * dx), yf1 = (float) (y1 + lpos * dy),
						xf2 = (float) (x1 + (lpos + ldif) * dx), yf2 = (float) (y1 + (lpos + ldif) * dy);
				if (li%2 == 0 && ldif >= 0d) {
					int col = lc;
					if(useColorGradient)
						col = JPlotMath.colorLerp(cs, ce, (lpos+0.5d*ldif)/l);
					linesh.addChild(new JLineShape(lw, col, xf1, yf1, xf2, yf2));
				}
				loff += ldif;
				switch (li) {
					case 0: if(loff >= lln) { loff -= lln; li = 1; } break;
					case 1: if(loff >= llf) { loff -= llf; li = 2; } break;
					case 2: if(loff >= lpn) { loff -= lpn; li = 3; } break;
					case 3: if(loff >= lpf) { loff -= lpf; li = 0; } break;
				}
				lpos += ldif;
			}
		}
	}
}
