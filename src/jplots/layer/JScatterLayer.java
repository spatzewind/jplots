package jplots.layer;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import jplots.JPlot;
import jplots.axes.JAxis;
import jplots.helper.GeometryTools;
import jplots.maths.JDPolygon;
import jplots.maths.JPlotMath;
import jplots.shapes.JEllipseShape;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JPolygonShape;
import jplots.shapes.JRectShape;
import processing.core.PGraphics;

public class JScatterLayer extends JPlotsLayer {

	private double[] xarrayx, yarrayy;
	private int col;
	private double lw;
	private String ls;
	private boolean points2polygon;

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
		if ("o".equalsIgnoreCase(ls) || "@".equals(ls)) {
			JPlotShape.fill(col);
		} else {
			JPlotShape.noFill();
		}
		JPlotShape.stroke(col);
		JPlotShape.strokeWeight((float) lw);
		JGroupShape xyShape = new JGroupShape();
		List<double[]> centers = new ArrayList<>();
		double cx = p[0]+0.5*p[2], cy = p[1]+0.5*p[3];
		double xs = p[2] / (Xax - Xin), ys = p[3] / (Yax - Yin);
		for (int i = 0; i < xarrayx.length; i++) {
			if (!Double.isFinite(xarrayx[i]) || !Double.isFinite(yarrayy[i]))
				continue;
			double[] xy = inputProj.fromPROJtoLATLON(ax.isXlogAxis() ? Math.log10(xarrayx[i]) : xarrayx[i],
					ax.isYlogAxis() ? Math.log10(yarrayy[i]) : yarrayy[i], false, true);
			if (ax.isGeoAxis())
				xy = ax.getGeoProjection().fromLATLONtoPROJ(xy[0], xy[1], false, true);
			double x1 = p[0] + xs * (invertAxisX ? Xax - xy[0] : xy[0] - Xin);
			double y1 = p[1] + ys * (invertAxisY ? xy[1] - Yin : Yax - xy[1]);
			if (x1 < p[0] || x1 > p[0] + p[2] || y1 < p[1] || y1 > p[1] + p[3])
				continue;
			centers.add(new double[] {x1, y1, (x1-cx)*(x1-cx)+(y1-cy)*(y1-cy)});
		}
		centers.sort(new Comparator<double[]>() {
			@Override
			public int compare(double[] o1, double[] o2) {
				return Double.compare(o1[2],o2[2]);
			}
		});
		if(points2polygon && ("o".equalsIgnoreCase(ls) || "@".equals(ls))) {
			List<JDPolygon> polys = new ArrayList<>();
			boolean q = "o".equalsIgnoreCase(ls);
			double r = q ? 6*lw : 2*lw;
			for(double[] c: centers)
				polys.add(q ? JDPolygon.circle(c[0], c[1], r) :
						new JDPolygon(new double[][] {{c[0]-r,c[1]-r},{c[0]+r,c[1]-r},{c[0]+r,c[1]+r},{c[0]-r,c[1]+r}})
						);
			centers.clear();
			int pl = polys.size()+1;
			boolean debug = ax.getPlot().isDebug();
			if(debug) System.out.println("[----+----|----+----|----+----|----+----|----+----|] ("+polys.size()+")");
			while(polys.size()<pl) {
				pl = polys.size();
				if(debug) System.out.print("["); int perc = 0;
				for(int i=pl-1; i>0; i-=2) {
					for(int j=Math.max(0,i-100); j<i; j++) {
						JDPolygon temp = GeometryTools.union(polys.get(i), polys.get(j), 0.0001d);
						if(temp==null) continue;
						polys.set(j, temp);
						polys.remove(i);
						break;
					}
					int crep = (pl+1+pl%2-i) * 50 / pl;
					if(debug && crep>perc)
						System.out.print("##################################################".substring(0, crep-perc));
					perc = crep;
				}
				if(debug) System.out.println("] ("+polys.size()+")");
			}
			for(JDPolygon poly: polys)
				xyShape.addChild(new JPolygonShape(poly, col, col, (float)lw, true, true));
		} else {
			for(double[] c: centers) {
				xyShape.addChild(getShape((float)c[0], (float)c[1], ls, col));
			}
		}
		s.addChild(xyShape);
	}
	
	public void setPoints2Polygon(boolean b) {
		points2polygon = b;
	}
	
	private JPlotShape getShape(float x, float y, String style, int col) {
		//size
		float r = 3*(float)lw;
		if("()".equals(style) || "o".equalsIgnoreCase(style)) r = 6*(float)lw;
		if("[]".equals(style) || "@".equals(style) || "x".equalsIgnoreCase(style)) r = 2*(float)lw;
		//shapes
		if("()".equals(style) || "(F)".equalsIgnoreCase(style) || "o".equalsIgnoreCase(style))
			return new JEllipseShape(x,y, r, r, !"()".equals(style));
		if("[]".equals(style) || "[F]".equalsIgnoreCase(style) || "@".equals(style))
			return new JRectShape(x-r, y-r, x+r, y+r, !"[]".equals(style));
		if("<>".equals(style))
			return new JLineShape((float)lw, col, x,y-r, x+r,y, x,y+r, x-r,y, x,y-r);
		if("<F>".equalsIgnoreCase(style) || "%".equals(style))
			return new JPolygonShape(new float[][] {{x,y-r},{x+r,y},{x,y+r},{x-r,y}}, true,true);
		if("v".equalsIgnoreCase(style))
			return new JPolygonShape(new float[][] {{x-0.866f*r,y-0.5f*r},{x+0.866f*r,y-0.5f*r},{x,y+r}}, !"v".equals(style),true);
		if("^".equals(style) || "^F".equalsIgnoreCase(style))
			return new JPolygonShape(new float[][] {{x,y-r},{x+0.866f*r,y+0.5f*r},{x-0.866f*r,y+0.5f*r}}, "^F".equalsIgnoreCase(style),true);
		JGroupShape gs = new JGroupShape();
		if("x".equalsIgnoreCase(ls)) {
			gs.addChild(new JLineShape((float)lw, col, x-r, y-r, x+r, y+r));
			gs.addChild(new JLineShape((float)lw, col, x-r, y+r, x+r, y-r));
			return gs;
		}
		if("+".equals(ls)) {
			gs.addChild(new JLineShape((float)lw, col, x-r, y, x+r, y));
			gs.addChild(new JLineShape((float)lw, col, x, y-r, x, y+r));
			return gs;
		}
		return gs;
	}
}
