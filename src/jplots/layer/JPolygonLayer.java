package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import jplots.JPlot;
import jplots.axes.JAxis;
import jplots.maths.AffineBuilder;
import jplots.maths.JDLine;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JPolygonShape;
import processing.core.PGraphics;

public class JPolygonLayer extends JPlotsLayer {
	
	private JDPolygon poly;
	private List<double[][]> coordsSave;

	public JPolygonLayer(JDPolygon polygon, int innerColour, int outerColour, double linesize) {
		super();
		pc = innerColour;
		lc = outerColour;
		lw = linesize;
		poly = polygon;
		coordsSave = new ArrayList<>();
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
		AffineBuilder affine = new AffineBuilder()
				.scale(invertAxisX?-1d:1d, invertAxisY?1d:-1d).translate(invertAxisX?Xax:-Xin, invertAxisY?-Yin:Yax)
				.scale(xs, ys).translate(p[0], p[1]);
		double[][] affmat = affine.getMatrix();
//		boolean debug = ax.getPlot().isDebug();
//		if(debug)
//			System.out.println("[DEBUG] JShapeLayer: begin reading shape file \""+connect.get("url")+"\"");
		JDPoint[] coords = poly.c.clone();
		for (int c = 0; c < coords.length; c++) {
			double[] xy = { coords[c].x,  coords[c].y };
			if (ax.isGeoAxis()) xy = ax.getGeoProjection().fromLATLONtoPROJ(xy[1], xy[0], true, false);
			if (ax.isXlogAxis()) xy[0] = Math.log10(coords[c].x);
			if (ax.isYlogAxis()) xy[1] = Math.log10(coords[c].y);
			coords[c] = new JDPoint(xy[0], xy[1]);
		}
		JPlotShape.fill(pc);
		JPlotShape.stroke(lc);
		JPlotShape.strokeWeight((float) lw);
		JDPolygon[] polyArr = null;
		if(ax.isGeoAxis())
			polyArr = ax.getGeoProjection().splitByMapBorder(new JDPolygon(coords)).toArray(new JDPolygon[0]);
		else
			polyArr = new JDPolygon[] { new JDPolygon(coords) };
		for(JDPolygon pa: polyArr)
			for(JDPolygon pl: pa.affine(affmat).intersectsAABB(p[0], p[0]+p[2], p[1], p[1]+p[3]))
				s.addChild(new JPolygonShape(pl.getCoords(), pc, lc, (float)lw, true, false));
		JDLine[] lineArr = null;
		if(ax.isGeoAxis())
			lineArr = ax.getGeoProjection().splitByMapBorder(new JDLine(coords)).toArray(new JDLine[0]);
		else
			lineArr = new JDLine[] { new JDLine(coords) };
		for(JDLine la: lineArr)
			for(JDLine l: la.affine(affmat).intersectsAABB(p[0], p[0]+p[2], p[1], p[1]+p[3]))
				s.addChild(new JLineShape((float)lw, lc, l.getCoordsAsFloats()));
	}
	
	public List<double[][]> getCoords() {
		return coordsSave;
	}
}
