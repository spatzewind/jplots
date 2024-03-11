package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;

import jplots.JPlot;
import jplots.axes.JAxis;
import jplots.axes.LogarithmicScale;
import jplots.helper.FileLoader;
import jplots.maths.AffineBuilder;
import jplots.maths.JDLine;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import jplots.shapes.JEllipseShape;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JPolygonShape;
import processing.core.PGraphics;

public class JShapesLayer extends JPlotsLayer {
	
	private String shapefileKey;
	private String shapeType;
	private List<double[][]> coordsSave;

	public JShapesLayer(String resource, String type) {
		super();
		pc = 0xff99ccff;
		lc = 0xff336699;
		shapefileKey = resource;
		shapeType = type;
		coordsSave = new ArrayList<>();
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		//TODO: deal with JMultiAxis
		boolean isXlog = (ax.getScaleX() instanceof LogarithmicScale);
		boolean isYlog = (ax.getScaleY() instanceof LogarithmicScale);
		double Xin = isXlog ? Math.log10(minX) : minX, Xax = isXlog ? Math.log10(maxX) : maxX;
		double Yin = isYlog ? Math.log10(minY) : minY, Yax = isYlog ? Math.log10(maxY) : maxY;
		double xs = p[2] / (Xax - Xin), ys = p[3] / (Yax - Yin);
		AffineBuilder affine = new AffineBuilder()
				.scale(invertAxisX?-1d:1d, invertAxisY?1d:-1d).translate(invertAxisX?Xax:-Xin, invertAxisY?-Yin:Yax)
				.scale(xs, ys).translate(p[0], p[1]);
		double[][] affmat = affine.getMatrix();
//		boolean debug = ax.getPlot().isDebug();
//		if(debug)
//			System.out.println("[DEBUG] JShapeLayer: begin reading shape file \""+connect.get("url")+"\"");
		int fidx = 0;
		for (Coordinate[] points : FileLoader.shapefiles.get(shapefileKey)) {
			int ppcc = pc, llcc = lc;
			if (!singleFillColour) ppcc = pcs[fidx % pcs.length];
			if (!singleLineColour) llcc = lcs[fidx % lcs.length];
			JDPoint[] coords = new JDPoint[points.length];
			for (int c = 0; c < points.length; c++) {
				double[] xy = { points[c].x,  points[c].y };
				if (ax.isGeoAxis()) xy = ax.getGeoProjection().fromLATLONtoPROJ(xy[1], xy[0], true, false);
				if (isXlog) xy[0] = Math.log10(points[c].x);
				if (isYlog) xy[1] = Math.log10(points[c].y);
				coords[c] = new JDPoint(xy[0], xy[1]);
			}
//			coordsSave.add(toCoords(coords));
			if (shapeType.equalsIgnoreCase("point")) {
				JPlotShape.noFill();
				JPlotShape.stroke(llcc);
				JPlotShape.strokeWeight((float) lw);
				for (int c = 0; c < coords.length; c++) {
					coords[c].affine(affmat);
					double x1 = coords[c].x, y1 = coords[c].y;
					if (x1 < p[0] || x1 >= p[0] + p[2] || y1 < p[1] || y1 >= p[1] + p[3])
						continue;
					s.addChild(new JEllipseShape((float) x1, (float) y1, 1f, 1f));
				}
			}
			if (shapeType.equalsIgnoreCase("line")) {
				JDLine[] lineArr = null;
				if(ax.isGeoAxis())
					lineArr = ax.getGeoProjection().splitByMapBorder(new JDLine(coords)).toArray(new JDLine[0]);
				else
					lineArr = new JDLine[] { new JDLine(coords) };
				for(JDLine la: lineArr)
					for(JDLine l: la.affine(affmat).intersectsAABB(p[0], p[0]+p[2], p[1], p[1]+p[3]))
						s.addChild(new JLineShape((float)lw, llcc, l.getCoordsAsFloats()));
			}
			if (shapeType.equalsIgnoreCase("polygon")) {
				JDPolygon[] polyArr = null;
				if(ax.isGeoAxis())
					polyArr = ax.getGeoProjection().splitByMapBorder(new JDPolygon(coords)).toArray(new JDPolygon[0]);
				else
					polyArr = new JDPolygon[] { new JDPolygon(coords) };
				for(JDPolygon pa: polyArr)
					for(JDPolygon pl: pa.affine(affmat).intersectsAABB(p[0], p[0]+p[2], p[1], p[1]+p[3]))
						s.addChild(new JPolygonShape(pl.getCoords(), ppcc, llcc, (float)lw, true, drawLines));
			}
			fidx++;
		}
	}
	
	public List<double[][]> getCoords() {
		return coordsSave;
	}
}
