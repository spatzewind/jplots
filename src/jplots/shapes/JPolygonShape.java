package jplots.shapes;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;

import jplots.JPlot;
import jplots.helper.GeometryTools;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import processing.core.PConstants;
import processing.core.PGraphics;
import processing.core.PShape;

public class JPolygonShape extends JPlotShape {

	private boolean isFilled, isStroked;
	private int inCol, outCol;
	private float sw;
	private float[] xx, yy;
	private boolean[] isInner;

	public JPolygonShape(float[]... coords) {
		this(coords, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight, JPlotShape.useFill, JPlotShape.useStroke);
	}

	public JPolygonShape(float[][] coords, boolean filled, boolean withOutline) {
		this(coords, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight, filled, withOutline);
	}

	public JPolygonShape(float[][] coords, int colour) {
		this(coords, colour, colour, JPlotShape.strokeWeight, JPlotShape.useFill, JPlotShape.useStroke);
	}

	public JPolygonShape(float[][] coords, int inner_colour, int outer_colour) {
		this(coords, inner_colour, outer_colour, JPlotShape.strokeWeight, true, true);
	}

	public JPolygonShape(float[][] coords, int inner_colour, int outer_colour, float stroke_weight) {
		this(coords, inner_colour, outer_colour, stroke_weight, true, true);
	}

	public JPolygonShape(float[][] coords, int colour, boolean filled, boolean withOutline) {
		this(coords, colour, colour, JPlotShape.strokeWeight, filled, withOutline);
	}

	public JPolygonShape(float[][] coords, int inner_colour, int outer_colour, float stroke_weight, boolean filled, boolean withOutline) {
		if (coords == null) {
			System.err.println("[JDPolygonShape] invalid geometry!");
			xx = null;
			yy = null;
			return;
		}
		int cc = coords.length;
		xx = new float[cc];
		yy = new float[cc];
		for (int c = 0; c < cc; c++) {
			xx[c] = coords[c][0];
			yy[c] = coords[c][1];
		}
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
		isStroked = withOutline;
		if(isStroked) checkInnerEdges();
	}

	public JPolygonShape(Geometry geom, int inner_colour, int outer_colour, float stroke_weight, boolean filled,
			boolean withOutline) {
		if (geom == null) {
			System.err.println("[JDPolygonShape] invalid geometry!");
			xx = null;
			yy = null;
			return;
		}
		int cc = geom.getCoordinates().length;
		xx = new float[cc];
		yy = new float[cc];
		for (int c = 0; c < cc; c++) {
			Coordinate coords = geom.getCoordinates()[c];
			xx[c] = (float) coords.x;
			yy[c] = (float) coords.y;
		}
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
		isStroked = withOutline;
		if(isStroked) checkInnerEdges();
	}

	public JPolygonShape(JDPolygon geom, int inner_colour, int outer_colour, float stroke_weight, boolean filled,
			boolean withOutline) {
		if (geom == null) {
			System.err.println("[JDPolygonShape] invalid geometry!");
			xx = null;
			yy = null;
			return;
		}
		int cc = geom.c.length;
		xx = new float[cc];
		yy = new float[cc];
		for (int c = 0; c < cc; c++) {
			JDPoint coords = geom.c[c];
			xx[c] = (float) coords.x;
			yy[c] = (float) coords.y;
		}
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
		isStroked = withOutline;
		if(isStroked) checkInnerEdges();
	}

	@Override
	public void draw(JPlot plot, PGraphics g) {
		float a = GeometryTools.area(xx, yy);
		if (Float.isNaN(a) || a < 0.0001d) {
//			System.err.println(a);
			return;
		}
//		System.out.println("[JPolyg.Shape] draw polygon shape ("+
//				(isFilled?"fill":"nofill")+","+
//				(isStroked?"stroke":"nostroke")+","+
//				+a+")");
//		for(int c=0; c<Math.min(10, xx.length); c++)
//			System.out.println("[JPolyg.Shape]     ["+xx[c]+", "+yy[c]+"]");
//		if(xx.length>10)
//			System.out.println("[JPolyg.Shape]     ... (and "+(xx.length-10)+" more)");
		if (isFilled) {
			PShape p = g.createShape();
			p.beginShape();
			if(isFilled) p.fill(inCol);
			else p.noFill();
			p.stroke(inCol); p.strokeWeight(1f);
//			String coords = "";
			for(int c=0; c<xx.length; c++) {
				p.vertex(xx[c], yy[c]);
//				coords += ", ["+xx[c]+","+yy[c]+"]";
			}
			p.endShape(PConstants.CLOSE);
//			System.out.println("draw polygon: {"+coords.substring(2)+"}");
			g.shape(p);
		}
		if (isStroked) {
			g.noFill();
			g.stroke(outCol);
			g.strokeWeight(sw);
			for(int c=0; c<xx.length; c++) {
				if(isInner[c]) continue; //skip inner edges
				int d = (c+1) % xx.length;
				g.line(xx[c],yy[c], xx[d],yy[d]);
			}
		}
	}

	private void checkInnerEdges() {
		isInner = new boolean[xx.length];
		for(int i=0; i<xx.length; i++)
			isInner[i] = false;
		for(int a=xx.length-1,b=0; b<xx.length; a=b++) {
			double dx = xx[b]-xx[a];
			double dy = yy[b]-yy[a];
			double dr = 1d / (dx*dx + dy*dy + 1.0e-20);
			dx *= dr;
			dy *= dr;
			for(int c=0; c+1<xx.length; c++) {
				int i = (b+c) % xx.length;
				int j = (i+1) % xx.length;
				double s = Math.abs(dx*(yy[i]-yy[a]) - dy*(xx[i]-xx[a]));
				if(s>0.0001d) continue; //they are not on the same line
				double t = Math.abs(dx*(yy[j]-yy[i]) - dy*(xx[j]-xx[i]));
				if(t>0.0001d) continue; //they are not parallel
				s = dx*(xx[i]-xx[a]) + dy*(yy[i]-yy[a]);
				t = dx*(xx[j]-xx[a]) + dy*(yy[j]-yy[a]);
				if((-0.0001d<s && s<1.0001d) || (-0.0001d<t && t<1.0001d)) {
					isInner[a] = true;
					isInner[i] = true;
				}
			}
		}
	}
}
