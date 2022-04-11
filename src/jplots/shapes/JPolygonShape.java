package jplots.shapes;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;

import jplots.JPlot;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import processing.core.PApplet;
import processing.core.PGraphics;

public class JPolygonShape extends JPlotShape {

	private boolean isFilled, isStroked;
	private int inCol, outCol;
	private float sw;
	private float[] xx, yy;
	
	public JPolygonShape(float[]... coords) {
		this(coords, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight, JPlotShape.useFill, JPlotShape.useStroke); }
	public JPolygonShape(float[][] coords, boolean filled, boolean withOutline) {
		this(coords, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight, filled, withOutline); }
	public JPolygonShape(float[][] coords, int colour) {
		this(coords, colour, colour, JPlotShape.strokeWeight, JPlotShape.useFill, JPlotShape.useStroke); }
	public JPolygonShape(float[][] coords, int inner_colour, int outer_colour) {
		this(coords, inner_colour, outer_colour, JPlotShape.strokeWeight, true, true); }
	public JPolygonShape(float[][] coords, int inner_colour, int outer_colour, float stroke_weight) {
		this(coords, inner_colour, outer_colour, stroke_weight, true, true); }
	public JPolygonShape(float[][] coords, int colour, boolean filled, boolean withOutline) {
		this(coords, colour, colour, JPlotShape.strokeWeight, filled, withOutline); }
	public JPolygonShape(float[][] coords, int inner_colour, int outer_colour, float stroke_weight, boolean filled, boolean withOutline) {
		if(coords==null) {
			System.err.println("[JDPolygonShape] invalid geometry!");
			xx = null;
			yy = null;
			return;
		}
		int cc = coords.length;
		xx = new float[cc];
		yy = new float[cc];
		for(int c=0; c<cc; c++) {
			xx[c] = coords[c][0];
			yy[c] = coords[c][1];
		}
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
		isStroked = withOutline;
	}
	public JPolygonShape(Geometry geom, int inner_colour, int outer_colour, float stroke_weight, boolean filled, boolean withOutline) {
		if(geom==null) {
			System.err.println("[JDPolygonShape] invalid geometry!");
			xx = null;
			yy = null;
			return;
		}
		int cc = geom.getCoordinates().length;
		xx = new float[cc];
		yy = new float[cc];
		for(int c=0; c<cc; c++) {
			Coordinate coords = geom.getCoordinates()[c];
			xx[c] = (float) coords.x;
			yy[c] = (float) coords.y;
		}
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
		isStroked = withOutline;
	}
	public JPolygonShape(JDPolygon geom, int inner_colour, int outer_colour, float stroke_weight, boolean filled, boolean withOutline) {
		if(geom==null) {
			System.err.println("[JDPolygonShape] invalid geometry!");
			xx = null;
			yy = null;
			return;
		}
		int cc = geom.c.length;
		xx = new float[cc];
		yy = new float[cc];
		for(int c=0; c<cc; c++) {
			JDPoint coords = geom.c[c];
			xx[c] = (float) coords.x;
			yy[c] = (float) coords.y;
		}
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
		isStroked = withOutline;
	}
	
	@Override
	public void draw(JPlot plot, PGraphics g) {
		if(xx==null || yy==null)
			return;
		float a = calcArea();
		if(Float.isNaN(a) || a<0.0001d)
			return;
//		System.out.println("[JPolyg.Shape] draw polygon shape ("+
//				(isFilled?"fill":"nofill")+","+
//				(isStroked?"stroke":"nostroke")+","+
//				+a+")");
//		for(int c=0; c<Math.min(10, xx.length); c++)
//			System.out.println("[JPolyg.Shape]     ["+xx[c]+", "+yy[c]+"]");
//		if(xx.length>10)
//			System.out.println("[JPolyg.Shape]     ... (and "+(xx.length-10)+" more)");
		if(isFilled) {
			g.fill(inCol);
		} else {
			g.noFill();
		}
		if(isStroked) {
			g.stroke(outCol);
			g.strokeWeight(sw);
		} else {
			g.noStroke();
		}
		g.beginShape();
		for(int c=0; c<xx.length; c++)
			g.vertex(xx[c], yy[c]);
		g.endShape(PApplet.CLOSE);
	}
	private float calcArea() {
		float a = 0f;
		for(int i=0,j=xx.length-1; i<xx.length; j=i++) {
			a += (xx[i]-xx[j]) * (yy[i]+yy[j]);
		}
		return 0.5f * a;
	}
}
