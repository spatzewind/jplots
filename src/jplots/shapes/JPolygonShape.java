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
}
