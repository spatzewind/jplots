package jplots.shapes;

import jplots.JPlot;
import processing.core.PGraphics;

public class JLineShape extends JPlotShape {
	
	private int col;
	private float xstart,ystart,xend,yend;
	private float sw;

	public JLineShape(float x1, float y1, float x2, float y2) {
		this(x1, y1, x2, y2, JPlotShape.strokeColour, JPlotShape.strokeWeight); }
	public JLineShape(float x1, float y1, float x2, float y2, float stroke_weight) {
		this(x1,y1,x2,y2,JPlotShape.strokeColour,stroke_weight); }
	public JLineShape(float x1, float y1, float x2, float y2, int colour) {
		this(x1, y1, x2, y2, colour, JPlotShape.strokeWeight); }
	public JLineShape(float x1, float y1, float x2, float y2, int colour, float stroke_weight) {
		xstart = x1;
		ystart = y1;
		xend = x2;
		yend = y2;
		col = colour;
		sw = stroke_weight;
	}

	@Override
	public void draw(JPlot plot, PGraphics g) {
		g.noFill(); g.stroke(col); g.strokeWeight(sw);
		g.line(xstart, ystart, xend, yend);
	}

}
