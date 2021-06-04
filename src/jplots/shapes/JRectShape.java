package jplots.shapes;

import jplots.JPlot;
import processing.core.PApplet;
import processing.core.PGraphics;

public class JRectShape extends JPlotShape {

	private boolean isFilled;
	private int inCol, outCol;
	private float xtl, ytl, xbr, ybr, sw;

	public JRectShape(float x1, float y1, float x2, float y2) {
		this(x1, y1, x2, y2, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight, JPlotShape.useFill); }
	public JRectShape(float x1, float y1, float x2, float y2, boolean filled) {
		this(x1, y1, x2, y2, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight, filled); }
	public JRectShape(float x1, float y1, float x2, float y2, int colour) {
		this(x1, y1, x2, y2, colour, colour, JPlotShape.strokeWeight, JPlotShape.useFill); }
	public JRectShape(float x1, float y1, float x2, float y2, int inner_colour, int outer_colour) {
		this(x1, y1, x2, y2, inner_colour, outer_colour, JPlotShape.strokeWeight, true); }
	public JRectShape(float x1, float y1, float x2, float y2, int inner_colour, int outer_colour, float stroke_weight) {
		this(x1, y1, x2, y2, inner_colour, outer_colour, stroke_weight, true); }
	public JRectShape(float x1, float y1, float x2, float y2, int colour, boolean filled) {
		this(x1, y1, x2, y2, colour, colour, JPlotShape.strokeWeight, filled); }
	public JRectShape(float x1, float y1, float x2, float y2, int inner_colour, int outer_colour, float stroke_weight, boolean filled) {
		xtl = Math.min(x1, x2);
		ytl = Math.min(y1, y2);
		xbr = Math.max(x1, x2);
		ybr = Math.max(y1, y2);
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
	}
	
	@Override
	public void draw(JPlot plot, PGraphics g) {
		if(isFilled) {
			g.fill(inCol);
		} else {
			g.noFill();
		}
		g.stroke(outCol);
		g.strokeWeight(sw);
		switch(g.rectMode) {
			case PApplet.CENTER:  g.rect(0.5f*(xtl+xbr), 0.5f*(ytl+ybr), xbr-xtl, ybr-ytl); break;
			case PApplet.CORNER:  g.rect(xtl, ytl, xbr-xtl, ybr-ytl); break;
			case PApplet.CORNERS: g.rect(xtl, ytl, xbr, ybr); break;
			default: break;
		}
	}
}
