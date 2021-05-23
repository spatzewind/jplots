package pplots.gfx;

import processing.core.PApplet;
import processing.core.PGraphics;

public class PEllipseShape extends PPlotShape {

	private boolean isFilled;
	private int inCol, outCol;
	private float xtl, ytl, xbr, ybr, sw;

	public PEllipseShape(float cx, float cy, float xw, float yh) {
		this(cx, cy, xw, yh, PPlotShape.fillColour, PPlotShape.strokeColour, PPlotShape.strokeWeight, PPlotShape.useFill); }
	public PEllipseShape(float cx, float cy, float xw, float yh, boolean filled) {
		this(cx, cy, xw, yh, PPlotShape.fillColour, PPlotShape.strokeColour, PPlotShape.strokeWeight, filled); }
	public PEllipseShape(float cx, float cy, float xw, float yh, int colour) {
		this(cx, cy, xw, yh, colour, colour, PPlotShape.strokeWeight, PPlotShape.useFill); }
	public PEllipseShape(float cx, float cy, float xw, float yh, int inner_colour, int outer_colour) {
		this(cx, cy, xw, yh, inner_colour, outer_colour, PPlotShape.strokeWeight, true); }
	public PEllipseShape(float cx, float cy, float xw, float yh, int inner_colour, int outer_colour, float stroke_weight) {
		this(cx, cy, xw, yh, inner_colour, outer_colour, stroke_weight, true); }
	public PEllipseShape(float cx, float cy, float xw, float yh, int colour, boolean filled) {
		this(cx, cy, xw, yh, colour, colour, PPlotShape.strokeWeight, filled); }
	public PEllipseShape(float cx, float cy, float xw, float yh, int inner_colour, int outer_colour, float stroke_weight, boolean filled) {
		xtl = cx-0.5f*Math.abs(xw);
		ytl = cy-0.5f*Math.abs(yh);
		xbr = cx+0.5f*Math.abs(xw);
		ybr = cy+0.5f*Math.abs(yh);
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
	}
	
	@Override
	public void draw(PGraphics g) {
		if(isFilled) {
			g.fill(inCol);
		} else {
			g.noFill();
		}
		g.stroke(outCol);
		g.strokeWeight(sw);
		switch(g.ellipseMode) {
			case PApplet.RADIUS:  g.ellipse(0.5f*(xtl+xbr), 0.5f*(ytl+ybr), 0.5f*(xbr-xtl), 0.5f*(ybr-ytl)); break;
			case PApplet.CENTER:  g.ellipse(0.5f*(xtl+xbr), 0.5f*(ytl+ybr), xbr-xtl, ybr-ytl); break;
			case PApplet.CORNER:  g.ellipse(xtl, ytl, xbr-xtl, ybr-ytl); break;
			case PApplet.CORNERS: g.ellipse(xtl, ytl, xbr, ybr); break;
			default: break;
		}
	}
}
