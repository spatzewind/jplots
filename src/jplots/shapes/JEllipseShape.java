package jplots.shapes;

import jplots.JPlot;
import processing.core.PConstants;
import processing.core.PGraphics;

public class JEllipseShape extends JPlotShape {

	private boolean isFilled;
	private int inCol, outCol;
	private float xtl, ytl, xbr, ybr, sw;

	public JEllipseShape(float cx, float cy, float xw, float yh) {
		this(cx, cy, xw, yh, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight,
				JPlotShape.useFill);
	}

	public JEllipseShape(float cx, float cy, float xw, float yh, boolean filled) {
		this(cx, cy, xw, yh, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight, filled);
	}

	public JEllipseShape(float cx, float cy, float xw, float yh, int colour) {
		this(cx, cy, xw, yh, colour, colour, JPlotShape.strokeWeight, JPlotShape.useFill);
	}

	public JEllipseShape(float cx, float cy, float xw, float yh, int inner_colour, int outer_colour) {
		this(cx, cy, xw, yh, inner_colour, outer_colour, JPlotShape.strokeWeight, true);
	}

	public JEllipseShape(float cx, float cy, float xw, float yh, int inner_colour, int outer_colour,
			float stroke_weight) {
		this(cx, cy, xw, yh, inner_colour, outer_colour, stroke_weight, true);
	}

	public JEllipseShape(float cx, float cy, float xw, float yh, int colour, boolean filled) {
		this(cx, cy, xw, yh, colour, colour, JPlotShape.strokeWeight, filled);
	}

	public JEllipseShape(float cx, float cy, float xw, float yh, int inner_colour, int outer_colour,
			float stroke_weight, boolean filled) {
		xtl = cx - 0.5f * Math.abs(xw);
		ytl = cy - 0.5f * Math.abs(yh);
		xbr = cx + 0.5f * Math.abs(xw);
		ybr = cy + 0.5f * Math.abs(yh);
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
	}

	@Override
	public void draw(JPlot plot, PGraphics g) {
		if (isFilled) {
			g.fill(inCol);
		} else {
			g.noFill();
		}
		g.stroke(outCol);
		g.strokeWeight(sw);
		switch (g.ellipseMode) {
		case PConstants.RADIUS:
			g.ellipse(0.5f * (xtl + xbr), 0.5f * (ytl + ybr), 0.5f * (xbr - xtl), 0.5f * (ybr - ytl));
			break;
		case PConstants.CENTER:
			g.ellipse(0.5f * (xtl + xbr), 0.5f * (ytl + ybr), xbr - xtl, ybr - ytl);
			break;
		case PConstants.CORNER:
			g.ellipse(xtl, ytl, xbr - xtl, ybr - ytl);
			break;
		case PConstants.CORNERS:
			g.ellipse(xtl, ytl, xbr, ybr);
			break;
		default:
			break;
		}
	}
	
	
	public float[] getCenter() {
		return new float[] {0.5f * (xtl + xbr), 0.5f * (ytl + ybr)};
	}
	public float[] getRadii() {
		return new float[] {0.5f * (xbr - xtl), 0.5f * (ybr - ytl)};
	}
}
