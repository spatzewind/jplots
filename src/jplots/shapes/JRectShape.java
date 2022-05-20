package jplots.shapes;

import jplots.JPlot;
import processing.core.PConstants;
import processing.core.PGraphics;

public class JRectShape extends JPlotShape {

	private boolean isFilled;
	private int inCol, outCol;
	private float xtl, ytl, xbr, ybr, sw, cr;

	public JRectShape(float x1, float y1, float x2, float y2) {
		this(x1, y1, x2, y2, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight,
				JPlotShape.useFill);
	}

	public JRectShape(float x1, float y1, float x2, float y2, boolean filled) {
		this(x1, y1, x2, y2, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight, filled);
	}

	public JRectShape(float x1, float y1, float x2, float y2, int colour) {
		this(x1, y1, x2, y2, colour, colour, JPlotShape.strokeWeight, JPlotShape.useFill);
	}

	public JRectShape(float x1, float y1, float x2, float y2, int inner_colour, int outer_colour) {
		this(x1, y1, x2, y2, inner_colour, outer_colour, JPlotShape.strokeWeight, true);
	}

	public JRectShape(float x1, float y1, float x2, float y2, int inner_colour, int outer_colour, float stroke_weight) {
		this(x1, y1, x2, y2, inner_colour, outer_colour, stroke_weight, true);
	}

	public JRectShape(float x1, float y1, float x2, float y2, int colour, boolean filled) {
		this(x1, y1, x2, y2, colour, colour, JPlotShape.strokeWeight, filled);
	}

	public JRectShape(float x1, float y1, float x2, float y2, int inner_colour, int outer_colour, float stroke_weight,
			boolean filled) {
		xtl = Math.min(x1, x2);
		ytl = Math.min(y1, y2);
		xbr = Math.max(x1, x2);
		ybr = Math.max(y1, y2);
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		cr = -1f;
		isFilled = filled;
	}

	public JRectShape(float x1, float y1, float x2, float y2, float corner, int inner_colour, int outer_colour,
			float stroke_weight, boolean filled) {
		xtl = Math.min(x1, x2);
		ytl = Math.min(y1, y2);
		xbr = Math.max(x1, x2);
		ybr = Math.max(y1, y2);
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
		cr = corner <= 0f ? -1f : corner;
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
		switch (g.rectMode) {
		case PConstants.CENTER:
			if (cr > 0f)
				g.rect(0.5f * (xtl + xbr), 0.5f * (ytl + ybr), xbr - xtl, ybr - ytl, cr);
			else
				g.rect(0.5f * (xtl + xbr), 0.5f * (ytl + ybr), xbr - xtl, ybr - ytl);
			break;
		case PConstants.CORNER:
			if (cr > 0f)
				g.rect(xtl, ytl, xbr - xtl, ybr - ytl, cr);
			else
				g.rect(xtl, ytl, xbr - xtl, ybr - ytl);
			break;
		case PConstants.CORNERS:
			if (cr > 0f)
				g.rect(xtl, ytl, xbr, ybr, cr);
			else
				g.rect(xtl, ytl, xbr, ybr);
			break;
		default:
			break;
		}
	}
}
