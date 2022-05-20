package jplots.shapes;

import jplots.JPlot;
import processing.core.PGraphics;

public class JTriangleShape extends JPlotShape {

	private boolean isFilled, isStroked;
	private int inCol, outCol;
	private float xa, ya, xb, yb, xc, yc, sw;

	public JTriangleShape(float x1, float y1, float x2, float y2, float x3, float y3) {
		this(x1, y1, x2, y2, x3, y3, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight,
				JPlotShape.useFill, JPlotShape.useStroke);
	}

	public JTriangleShape(float x1, float y1, float x2, float y2, float x3, float y3, boolean filled,
			boolean withOutline) {
		this(x1, y1, x2, y2, x3, y3, JPlotShape.fillColour, JPlotShape.strokeColour, JPlotShape.strokeWeight, filled,
				withOutline);
	}

	public JTriangleShape(float x1, float y1, float x2, float y2, float x3, float y3, int colour) {
		this(x1, y1, x2, y2, x3, y3, colour, colour, JPlotShape.strokeWeight, JPlotShape.useFill, JPlotShape.useStroke);
	}

	public JTriangleShape(float x1, float y1, float x2, float y2, float x3, float y3, int inner_colour,
			int outer_colour) {
		this(x1, y1, x2, y2, x3, y3, inner_colour, outer_colour, JPlotShape.strokeWeight, true, true);
	}

	public JTriangleShape(float x1, float y1, float x2, float y2, float x3, float y3, int inner_colour,
			int outer_colour, float stroke_weight) {
		this(x1, y1, x2, y2, x3, y3, inner_colour, outer_colour, stroke_weight, true, true);
	}

	public JTriangleShape(float x1, float y1, float x2, float y2, float x3, float y3, int colour, boolean filled,
			boolean withOutline) {
		this(x1, y1, x2, y2, x3, y3, colour, colour, JPlotShape.strokeWeight, filled, withOutline);
	}

	public JTriangleShape(float x1, float y1, float x2, float y2, float x3, float y3, int inner_colour,
			int outer_colour, float stroke_weight, boolean filled, boolean withOutline) {
		xa = x1;
		ya = y1;
		xb = x2;
		yb = y2;
		xc = x3;
		yc = y3;
		inCol = inner_colour;
		outCol = outer_colour;
		sw = stroke_weight;
		isFilled = filled;
		isStroked = withOutline;
	}

	@Override
	public void draw(JPlot plot, PGraphics g) {
		if (isFilled) {
			g.fill(inCol);
		} else {
			g.noFill();
		}
		if (isStroked) {
			g.stroke(outCol);
			g.strokeWeight(sw);
		} else {
			g.noStroke();
		}
		g.triangle(xa, ya, xb, yb, xc, yc);
	}
}
