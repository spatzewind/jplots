package jplots.shapes;

import jplots.JPlot;
import processing.core.PGraphics;

public class JLineShape extends JPlotShape {

	private int col;
	private float sw;
	private float[] coords;
	
//	public JLineShape(float... xy) {
//		this(JPlotShape.strokeWeight, JPlotShape.strokeColour, xy);
//	}
	public JLineShape(float stroke_weight, int colour, float... xy) {
		col = colour;
		sw = stroke_weight;
		coords = xy;
	}

	@Override
	public void draw(JPlot plot, PGraphics g) {
		g.stroke(col);
		g.strokeWeight(sw);
		for(int i=0; i+3<coords.length; i+=2) {
			if(Float.isNaN(coords[i]) || Float.isNaN(coords[i+1]) || Float.isNaN(coords[i+2]) || Float.isNaN(coords[i+3]))
				continue;
			g.line(coords[i], coords[i+1], coords[i+2], coords[i+3]);
		}
	}
}
