package jplots.shapes;

import jplots.JPlot;
import processing.core.PApplet;
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
		g.noFill();
		g.stroke(col);
		g.strokeWeight(sw);
//		if(coords.length==4) {
//			g.line(coords[0], coords[1], coords[2], coords[3]);
//		} else {
			g.beginShape(PApplet.LINE_STRIP);
			for(int i=0; i+1<coords.length; i+=2) {
				g.vertex(coords[i],coords[i+1]);
			}
			g.endShape();
//		}
	}

}
