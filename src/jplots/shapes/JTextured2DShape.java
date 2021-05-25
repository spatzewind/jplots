package jplots.shapes;

import java.util.ArrayList;
import java.util.List;

import processing.core.PApplet;
import processing.core.PGraphics;
import processing.core.PImage;
import processing.core.PShape;

public class JTextured2DShape extends JPlotShape {

	private List<double[]> vectors;
	private PImage texture;

	public JTextured2DShape(PImage image) {
		texture = image;
		vectors = new ArrayList<double[]>();
	}
	
	public void addVector(double x, double y, double u, double v) {
		vectors.add(new double[] {x,y,u,v});
	}
	
	@Override
	public void draw(PGraphics g) {
		g.noFill(); g.noStroke();
		g.beginShape();
		g.texture(texture);
		for(int v=0; v<vectors.size(); v++) {
			double[] vv = vectors.get(v);
			g.vertex((float)vv[0], (float)vv[1], (float)vv[2], (float)vv[3]);
		}
		g.endShape(PShape.CLOSE);
	}
}
