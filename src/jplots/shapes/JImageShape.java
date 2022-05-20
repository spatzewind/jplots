package jplots.shapes;

import jplots.JPlot;
import processing.core.PGraphics;
import processing.core.PImage;

public class JImageShape extends JPlotShape {

	private double tx, ty, tw, th;
	private PImage texture;

	public JImageShape(PImage image, double x, double y, double w, double h) {
		texture = image;
		tx = x;
		ty = y;
		tw = w;
		th = h;
	}

	@Override
	public void draw(JPlot plot, PGraphics g) {
		g.image(texture, (float) tx, (float) ty, (float) tw, (float) th);
	}
}
