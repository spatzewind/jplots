package pplots.layer;

import pplots.PPlot;
import processing.core.PGraphics;
import processing.core.PShape;

public abstract class PLayer {
	public abstract void createRasterImg(PPlot plot, PGraphics g);
	public abstract void createVectorImg(PPlot plot, PShape s, int x, int y, int w, int h);
	public abstract void setRange(double xmin,double xmax,double ymin,double ymax);
	public abstract double[] getRange();
}
