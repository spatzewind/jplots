package pplots.layer;

import pplots.PMath;
import pplots.PPlot;
import processing.core.PApplet;
import processing.core.PGraphics;
import processing.core.PShape;

public class PXYLayer extends PLayer {
	
	private double[] xarrayx, yarrayy;
	private int col;
	private double lw;
	private String ls;
	private double minX,maxX,minY,maxY;
	
	public PXYLayer(float[] x, float[] y, int colour, float linewidth, String linestyle) {
		xarrayx = new double[x.length];
		for(int i=0; i<x.length; i++)
			xarrayx[i] = x[i];
		yarrayy = new double[y.length];
		for(int i=0; i<y.length; i++)
			yarrayy[i] = y[i];
		minX = PMath.dmin(xarrayx);
		maxX = PMath.dmax(xarrayx);
		minY = PMath.dmin(yarrayy);
		maxY = PMath.dmax(yarrayy);
		col = colour;
		lw = linewidth;
		ls = linestyle;
	}
	public PXYLayer(double[] x, double[] y, int colour, double linewidth, String linestyle) {
		xarrayx = x;
		yarrayy = y;
		minX = PMath.dmin(xarrayx);
		maxX = PMath.dmax(xarrayx);
		minY = PMath.dmin(yarrayy);
		maxY = PMath.dmax(yarrayy);
		col = colour;
		lw = linewidth;
		ls = linestyle;
	}
	
	
	@Override
	public void setRange(double xmin, double xmax, double ymin, double ymax) {
		minX = xmin; maxX = xmax;
		minY = ymin; maxY = ymax;
	}
	
	@Override
	public double[] getRange() {
		return new double[] {minX,maxX,minY,maxY};
	}

	@Override
	public void createRasterImg(PPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(PPlot plot, PShape s, int x, int y, int w, int h) {
		int shapeType = PApplet.LINES;
		PShape xyShape = plot.getApplet().createShape(PShape.GEOMETRY);
		xyShape.beginShape(shapeType);
		xyShape.noFill(); xyShape.stroke(col); xyShape.strokeWeight((float)lw);
		double lln=1d, llf=0d, lpn=0d, lpf=0d, loff = 0d;
		if("-".equals(ls)) { lln=1000*lw; llf=0; lpn=0; lpf=0; }
		if(".".equals(ls)) { lln=0; llf=0; lpn=1*lw; lpf=3*lw; }
		if(",".equals(ls)) { lln=12*lw; llf=3*lw; lpn=0; lpf=0; }
		if(";".equals(ls)) { lln=8*lw; llf=3*lw; lpn=1*lw; lpf=3*lw; }
		double xs = w/(maxX-minX), ys = h/(maxY-minY);
		int li = 0;
		for(int i=0; i+1<xarrayx.length; i++)
			if(Double.isFinite(xarrayx[i]) && Double.isFinite(xarrayx[i+1]) &&
					Double.isFinite(yarrayy[i]) && Double.isFinite(yarrayy[i+1])) {
				double x1 = x+xs*(xarrayx[i  ]-minX);
				double x2 = x+xs*(xarrayx[i+1]-minX);
				double y1 = y+ys*(maxY-yarrayy[i  ]);
				double y2 = y+ys*(maxY-yarrayy[i+1]);
				double dx = x2-x1, dy = y2-y1;
				double l = Math.sqrt(dx*dx+dy*dy);
				dx /= l; dy /= l;
				double lpos = 0d, ldif = 0d;
				//TODO cutoffs, if line outside of PAxis range
				while(lpos<l) {
					switch(li) {
						case 0: if(lln==0d) break; ldif = Math.min(l-lpos, lln-loff); break;
						case 1: if(llf==0d) break; ldif = Math.min(l-lpos, llf-loff); break;
						case 2: if(lpn==0d) break; ldif = Math.min(l-lpos, lpn-loff); break;
						case 3: if(lpf==0d) break; ldif = Math.min(l-lpos, lpf-loff); break;
					}
					if(li%2==0) {
						xyShape.vertex((float)(x1+lpos*dx), (float)(y1+lpos*dy));
						lpos += ldif;
						xyShape.vertex((float)(x1+lpos*dx), (float)(y1+lpos*dy));
					}
					loff += ldif;
					switch(li) {
						case 0: if(loff>=lln) loff -= lln; break;
						case 1: if(loff>=llf) loff -= llf; break;
						case 2: if(loff>=lpn) loff -= lpn; break;
						case 3: if(loff>=lpf) loff -= lpf; break;
					}
					li = (li+1)&3;
				}
			} else {
				loff = 0d;
				li = 0;
			}
		xyShape.endShape();
		s.addChild(xyShape);
	}

}
