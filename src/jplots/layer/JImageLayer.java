package jplots.layer;

import jplots.JAxis;
import jplots.JPlot;
import jplots.shapes.JGroupShape;
import jplots.shapes.JImageShape;
import jplots.transform.JProjection;
import jplots.transform.PlateCarreeJProjection;
import processing.core.PApplet;
import processing.core.PGraphics;
import processing.core.PImage;

public class JImageLayer extends JPlotsLayer {
	
	private double[] srcExt;
	private PImage srcImg;
	private JProjection srcProj;

	public JImageLayer(PImage image) {
		this(image, new PlateCarreeJProjection(true), -180, 90, 180, -90); }
	public JImageLayer(PImage image, JProjection proj, double[] extent) {
		this(image, proj, extent[0], extent[1], extent[2], extent[3]); }
	public JImageLayer(PImage image, JProjection proj, double x_lefttop, double y_lefttop, double x_rightbottom, double y_rightbottom) {
		srcImg = image;
		srcExt = new double[] {x_lefttop, y_lefttop, x_rightbottom, y_rightbottom};
		srcProj = proj;
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		double Xin = ax.isXlogAxis()?Math.log10(minX):minX, Xax = ax.isXlogAxis()?Math.log10(maxX):maxX;
		double Yin = ax.isYlogAxis()?Math.log10(minY):minY, Yax = ax.isYlogAxis()?Math.log10(maxY):maxY;
		double xs = p[2]/(Xax-Xin), ys = p[3]/(Yax-Yin);
		double us = srcImg.width/(srcExt[2]-srcExt[0]), vs = srcImg.height/(srcExt[3]-srcExt[1]);
		if(img==null) {
			img = ax.getPlot().getApplet().createImage(p[2], p[3], PApplet.ARGB);
		} else if(img.width!=p[2] || img.height!=p[3]) {
			img = ax.getPlot().getApplet().createImage(p[2], p[3], PApplet.ARGB);
		}
		img.loadPixels();
		srcImg.loadPixels();
		for(int j=0; j<p[3]; j++) {
			for(int i=0; i<p[2]; i++) {
				int idx = j*p[2]+i;
				double[] xy = { invertAxisX ? Xax-i/xs : Xin+i/xs,
								invertAxisY ? Yin+j/ys : Yax-j/ys };
				if(ax.isXlogAxis()) xy[0] = Math.pow(10d, xy[0]);
				if(ax.isYlogAxis()) xy[1] = Math.pow(10d, xy[1]);
				xy = ax.getGeoProjection().fromPROJtoLATLON(xy[0], xy[1], false);
				xy = srcProj.fromLATLONtoPROJ(xy[0], xy[1], false);
				if(Double.isNaN(xy[0]) || Double.isNaN(xy[1])) {
					img.pixels[idx] = 0x00ffffff;
					continue;
				}
				double u = us * (xy[0]-srcExt[0]), v = vs * (xy[1]-srcExt[1]);
				if(u<0d || v<0d) {
					img.pixels[idx] = 0x00ffffff;
					continue;
				}
				int ui = (int) u, vi = (int) v;
				if(ui>=srcImg.width || vi>=srcImg.height) {
					img.pixels[idx] = 0x00ffffff;
					continue;
				}
				img.pixels[idx] = srcImg.pixels[vi*srcImg.width+ui];
			}
		}
		img.updatePixels();
		s.addChild(new JImageShape(img, p[0], p[1], p[2], p[3]));
	}

}
