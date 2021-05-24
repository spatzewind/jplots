package pplots.layer;

import pplots.PAxis;
import pplots.PPlot;
import pplots.shapes.PGroupShape;
import pplots.shapes.PTextured2DShape;
import pplots.transform.PPlateCarreeProjection;
import pplots.transform.PProjection;
import processing.core.PApplet;
import processing.core.PGraphics;
import processing.core.PImage;

public class PImageLayer extends PLayer {
	
	private double[] srcExt;
	private PImage srcImg;
	private PProjection srcProj;

	public PImageLayer(PImage image) {
		this(image, new PPlateCarreeProjection(true), -180, 90, 180, -90); }
	public PImageLayer(PImage image, PProjection proj, double[] extent) {
		this(image, proj, extent[0], extent[1], extent[2], extent[3]); }
	public PImageLayer(PImage image, PProjection proj, double x_lefttop, double y_lefttop, double x_rightbottom, double y_rightbottom) {
		srcImg = image;
		srcExt = new double[] {x_lefttop, y_lefttop, x_rightbottom, y_rightbottom};
		srcProj = proj;
	}

	@Override
	public void createRasterImg(PPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(PAxis ax, int layernum, PGroupShape s) {
		int[] p = ax.getSize();
		double xs = p[2]/(maxX-minX), ys = p[3]/(maxY-minY);
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
				double[] xy = ax.getGeoProjection().fromPROJtoLATLON(minX+i/xs, maxY-j/ys, false);
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
		PTextured2DShape texs = new PTextured2DShape(img);
		texs.addVector(p[0],      p[1],      0,    0   );
		texs.addVector(p[0]+p[2], p[1],      p[2], 0   );
		texs.addVector(p[0]+p[2], p[1]+p[3], p[2], p[3]);
		texs.addVector(p[0],      p[1]+p[3], 0,    p[3]);
		s.addChild(texs);
	}

}
