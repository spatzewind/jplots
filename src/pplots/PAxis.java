package pplots;

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import pplots.layer.PImageLayer;
import pplots.layer.PLayer;
import pplots.layer.PScatterLayer;
import pplots.layer.PXYLayer;
import pplots.shapes.PGroupShape;
import pplots.shapes.PLineShape;
import pplots.shapes.PPlotShape;
import pplots.shapes.PTextShape;
import pplots.transform.PIdentityProjection;
import pplots.transform.PProjection;
import pplots.transform.PRectangleProjection;
import processing.core.PApplet;
import processing.core.PFont;
import processing.core.PImage;

public class PAxis {
	
	private static Map<String, PImage> loadedPreDefImgs = new HashMap<>();

	private PPlot pplot;
	private boolean xRangeFix,yRangeFix, isGeoAxis;
	private boolean xAxOn, yAxOn, xGrdOn, yGrdOn;
	private int px, py, pw, ph;
	private double minX,maxX,minY,maxY;
	private double txtsize;
	private List<PLayer> layers;
	private PFont pfont;
	private PProjection projection;
	
	PAxis(PPlot plot, int pos_x, int pos_y, int width, int height) {
		pplot = plot;
		px = pos_x;
		py = pos_y;
		pw = width;
		ph = height;
		layers = new ArrayList<PLayer>();
		if(plot.isDebug())
			System.out.println("[DEBUG] created PAxis-object: x/y="+pos_x+"/"+pos_y+" w/h="+width+"/"+height);
		defaults();
	}
	
	private void defaults() {
		xRangeFix = false;
		xAxOn = true;
		yAxOn = true;
		xGrdOn = false;
		yGrdOn = false;
		minX = -1d;
		maxX =  1d;
		minY = -1d;
		maxY =  1d;
		txtsize = 300d*10d/72d;
		projection = new PIdentityProjection();
	}


	//************************************
	//**** PUBLIC ************************
	//************************************

	//....
	public void contour(double[] x, double[] y, double[][] z) {
		this.contour(x, y, z, 10, (Object[])null); }
	public void contour(double[] x, double[] y, double[][] z, int levels, Object... params) {
		
	}

	public void plot(float[] x, float[] y) {
		this.plot(x, y, 0xff000000, 3f, "-", (Object)null); }
	public void plot(float[] x, float[] y, int colour, float linewidth, String linestyle, Object... params) {
		PLayer xyl = new PXYLayer(x, y, colour, linewidth, linestyle); layers.add(xyl);
		readLineParams(xyl, params); updateRange(xyl); }
	public void plot(double[] x, double[] y) {
		this.plot(x, y, 0xff000000, 3d, "-", (Object)null); }
	public void plot(double[] x, double[] y, int colour, double linewidth, String linestyle, Object... params) {
		PLayer xyl = new PXYLayer(x, y, colour, linewidth, linestyle); layers.add(xyl);
		readLineParams(xyl, params); updateRange(xyl); }

	public void scatter(float[] x, float[] y) {
		this.scatter(x, y, 0xff000000, 1f, "c", (Object)null); }
	public void scatter(float[] x, float[] y, int colour, float iconsize, String symbol, Object... params) {
		PLayer scl = new PScatterLayer(x, y, colour, iconsize, symbol); layers.add(scl);
		readLineParams(scl, params); updateRange(scl); }
	public void scatter(double[] x, double[] y) {
		this.scatter(x, y, 0xff000000, 1d, "c", (Object)null); }
	public void scatter(double[] x, double[] y, int colour, double iconsize, String symbol, Object... params) {
		PLayer scl = new PScatterLayer(x, y, colour, iconsize, symbol); layers.add(scl);
		readLineParams(scl, params); updateRange(scl); }

	/**
	 * predefined images are used as background images in plot, especially with geographical projections
	 * <p>
	 *     earth&nbsp; -- NASA image of earths surface<br>
	 *     earth2 -- ETOPO-like image of earths surface<br>
	 *     mercury,venus,mars,jupiter,saturn,uranus,neptune -- NASA images of surfaces of planets in our solar system<br>
	 *     hsbpalette -- rainbow palette of all possible hsb-colors
	 * </p>
	 * 
	 * @param predefined_image one name of above list of predefined images
	 */
	public void predefImgShow(String predefined_image) {
		PLayer iml;
		if(loadedPreDefImgs.containsKey(predefined_image)) {
			iml = new PImageLayer(loadedPreDefImgs.get(predefined_image));
		} else {
			iml = new PImageLayer(loadPreDefImg(pplot.getApplet(), predefined_image));
		}
		layers.add(0, iml);
	}
	public void imgShow(PImage img) {
		PLayer iml = new PImageLayer(img); layers.add(iml);
		//updateRange(iml);
	}

	//....
	public PAxis setPositionAndSize(int pos_x, int pos_y, int width, int height) {
		px = pos_x;
		py = pos_y;
		pw = width;
		ph = height;
		if(pplot.isDebug())
			System.out.println("[DEBUG] resize PAxis-object: x/y="+pos_x+"/"+pos_y+" w/h="+width+"/"+height);
		return this;
	}
	public PAxis setXRange(float xmin, float xmax) { minX = xmin; maxX = xmax; xRangeFix = true; return this; }
	public PAxis setXRange(double xmin, double xmax) { minX = xmin; maxX = xmax; xRangeFix = true; return this; }
	public PAxis setYRange(float ymin, float ymax) { minY = ymin; maxY = ymax; yRangeFix = true; return this; }
	public PAxis setYRange(double ymin, double ymax) { minY = ymin; maxY = ymax; yRangeFix = true; return this; }
	public PAxis setRange(float xmin, float xmax, float ymin, float ymax) {
		minX = xmin; maxX = xmax; minY = ymin; maxY = ymax; xRangeFix = true; yRangeFix = true; return this; }
	public PAxis setRange(double xmin, double xmax, double ymin, double ymax) {
		minX = xmin; maxX = xmax; minY = ymin; maxY = ymax; xRangeFix = true; yRangeFix = true; return this; }
	public void setAxis(String axis, String which, boolean onoff) {
		boolean setX=false, setY=false;
		boolean setAx=false, setGrd=false;
		if("both".equals(axis.toLowerCase())) { setX=true; setY=true; }
		if("x".equals(axis.toLowerCase())) setX=true;
		if("y".equals(axis.toLowerCase())) setY=true;
		if("both".equals(which.toLowerCase())) { setAx=true; setGrd=true; }
		if("a".equals(which.toLowerCase()) || "axis".equals(which.toLowerCase())) setAx=true;
		if("g".equals(which.toLowerCase()) || "grid".equals(which.toLowerCase())) setGrd=true;
		if(setAx) { if(setX) xAxOn=onoff; if(setY) yAxOn=onoff; }
		if(setGrd) { if(setX) xGrdOn=onoff; if(setY) yGrdOn=onoff; }
	}
	public void setGrid() {
		setGrid("both", true); }
	public void setGrid(String axis, boolean onoff) {
		if(isGeoAxis) {
			
		} else {
			boolean setX=false, setY=false;
			if("both".equals(axis.toLowerCase())) { setX=true; setY=true; }
			if("x".equals(axis.toLowerCase())) setX=true;
			if("y".equals(axis.toLowerCase())) setY=true;
			if(setX) xGrdOn=onoff;
			if(setY) yGrdOn=onoff;
		}
	}
	public PAxis setGeoProjection(PProjection proj) {
		projection = proj;
		isGeoAxis = true;
		double[] rr = projection.defaultMapExtend();
		setRange(rr[0], rr[1], rr[2], rr[3]);
		xAxOn = false;  yAxOn = false;
		xGrdOn = false; yGrdOn = false;
		pplot.redraw(true);
		return this;
	}

	//....
	public void setFont(PFont font) {
		pfont = font; }
	public void setTextSize(double ts) {
		txtsize = PPlot.dpi*ts/72d; }


	//************************************
	//**** GETTER ************************
	//************************************
	
	public PPlot getPlot() { return pplot; }
	public int[] getSize() { return new int[] {px,py,pw,ph}; }
	public double[] getRange() { return new double[] {minX,maxX,minY,maxY}; }
	public boolean isGeoAxis() { return isGeoAxis; }
	public PProjection getGeoProjection() { return projection; }
	public PLayer getLayer(int layernum) { return layers.get(layernum); }


	//************************************
	//**** PACKAGE PRIVATE ***************
	//************************************
	
	PGroupShape createPlot(PApplet applet, int w, int h) {
		if(isGeoAxis) {
			double r = 0.5d * Math.max((maxX-minX)/pw, (maxY-minY)/ph);
			double xm = 0.5d*(minX+maxX);
			double ym = 0.5d*(minY+maxY);
			minX = xm - r*pw; maxX = xm + r*pw;
			minY = ym - r*pw; maxY = ym + r*pw;
		}
		if(pplot.isDebug())
			System.out.println("[DEBUG] PAxis-object: min/max={x:"+minX+"/"+maxX+", y:"+minY+"/"+maxY+
				"} with "+layers.size()+" layers");
		PGroupShape graph = new PGroupShape();
		if(isGeoAxis) {
			projection.drawBorder(this, graph);
		} else {
			if(xAxOn || xGrdOn) graph.addChild(createXAxis());
			if(yAxOn || yGrdOn) graph.addChild(createYAxis());
			if(xAxOn || yAxOn) {
				if(xAxOn) {
					graph.addChild(new PLineShape(px, py,    px+pw, py   ));
					graph.addChild(new PLineShape(px, py+ph, px+pw, py+ph));
				}
				if(yAxOn) {
					graph.addChild(new PLineShape(px,    py, px,    py+ph));
					graph.addChild(new PLineShape(px+pw, py, px+pw, py+ph));
				}
			}
		}
		for(int l=0; l<layers.size(); l++) {
			PLayer layer = layers.get(l);
			layer.setRange(minX,maxX,minY,maxY);
			layer.createVectorImg(this, l, graph);
		}
		return graph;
	}
	PProjection getProjection() {
		return projection; }

	//************************************
	//**** PRIVATE ***********************
	//************************************
	
	private void readLineParams(PLayer layer, Object... params) {
		if(params==null)
			return;
		int o=0;
		while(o<params.length) {
			if(params[o] instanceof String) {
				String p = ((String) params[o]).toLowerCase();
				boolean isunread = true;
				if(isunread && "transform".equals(p) && o+1<params.length) {
					layer.setSourceProjection((PProjection) params[o+1]); o++; isunread=false; }
				if(isunread && "anglemode".equals(p) && o+1<params.length) {
					layer.angleMode((String) params[o+1]); o++; isunread=false; }
			} else {
				System.err.println("[ERROR] Cannot interprete param "+o+": "+params[o]);
			}
			o++;
		}
	}
	
	private void updateRange(PLayer layer) {
		double[] r = layer.getRange();
		if(layers.size()==1) {
			if(!xRangeFix) {
				minX = r[0];
				maxX = r[1];
			}
			if(!yRangeFix) {
				minY = r[2];
				maxY = r[3];
			}
		} else {
			if(!xRangeFix) {
				if(r[0]<minX) minX = r[0];
				if(r[1]>maxX) maxX = r[1];
			}
			if(!yRangeFix) {
				if(r[2]<minY) minY = r[2];
				if(r[3]>maxY) maxY = r[3];
			}
		}
	}
	private PGroupShape createXAxis() {
		PGroupShape axisgrid = new PGroupShape();
		double[] ticks = PPlotMath.optimalLinearTicks(minX, maxX);
		double[] tcpos = PPlotMath.dlerp(ticks,minX,maxX,px,px+pw);
		if(pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for(int t=2; t<ticks.length; t++) {
				tickStr += ", "+PApplet.nf((float)ticks[t],0,2);
				posStr  += ", "+PApplet.nf((float)tcpos[t],0,2);
			}
			System.out.println("[DEBUG] PAxis-object: Xtickfactors={p10: "+ticks[0]+", f: "+ticks[1]+"}");
			System.out.println("[DEBUG] PAxis-object: Xtickval={"+tickStr.substring(2)+"}");
			System.out.println("[DEBUG] PAxis-object: Xtickpos={"+posStr.substring(2)+"}");
		}
		if(xGrdOn) {
			PPlotShape.stroke(0xff999999); PPlotShape.strokeWeight(2f);
			for(int t=2; t<ticks.length; t++)
				if(ticks[t]>=Math.min(minX, maxX) && ticks[t]<=Math.max(minX, maxX))
					axisgrid.addChild(new PLineShape((float)tcpos[t],py,(float)tcpos[t],py+ph));
		}
		if(xAxOn) {
			PPlotShape.stroke(0xff000000); PPlotShape.strokeWeight(2f);
			double vf = 1d/(ticks[0]);
			int decimal = (int) (1000d*ticks[1]+0.5d);
			decimal = decimal%100==0 ? 1 : decimal%10==0 ? 2 : 3;
			for(int t=2; t<ticks.length; t++)
				if(ticks[t]>=Math.min(minX, maxX) && ticks[t]<=Math.max(minX, maxX)) {
					axisgrid.addChild(new PLineShape((float)tcpos[t],py+ph,(float)tcpos[t],py+1.1f*ph));
					//axisgrid.addChild(ap.createShape(PShape.TEXT, "H", (float)tcpos[t],py-0.1f*ph,(float)tcpos[t],py));
					axisgrid.addChild(new PTextShape(PApplet.nf((float)(ticks[t]*vf),0,decimal), (float)tcpos[t], py+1.11f*ph, (float)txtsize, PApplet.CENTER, PApplet.TOP, 0xff000000));
				}
		}
		return axisgrid;
	}
	private PGroupShape createYAxis() {
		PGroupShape axisgrid = new PGroupShape();
		double[] ticks = PPlotMath.optimalLinearTicks(minY, maxY);
		double[] tcpos = PPlotMath.dlerp(ticks,minY,maxY,py+ph,py);
		if(pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for(int t=2; t<ticks.length; t++) {
				tickStr += ", "+PApplet.nf((float)ticks[t],0,2);
				posStr  += ", "+PApplet.nf((float)tcpos[t],0,2);
			}
			System.out.println("[DEBUG] PAxis-object: Ytickfactors={p10: "+ticks[0]+", f: "+ticks[1]+"}");
			System.out.println("[DEBUG] PAxis-object: Ytickval={"+tickStr.substring(2)+"}");
			System.out.println("[DEBUG] PAxis-object: Ytickpos={"+posStr.substring(2)+"}");
		}
		if(yGrdOn) {
			PPlotShape.stroke(0xff999999); PPlotShape.strokeWeight(2f);
			for(int t=2; t<ticks.length; t++)
				if(ticks[t]>=Math.min(minY, maxY) && ticks[t]<=Math.max(minY, maxY))
					axisgrid.addChild(new PLineShape(px,(float)tcpos[t],px+pw,(float)tcpos[t]));
		}
		if(yAxOn) {
			PPlotShape.stroke(0xff000000); PPlotShape.strokeWeight(2f);
			double vf = 1d/(ticks[0]);
			int decimal = (int) (1000d*ticks[1]+0.5d);
			decimal = decimal%100==0 ? 1 : decimal%10==0 ? 2 : 3;
			for(int t=2; t<ticks.length; t++)
				if(ticks[t]>=Math.min(minY, maxY) && ticks[t]<=Math.max(minY, maxY)) {
					axisgrid.addChild(new PLineShape(px-0.1f*pw,(float)tcpos[t],px,(float)tcpos[t]));
					//axisgrid.addChild(ap.createShape(PShape.TEXT, "H", (float)tcpos[t],py-0.1f*ph,(float)tcpos[t],py));
					axisgrid.addChild(new PTextShape(PApplet.nf((float)(ticks[t]*vf),0,decimal), px-0.11f*pw, (float)tcpos[t], (float)txtsize, PApplet.RIGHT, PApplet.CENTER, 0xff000000));
				}
		}
		return axisgrid;
	}

	private static PImage loadPreDefImg(PApplet applet, String name) {
		BufferedImage bimg;
		try {
			bimg = ImageIO.read(PPlot.class.getResourceAsStream("/data/"+name+".png"));
			PImage limg = applet.createImage(bimg.getWidth(), bimg.getHeight(), PApplet.ARGB);
			bimg.getRGB(0, 0, bimg.getWidth(), bimg.getHeight(), limg.pixels, 0, bimg.getWidth());
			loadedPreDefImgs.put(name, limg);
			return limg;
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.err.println("The image \""+name+".png\" is missing in the jar or inaccessible. Contact the author of ##library.name## ##library.prettyVersion## if the image should be there.");
		return null;
	}
}
