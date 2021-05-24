package pplots;

import pplots.shapes.PGroupShape;
import pplots.transform.PProjection;
import processing.core.*;

/**
 * This is a template class and can be used to start a new processing Library.
 * Make sure you rename this class as well as the name of the example package 'template' 
 * to your own Library naming convention.
 * 
 * (the tag example followed by the name of an example included in folder 'examples' will
 * automatically include the example in the javadoc.)
 *
 * @example SimpleXY
 */

public class PPlot {
	
	public static double dpi = 300d;

	private PApplet myParent;
	private PGraphics plotImg;
	private PGroupShape plotShp;

	private boolean img_is_created, useDebug;
	private int width, height, lastAxisNum;
	private PAxis[] axes;


	//************************************
	//**** CONSTRUCTOR *******************
	//************************************

	/**
	 * a Constructor, usually called in the setup() method in your sketch to
	 * initialize and start the Library.
	 * 
	 * @example Hello
	 * @param applet the parent PApplet
	 */
	public PPlot(PApplet applet) {
		this.myParent = applet;
		this.width = 0;
		this.height = 0;
		this.img_is_created = false;
		this.useDebug = false;
		this.welcome();
	}


	//************************************
	//**** STATIC ************************
	//************************************

	/**
	 * return the version of the Library.
	 * 
	 * @return String
	 */
	public static String version() {
		return PConstants.VERSION;
	}


	//************************************
	//**** PURE PPLOT ********************
	//************************************
	
	/**
	 * starts the new figure/plot with specified width and height in inch
	 * 
	 * @return the PPlot object
	 */
	public PPlot figure() { return subplots(1.7066666667d, 1.7066666667d, 1, 1); }
	/**
	 * starts the new figure/plot with specified width and height in inch
	 * 
	 * @param width  width of figure in inch
	 * @param height height of figure in inch
	 * @return the PPlot object
	 */
	public PPlot figure(double width, double height) { return subplots(width, height, 1, 1); }

	/**
	 * starts the new figure/plot
	 * 
	 * @param nrows  number of rows for array of diagrams
	 * @param ncols  number of columns for array of diagrams
	 * @return the PPlot object
	 */
	public PPlot subplots(int nrows, int ncols) { return subplots(1.7066666667d, 1.7066666667d, nrows, ncols); }
	/**
	 * starts the new figure/plot with specified width and height in inch
	 * 
	 * @param width  width of figure in inch
	 * @param height height of figure in inch
	 * @param nrows  number of rows for array of diagrams
	 * @param ncols  number of columns for array of diagrams
	 * @return the PPlot object
	 */
	public PPlot subplots(double width, double height, int nrows, int ncols) {
		this.width  = (int) (dpi * width);
		this.height = (int) (dpi * height);
		this.axes   = new PAxis[nrows*ncols];
		if(useDebug)
			System.out.println(
				"[DEBUG] PPlot: width/height: "+this.width+"/"+this.height+"px ("+width+"/"+height+"inch)\n"+
				"               axes:         "+ncols+"x"+nrows);
		double ws = 1d / (0.25d+(ncols-1)*1.5d+1.25d);
		double hs = 1d / (0.25d+(nrows-1)*1.5d+1.25d);
		int pawid = (int) (dpi * width * ws);
		int pahei = (int) (dpi * height * hs);
		for(int r=0; r<nrows; r++) {
			int py = (int) (dpi * height * (0.25d+1.5d*r) * hs);
			for(int c=0; c<ncols; c++) {
				int px = (int) (dpi * width * (0.25d+1.5d*c) * ws);
				this.axes[r*ncols+c] = new PAxis(this, px, py, pawid, pahei);
			}
		}
		this.img_is_created = false;
		return this;
	}
	/**
	 * add new PAxis object as a subplot in the PPlot object
	 * 
	 * @param pos_x  left horizontal position relative to figure width
	 * @param pos_y  top vertical position relative to figure height
	 * @param width  width relative to the figure width
	 * @param height height relative to the figure height
	 * @return the PPlot object
	 */
	public PPlot addSubplot(double pos_x, double pos_y, double width, double height) {
		PAxis[] tempa = new PAxis[axes.length+1];
		for(int a=0; a<axes.length; a++)
			tempa[a] = axes[a];
		tempa[axes.length] = new PAxis(this,
									   (int)(pos_x*this.width),(int)(pos_y*this.height),
									   (int)(width*this.width+0.5d),(int)(height*this.height+0.5d));
		axes = new PAxis[tempa.length];
		for(int a=0; a<axes.length; a++)
			axes[a] = tempa[a];
		lastAxisNum = axes.length-1;
		if(useDebug)
			System.out.println("[DEBUG] addSubplot: "+axes.length+" PAxis-objects");
		return this;
	}

	public PPlot setGeoProjection(PProjection proj) {
		gca().setGeoProjection(proj);
		return this;
	}
	
	/**
	 * creates image of the PPlot -- the figure
	 */
	public PPlot createImage() {
		if(useDebug)
			System.out.println("[DEBUG] Start image-creation ...");
		if(plotImg==null) {
			plotImg = myParent.createGraphics(width, height, PApplet.P2D);
		} else if(plotImg.width!=width || plotImg.height!=height) {
			plotImg = myParent.createGraphics(width, height, PApplet.P2D);
		}
		if(useDebug)
			System.out.println("[DEBUG]   - create PGraphics-object "+plotImg.width+"x"+plotImg.height+"px");
		plotShp = new PGroupShape();
		if(useDebug)
			System.out.println("[DEBUG]   - add "+axes.length+" subplots to plotting queue");
		for(PAxis ax: axes)
			plotShp.addChild(ax.createPlot(myParent,width,height));
		plotImg.beginDraw();
		plotImg.clear();
		if(useDebug)
			System.out.println("[DEBUG]   - begin drawing ...");
		plotShp.draw(plotImg);
		plotImg.endDraw();
		if(useDebug)
			System.out.println("[DEBUG]   - end drawing");
		img_is_created = true;
		return this;
	}

	/**
	 * gives the figure<br>
	 * the figure will be created automatically, if it isn't allready
	 * 
	 * @return the figure as PImage
	 */
	public PImage show() {
		if(!img_is_created) {
			createImage();
			if(useDebug)
				System.out.println("[DEBUG] start image creation, because it ws not done before or changes happend.");
		}
		if(useDebug)
			System.out.println("[DEBUG] check: \"plotImg\"="+plotImg);
		return plotImg;
	}
	public PGraphics getGraphic() {
		return plotImg;
	}
	public PPlot redraw(boolean redraw) {
		img_is_created = redraw?false:img_is_created;
		return this;
	}
	
	//....
	/**
	 * set debugmode on/off
	 * 
	 * @param bdebug use true for debugmode on
	 */
	public void debug(boolean bdebug) {
		useDebug = bdebug; }
	public boolean isDebug() {
		return useDebug; }


	//************************************
	//**** GETTER ************************
	//************************************
	
	public PApplet getApplet() { return myParent; }
	
	public PAxis[] ga() { return axes; }
	public PAxis gca() { return axes[lastAxisNum]; }
	public PAxis ga(int a) {
		if(a<-axes.length || a>axes.length-1)
			throw new IndexOutOfBoundsException("Number out of range of axis count: "+a+" <> {0..."+axes.length+"}");
		lastAxisNum = a + (a<0 ? axes.length : 0);
		return axes[lastAxisNum];
	}


	//************************************
	//**** PRIVATE ***********************
	//************************************

	private void welcome() {
		System.out.println("##library.name## ##library.prettyVersion## by ##author##");
	}


	//************************************
	//**** PASS ON ***********************
	//************************************
	
	public void plot(float[] x, float[] y) {
		gca().plot(x, y); }
	public void plot(float[] x, float[] y, int colour, float linewidth, String linestyle, Object... params) {
		gca().plot(x,y,colour,linewidth,linestyle,params); }
	public void plot(double[] x, double[] y) {
		gca().plot(x,y); }
	public void plot(double[] x, double[] y, int colour, double linewidth, String linestyle, Object... params) {
		gca().plot(x,y,colour,linewidth,linestyle,params); }
	public void scatter(float[] x, float[] y) {
		gca().scatter(x,y); }
	public void scatter(float[] x, float[] y, int colour, float iconsize, String symbol, Object... params) {
		gca().scatter(x,y,colour,iconsize,symbol,params); }
	public void scatter(double[] x, double[] y) {
		gca().scatter(x,y); }
	public void scatter(double[] x, double[] y, int colour, double iconsize, String symbol, Object... params) {
		gca().scatter(x,y,colour,iconsize,symbol,params); }
	public void contour(double[] x, double[] y, double[][] z) {
		gca().contour(x,y,z); }
	public void contour(double[] x, double[] y, double[][] z, int levels, Object... params) {
		gca().contour(x, y, z, levels, params); }

	public void predefImgShow(String predefined_images) {
		gca().predefImgShow(predefined_images);
	}
	public void imgShow(PImage img) {
		gca().imgShow(img);
	}

	public void setRange(float xmin, float xmax, float ymin, float ymax) {
		gca().setRange(xmin, xmax, ymin, ymax);
	}
	public void setRange(double xmin, double xmax, double ymin, double ymax) {
		gca().setRange(xmin, xmax, ymin, ymax);
	}
	public void setFont(PFont font) {
		gca().setFont(font); }
}

