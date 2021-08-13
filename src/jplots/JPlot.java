package jplots;

import jplots.color.JColourbar;
import jplots.shapes.JPlotShape;
import jplots.transform.JProjection;
import processing.core.*;
import processing.data.IntList;

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

public class JPlot {

	/** the version of the Library. */
	public final static String VERSION = "##library.prettyVersion##";

	private static boolean hasWelcomed = false;
	public static double dpi = 300d;

	private PApplet myParent;
	private PGraphics plotImg;

	private boolean img_is_created, useDebug;
	private int width, height, lastAxisNum;
	private int n_cols, n_rows;
	private JAxis[] axes;


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
	public JPlot(PApplet applet) {
		this.myParent = applet;
		this.width = 0;
		this.height = 0;
		this.img_is_created = false;
		this.useDebug = false;
		//this.plotImg = applet.createGraphics(16, 16);
		if(!hasWelcomed)
			this.welcome();
	}


	//************************************
	//**** STATIC ************************
	//************************************


	//************************************
	//**** PURE PPLOT ********************
	//************************************
	
	/**
	 * starts the new figure/plot width standard width and height about 1.7067inch x 1.7067inch and 300dpi
	 * 
	 * @return the JPlot object
	 */
	public JPlot figure() { return subplots(1.7066666667d, 1.7066666667d, 1, 1); }
	/**
	 * starts the new figure/plot with specified width and height in inch
	 * 
	 * @param width  width of figure in inch
	 * @param height height of figure in inch
	 * @return the JPlot object
	 */
	public JPlot figure(double width, double height) { return subplots(width, height, 1, 1); }
	/**
	 * starts the new figure/plot width standard width and height about 1.7067inch x 1.7067inch and 300dpi
	 * 
	 * @param nrows  number of rows for array of diagrams
	 * @param ncols  number of columns for array of diagrams
	 * @return the JPlot object
	 */
	public JPlot subplots(int nrows, int ncols) { return subplots(1.7066666667d, 1.7066666667d, nrows, ncols); }
	/**
	 * starts the new figure/plot with specified width and height in inch
	 * 
	 * @param width  width of figure in inch
	 * @param height height of figure in inch
	 * @param nrows  number of rows for array of diagrams
	 * @param ncols  number of columns for array of diagrams
	 * @return the JPlot object
	 */
	public JPlot subplots(double width, double height, int nrows, int ncols) {
		this.width  = (int) (dpi * width);
		this.height = (int) (dpi * height);
		this.axes   = new JAxis[nrows*ncols];
		if(useDebug)
			System.out.println(
				"[DEBUG] PPlot: width/height: "+this.width+"/"+this.height+"px ("+width+"/"+height+"inch)\n"+
				"               axes:         "+ncols+"x"+nrows);
		double ws = 1d / (0.25d+(ncols-1)*1.5d+1.25d);
		double hs = 1d / (0.25d+(nrows-1)*1.5d+1.25d);
		int pawid = (int) (dpi * width * ws);
		int pahei = (int) (dpi * height * hs);
		n_rows = nrows;
		n_cols = ncols;
		for(int r=0; r<nrows; r++) {
			int py = (int) (dpi * height * (0.25d+1.5d*r) * hs);
			for(int c=0; c<ncols; c++) {
				int px = (int) (dpi * width * (0.25d+1.5d*c) * ws);
				this.axes[r*ncols+c] = new JAxis(this, px, py, pawid, pahei);
			}
		}
		this.img_is_created = false;
		return this;
	}
	/**
	 * add new JAxis object as a subplot in the JPlot object
	 * 
	 * @param pos_x  left horizontal position relative to figure width
	 * @param pos_y  top vertical position relative to figure height
	 * @param width  width relative to the figure width
	 * @param height height relative to the figure height
	 * @return the JPlot object
	 */
	public JPlot addSubplot(double pos_x, double pos_y, double width, double height) {
		JAxis[] tempa = new JAxis[axes.length+1];
		for(int a=0; a<axes.length; a++)
			tempa[a] = axes[a];
		tempa[axes.length] = new JAxis(this,
									   (int)(pos_x*this.width),(int)(pos_y*this.height),
									   (int)(width*this.width+0.5d),(int)(height*this.height+0.5d));
		axes = new JAxis[tempa.length];
		for(int a=0; a<axes.length; a++)
			axes[a] = tempa[a];
		lastAxisNum = axes.length-1;
		if(useDebug)
			System.out.println("[DEBUG] addSubplot: "+axes.length+" PAxis-objects");
		return this;
	}
	/**
	 * add new JAxis object as a subplot in the JPlot object
	 * 
	 * @param axis new JAxis object
	 * @return the JPlot object
	 */
	public JPlot addSubplot(JAxis axis) {
		JAxis[] tempa = new JAxis[axes.length+1];
		for(int a=0; a<axes.length; a++)
			tempa[a] = axes[a];
		tempa[axes.length] = axis;
		axes = new JAxis[tempa.length];
		for(int a=0; a<axes.length; a++)
			axes[a] = tempa[a];
		lastAxisNum = axes.length-1;
		if(useDebug)
			System.out.println("[DEBUG] addSubplot: "+axes.length+" PAxis-objects");
		return this;
	}
	/**
	 * removes subplot with indices a
	 * 
	 * @param a  indicex of subplot
	 * @return this JPlot object
	 */
	public JPlot removeSubplots(int... a) {
		IntList il = new IntList(a);
		for(int ai=il.size()-1; ai>=0; ai--) {
			boolean isDublicate = false;
			for(int aj=0; aj<ai && !isDublicate; aj++)
				if(il.get(ai)==il.get(aj))
					isDublicate = true;
			if(isDublicate || il.get(ai)>=n_rows*n_cols || il.get(ai)<0)
				il.remove(ai);
		}
		il.sort();
		if(useDebug)
			System.out.println(il.toString());
		int[] srcIdx = new int[axes.length-il.size()];
		for(int ai=0; ai<srcIdx.length; ai++)
			srcIdx[ai] = ai;
		for(int aj: il)
			for(int ai=0; ai<srcIdx.length; ai++)
				if(srcIdx[ai]>=aj)
					srcIdx[ai]++;
		JAxis[] tempa = new JAxis[srcIdx.length];
		for(int ai=0; ai<tempa.length; ai++) {
			tempa[ai] = axes[srcIdx[ai]];
		}
		axes = new JAxis[tempa.length];
		for(int ai=0; ai<axes.length; ai++)
			axes[ai] = tempa[ai];
		lastAxisNum = axes.length-1;
		if(useDebug)
			System.out.println("[DEBUG] addSubplot: "+axes.length+" PAxis-objects");
		return this;
	}
	/**
	 * set spacing between different subplots in the plot
	 * 
	 * @param wspacing horizontal spacing as fraction of width of subplots
	 * @param hspacing vertical spacing as fraction of height of subplots
	 * @return this JPlot object
	 */
	public JPlot setSpacing(double wspacing, double hspacing) {
		if(axes.length!=n_rows*n_cols) {
			System.err.println("JAxis-configuration altered, cannot reorder images!");
			return this;
		}
		double ws = 1d / (0.5d*wspacing+(n_cols-1)*(1.0d+wspacing)+1.0d+0.5d*wspacing);
		double hs = 1d / (0.5d*hspacing+(n_rows-1)*(1.0d+hspacing)+1.0d+0.5d*hspacing);
		int pawid = (int) (width * ws+0.1d);
		int pahei = (int) (height * hs+0.1d);
		for(int r=0; r<n_rows; r++) {
			int py = (int) (height * (0.5d*hspacing+(1.0d+hspacing)*r) * hs + 0.1d);
			for(int c=0; c<n_cols; c++) {
				int px = (int) (width * (0.5d*wspacing+(1.0d+wspacing)*c) * ws + 0.1d);
				this.axes[r*n_cols+c].setPositionAndSize(px,py,pawid,pahei);
			}
		}
		return this;
	}
	
	public JPlot setGeoProjection(JProjection proj) {
		gca().setGeoProjection(proj);
		return this; }
	public JPlot setGeoProjection(JProjection proj, int axis_num) {
		ga(axis_num).setGeoProjection(proj);
		return this; }
	
	/**
	 * creates image of the PPlot -- the figure
	 */
	private JPlot createImage() {
		if(useDebug)
			System.out.println("[DEBUG] Start image-creation ...");
		if(plotImg==null) {
			plotImg = myParent.createGraphics(width, height, PApplet.P2D);
		} else if(plotImg.width!=width || plotImg.height!=height) {
			plotImg = myParent.createGraphics(width, height, PApplet.P2D);
		}
		if(useDebug) {
			System.out.println("[DEBUG]   - create PGraphics-object "+plotImg.width+"x"+plotImg.height+"px");
			System.out.println("[DEBUG]   - begin drawing");
		}
		plotImg.beginDraw();
		plotImg.clear();
		plotImg.textSize(200);
		createImage(plotImg);
		if(useDebug)
			System.out.println("[DEBUG]   - end drawing");
		plotImg.endDraw();
		img_is_created = true;
		return this;
	}
	private JPlot createImage(PGraphics g) {
		if(useDebug)
			System.out.println("[DEBUG]   - add "+axes.length+" subplots to plotting queue");
		for(JAxis ax: axes) {
			if(ax instanceof JColourbar)
				continue;
			JPlotShape plotShp = ax.createPlot(myParent,width,height);
			if(ax.getFont()!=null)
				g.textFont(ax.getFont());
			plotShp.draw(this, g);
		}
		for(JAxis ax: axes) {
			if(ax instanceof JColourbar) {
				JPlotShape plotShp = ax.createPlot(myParent,width,height);
				if(ax.getFont()!=null)
					g.textFont(ax.getFont());
				plotShp.draw(this, g);
			}
		}
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
//		if(useDebug)
//			System.out.println("[DEBUG] check: \"plotImg\"="+plotImg);
		return plotImg;
	}
	public void draw(PGraphics g) {
		createImage(g);
	}
	public PGraphics getGraphic() {
		return plotImg;
	}
	public JPlot redraw(boolean redraw) {
		img_is_created = redraw?false:img_is_created;
		return this;
	}

	/**
	 * removes all plotting infos
	 * also all configurations except number and position of axes will be reseted
	 */
	public JPlot clear() {
		int cb_count = 0;
		for(JAxis ax: axes) {
			if(ax instanceof JColourbar)
				cb_count++;
			ax.clear();
		}
		if(cb_count>0) {
			JAxis[] tempa = new JAxis[axes.length-cb_count];
			cb_count = 0;
			for(int a=0; a<axes.length; a++) {
				if(axes[a] instanceof JColourbar) {
					cb_count++;
					continue;
				}
				tempa[a-cb_count] = axes[a];
			}
			axes = new JAxis[tempa.length];
			for(int a=0; a<axes.length; a++)
				axes[a] = tempa[a];
			lastAxisNum = axes.length-1;
		}
		redraw(true);
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
	public int[] getSize() { return new int[] {width, height}; }
	
	public JAxis[] ga() { return axes; }
	public JAxis gca() { return axes[lastAxisNum]; }
	public JAxis ga(int a) {
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
		hasWelcomed = true;
	}


	//************************************
	//**** PASS ON ***********************
	//************************************
	
	public void hline(float y) {
		gca().axhline(y); }
	public void hline(float y, int colour, float linewidth, String linestyle) {
		gca().axhline(y, colour, linewidth, linestyle); }
	public void hline(double y) {
		gca().axhline(y); }
	public void hline(double y, int colour, double linewidth, String linestyle) {
		gca().axhline(y, colour, linewidth, linestyle); }
	public void vline(float x) {
		gca().axvline(x); }
	public void vline(float x, int colour, float linewidth, String linestyle) {
		gca().axvline(x, colour, linewidth, linestyle); }
	public void vline(double x) {
		gca().axvline(x); }
	public void vline(double x, int colour, double linewidth, String linestyle) {
		gca().axvline(x, colour, linewidth, linestyle); }
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
	public void contour(float[] x, float[] y, float[][] z) {
		gca().contour(x,y,z); }
	public void contour(float[] x, float[] y, float[][] z, int levels, Object... params) {
		gca().contour(x, y, z, levels, params); }
	public void contour(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		gca().contour(x, y, z, levels, params); }
	public void contour(float[][] x, float[][] y, float[][] z) {
		gca().contour(x,y,z); }
	public void contour(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		gca().contour(x, y, z, levels, params); }
	public void contour(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		gca().contour(x, y, z, levels, params); }
	public void contour(double[] x, double[] y, double[][] z) {
		gca().contour(x,y,z); }
	public void contour(double[] x, double[] y, double[][] z, int levels, Object... params) {
		gca().contour(x, y, z, levels, params); }
	public void contour(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		gca().contour(x, y, z, levels, params); }
	public void contour(double[][] x, double[][] y, double[][] z) {
		gca().contour(x,y,z); }
	public void contour(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		gca().contour(x, y, z, levels, params); }
	public void contour(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		gca().contour(x, y, z, levels, params); }
	public void contourf(float[] x, float[] y, float[][] z) {
		gca().contourf(x,y,z); }
	public void contourf(float[] x, float[] y, float[][] z, int levels, Object... params) {
		gca().contourf(x, y, z, levels, params); }
	public void contourf(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		gca().contourf(x, y, z, levels, params); }
	public void contourf(float[][] x, float[][] y, float[][] z) {
		gca().contourf(x,y,z); }
	public void contourf(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		gca().contourf(x, y, z, levels, params); }
	public void contourf(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		gca().contourf(x, y, z, levels, params); }
	public void contourf(double[] x, double[] y, double[][] z) {
		gca().contourf(x,y,z); }
	public void contourf(double[] x, double[] y, double[][] z, int levels, Object... params) {
		gca().contourf(x, y, z, levels, params); }
	public void contourf(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		gca().contourf(x, y, z, levels, params); }
	public void contourf(double[][] x, double[][] y, double[][] z) {
		gca().contourf(x,y,z); }
	public void contourf(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		gca().contourf(x, y, z, levels, params); }
	public void contourf(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		gca().contourf(x, y, z, levels, params); }
	public void contourp(float[] x, float[] y, float[][] z) {
		gca().contourp(x,y,z); }
	public void contourp(float[] x, float[] y, float[][] z, int levels, Object... params) {
		gca().contourp(x, y, z, levels, params); }
	public void contourp(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		gca().contourp(x, y, z, levels, params); }
	public void contourp(float[][] x, float[][] y, float[][] z) {
		gca().contourp(x,y,z); }
	public void contourp(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		gca().contourp(x, y, z, levels, params); }
	public void contourp(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		gca().contourp(x, y, z, levels, params); }
	public void contourp(double[] x, double[] y, double[][] z) {
		gca().contourp(x,y,z); }
	public void contourp(double[] x, double[] y, double[][] z, int levels, Object... params) {
		gca().contourp(x, y, z, levels, params); }
	public void contourp(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		gca().contourp(x, y, z, levels, params); }
	public void contourp(double[][] x, double[][] y, double[][] z) {
		gca().contourp(x,y,z); }
	public void contourp(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		gca().contourp(x, y, z, levels, params); }
	public void contourp(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		gca().contourp(x, y, z, levels, params); }

	public void coastLines() {
		gca().coastLines(); }
	public void coastLines(int resolution) {
		gca().coastLines(resolution); }

	public void colourbar() {
		colourbar(gca(), ""); }
	public void colourbar(JAxis axis) {
		colourbar(axis, ""); }
	public void colourbar(String name) {
		colourbar(gca(), name); }
	public void colourbar(JAxis axis, String name) {
		addSubplot(new JColourbar(axis, name)); }

	public void predefImgShow(String predefined_images) {
		gca().predefImgShow(predefined_images);
	}
	public void imgShow(PImage img) {
		gca().imgShow(img);
	}

	public void setXRange(float xmin, float xmax) {
		gca().setXRange(xmin, xmax); }
	public void setXRange(double xmin, double xmax) {
		gca().setXRange(xmin, xmax); }
	public void setYRange(float ymin, float ymax) {
		gca().setYRange(ymin, ymax); }
	public void setYRange(double ymin, double ymax) {
		gca().setYRange(ymin, ymax); }
	public void setRange(float xmin, float xmax, float ymin, float ymax) {
		gca().setRange(xmin, xmax, ymin, ymax);
	}
	public void setRange(double xmin, double xmax, double ymin, double ymax) {
		gca().setRange(xmin, xmax, ymin, ymax);
	}
	public void setFont(PFont font) {
		gca().setFont(font); }
	public void setXTitle(String xtitle) {
		gca().setXTitle(xtitle); }
	public void setYTitle(String ytitle) {
		gca().setYTitle(ytitle); }
	public void setTitle(String _title) {
		gca().setTitle(_title); }
}

