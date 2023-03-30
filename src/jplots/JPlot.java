package jplots;

import org.opengis.referencing.crs.CoordinateReferenceSystem;

import jplots.axes.JSingleAxis;
import jplots.axes.JAxis;
import jplots.axes.JColourbar;
import jplots.axes.JMultiAxis;
import jplots.maths.JDPolygon;
import jplots.shapes.JPlotShape;
import jplots.transform.JProjection;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PFont;
import processing.core.PGraphics;
import processing.core.PImage;
import processing.data.IntList;

/**
 * This is a template class and can be used to start a new processing Library.
 * Make sure you rename this class as well as the name of the example package
 * 'template' to your own Library naming convention.
 *
 * (the tag example followed by the name of an example included in folder
 * 'examples' will automatically include the example in the javadoc.)
 *
 * @example SimpleXY
 */

public class JPlot {
	
	/** the version of the Library. */
	public final static String VERSION = "##library.prettyVersion##";
	
	public static double dpi = 300d;
	public static boolean supportLatex = false;
	
	static {
		System.out.println("##library.name## ##library.prettyVersion## by ##author##");
		try {
			Class.forName("latex.PTeX");
			supportLatex = true;
		} catch (ClassNotFoundException cnfe) {
			supportLatex = false;
		}
	}
	
	private PApplet myParent;
	private PGraphics plotImg;

	private boolean img_is_created, useDebug;
	private int width, height, lastAxisNum;
	private int n_cols, n_rows, currAxNum;
	private JAxis[] axes;

	// ************************************
	// **** CONSTRUCTOR *******************
	// ************************************

	/**
	 * a Constructor, usually called in the setup() method in your sketch to
	 * initialize and start the Library.
	 *
	 * @example SimpleXY
	 * @param applet the parent PApplet
	 */
	public JPlot(PApplet applet) {
		this.myParent = applet;
		this.width = 0;
		this.height = 0;
		this.img_is_created = false;
		this.useDebug = false;
		this.currAxNum = -1;
		// this.plotImg = applet.createGraphics(16, 16);
		System.out.println("Support LATEX: "+supportLatex);
	}

	// ************************************
	// **** STATIC ************************
	// ************************************
	
	public void activateLatex(boolean b) {
		if(b) {
			try {
				Class.forName("latex.PTeX");
				supportLatex = true;
				System.out.println("Now all titles and anotation are rendered with Latex.");
			} catch (ClassNotFoundException cnfe) {
				supportLatex = false;
			}
		} else {
			supportLatex = false;
			System.out.println("Rendering with LaTex is deactivated.");
		}
	}
	
	// ************************************
	// **** PURE PPLOT ********************
	// ************************************

	/**
	 * starts the new figure/plot width standard width and height about 1.7067inch x
	 * 1.7067inch and 300dpi
	 *
	 * @return the JPlot object
	 */
	public JPlot figure() {
		return subplots(1.7066666667d, 1.7066666667d, 1, 1);
	}
	/**
	 * starts the new figure/plot width standard width and height about 1.7067inch x
	 * 1.7067inch and 300dpi
	 *
	 * @param multi set axis as a JMultiAxis or normal JAxis.
	 * @return the JPlot object
	 */
	public JPlot figure(boolean multi, int Xparts, int Yparts) {
		Object[] arg = null;
		if(multi) arg = new Object[] {"multi", Xparts, Yparts};
		return subplots(1.7066666667d, 1.7066666667d, 1, 1, arg);
	}
	/**
	 * starts the new figure/plot with specified width and height in inch
	 *
	 * @param width  width of figure in inch
	 * @param height height of figure in inch
	 * @return the JPlot object
	 */
	public JPlot figure(double width, double height) {
		return subplots(width, height, 1, 1);
	}
	/**
	 * starts the new figure/plot with specified width and height in inch
	 *
	 * @param width  width of figure in inch
	 * @param height height of figure in inch
	 * @param multi  set axis as a JMultiAxis or normal JAxis
	 * @return the JPlot object
	 */
	public JPlot figure(double width, double height, boolean multi, int Xparts, int Yparts) {
		Object[] arg = null;
		if(multi) arg = new Object[] {"multi", Xparts, Yparts};
		return subplots(width, height, 1, 1, arg);
	}
	
	/**
	 * starts the new figure/plot width standard width and height about 1.7067inch x
	 * 1.7067inch and 300dpi
	 *
	 * @param nrows number of rows for array of diagrams
	 * @param ncols number of columns for array of diagrams
	 * @return the JPlot object
	 */
	public JPlot subplots(int nrows, int ncols) {
		return subplots(1.7066666667d, 1.7066666667d, nrows, ncols);
	}
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
		return subplots(width, height, nrows, ncols, (Object) null);
	}
	/**
	 * starts the new figure/plot with specified width and height in inch
	 *
	 * @param width  width of figure in inch
	 * @param height height of figure in inch
	 * @param nrows  number of rows for array of diagrams
	 * @param ncols  number of columns for array of diagrams
	 * @param kwargs additional keyword arguments
	 * @return the JPlot object
	 */
	public JPlot subplots(double width, double height, int nrows, int ncols, Object... kwargs) {
		this.width = (int) (dpi * width);
		this.height = (int) (dpi * height);
		this.axes = new JAxis[nrows * ncols];
		if (useDebug)
			System.out.println("[DEBUG] PPlot: width/height: " + this.width + "/" + this.height + "px (" + width + "/"
					+ height + "inch)\n" + "               axes:         " + ncols + "x" + nrows);
		boolean isMulti = false;
		int xparts=0, yparts=0;
		boolean isShareX = false;
		boolean isShareY = false;
		if(kwargs!=null) {
			for(int o=0; o<kwargs.length; o++) {
				if(!(kwargs[o] instanceof String))
					continue;
				String s = ((String)kwargs[o]).toLowerCase();
				if(s.equals("sharex")) isShareX = true;
				if(s.equals("sharey")) isShareY = true;
				if(s.equals("multi")) { isMulti = true; xparts = (int)kwargs[++o]; yparts = (int)kwargs[++o]; }
			}
		}
		if(isShareX||isShareY) {
			if(isMulti) System.err.println("Cannot set JMultiAxis with shared x-/y-axes.");
			isMulti = false;
		}
		double ws = 1d / (0.25d + (ncols - 1) * 1.5d + 1.25d);
		double hs = 1d / (0.25d + (nrows - 1) * 1.5d + 1.25d);
		int pawid = (int) (dpi * width * ws);
		int pahei = (int) (dpi * height * hs);
		n_rows = nrows;
		n_cols = ncols;
		for (int r = 0; r < nrows; r++) {
			int py = (int) (dpi * height * (0.25d + 1.5d * r) * hs);
			for (int c = 0; c < ncols; c++) {
				int px = (int) (dpi * width * (0.25d + 1.5d * c) * ws);
				if(isMulti)
					this.axes[r * ncols + c] = new JMultiAxis(this, px, py, pawid, pahei, xparts, yparts);
				else
					this.axes[r * ncols + c] = new JSingleAxis(this, px, py, pawid, pahei);
			}
		}
		this.img_is_created = false;
		if(isShareX)
			for (int col = 0; col < n_cols; col++)
				for (int row = 1; row < n_rows; row++)
					((JSingleAxis)this.axes[col]).addSharedAxis('x', (JSingleAxis)this.axes[row * n_cols + col]);
		if(isShareY)
			for (int row = 0; row < n_rows; row++)
				for (int col = 1; col < n_cols; col++)
					((JSingleAxis)this.axes[row * n_cols]).addSharedAxis('y', (JSingleAxis)this.axes[row * n_cols + col]);
		this.currAxNum = nrows*ncols-1;
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
		JAxis[] tempa = new JAxis[axes.length + 1];
		for (int a = 0; a < axes.length; a++)
			tempa[a] = axes[a];
		tempa[axes.length] = new JSingleAxis(this, (int) (pos_x * this.width), (int) (pos_y * this.height),
				(int) (width * this.width + 0.5d), (int) (height * this.height + 0.5d));
		axes = new JAxis[tempa.length];
		for (int a = 0; a < axes.length; a++)
			axes[a] = tempa[a];
		lastAxisNum = axes.length - 1;
		currAxNum = lastAxisNum;
		if (useDebug)
			System.out.println("[DEBUG] addSubplot: " + axes.length + " PAxis-objects");
		return this;
	}
	/**
	 * add new JAxis object as a subplot in the JPlot object
	 *
	 * @param axis new JAxis object
	 * @return the JPlot object
	 */
	public JPlot addSubplot(JAxis axis) {
		JAxis[] tempa = new JAxis[axes.length + 1];
		for (int a = 0; a < axes.length; a++)
			tempa[a] = axes[a];
		tempa[axes.length] = axis;
		axes = new JAxis[tempa.length];
		for (int a = 0; a < axes.length; a++)
			axes[a] = tempa[a];
		lastAxisNum = axes.length - 1;
		if(!(axes[lastAxisNum] instanceof JColourbar))
			currAxNum = lastAxisNum;
		if (useDebug)
			System.out.println("[DEBUG] addSubplot: " + axes.length + " PAxis-objects");
		return this;
	}

	/**
	 * removes subplot with indices a
	 *
	 * @param a indicex of subplot
	 * @return this JPlot object
	 */
	public JPlot removeSubplots(int... a) {
		IntList il = new IntList(a);
		for (int ai = il.size() - 1; ai >= 0; ai--) {
			boolean isDublicate = false;
			for (int aj = 0; aj < ai && !isDublicate; aj++)
				if (il.get(ai) == il.get(aj))
					isDublicate = true;
			if (isDublicate || il.get(ai) >= n_rows * n_cols || il.get(ai) < 0)
				il.remove(ai);
		}
		il.sort();
		if (useDebug)
			System.out.println(il.toString());
		int[] srcIdx = new int[axes.length - il.size()];
		for (int ai = 0; ai < srcIdx.length; ai++)
			srcIdx[ai] = ai;
		for (int aj : il)
			for (int ai = 0; ai < srcIdx.length; ai++)
				if (srcIdx[ai] >= aj)
					srcIdx[ai]++;
		JAxis[] tempa = new JAxis[srcIdx.length];
		for (int ai = 0; ai < tempa.length; ai++) {
			tempa[ai] = axes[srcIdx[ai]];
		}
		axes = new JAxis[tempa.length];
		for (int ai = 0; ai < axes.length; ai++)
			axes[ai] = tempa[ai];
		lastAxisNum = axes.length - 1;
		if (useDebug)
			System.out.println("[DEBUG] addSubplot: " + axes.length + " PAxis-objects");
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
		return setSpacing(wspacing, 0.25d, hspacing, 0.25d);
	}

	/**
	 * set spacing between different subplots in the plot
	 *
	 * @param wspacing horizontal spacing as fraction of width of subplots
	 * @param wbounds  width of space to canvas boundary as fraction of subplot
	 *                 width
	 * @param hspacing vertical spacing as fraction of height of subplots
	 * @param hbounds  height of space to canvas boundary as fraction of subplot
	 *                 height
	 * @return this JPlot object
	 */
	public JPlot setSpacing(double wspacing, double wbounds, double hspacing, double hbounds) {
		if (axes.length != n_rows * n_cols) {
			System.err.println("JAxis-configuration altered, cannot reorder images!");
			return this;
		}
		double ww = wspacing * (n_cols - 1) + wbounds * 2 + n_cols;
		int awid = (int) (1d * width / ww + 0.5d), ws = (int) ((1d + wspacing) * width / ww + 0.2d),
				wo = (width - ws * (n_cols - 1) - awid) / 2;
		double hh = hspacing * (n_rows - 1) + hbounds * 2 + n_rows;
		int ahei = (int) (1d * height / hh + 0.5d), hs = (int) ((1d + hspacing) * height / hh + 0.2d),
				ho = (height - hs * (n_rows - 1) - ahei) / 2;
		for (int r = 0; r < n_rows; r++) {
			int py = ho + r * hs;
			for (int c = 0; c < n_cols; c++) {
				int px = wo + c * ws;
				this.axes[r * n_cols + c].setPositionAndSize(px, py, awid, ahei);
			}
		}
		return this;
	}

	public JPlot setGeoProjection(JProjection proj) {
		gcja().setGeoProjection(proj);
		return this;
	}

	public JPlot setGeoProjection(JProjection proj, int axis_num) {
		if(ga(axis_num) instanceof JSingleAxis)
			((JSingleAxis)ga(axis_num)).setGeoProjection(proj);
		else
			System.err.println("This axis-object cannot have a GeoProjection.");
		return this;
	}
	
	/**
	 * creates image of the PPlot -- the figure
	 */
	private JPlot createImage() {
		if (useDebug)
			System.out.println("[DEBUG] Start image-creation ...");
		if (plotImg == null) {
			plotImg = myParent.createGraphics(width, height, PConstants.P2D);
			plotImg.beginDraw();
			plotImg.clear();
			plotImg.endDraw();
		} else if (plotImg.width != width || plotImg.height != height) {
			plotImg = myParent.createGraphics(width, height, PConstants.P2D);
			plotImg.beginDraw();
			plotImg.clear();
			plotImg.endDraw();
		}
		if (useDebug) {
			System.out.println("[DEBUG]   - create PGraphics-object " + plotImg.width + "x" + plotImg.height + "px");
			System.out.println("[DEBUG]   - begin drawing");
		}
		onlyPlotAxes();
		plotImg.beginDraw();
		plotImg.clear();
		createImage(plotImg);
		if (useDebug)
			System.out.println("[DEBUG]   - end drawing");
		plotImg.endDraw();
		img_is_created = true;
		return this;
	}

	private void onlyPlotAxes() {
		boolean tempDebug = isDebug();
		debug(false);
		plotImg.beginDraw();
		plotImg.clear();
		plotImg.textSize(200);
		for (JAxis ax : axes) {
			if (ax instanceof JColourbar)
				continue;
			JPlotShape plotShp = ax.createPlotOnlyAxes(myParent, width, height);
			if (ax.getFont() != null)
				plotImg.textFont(ax.getFont());
			plotShp.draw(this, plotImg);
		}
		for (JAxis ax : axes) {
			if (ax instanceof JColourbar) {
				JPlotShape plotShp = ax.createPlotOnlyAxes(myParent, width, height);
				if (ax.getFont() != null)
					plotImg.textFont(ax.getFont());
				plotShp.draw(this, plotImg);
			}
		}
		plotImg.endDraw();
		debug(tempDebug);
	}
	private JPlot createImage(PGraphics g) {
		if (useDebug)
			System.out.println("[DEBUG]   - add " + axes.length + " subplots to plotting queue");
		g.textSize(200);
//		PrintWriter pw = myParent.createWriter("./shapestack_"+getDate()+".txt");
		for (JAxis ax : axes) {
			if (ax instanceof JColourbar)
				continue;
//			pw.println("\nAxes: "+ax);
			JPlotShape plotShp = ax.createPlot(myParent, width, height);
//			plotShp.printStack(pw,"");
			if (ax.getFont() != null)
				g.textFont(ax.getFont());
			plotShp.draw(this, g);
		}
		for (JAxis ax : axes) {
			if (ax instanceof JColourbar) {
//				pw.println("\nAxes: "+ax);
				JPlotShape plotShp = ax.createPlot(myParent, width, height);
//				plotShp.printStack(pw,"");
				if (ax.getFont() != null)
					g.textFont(ax.getFont());
				plotShp.draw(this, g);
			}
		}
//		pw.flush();
//		pw.close();
		return this;
	}

	/**
	 * gives the figure<br>
	 * the figure will be created automatically, if it isn't allready
	 *
	 * @return the figure as PImage
	 */
	public PImage show() {
		if (!img_is_created) {
			if (useDebug)
				System.out.println("[DEBUG] start image creation, because it ws not done before or changes happend.");
			createImage();
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
		img_is_created = redraw ? false : img_is_created;
		return this;
	}

	public boolean hasBeenDrawn() {
		return img_is_created;
	}

	/**
	 * removes all plotting infos also all configurations except number and position
	 * of axes will be reseted
	 *
	 * @return this JPlot object
	 */
	public JPlot clear() {
		int cb_count = 0;
		for (JAxis ax : axes) {
			if (ax instanceof JColourbar)
				cb_count++;
			ax.clear();
		}
		if (cb_count > 0) {
			JAxis[] tempa = new JAxis[axes.length - cb_count];
			cb_count = 0;
			for (int a = 0; a < axes.length; a++) {
				if (axes[a] instanceof JColourbar) {
					cb_count++;
					continue;
				}
				tempa[a - cb_count] = axes[a];
			}
			axes = new JAxis[tempa.length];
			for (int a = 0; a < axes.length; a++)
				axes[a] = tempa[a];
			lastAxisNum = axes.length - 1;
		}
		redraw(true);
		return this;
	}

	// ....
	/**
	 * set debugmode on/off
	 *
	 * @param bdebug use true for debugmode on
	 */
	public void debug(boolean bdebug) {
		useDebug = bdebug;
	}

	public boolean isDebug() {
		return useDebug;
	}

	// ************************************
	// **** GETTER ************************
	// ************************************

	public PApplet getApplet() {
		return myParent;
	}

	public int[] getSize() {
		return new int[] { width, height };
	}

	public JAxis[] ga() {
		return axes;
	}

	public JAxis gca() {
		if(currAxNum<0) return null;
		return axes[currAxNum];
	}
	private JSingleAxis gcja() {
		if(currAxNum<0) return null;
		return (JSingleAxis) axes[currAxNum];
	}
	private JMultiAxis gcma() {
		if(currAxNum<0) return null;
		return (JMultiAxis) axes[currAxNum];
	}

	public JAxis ga(int a) {
		if (a < -axes.length || a > axes.length-1)
			throw new IndexOutOfBoundsException(
					"Number out of range of axis count: " + a + " <> {0..." + axes.length + "}");
		currAxNum = a + (a < 0 ? axes.length : 0);
		return axes[currAxNum];
	}

	// ************************************
	// **** PRIVATE ***********************
	// ************************************
	
	public static String getDate() {
		String ds = ""+
				PApplet.nf(PApplet.year(),4)+
				PApplet.nf(PApplet.month(),2)+
				PApplet.nf(PApplet.day(),2)+
				"-"+
				PApplet.nf(PApplet.hour(),2)+
				PApplet.nf(PApplet.minute(),2)+
				PApplet.nf(PApplet.second(),2);
		return ds;
	}
	
	public int getNumColumns() {
		return n_cols;
	}
	public int getNumRows() {
		return n_rows;
	}
	
	// ************************************
	// **** PASS ON ***********************
	// ************************************

	public void hline(float y) {
		gcja().axhline(y);
	}
	public void hline(float y, int colour, float linewidth, String linestyle) {
		gcja().axhline(y, colour, linewidth, linestyle);
	}
	public void hline(double y) {
		gcja().axhline(y);
	}
	public void hline(double y, int colour, double linewidth, String linestyle) {
		gcja().axhline(y, colour, linewidth, linestyle);
	}

	public void vline(float x) {
		gcja().axvline(x);
	}
	public void vline(float x, int colour, float linewidth, String linestyle) {
		gcja().axvline(x, colour, linewidth, linestyle);
	}
	public void vline(double x) {
		gcja().axvline(x);
	}
	public void vline(double x, int colour, double linewidth, String linestyle) {
		gcja().axvline(x, colour, linewidth, linestyle);
	}
	
	public void plot(float[] x, float[] y) {
		gcja().plot(x, y);
	}
	public void plot(float[] x, float[] y, int colour, float linewidth, String linestyle, Object... params) {
		gcja().plot(x, y, colour, linewidth, linestyle, params);
	}
	public void plot(double[] x, double[] y) {
		gcja().plot(x, y);
	}
	public void plot(double[] x, double[] y, int colour, double linewidth, String linestyle, Object... params) {
		gcja().plot(x, y, colour, linewidth, linestyle, params);
	}
	public void plot(int subaxis, float[] x, float[] y) {
		gcma().plot(subaxis, x, y);
	}
	public void plot(int subaxis, float[] x, float[] y, int colour, float linewidth, String linestyle, Object... params) {
		gcma().plot(subaxis, x, y, colour, linewidth, linestyle, params);
	}
	public void plot(int subaxis, double[] x, double[] y) {
		gcma().plot(subaxis, x, y);
	}
	public void plot(int subaxis, double[] x, double[] y, int colour, double linewidth, String linestyle, Object... params) {
		gcma().plot(subaxis, x, y, colour, linewidth, linestyle, params);
	}
	
	public void scatter(float[] x, float[] y) {
		gcja().scatter(x, y);
	}
	public void scatter(float[] x, float[] y, int colour, float iconsize, String symbol, Object... params) {
		gcja().scatter(x, y, colour, iconsize, symbol, params);
	}
	public void scatter(double[] x, double[] y) {
		gcja().scatter(x, y);
	}
	public void scatter(double[] x, double[] y, int colour, double iconsize, String symbol, Object... params) {
		gcja().scatter(x, y, colour, iconsize, symbol, params);
	}
	public void scatter(int subaxis, float[] x, float[] y) {
		gcma().scatter(subaxis, x, y);
	}
	public void scatter(int subaxis, float[] x, float[] y, int colour, float iconsize, String symbol, Object... params) {
		gcma().scatter(subaxis, x, y, colour, iconsize, symbol, params);
	}
	public void scatter(int subaxis, double[] x, double[] y) {
		gcma().scatter(subaxis, x, y);
	}
	public void scatter(int subaxis, double[] x, double[] y, int colour, double iconsize, String symbol, Object... params) {
		gcma().scatter(subaxis, x, y, colour, iconsize, symbol, params);
	}
	
	public void contour(float[] x, float[] y, float[][] z) {
		gcja().contour(x, y, z);
	}
	public void contour(float[] x, float[] y, float[][] z, int levels, Object... params) {
		gcja().contour(x, y, z, levels, params);
	}
	public void contour(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		gcja().contour(x, y, z, levels, params);
	}
	public void contour(float[][] x, float[][] y, float[][] z) {
		gcja().contour(x, y, z);
	}
	public void contour(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		gcja().contour(x, y, z, levels, params);
	}
	public void contour(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		gcja().contour(x, y, z, levels, params);
	}
	public void contour(int subaxis, float[] x, float[] y, float[][] z) {
		gcma().contour(subaxis, x, y, z);
	}
	public void contour(int subaxis, float[] x, float[] y, float[][] z, int levels, Object... params) {
		gcma().contour(subaxis, x, y, z, levels, params);
	}
	public void contour(int subaxis, float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		gcma().contour(subaxis, x, y, z, levels, params);
	}
	public void contour(int subaxis, float[][] x, float[][] y, float[][] z) {
		gcma().contour(subaxis, x, y, z);
	}
	public void contour(int subaxis, float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		gcma().contour(subaxis, x, y, z, levels, params);
	}
	public void contour(int subaxis, float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		gcma().contour(subaxis, x, y, z, levels, params);
	}
	public void contour(double[] x, double[] y, double[][] z) {
		gcja().contour(x, y, z);
	}
	public void contour(double[] x, double[] y, double[][] z, int levels, Object... params) {
		gcja().contour(x, y, z, levels, params);
	}
	public void contour(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		gcja().contour(x, y, z, levels, params);
	}
	public void contour(double[][] x, double[][] y, double[][] z) {
		gcja().contour(x, y, z);
	}
	public void contour(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		gcja().contour(x, y, z, levels, params);
	}
	public void contour(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		gcja().contour(x, y, z, levels, params);
	}
	public void contour(int subaxis, double[] x, double[] y, double[][] z) {
		gcma().contour(subaxis, x, y, z);
	}
	public void contour(int subaxis, double[] x, double[] y, double[][] z, int levels, Object... params) {
		gcma().contour(subaxis, x, y, z, levels, params);
	}
	public void contour(int subaxis, double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		gcma().contour(subaxis, x, y, z, levels, params);
	}
	public void contour(int subaxis, double[][] x, double[][] y, double[][] z) {
		gcma().contour(subaxis, x, y, z);
	}
	public void contour(int subaxis, double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		gcma().contour(subaxis, x, y, z, levels, params);
	}
	public void contour(int subaxis, double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		gcma().contour(subaxis, x, y, z, levels, params);
	}
	
	public void contourf(float[] x, float[] y, float[][] z) {
		gcja().contourf(x, y, z);
	}
	public void contourf(float[] x, float[] y, float[][] z, int levels, Object... params) {
		gcja().contourf(x, y, z, levels, params);
	}
	public void contourf(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		gcja().contourf(x, y, z, levels, params);
	}
	public void contourf(float[][] x, float[][] y, float[][] z) {
		gcja().contourf(x, y, z);
	}
	public void contourf(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		gcja().contourf(x, y, z, levels, params);
	}
	public void contourf(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		gcja().contourf(x, y, z, levels, params);
	}
	public void contourf(int subaxis, float[] x, float[] y, float[][] z) {
		gcma().contourf(subaxis, x, y, z);
	}
	public void contourf(int subaxis, float[] x, float[] y, float[][] z, int levels, Object... params) {
		gcma().contourf(subaxis, x, y, z, levels, params);
	}
	public void contourf(int subaxis, float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		gcma().contourf(subaxis, x, y, z, levels, params);
	}
	public void contourf(int subaxis, float[][] x, float[][] y, float[][] z) {
		gcma().contourf(subaxis, x, y, z);
	}
	public void contourf(int subaxis, float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		gcma().contourf(subaxis, x, y, z, levels, params);
	}
	public void contourf(int subaxis, float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		gcma().contourf(subaxis, x, y, z, levels, params);
	}
	public void contourf(double[] x, double[] y, double[][] z) {
		gcja().contourf(x, y, z);
	}
	public void contourf(double[] x, double[] y, double[][] z, int levels, Object... params) {
		gcja().contourf(x, y, z, levels, params);
	}
	public void contourf(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		gcja().contourf(x, y, z, levels, params);
	}
	public void contourf(double[][] x, double[][] y, double[][] z) {
		gcja().contourf(x, y, z);
	}
	public void contourf(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		gcja().contourf(x, y, z, levels, params);
	}
	public void contourf(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		gcja().contourf(x, y, z, levels, params);
	}
	public void contourf(int subaxis, double[] x, double[] y, double[][] z) {
		gcma().contourf(subaxis, x, y, z);
	}
	public void contourf(int subaxis, double[] x, double[] y, double[][] z, int levels, Object... params) {
		gcma().contourf(subaxis, x, y, z, levels, params);
	}
	public void contourf(int subaxis, double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		gcma().contourf(subaxis, x, y, z, levels, params);
	}
	public void contourf(int subaxis, double[][] x, double[][] y, double[][] z) {
		gcma().contourf(subaxis, x, y, z);
	}
	public void contourf(int subaxis, double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		gcma().contourf(subaxis, x, y, z, levels, params);
	}
	public void contourf(int subaxis, double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		gcma().contourf(subaxis, x, y, z, levels, params);
	}
	
	public void contourp(float[] x, float[] y, float[][] z) {
		gcja().contourp(x, y, z);
	}
	public void contourp(float[] x, float[] y, float[][] z, int levels, Object... params) {
		gcja().contourp(x, y, z, levels, params);
	}
	public void contourp(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		gcja().contourp(x, y, z, levels, params);
	}
	public void contourp(float[][] x, float[][] y, float[][] z) {
		gcja().contourp(x, y, z);
	}
	public void contourp(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		gcja().contourp(x, y, z, levels, params);
	}
	public void contourp(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		gcja().contourp(x, y, z, levels, params);
	}
	public void contourp(int subaxis, float[] x, float[] y, float[][] z) {
		gcma().contourp(subaxis, x, y, z);
	}
	public void contourp(int subaxis, float[] x, float[] y, float[][] z, int levels, Object... params) {
		gcma().contourp(subaxis, x, y, z, levels, params);
	}
	public void contourp(int subaxis, float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		gcma().contourp(subaxis, x, y, z, levels, params);
	}
	public void contourp(int subaxis, float[][] x, float[][] y, float[][] z) {
		gcma().contourp(subaxis, x, y, z);
	}
	public void contourp(int subaxis, float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		gcma().contourp(subaxis, x, y, z, levels, params);
	}
	public void contourp(int subaxis, float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		gcma().contourp(subaxis, x, y, z, levels, params);
	}
	public void contourp(double[] x, double[] y, double[][] z) {
		gcja().contourp(x, y, z);
	}
	public void contourp(double[] x, double[] y, double[][] z, int levels, Object... params) {
		gcja().contourp(x, y, z, levels, params);
	}
	public void contourp(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		gcja().contourp(x, y, z, levels, params);
	}
	public void contourp(double[][] x, double[][] y, double[][] z) {
		gcja().contourp(x, y, z);
	}
	public void contourp(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		gcja().contourp(x, y, z, levels, params);
	}
	public void contourp(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		gcja().contourp(x, y, z, levels, params);
	}
	public void contourp(int subaxis, double[] x, double[] y, double[][] z) {
		gcma().contourp(subaxis, x, y, z);
	}
	public void contourp(int subaxis, double[] x, double[] y, double[][] z, int levels, Object... params) {
		gcma().contourp(subaxis, x, y, z, levels, params);
	}
	public void contourp(int subaxis, double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		gcma().contourp(subaxis, x, y, z, levels, params);
	}
	public void contourp(int subaxis, double[][] x, double[][] y, double[][] z) {
		gcma().contourp(subaxis, x, y, z);
	}
	public void contourp(int subaxis, double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		gcma().contourp(subaxis, x, y, z, levels, params);
	}
	public void contourp(int subaxis, double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		gcma().contourp(subaxis, x, y, z, levels, params);
	}
	
	public void hatch(float[] x, float[] y, float[][] z, float lower, float upper, String pattern) {
		gcja().hatch(x,y,z, lower,upper, pattern, (Object)null);
	}
	public void hatch(float[] x, float[] y, float[][] z, float lower, float upper, String pattern, Object... params) {
		gcja().hatch(x, y, z, lower, upper, pattern, params);
	}
	public void hatch(float[][] x, float[][] y, float[][] z, float lower, float upper, String pattern) {
		gcja().hatch(x,y,z, lower,upper, pattern, (Object)null);
	}
	public void hatch(float[][] x, float[][] y, float[][] z, float lower, float upper, String pattern, Object... params) {
		gcja().hatch(x, y, z, lower, upper, pattern, params);
	}
	public void hatch(int subaxis, float[] x, float[] y, float[][] z, float lower, float upper, String pattern) {
		gcma().hatch(subaxis, x,y,z, lower,upper, pattern, (Object)null);
	}
	public void hatch(int subaxis, float[] x, float[] y, float[][] z, float lower, float upper, String pattern, Object... params) {
		gcma().hatch(subaxis, x, y, z, lower, upper, pattern, params);
	}
	public void hatch(int subaxis, float[][] x, float[][] y, float[][] z, float lower, float upper, String pattern) {
		gcma().hatch(subaxis, x,y,z, lower,upper, pattern, (Object)null);
	}
	public void hatch(int subaxis, float[][] x, float[][] y, float[][] z, float lower, float upper, String pattern, Object... params) {
		gcma().hatch(subaxis, x, y, z, lower, upper, pattern, params);
	}
	public void hatch(double[] x, double[] y, double[][] z, double lower, double upper, String pattern) {
		gcja().hatch(x, y, z, lower, upper, pattern, (Object)null);
	}
	public void hatch(double[] x, double[] y, double[][] z, double lower, double upper, String pattern, Object... params) {
		gcja().hatch(x, y, z, lower, upper, pattern, params);
	}
	public void hatch(double[][] x, double[][] y, double[][] z, double lower, double upper, String pattern) {
		gcja().hatch(x, y, z, lower, upper, pattern, (Object)null);
	}
	public void hatch(double[][] x, double[][] y, double[][] z, double lower, double upper, String pattern, Object... params) {
		gcja().hatch(x, y, z, lower, upper, pattern, params);
	}
	public void hatch(int subaxis, double[] x, double[] y, double[][] z, double lower, double upper, String pattern) {
		gcma().hatch(subaxis, x, y, z, lower, upper, pattern, (Object)null);
	}
	public void hatch(int subaxis, double[] x, double[] y, double[][] z, double lower, double upper, String pattern, Object... params) {
		gcma().hatch(subaxis, x, y, z, lower, upper, pattern, params);
	}
	public void hatch(int subaxis, double[][] x, double[][] y, double[][] z, double lower, double upper, String pattern) {
		gcma().hatch(subaxis, x, y, z, lower, upper, pattern, (Object)null);
	}
	public void hatch(int subaxis, double[][] x, double[][] y, double[][] z, double lower, double upper, String pattern, Object... params) {
		gcma().hatch(subaxis, x, y, z, lower, upper, pattern, params);
	}
	
	public void pcolour(float[] x, float[] y, float[][] z, float cmin, float cmax, Object... params) {
		gcja().pcolour(x, y, z, cmin, cmax, params);
	}
	public void pcolour(float[][] x, float[][] y, float[][] z, float cmin, float cmax, Object... params) {
		gcja().pcolour(x, y, z, cmin, cmax, params);
	}
	public void pcolour(double[] x, double[] y, double[][] z, double cmin, double cmax, Object... params) {
		gcja().pcolour(x, y, z, cmin, cmax, params);
	}
	public void pcolour(double[][] x, double[][] y, double[][] z, double cmin, double cmax, Object... params) {
		gcja().pcolour(x, y, z, cmin, cmax, params);
	}
	public void pcolour(int subaxis, float[] x, float[] y, float[][] z, float cmin, float cmax, Object... params) {
		gcma().pcolour(subaxis, x, y, z, cmin, cmax, params);
	}
	public void pcolour(int subaxis, float[][] x, float[][] y, float[][] z, float cmin, float cmax, Object... params) {
		gcma().pcolour(subaxis, x, y, z, cmin, cmax, params);
	}
	public void pcolour(int subaxis, double[] x, double[] y, double[][] z, double cmin, double cmax, Object... params) {
		gcma().pcolour(subaxis, x, y, z, cmin, cmax, params);
	}
	public void pcolour(int subaxis, double[][] x, double[][] y, double[][] z, double cmin, double cmax, Object... params) {
		gcma().pcolour(subaxis, x, y, z, cmin, cmax, params);
	}
	
	public void annotate(double x, double y, String text) {
		gcja().annotate(x, y, text, new Object[0]);
	}
	public void annotate(double x, double y, String text, Object... params) {
		gcja().annotate(x, y, text, params);
	}
	
	public void coastLines() {
		gcja().coastLines();
	}
	public void coastLines(int resolution) {
		gcja().coastLines(resolution);
	}
	
	public void land() {
		gcja().land(0xff676767, 0xff000000);
	}
	public void land(int land_colour, int coast_colour) {
		gcja().land(land_colour, coast_colour);
	}
	
	public void showShapefile(String path_to_shapefile, String shapeType) {
		gcja().showShapefile(path_to_shapefile, shapeType);
	}
	public void showShapefile(String path_to_shapefile, String shapeType, CoordinateReferenceSystem user_crs,
			Object... params) {
		gcja().showShapefile(path_to_shapefile, shapeType, user_crs, params);
	}
	public void showShapefile(String path_to_shapefile, String shapeType, int user_epsg_code, Object... params) {
		gcja().showShapefile(path_to_shapefile, shapeType, user_epsg_code, params);
	}
	
	public JColourbar colourbar() {
		return colourbar(gca(), "", "meither");
	}
	public JColourbar colourbar(String name) {
		return colourbar(gca(), name, "neither");
	}
	public JColourbar colourbar(String name, String extent) {
		return colourbar(gca(), name, extent);
	}
	public JColourbar colourbar(JAxis axis) {
		return colourbar(axis, "", "neither");
	}
	public JColourbar colourbar(JAxis axis, String name) {
		return colourbar(axis, name, "neither");
	}
	public JColourbar colourbar(JAxis axis, String name, String extent) {
		JColourbar cb = new JColourbar(axis, name);
		cb.setExtent(extent);
		addSubplot(cb);
		return cb;
	}
	
	/**
	 * adds a legend to the current JAxis object
	 * 
	 * @example SimpleLegend
	 */
	public void legend() {
		gcja().legend(PConstants.RIGHT, PConstants.TOP, 1d);
	}
	public void legend(double rts) {
		gcja().legend(PConstants.RIGHT, PConstants.TOP, rts);
	}
	public void legend(int left_right, int top_bottom) {
		gcja().legend(left_right, top_bottom, 1d);
	}
	public void legend(int left_right, int top_bottom, double rts) {
		gcja().legend(left_right, top_bottom, rts);
	}
	public void legend(int subaxis) {
		gcma().legend(subaxis, PConstants.RIGHT, PConstants.TOP, 1d);
	}
	public void legend(int subaxis, double rts) {
		gcma().legend(subaxis, PConstants.RIGHT, PConstants.TOP, rts);
	}
	public void legend(int subaxis, int left_right, int top_bottom) {
		gcma().legend(subaxis, left_right, top_bottom, 1d);
	}
	public void legend(int subaxis, int left_right, int top_bottom, double rts) {
		gcma().legend(subaxis, left_right, top_bottom, rts);
	}
	
	public void addText(double x, double y, String text) {
		gcja().addText(x, y, text, 1.0d, 0xff000000, PConstants.LEFT, PConstants.BOTTOM, 0d, null);
	}
	public void addText(double x, double y, String text, double textsize, int colour, String style) {
		gcja().addText(x, y, text, textsize, colour, PConstants.LEFT, PConstants.BOTTOM, 0d, style);
	}
	public void addText(double x, double y, String text, double textsize, int colour, int alignx, int aligny, String style) {
		gcja().addText(x, y, text, textsize, colour, alignx, aligny, 0d, style);
	}
	public void addText(double x, double y, String text, double textsize, int colour, int alignx, int aligny, double rotation, String style) {
		gcja().addText(x, y, text, textsize, colour, alignx, aligny, rotation, style);
	}
	public void addPolygon(JDPolygon poly, int inn_colour, int out_colour, double linewidth) {
		gcja().addPolygon(poly, inn_colour, out_colour, linewidth);
	}

	public void predefImgShow(String predefined_images) {
		gcja().predefImgShow(predefined_images);
	}

	public void imgShow(PImage img) {
		gcja().imgShow(img);
	}

	public void setXRange(float xmin, float xmax) {
		gcja().setXRange(xmin, xmax);
	}
	public void setXRange(double xmin, double xmax) {
		gcja().setXRange(xmin, xmax);
	}
	public void setYRange(float ymin, float ymax) {
		gcja().setYRange(ymin, ymax);
	}
	public void setYRange(double ymin, double ymax) {
		gcja().setYRange(ymin, ymax);
	}
	public void setRange(float xmin, float xmax, float ymin, float ymax) {
		gcja().setRange(xmin, xmax, ymin, ymax);
	}
	public void setRange(double xmin, double xmax, double ymin, double ymax) {
		gcja().setRange(xmin, xmax, ymin, ymax);
	}
	
	public void setFont(PFont font) {
		gcja().setFont(font);
	}
	
	public void setXTitle(String xtitle) {
		gcja().setXTitle(xtitle,"");
	}
	public void setXTitle(String xtitle, String xunit) {
		gcja().setXTitle(xtitle, xunit);
	}
	public void setYTitle(String ytitle) {
		gcja().setYTitle(ytitle,"");
	}
	public void setYTitle(String ytitle, String yunit) {
		gcja().setYTitle(ytitle, yunit);
	}
	public void setXTitle(int x_section, String xtitle) {
		gcma().setXTitle(x_section, xtitle,"");
	}
	public void setXTitle(int x_section, String xtitle, String xunit) {
		gcma().setXTitle(x_section, xtitle, xunit);
	}
	public void setYTitle(int y_section, String ytitle) {
		gcma().setYTitle(y_section, ytitle,"");
	}
	public void setYTitle(int y_section, String ytitle, String yunit) {
		gcma().setYTitle(y_section, ytitle, yunit);
	}
	public void setTitle(String _title) {
		gca().setTitle(_title);
	}
	
	public void setLogarithmicAxis(char axis) {
		gcja().setLogarithmicAxis(axis);
	}
	public void setAsTimeAxis(char axis, String unit) {
		gcja().setAsTimeAxis(axis, unit);
	}
	public void setAsTimeAxis(char axis, String unit, String calendar) {
		gcja().setAsTimeAxis(axis, unit, calendar);
	}
	public void setAsTimeAxis(char axis, String unit, String calendar, String format) {
		gcja().setAsTimeAxis(axis, unit, calendar, format);
	}
	
	public void showXAxisSide(String side) { gca().showXAxisSide(side); }
	public void showXAxisSide(int side)    { gca().showXAxisSide(side); }
	public void showYAxisSide(String side) { gca().showYAxisSide(side); }
	public void showYAxisSide(int side)    { gca().showYAxisSide(side); }
}
