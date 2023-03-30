package jplots.axes;

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import org.geotools.referencing.CRS;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import jplots.JPlot;
import jplots.colour.JColourtable;
import jplots.helper.FileLoader;
import jplots.layer.JContourLayer;
import jplots.layer.JContourLayer2D;
import jplots.layer.JHatchLayer;
import jplots.layer.JHatchLayer2D;
import jplots.layer.JImageLayer;
import jplots.layer.JLegend;
import jplots.layer.JLineLayer;
import jplots.layer.JPColourLayer;
import jplots.layer.JPlotsLayer;
import jplots.layer.JPolygonLayer;
import jplots.layer.JScatterLayer;
import jplots.layer.JShapesLayer;
import jplots.layer.JTextLayer;
import jplots.layer.JXYLayer;
import jplots.maths.JDPolygon;
import jplots.maths.JPlotMath;
import jplots.maths.JPlotMath.DateTime;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLatexShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JTextShape;
import jplots.transform.IdentityJProjection;
import jplots.transform.JProjection;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PImage;

public class JSingleAxis extends JAxis {

	private static CoordinateReferenceSystem coast_crs = null;
	private static Map<String, PImage> loadedPreDefImgs = new HashMap<>();
	static {
		try {
			coast_crs = CRS.parseWKT("GEOGCS[\"WGS 84\"," + "  DATUM[\"World Geodetic System 1984\","
					+ "    SPHEROID[\"WGS 84\", 6378137.0, 298.257223563, AUTHORITY[\"EPSG\",\"7030\"]],"
					+ "    AUTHORITY[\"EPSG\",\"6326\"]],"
					+ "  PRIMEM[\"Greenwich\", 0.0, AUTHORITY[\"EPSG\",\"8901\"]],"
					+ "  UNIT[\"degree\", 0.017453292519943295]," + "  AXIS[\"Geodetic longitude\", EAST],"
					+ "  AXIS[\"Geodetic latitude\", NORTH]," + "  AUTHORITY[\"EPSG\",\"4326\"]]");
		} catch (FactoryException fe) {
			fe.printStackTrace();
		}
	}
	
	protected boolean xRangeFix, yRangeFix;
	protected boolean isGeoAxis, xTim, yTim, xLog, yLog, xAxInv, yAxInv;
	protected double minX, maxX, minY, maxY;
	protected String titleX, titleY, unitX, unitY;
	protected String xTimUnit, yTimUnit, xTimCal, yTimCal, xTimFormat, yTimFormat;
	protected List<JPlotsLayer> layers;
	protected JSingleAxis[] shareXaxis, shareYaxis;
	
	public JSingleAxis(JPlot plot, int pos_x, int pos_y, int width, int height) {
		super(plot, pos_x, pos_y, width, height);
		layers = new ArrayList<>();
		shareXaxis = new JSingleAxis[] {};
		shareYaxis = new JSingleAxis[] {};
		if (plot.isDebug())
			System.out.println("[DEBUG] created PAxis-object: x/y=" + px + "/" + py + " w/h=" + pw + "/" + ph);
		defaults();
	}
	
	public JSingleAxis(JSingleAxis src_axis) {
		super(src_axis);
		layers = new ArrayList<>();
		shareXaxis = new JSingleAxis[] {};
		shareYaxis = new JSingleAxis[] {};
		defaults();
	}
	
	@Override
	protected void defaults() {
		super.defaults();
		xRangeFix = false;
		yRangeFix = false;
		xAxInv = false;
		yAxInv = false;
		minX = -1d;
		maxX = 1d;
		minY = -1d;
		maxY = 1d;
		xLog = false;
		yLog = false;
		xTim = false;
		yTim = false;
		xTimUnit = null;
		yTimUnit = null;
		xTimCal = null;
		yTimCal = null;
		xTimFormat = "dd.mm.yyyy";
		yTimFormat = "dd.mm.yyyy";
		titleX = ""; unitX = "";
		titleY = ""; unitY = "";
		layers.clear();
	}
	
	// ************************************
	// **** PUBLIC ************************
	// ************************************
	
	// ....
	@Override
	public void printInfo() {
		System.out.print(
				 this.getClass().getSimpleName()+":\n"
				+"    pos/size: "+px+"|"+py+" / "+pw+"|"+ph+"\n"
				+"    font/txtsize: "+pfont+"/"+txtsize+"\n"
				+"    x-Axis: log="+isXlogAxis()+" draw="+xAxOn+" ticks="+xTkOn+" grid="+xGrdOn+"\n"
				+"            label="+(titleX!=null?"\""+titleX+"\"":null)+" unit="+(unitX!=null?"\""+unitX+"\"":null)+"\n"
				+"    y-Axis: log="+isYlogAxis()+" draw="+yAxOn+" ticks="+yTkOn+" grid="+yGrdOn+"\n"
				+"            label="+(titleY!=null?"\""+titleY+"\"":null)+" unit="+(unitY!=null?"\""+unitY+"\"":null)+"\n"
				+"    proj:   geo="+isGeoAxis()+" which="+projection+"\n"
				+"    title:  "+(titleP!=null?"\""+titleP+"\"":null)+"\n"
			);
		}
	
	public void contour(float[] x, float[] y, float[][] z) {
		this.contour(x, y, z, 10, (Object[]) null);
	}
	public void contour(float[] x, float[] y, float[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, null);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contour(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), 0, levels);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contour(float[][] x, float[][] y, float[][] z) {
		this.contour(x, y, z, 10, (Object[]) null);
	}
	public void contour(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, null);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contour(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), 0, levels);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contour(double[] x, double[] y, double[][] z) {
		this.contour(x, y, z, 10, (Object[]) null);
	}
	public void contour(double[] x, double[] y, double[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, null);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contour(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), 0, levels);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contour(double[][] x, double[][] y, double[][] z) {
		this.contour(x, y, z, 10, (Object[]) null);
	}
	public void contour(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, null);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contour(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), 0, levels);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	
	public void contourf(float[] x, float[] y, float[][] z) {
		this.contourf(x, y, z, 10, (Object[]) null);
	}
	public void contourf(float[] x, float[] y, float[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels,
				JColourtable.pctables.get("default"), 2.0f, false, true, false);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourf(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0f, false, true,
				false);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourf(float[][] x, float[][] y, float[][] z) {
		this.contourf(x, y, z, 10, (Object[]) null);
	}
	public void contourf(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels,
				JColourtable.pctables.get("default"), 2.0f, false, true, false);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourf(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0f, false, true,
				false);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourf(double[] x, double[] y, double[][] z) {
		this.contourf(x, y, z, 10, (Object[]) null);
	}
	public void contourf(double[] x, double[] y, double[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels,
				JColourtable.pctables.get("default"), 2.0d, false, true, false);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourf(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0d, false, true,
				false);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourf(double[][] x, double[][] y, double[][] z) {
		this.contourf(x, y, z, 10, (Object[]) null);
	}
	public void contourf(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels,
				JColourtable.pctables.get("default"), 2.0d, false, true, false);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourf(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0d, false, true,
				false);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	
	public void contourp(float[] x, float[] y, float[][] z) {
		this.contourf(x, y, z, 10, (Object[]) null);
	}
	public void contourp(float[] x, float[] y, float[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, null);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourp(float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z),JPlotMath.fmax(z),0, levels);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourp(float[][] x, float[][] y, float[][] z) {
		this.contourf(x, y, z, 10, (Object[]) null);
	}
	public void contourp(float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, null);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourp(float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z),JPlotMath.fmax(z),0, levels);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourp(double[] x, double[] y, double[][] z) {
		this.contourp(x, y, z, 10, (Object[]) null);
	}
	public void contourp(double[] x, double[] y, double[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, null);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourp(double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z),JPlotMath.dmax(z), 0, levels);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourp(double[][] x, double[][] y, double[][] z) {
		this.contourp(x, y, z, 10, (Object[]) null);
	}
	public void contourp(double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, null);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	public void contourp(double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z),JPlotMath.dmax(z),0, levels);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers.add(cnl);
		readParams(cnl, params);
		updateRange(cnl);
	}
	
	public void pcolour(float[] x, float[] y, float[][] z, float cmin, float cmax, Object... params) {
		JPColourLayer pcl = new JPColourLayer(x, y, z, cmin, cmax);
		pcl.setColourtable(JColourtable.pctables.get("default"));
		layers.add(pcl);
		readParams(pcl, params);
		updateRange(pcl);
	}
	public void pcolour(float[][] x, float[][] y, float[][] z, float cmin, float cmax, Object... params) {
		JPColourLayer pcl = new JPColourLayer(x, y, z, cmin, cmax);
		pcl.setColourtable(JColourtable.pctables.get("default"));
		layers.add(pcl);
		readParams(pcl, params);
		updateRange(pcl);
	}
	public void pcolour(double[] x, double[] y, double[][] z, double cmin, double cmax, Object... params) {
		JPColourLayer pcl = new JPColourLayer(x, y, z, cmin, cmax);
		pcl.setColourtable(JColourtable.pctables.get("default"));
		layers.add(pcl);
		readParams(pcl, params);
		updateRange(pcl);
	}
	public void pcolour(double[][] x, double[][] y, double[][] z, double cmin, double cmax, Object... params) {
		JPColourLayer pcl = new JPColourLayer(x, y, z, cmin, cmax);
		pcl.setColourtable(JColourtable.pctables.get("default"));
		layers.add(pcl);
		readParams(pcl, params);
		updateRange(pcl);
	}
	
	public void plot(float[] x, float[] y) {
		this.plot(x, y, 0xff000000, 3f, "-", (Object) null);
	}
	public void plot(float[] x, float[] y, int colour, float linewidth, String linestyle, Object... params) {
		JPlotsLayer xyl = new JXYLayer(x, y, colour, linewidth, linestyle);
		layers.add(xyl);
		readParams(xyl, params);
		updateRange(xyl);
	}
	public void plot(double[] x, double[] y) {
		this.plot(x, y, 0xff000000, 3d, "-", (Object) null);
	}
	public void plot(double[] x, double[] y, int colour, double linewidth, String linestyle, Object... params) {
		JPlotsLayer xyl = new JXYLayer(x, y, colour, linewidth, linestyle);
		layers.add(xyl);
		readParams(xyl, params);
		updateRange(xyl);
	}
	
	public void scatter(float[] x, float[] y) {
		this.scatter(x, y, 0xff000000, 1f, "c", (Object) null);
	}
	public void scatter(float[] x, float[] y, int colour, float iconsize, String symbol, Object... params) {
		JPlotsLayer scl = new JScatterLayer(x, y, colour, iconsize, symbol);
		layers.add(scl);
		readParams(scl, params);
		updateRange(scl);
	}
	public void scatter(double[] x, double[] y) {
		this.scatter(x, y, 0xff000000, 1d, "c", (Object) null);
	}
	public void scatter(double[] x, double[] y, int colour, double iconsize, String symbol, Object... params) {
		JPlotsLayer scl = new JScatterLayer(x, y, colour, iconsize, symbol);
		layers.add(scl);
		readParams(scl, params);
		updateRange(scl);
	}
	
	public void hatch(float[] x, float[] y, float[][] z, float lower, float upper, String pattern) {
		this.hatch(x,y,z, lower,upper, pattern, (Object)null);
	}
	public void hatch(float[] x, float[] y, float[][] z, float lower, float upper, String pattern, Object... params) {
		JPlotsLayer hl = new JHatchLayer(x, y, null, null, z, 1f, lower, upper, pattern);
		layers.add(hl);
		readParams(hl, params);
		updateRange(hl);
	}
	public void hatch(float[][] x, float[][] y, float[][] z, float lower, float upper, String pattern) {
		this.hatch(x,y,z, lower,upper, pattern, (Object)null);
	}
	public void hatch(float[][] x, float[][] y, float[][] z, float lower, float upper, String pattern, Object... params) {
		JPlotsLayer hl = new JHatchLayer(null, null, x, y, z, 1f, lower, upper, pattern);
		layers.add(hl);
		readParams(hl, params);
		updateRange(hl);
	}
	public void hatch(double[] x, double[] y, double[][] z, double lower, double upper, String pattern) {
		this.hatch(x, y, z, lower, upper, pattern, (Object)null);
	}
	public void hatch(double[] x, double[] y, double[][] z, double lower, double upper, String pattern, Object... params) {
		JPlotsLayer hl = new JHatchLayer(x, y, null, null, z, 1d, lower, upper, pattern);
		layers.add(hl);
		readParams(hl, params);
		updateRange(hl);
	}
	public void hatch(double[][] x, double[][] y, double[][] z, double lower, double upper, String pattern) {
		this.hatch(x, y, z, lower, upper, pattern, (Object)null);
	}
	public void hatch(double[][] x, double[][] y, double[][] z, double lower, double upper, String pattern, Object... params) {
		JPlotsLayer hl = new JHatchLayer(null, null, x, y, z, 1d, lower, upper, pattern);
		layers.add(hl);
		readParams(hl, params);
		updateRange(hl);
	}
	
	public void axhline(double y) {
		axhline(y, 0xff000000, 3f, "-");
	}
	public void axhline(double y, int colour, double linewidth, String linestyle) {
		JPlotsLayer xyl = new JLineLayer(y, 'h', colour, linewidth, linestyle);
		layers.add(xyl);
		updateRange(xyl, "y");
	}
	public void axvline(double x) {
		axvline(x, 0xff000000, 3f, "-");
	}
	public void axvline(double x, int colour, double linewidth, String linestyle) {
		JPlotsLayer xyl = new JLineLayer(x, 'v', colour, linewidth, linestyle);
		layers.add(xyl);
		updateRange(xyl, "x");
	}
	
	public void addText(double x, double y, String text) {
		JTextLayer tl = new JTextLayer(false, text, x, y, 1.0d, 0xff000000, PConstants.LEFT, PConstants.BOTTOM, 0d, null);
		layers.add(tl);
	}
	public void addText(double x, double y, String text, double textsize, int colour, String style) {
		JTextLayer tl = new JTextLayer(false, text, x, y, textsize, colour, PConstants.LEFT, PConstants.BOTTOM, 0d, style);
		layers.add(tl);
	}
	public void addText(double x, double y, String text, double textsize, int colour, int alignx, int aligny, String style) {
		JTextLayer tl = new JTextLayer(false, text, x, y, textsize, colour, alignx, aligny, 0d, style);
		layers.add(tl);
	}
	public void addText(double x, double y, String text, double textsize, int colour, int alignx, int aligny, double rotation, String style) {
		JTextLayer tl = new JTextLayer(false, text, x, y, textsize, colour, alignx, aligny, rotation, style);
		layers.add(tl);
	}
	public void addPolygon(JDPolygon poly, int inn_colour, int out_colour, double linewidth) {
		layers.add(new JPolygonLayer(poly, inn_colour, out_colour, linewidth));
	}
	
	public void annotate(double x, double y, String text) {
		annotate(x, y, text, new Object[0]);
	}
	public void annotate(double x, double y, String text, Object... params) {
		JTextLayer tl = new JTextLayer(true, text, x, y, 1d, 0xff000000, PConstants.LEFT, PConstants.BOTTOM, 0d, null);
		layers.add(tl); readParams(tl, params);
	}
	
	public void legend() {
		JPlotsLayer lgl = new JLegend(this, PConstants.RIGHT, PConstants.TOP, false, 1d);
		layers.add(lgl);
	}
	public void legend(double rts) {
		JPlotsLayer lgl = new JLegend(this, PConstants.RIGHT, PConstants.TOP, false, rts);
		layers.add(lgl);
	}
	public void legend(int left_right, int top_bottom) {
		JPlotsLayer lgl = new JLegend(this, left_right, top_bottom, false, 1d);
		layers.add(lgl);
	}
	public void legend(int left_right, int top_bottom, double rts) {
		JPlotsLayer lgl = new JLegend(this, left_right, top_bottom, false, rts);
		layers.add(lgl);
	}
	
	public void coastLines() {
		coastLines(110);
	}
	public void coastLines(int resolution) {
		JPlotsLayer shl = new JShapesLayer(
				FileLoader.loadResourceShapeFile("/data/ne_" + resolution + "m_coastline", coast_crs), "line");
		layers.add(shl);
		shl.setLineColour(0xff000000);
	}
	
	public void land() {
		land(0xff676767, 0xff000000);
	}
	public void land(int land_colour, int coast_colour) {
		// JPlotsLayer shl = new
		// JShapesLayer(FileLoader.loadResourceShapeFile("/data/simplified_land_polygons"),
		// "polygon", 3857);
		JPlotsLayer shl = new JShapesLayer(FileLoader.loadResourceShapeFile("/data/ne_110m_land", 4326), "polygon");
		layers.add(shl);
		shl.setFillColour(land_colour);
		shl.setLineColour(coast_colour);
	}
	
	public void showShapefile(String path_to_shapefile, String shapeType) {
		showShapefile(path_to_shapefile, shapeType, null);
	}
	public void showShapefile(String path_to_shapefile, String shapeType, CoordinateReferenceSystem user_crs,
			Object... params) {
//		int i = path_to_shapefile.lastIndexOf(".");
//		if(i<0) i = path_to_shapefile.length();
//		String path = path_to_shapefile.substring(0, i);
		JPlotsLayer shl = new JShapesLayer(FileLoader.loadResourceShapeFile(path_to_shapefile, user_crs), shapeType);
		layers.add(shl);
		readParams(shl, params);
	}
	public void showShapefile(String path_to_shapefile, String shapeType, int user_epsg_code, Object... params) {
//		int i = path_to_shapefile.lastIndexOf(".");
//		if(i<0) i = path_to_shapefile.length();
//		String path = path_to_shapefile.substring(0, i);
		JPlotsLayer shl = new JShapesLayer(FileLoader.loadResourceShapeFile(path_to_shapefile, user_epsg_code),
				shapeType);
		layers.add(shl);
		readParams(shl, params);
	}
	
	/**
	 * predefined images are used as background images in plot, especially with
	 * geographical projections
	 * <p>
	 * earth&nbsp; -- NASA image of earths surface<br>
	 * earth2 -- ETOPO-like image of earths surface<br>
	 * etopo1 -- ETOPO image from <a href=
	 * "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/image/">NOAA</a><br>
	 * mercury,venus,mars,jupiter,saturn,uranus,neptune -- NASA images of surfaces
	 * of planets in our solar system<br>
	 * hsbpalette -- rainbow palette of all possible hsb-colors
	 * </p>
	 *
	 * @param predefined_image one name of above list of predefined images
	 */
	public void predefImgShow(String predefined_image) {
		JPlotsLayer iml;
		if (loadedPreDefImgs.containsKey(predefined_image)) {
			iml = new JImageLayer(loadedPreDefImgs.get(predefined_image));
		} else {
			iml = new JImageLayer(loadPreDefImg(pplot.getApplet(), predefined_image));
		}
		layers.add(0, iml);
	}
	
	public void imgShow(PImage img) {
		JPlotsLayer iml = new JImageLayer(img);
		layers.add(iml);
		// updateRange(iml);
	}
	
	// ....
	public void setXTitle(String xtitle) {
		titleX = "";
		if(xtitle!=null) titleX = xtitle;
		unitX  = "";
	}
	public void setXTitle(String xtitle, String xunit) {
		titleX = "";
		if(xtitle!=null) titleX = xtitle;
		unitX  = "";
		if(xunit!=null) unitX = xunit;
	}
	public void setYTitle(String ytitle) {
		titleY = "";
		if(ytitle!=null) titleY = ytitle;
		unitY  = "";
	}
	public void setYTitle(String ytitle, String yunit) {
		titleY = "";
		if(ytitle!=null) titleY = ytitle;
		unitY  = "";
		if(yunit!=null) unitY = yunit;
	}
	
	public JSingleAxis setXRange(double xmin, double xmax) {
		setXRange(xmin, xmax, true);
		return this;
	}
	private void setXRange(double xmin, double xmax, boolean notify) {
		minX = Math.min(xmin,xmax);
		maxX = Math.max(xmin,xmax);
		xRangeFix = true;
		if (notify)
			for (JSingleAxis a : shareXaxis)
				if (!a.equals(this))
					a.setXRange(xmin, xmax, false);
	}
	public JSingleAxis setYRange(double ymin, double ymax) {
		setYRange(ymin, ymax, true);
		return this;
	}
	private void setYRange(double ymin, double ymax, boolean notify) {
		minY = Math.min(ymin,ymax);
		maxY = Math.max(ymin,ymax);
		yRangeFix = true;
		if (notify)
			for (JSingleAxis a : shareYaxis)
				if (!a.equals(this))
					a.setYRange(ymin, ymax, false);
	}
	public JSingleAxis setRange(double xmin, double xmax, double ymin, double ymax) {
		setXRange(xmin, xmax);
		setYRange(ymin, ymax);
		return this;
	}
	
	@Override
	public void setGeoProjection(JProjection proj) {
		projection = proj;
		isGeoAxis = !(proj instanceof IdentityJProjection);
		double[] rr = projection.defaultMapExtend();
		setRange(rr[0], rr[1], rr[2], rr[3]);
		xAxOn = false;
		yAxOn = false;
		xGrdOn = false;
		yGrdOn = false;
		pplot.redraw(true);
	}
	
	public void setLogarithmicAxis(char axis) {
		switch (axis) {
		case 'b':
			setLogarithmicAxis('x');
			setLogarithmicAxis('y');
			break;
		case 'x':
			xLog = true;
			xTim = false;
			break;
		case 'y':
			yLog = true;
			yTim = false;
			break;
		default:
			System.err.println("Unknown parameter '" + axis + "' for axis in <setLogarithmicAxis(axis)>.");
			break;
		}
	}
	public void setAsTimeAxis(char axis, String unit) {
		setAsTimeAxis(axis, unit, "gregorian", "dd.mm.yyyy");
	}
	public void setAsTimeAxis(char axis, String unit, String calendar) {
		setAsTimeAxis(axis, unit, calendar, "dd.mm.yyyy");
	}
	public void setAsTimeAxis(char axis, String unit, String calendar, String format) {
		switch (axis) {
		case 'b':
			setAsTimeAxis('x', unit, calendar, format);
			setAsTimeAxis('y', unit, calendar, format);
			break;
		case 'x':
			xLog = false;
			xTim = true;
			xTimUnit = unit;
			xTimCal = calendar;
			if (format != null)
				xTimFormat = format;
			break;
		case 'y':
			yLog = false;
			yTim = true;
			yTimUnit = unit;
			yTimCal = calendar;
			if (format != null)
				yTimFormat = format;
			break;
		default:
			System.err.println("Unknown parameter '" + axis + "' for axis in <setLogarithmicAxis(axis)>.");
			break;
		}
	}
	
	public JSingleAxis addSharedAxis(char which, JSingleAxis new_axis) {
		addSharedAxis(which, new_axis, false);
		return this;
	}
	
	// ************************************
	// **** GETTER ************************
	// ************************************
	
	/**
	 * gives a new JAxis instance with same position and dimensions
	 * 
	 * @return new JAxis object
	 */
	@Override
	public JAxis copy() {
		return new JSingleAxis(this);
	}
	
	@Override
	public double[] getRange() {
		return new double[] { minX, maxX, minY, maxY };
	}
	public int[] getSize() {
		return new int[] { px, py, pw, ph };
	}
	
	public boolean isXlogAxis() {
		return xLog;
	}
	public boolean isYlogAxis() {
		return yLog;
	}
	public boolean isGeoAxis() {
		return isGeoAxis;
	}
	
	public JPlotsLayer getLayer(int layernum) {
		return layers.get(layernum);
	}

	@Override
	public List<JPlotsLayer> getLayers() {
		return layers;
	}
	
	public void addLayer(JPlotsLayer layer, Object... params) {
		layers.add(layer);
		readParams(layer, params);
		updateRange(layer);
	}
	
	// ************************************
	// **** PACKAGE PRIVATE ***************
	// ************************************
	
	@Override
	public JGroupShape createPlot(PApplet applet, int w, int h) {
		if (isGeoAxis) {
			double r = 0.5d * Math.max((maxX - minX) / pw, (maxY - minY) / ph);
			double xm = 0.5d * (minX + maxX);
			double ym = 0.5d * (minY + maxY);
			minX = xm - r * pw;
			maxX = xm + r * pw;
			minY = ym - r * ph;
			maxY = ym + r * ph;
			xLog = false;
			yLog = false;
			xTim = false;
			yTim = false;
		} else if (xLog || yLog) {
			if (xLog && (minX < 0d || maxX < 0d)) {
				xLog = false;
				System.err.println("found negative values in X-range [x=" + minX + " ... " + maxX + "]");
			}
			if (yLog && (minY < 0d || maxY < 0d)) {
				yLog = false;
				System.err.println("found negative values in Y-range [y=" + minY + " ... " + maxY + "]");
			}
		}
		if (pplot.isDebug())
			System.out.println("[DEBUG] JAxis-object: min/max={x:" + minX + "/" + maxX + ", y:" + minY + "/" + maxY
					+ "} with " + layers.size() + " layer" + (layers.size() > 1 ? "s" : ""));
		JGroupShape graph = new JGroupShape();
		for (int l = 0; l < layers.size(); l++) {
			JPlotsLayer layer = layers.get(l);
			if (layer instanceof JLegend)
				continue;
			layer.setRange(minX, maxX, minY, maxY);
			if (xAxInv || yAxInv)
				layer.invert(xAxInv ? (yAxInv ? "both" : "x") : "y", true);
			layer.createVectorImg(this, l, graph);
		}
		for (int l = 0; l < layers.size(); l++) {
			JPlotsLayer layer = layers.get(l);
			if (layer instanceof JLegend) {
				layer.setRange(minX, maxX, minY, maxY);
				if (xAxInv || yAxInv)
					layer.invert(xAxInv ? (yAxInv ? "both" : "x") : "y", true);
				layer.createVectorImg(this, l, graph);
			}
		}
		if (isGeoAxis) {
			projection.addGrid(this, graph);
			projection.drawBorder(this, graph);
		} else {
			if (xAxOn || xGrdOn)
				graph.addChild(createXAxis());
			if (yAxOn || yGrdOn)
				graph.addChild(createYAxis());
			if (xAxOn) {
				graph.addChild(new JLineShape(3f, 0xff000000, px, (float)py, px + pw, py));
				graph.addChild(new JLineShape(3f, 0xff000000, px, (float)py + ph, px + pw, py + ph));
			}
			if (yAxOn) {
				graph.addChild(new JLineShape(3f, 0xff000000, px, (float)py, px, py + ph));
				graph.addChild(new JLineShape(3f, 0xff000000, px + pw, (float)py, px + pw, py + ph));
			}
		}
		if (titleP.length() > 0) {
			if (pplot.isDebug())
				System.out.println("[DEBUG] JAxis: add title \"" + titleP + "\" to graphic.");
			if(JPlot.supportLatex)
				graph.addChild(new JLatexShape(titleP, px+0.5f*pw, py-0.04f*ph, (float)(1.3d*txtsize), PConstants.CENTER, PConstants.BOTTOM, 0xff000000, 0f, null));
			else
				graph.addChild(new JTextShape(titleP, px+0.5f*pw, py-0.04f*ph, (float)(1.3d*txtsize), PConstants.CENTER, PConstants.BOTTOM, 0xff000000, 0f, null));
		}
		graph.addChild(new JLineShape(0f, 0x00999999, 0f,0f, 1f,1f));
		return graph;
	}
	@Override
	public JPlotShape createPlotOnlyAxes(PApplet applet, int w, int h) {
		if (isGeoAxis) {
			double r = 0.5d * Math.max((maxX - minX) / pw, (maxY - minY) / ph);
			double xm = 0.5d * (minX + maxX);
			double ym = 0.5d * (minY + maxY);
			minX = xm - r * pw;
			maxX = xm + r * pw;
			minY = ym - r * ph;
			maxY = ym + r * ph;
			xLog = false;
			yLog = false;
			xTim = false;
			yTim = false;
		} else if (xLog || yLog) {
			if (xLog && (minX < 0d || maxX < 0d)) {
				xLog = false;
//				System.err.println("found negative values in X-range [x=" + minX + " ... " + maxX + "]");
			}
			if (yLog && (minY < 0d || maxY < 0d)) {
				yLog = false;
//				System.err.println("found negative values in Y-range [y=" + minY + " ... " + maxY + "]");
			}
		}
		JGroupShape graph = new JGroupShape();
		for(JPlotsLayer layer: layers) {
			if(layer instanceof JLegend) layer.createVectorImg(this, 0, graph);
			if(layer instanceof JTextLayer) layer.createVectorImg(this, 0, graph);
		}
		if (isGeoAxis) {
			projection.addGrid(this, graph);
			projection.drawBorder(this, graph);
		} else {
			if (xAxOn || xGrdOn)
				graph.addChild(createXAxis());
			if (yAxOn || yGrdOn)
				graph.addChild(createYAxis());
			if (xAxOn) {
				graph.addChild(new JLineShape(3f, 0xff000000, px, (float)py, px + pw, py));
				graph.addChild(new JLineShape(3f, 0xff000000, px, (float)py + ph, px + pw, py + ph));
			}
			if (yAxOn) {
				graph.addChild(new JLineShape(3f, 0xff000000, px, (float)py, px, py + ph));
				graph.addChild(new JLineShape(3f, 0xff000000, px + pw, (float)py, px + pw, py + ph));
			}
		}
		if (titleP.length() > 0) {
			if (pplot.isDebug())
				System.out.println("[DEBUG] JAxis: add title \"" + titleP + "\" to graphic.");
			if(JPlot.supportLatex) {
				graph.addChild(new JLatexShape(titleP, px + 0.5f * pw, py - 0.04f * ph, (float) (1.3d * txtsize), PConstants.CENTER,
						PConstants.BOTTOM, 0xff000000, 0f, null));
			} else {
				graph.addChild(new JTextShape(titleP, px + 0.5f * pw, py - 0.04f * ph, (float) (1.3d * txtsize), PConstants.CENTER,
						PConstants.BOTTOM, 0xff000000, 0f, null));
			}
		}
		graph.addChild(new JLineShape(0f, 0x00999999, 0f,0f, 1f,1f));
		return graph;
	}
	
	protected void addSharedAxis(char which, JSingleAxis new_axis, boolean notify) {
		if (new_axis.equals(this))
			return;
		switch (which) {
		case 'x':
			JSingleAxis[] tempx = new JSingleAxis[shareXaxis.length];
			for (int a = 0; a < tempx.length; a++)
				tempx[a] = shareXaxis[a];
			shareXaxis = new JSingleAxis[tempx.length + 1];
			for (int a = 0; a < tempx.length; a++)
				shareXaxis[a] = tempx[a];
			shareXaxis[tempx.length] = new_axis;
			if (notify)
				for (int a = 0; a < tempx.length; a++)
					shareXaxis[a].addSharedAxis('x', new_axis, false);
			break;
		case 'y':
			JSingleAxis[] tempy = new JSingleAxis[shareYaxis.length];
			for (int a = 0; a < tempy.length; a++)
				tempy[a] = shareYaxis[a];
			shareYaxis = new JSingleAxis[tempy.length + 1];
			for (int a = 0; a < tempy.length; a++)
				shareYaxis[a] = tempy[a];
			shareYaxis[tempy.length] = new_axis;
			if (notify)
				for (int a = 0; a < tempy.length; a++)
					shareYaxis[a].addSharedAxis('y', new_axis, false);
			break;
		default:
			break;
		}
	}

	// ************************************
	// **** PRIVATE ***********************
	// ************************************
	
	private void readParams(JPlotsLayer layer, Object... params) {
		if (params == null)
			return;
		int o = 0;
		while (o < params.length) {
			if (params[o] instanceof String) {
				String p = ((String) params[o]).toLowerCase();
				boolean isunread = true;
				if (isunread && ("tf".equals(p) || "transform".equals(p)) && o + 1 < params.length) {
					layer.setSourceProjection((JProjection) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("am".equals(p) || "anglemode".equals(p)) && o + 1 < params.length) {
					layer.angleMode((String) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("l".equals(p) || "lines".equals(p)) && o + 1 < params.length) {
					layer.lines((boolean) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("ls".equals(p) || "linestyle".equals(p)) && o + 1 < params.length) {
					if(params[o+1] instanceof String[] && (layer instanceof JContourLayer || layer instanceof JContourLayer2D)) {
						if(layer instanceof JContourLayer)   ((JContourLayer)   layer).setStyles((String[]) params[o+1]);
						if(layer instanceof JContourLayer2D) ((JContourLayer2D) layer).setStyles((String[]) params[o+1]);
					} else
						layer.setStyle((String) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("lw".equals(p) || "linewidth".equals(p)) && o + 1 < params.length) {
					layer.setLineWidth(params[o + 1] instanceof Float ? (float) params[o + 1] : (double) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("lc".equals(p) || "linecolor".equals(p) || "linecolour".equals(p))
						&& o + 1 < params.length) {
					if (params[o + 1] instanceof Integer) {
						layer.setLineColour((int) params[o + 1]);
						o++;
						isunread = false;
					} else if (params[o + 1] instanceof int[]) {
						layer.setLineColours((int[]) params[o + 1]);
						o++;
						isunread = false;
					}
				}
				if (isunread && ("fc".equals(p) || "fillcolor".equals(p) || "fillcolour".equals(p))
						&& o + 1 < params.length) {
					if (params[o + 1] instanceof Integer) {
						layer.setFillColour((int) params[o + 1]);
						o++;
						isunread = false;
					} else if (params[o + 1] instanceof int[]) {
						layer.setFillColours((int[]) params[o + 1]);
						o++;
						isunread = false;
					}
				}
				if (isunread && ("fs".equals(p) || "ts".equals(p) || "fontsize".equals(p) || "textsize".equals(p))
						&& o + 1 < params.length) {
					layer.setLineWidth((double) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("lb".equals(p) || "label".equals(p))) {
					if (params[o + 1] instanceof String) {
						layer.setLabel((String) params[o + 1]);
						o++;
						isunread = false;
					}
				}
				if (isunread && ("ct".equals(p) || "colortable".equals(p) || "colourtable".equals(p))
						&& o + 1 < params.length) {
					layer.setColourtable((JColourtable) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("z".equals(p)) && o + 1 < params.length) {
					layer.addParallelArray(params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("zrange".equals(p)) && o + 2 < params.length) {
					layer.setZRange((double)params[o+1], (double)params[o+2]);
					o+=2;
					isunread = false;
				}
				if (isunread && ("density".equals(p) || "pd".equals(p)) && o + 1 < params.length) {
					if(layer instanceof JHatchLayer2D)
						((JHatchLayer2D) layer).setDensity((double)params[o+1]);
					o+=1;
					isunread = false;
				}
				if (isunread && "invertxaxis".equals(p)) {
					layer.invert("x", true);
					xAxInv = true;
					isunread = false;
				}
				if (isunread && "invertyaxis".equals(p)) {
					layer.invert("y", true);
					yAxInv = true;
					isunread = false;
				}
			} else {
				System.err.println("[ERROR] Cannot interprete param " + o + ": " + params[o]);
			}
			o++;
		}
	}

	private void updateRange(JPlotsLayer layer) {
		if(isGeoAxis) return;
		double[] r = layer.getRange();
		double xmin = r[0], xmax = r[1], ymin = r[2], ymax = r[3];
		if (r[1] - r[0] < 1.e-20) {
			double xm = 0.5d * (r[0] + r[1]);
			double xr = Math.max(1.0e-10d, Math.abs(xm) * 1.0e-10d);
			xmin = xm - xr;
			xmax = xm + xr;
		}
		if (r[3] - r[2] < 1.e-20) {
			double ym = 0.5d * (r[2] + r[3]);
			double yr = Math.max(1.0e-10d, Math.abs(ym) * 1.e-10d);
			ymin = ym - yr;
			ymax = ym + yr;
		}
		if (layers.size() == 1) {
			if (!xRangeFix) {
				minX = xmin;
				maxX = xmax;
			}
			if (!yRangeFix) {
				minY = ymin;
				maxY = ymax;
			}
		} else {
			if (!xRangeFix) {
				if (xmin < minX)
					minX = xmin;
				if (xmax > maxX)
					maxX = xmax;
			}
			if (!yRangeFix) {
				if (ymin < minY)
					minY = ymin;
				if (ymax > maxY)
					maxY = ymax;
			}
		}
	}

	private void updateRange(JPlotsLayer layer, String axis) {
		if(isGeoAxis) return;
		double[] r = layer.getRange();
		double xmin = r[0], xmax = r[1], ymin = r[2], ymax = r[3];
		if (r[1] - r[0] < 1.e-20) {
			double xm = 0.5d * (r[0] + r[1]);
			double xr = Math.max(1.0e-10d, Math.abs(xm) * 1.0e-10d);
			xmin = xm - xr;
			xmax = xm + xr;
		}
		if (r[3] - r[2] < 1.e-20) {
			double ym = 0.5d * (r[2] + r[3]);
			double yr = Math.max(1.0e-10d, Math.abs(ym) * 1.e-10d);
			ymin = ym - yr;
			ymax = ym + yr;
		}
		if (layers.size() == 1) {
			if (!xRangeFix && axis.equals("x")) {
				minX = xmin;
				maxX = xmax;
			}
			if (!yRangeFix && axis.equals("y")) {
				minY = ymin;
				maxY = ymax;
			}
		} else {
			if (!xRangeFix && axis.equals("x")) {
				if (xmin < minX)
					minX = xmin;
				if (xmax > maxX)
					maxX = xmax;
			}
			if (!yRangeFix && axis.equals("y")) {
				if (ymin < minY)
					minY = ymin;
				if (ymax > maxY)
					maxY = ymax;
			}
		}
	}

	private JGroupShape createXAxis() {
		double Xin = xLog ? Math.log10(minX) : minX;
		double Xax = xLog ? Math.log10(maxX) : maxX;
		JGroupShape axisgrid = new JGroupShape();
		// first estimate of ticks
		double[] oticks = null;
		String[] otickmark = null;
		if (xLog) {
			oticks = JPlotMath.optimalLogarithmicTicks(minX, maxX);
			otickmark = new String[oticks.length];
			for (int t = 2; t < otickmark.length; t++) {
				otickmark[t] = oticks[t] + "";
				oticks[t] = Math.log10(oticks[t]);
			}
		} else if (xTim) {
			oticks = JPlotMath.optimalTimeTicks(minX, maxX, xTimUnit, xTimCal);
			otickmark = new String[oticks.length];
			for (int t = 2; t < otickmark.length; t++)
				otickmark[t] = DateTime.fromDouble(oticks[t], xTimUnit, xTimCal).format(xTimFormat, xTimCal);
		} else {
			oticks = JPlotMath.optimalLinearTicks(minX, maxX);
			double vf = 1d / (oticks[0]);
			int decimal = (int) (1000d * oticks[1] + 0.5d);
			decimal = decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
			otickmark = new String[oticks.length];
			for (int t = 2; t < otickmark.length; t++)
				otickmark[t] = PApplet.nf((float) (oticks[t] * vf), 0, decimal).replace(",", ".");
		}
		otickmark[0] = "";
		otickmark[1] = "";
		double tmlen = 0d;
		pplot.getGraphic().textSize(200);
		pplot.getGraphic().textAlign(PConstants.LEFT, PConstants.TOP);
		// create tickmark strings and calc mean tickmark text width
		for (int t = 2; t < oticks.length; t++) {
			tmlen += pplot.getGraphic().textWidth(otickmark[t]) / 200f;
		}
		tmlen *= this.txtsize / (oticks.length - 2);
		// with upper bound of number of ticks
		int tickcount = Math.max(2, (int) (pw / (1.2d * tmlen) + 0.99999999d));
		if (pplot.isDebug())
			System.out.println("[DEBUG] JAxis-object: tmlen=" + tmlen + " -> tickcount approx. " + tickcount);
		// create new ticks
		double[] ticks = null;
		String[] tickmark = null;
		String tickmarkFactor = "";
		if (xLog) {
			ticks = JPlotMath.optimalLogarithmicTicks(minX, maxX, tickcount);
			tickmark = new String[ticks.length];
			for (int t = 2; t < tickmark.length; t++) {
				tickmark[t] = ticks[t] + "";
				ticks[t] = Math.log10(ticks[t]);
			}
		} else if (xTim) {
			ticks = JPlotMath.optimalTimeTicks(minX, maxX, xTimUnit, xTimCal, tickcount);
			tickmark = new String[ticks.length];
			for (int t = 2; t < tickmark.length; t++)
				tickmark[t] = DateTime.fromDouble(ticks[t], xTimUnit, xTimCal).format(xTimFormat, xTimCal);
		} else {
			ticks = JPlotMath.optimalLinearTicks(minX, maxX, tickcount);
			double vf = 1d / (ticks[0]);
			int decimal = (int) (1000d * ticks[1] + 0.5d);
			decimal = decimal % 1000 == 0 ? 0 : decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
			tickmark = new String[ticks.length];
			for (int t = 2; t < tickmark.length; t++) {
				if(decimal==0) tickmark[t] = ""+(int)(ticks[t]*vf+0.0005d-(ticks[t]*vf<0d?1:0));
				else tickmark[t] = PApplet.nf((float) (ticks[t] * vf), 0, decimal).replace(",", ".");
			}
			double lvf = Math.log10(ticks[0]);
			if(Math.abs(lvf)>2.9d) {
				int ivf = (int) (lvf+0.5d) - (lvf<0d ? -1 : 0);
				tickmarkFactor = JPlot.supportLatex ? "10$^{"+ivf+"}$" : "10^"+ivf;
			}
		}
		tickmark[0] = "";
		tickmark[1] = "";
		double[] tcpos = JPlotMath.map(ticks, Xin, Xax, px, px + pw);
		if (xAxInv)
			for (int t = 0; t < tcpos.length; t++)
				tcpos[t] = 2 * px + pw - tcpos[t];
		if (pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for (int t = 2; t < ticks.length; t++) {
				tickStr += ", " + tickmark[t];
				posStr += ", " + PApplet.nf((float) tcpos[t], 0, 2).replace(",", ".");
			}
			System.out.println("[DEBUG] JAxis-object: Xtickfactors={p10: " + ticks[0] + ", f: " + ticks[1] + "}");
			System.out.println("[DEBUG] JAxis-object: Xtickval={" + tickStr.substring(2) + "}");
			System.out.println("[DEBUG] JAxis-object: Xtickpos={" + posStr.substring(2) + "}");
		}
		if (xGrdOn) {
			for (int t = 2; t < ticks.length; t++)
				if (ticks[t] >= Math.min(Xin, Xax) && ticks[t] <= Math.max(Xin, Xax))
					axisgrid.addChild(new JLineShape(2f, 0xff999999, (float) tcpos[t], py, (float) tcpos[t], py + ph));
		}
		if (xAxOn) {
			if (titleX.length() > 0 || unitX.length() > 0) {
				if (pplot.isDebug())
					System.out.println("[DEBUG] JAxi-object: add x-axis title \""+titleX+"\" and unit \""+unitX+"\" with text size "+txtsize);
				String txtemp = ""+titleX;
				if(unitX.length()>0 || tickmarkFactor.length()>0) txtemp += (titleX.length()>0?" ":"")+"[";
				if(tickmarkFactor.length()>0) txtemp += tickmarkFactor;
				if(unitX.length()>0) txtemp += (tickmarkFactor.length()>0?" ":"")+unitX;
				if(unitX.length()>0 || tickmarkFactor.length()>0) txtemp += "]";
				//TODO set txtsize back to 1.1d*txtsize
				if(JPlot.supportLatex) {
					axisgrid.addChild(new JLatexShape(txtemp, px + 0.5f * pw, py + 1.04f * ph + (float) txtsize,
							(float) (1.1d * txtsize), PConstants.CENTER, PConstants.TOP, 0xff000000, 0, null));
				} else {
					axisgrid.addChild(new JTextShape(txtemp, px + 0.5f * pw, py + 1.04f * ph + (float) txtsize,
							(float) (1.1d * txtsize), PConstants.CENTER, PConstants.TOP, 0xff000000, 0, null));
				}
			}
			if (xTkOn) {
				for (int t = 2; t < ticks.length; t++)
					if (ticks[t] >= Math.min(Xin, Xax) && ticks[t] <= Math.max(Xin, Xax)) {
						axisgrid.addChild(new JLineShape(2f, 0xff000000, (float) tcpos[t], py + ph, (float) tcpos[t], py + 1.02f * ph));
						// axisgrid.addChild(ap.createShape(PShape.TEXT, "H",
						// (float)tcpos[t],py-0.1f*ph,(float)tcpos[t],py));
						axisgrid.addChild(new JTextShape(tickmark[t], (float) tcpos[t], py + 1.03f * ph,
								(float) txtsize, PConstants.CENTER, PConstants.TOP, 0xff000000, 0, null));
					}
				if(titleX.length()==0 && unitX.length()==0 && tickmarkFactor.length()>0) {
					if(JPlot.supportLatex) {
						axisgrid.addChild(new JLatexShape(tickmarkFactor, px+pw, py+ph, (float) txtsize,
								PConstants.LEFT, PConstants.CENTER, 0xff000000, 0, null));
					} else {
						axisgrid.addChild(new JTextShape(tickmarkFactor, px+pw, py+ph, (float) txtsize,
								PConstants.LEFT, PConstants.CENTER, 0xff000000, 0, null));
					}
				}
			}
		}
		return axisgrid;
	}

	private JGroupShape createYAxis() {
		double Yin = yLog ? Math.log10(minY) : minY;
		double Yax = yLog ? Math.log10(maxY) : maxY;
		JGroupShape axisgrid = new JGroupShape();
		double[] ticks = null;
		String[] tickmark = null;
		String   tickmarkFactor = "";
		if (yLog) {
			ticks = JPlotMath.optimalLogarithmicTicks(minY, maxY);
			tickmark = new String[ticks.length];
			for (int t = 2; t < ticks.length; t++) {
				tickmark[t] = "" + ticks[t] + "";
				ticks[t] = Math.log10(ticks[t]);
			}
		} else if (yTim) {
			ticks = JPlotMath.optimalTimeTicks(minY, maxY, yTimUnit, yTimCal);
			tickmark = new String[ticks.length];
			for (int t = 2; t < ticks.length; t++)
				tickmark[t] = DateTime.fromDouble(ticks[t], yTimUnit, yTimCal).format(yTimFormat, yTimCal);
		} else {
			ticks = JPlotMath.optimalLinearTicks(minY, maxY);
			double vf = 1d / (ticks[0]);
			int decimal = (int) (1000d * ticks[1] + 0.5d);
			decimal = decimal % 1000 == 0 ? 0 : decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
			tickmark = new String[ticks.length];
			for (int t = 2; t < ticks.length; t++) {
				if(decimal==0) tickmark[t] = ""+(int)(ticks[t]*vf+0.0005d - (ticks[t]*vf<0d?1:0));
				tickmark[t] = PApplet.nf((float) (ticks[t] * vf), 0, decimal).replace(",", ".");
			}
			double lvf = Math.log10(ticks[0]);
			if(Math.abs(lvf)>2.9d) {
				int ivf = (int) lvf - (lvf<0d ? 1 : 0);
				tickmarkFactor = JPlot.supportLatex ? "10$^{"+ivf+"}$" : "10^"+ivf;
			}
		}
		double[] tcpos = JPlotMath.map(ticks, Yin, Yax, py + ph, py);
		if (yAxInv)
			for (int t = 0; t < tcpos.length; t++)
				tcpos[t] = 2 * py + ph - tcpos[t];
		if (pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for (int t = 2; t < ticks.length; t++) {
				tickStr += ", " + tickmark[t];
				posStr += ", " + PApplet.nf((float) tcpos[t], 0, 2).replace(",", ".");
			}
			System.out.println("[DEBUG] JAxis-object: Ytickfactors={p10: " + ticks[0] + ", f: " + ticks[1] + "}");
			System.out.println("[DEBUG] JAxis-object: Ytickval={" + tickStr.substring(2) + "}");
			System.out.println("[DEBUG] JAxis-object: Ytickpos={" + posStr.substring(2) + "}");
		}
		if (yGrdOn) {
			for (int t = 2; t < ticks.length; t++)
				if (ticks[t] >= Math.min(Yin, Yax) && ticks[t] <= Math.max(Yin, Yax))
					axisgrid.addChild(new JLineShape(2f, 0xff999999, px, (float) tcpos[t], px + pw, (float) tcpos[t]));
		}
		if (yAxOn) {
			float tw = 0f;
			if (yTkOn) {
				for (int t = 2; t < ticks.length; t++)
					if (ticks[t] >= Math.min(Yin, Yax) && ticks[t] <= Math.max(Yin, Yax)) {
						axisgrid.addChild(new JLineShape(2f, 0xff000000, px - 0.02f * pw, (float) tcpos[t], px, (float) tcpos[t]));
						tw = Math.max(tw, (float) txtsize * pplot.getGraphic().textWidth(tickmark[t])
								/ pplot.getGraphic().textSize);
						axisgrid.addChild(new JTextShape(tickmark[t], px - 0.03f * pw, (float) tcpos[t],
								(float) txtsize, PConstants.RIGHT, PConstants.CENTER, 0xff000000, 0, null));
					}
				if(titleY.length()==0 && tickmarkFactor.length()>0) {
					if(JPlot.supportLatex) {
						axisgrid.addChild(new JLatexShape(tickmarkFactor, px, py, (float) txtsize,
								PConstants.CENTER, PConstants.BOTTOM, 0xff000000, 0, null));
					} else {
						axisgrid.addChild(new JTextShape(tickmarkFactor, px, py, (float) txtsize,
								PConstants.CENTER, PConstants.BOTTOM, 0xff000000, 0, null));
					}
				}
			}
			if (titleY.length() > 0 || unitY.length() > 0) {
				if (pplot.isDebug())
					System.out.println(
							"[DEBUG] JAxi-object: add y-axis title \""+titleY+"\" and unit \""+unitY+"\" with text size "+txtsize);
				String tytemp = ""+titleY;
				if(unitY.length()>0 || tickmarkFactor.length()>0) tytemp += (titleY.length()>0?" ":"")+"[";
				if(tickmarkFactor.length() > 0) tytemp += tickmarkFactor;
				if(unitY.length()>0) tytemp += (tickmarkFactor.length()>0?" ":"")+unitY;
				if(unitY.length()>0 || tickmarkFactor.length()>0) tytemp += "]";
				if(JPlot.supportLatex) {
					axisgrid.addChild(new JLatexShape(tytemp, px - 0.03f * pw - tw, py + 0.5f * ph, (float) (1.1d * txtsize),
							PConstants.CENTER, PConstants.BOTTOM, 0xff000000, JPlotShape.ROTATE_COUNTERCLOCKWISE, null));
				} else {
					axisgrid.addChild(new JTextShape(tytemp, px - 0.03f * pw - tw, py + 0.5f * ph, (float) (1.1d * txtsize),
							PConstants.CENTER, PConstants.BOTTOM, 0xff000000, JPlotShape.ROTATE_COUNTERCLOCKWISE, null));
				}
			}
		}
		return axisgrid;
	}
	
	private static PImage loadPreDefImg(PApplet applet, String name) {
		BufferedImage bimg;
		try {
			bimg = ImageIO.read(JPlot.class.getResourceAsStream("/data/" + name + ".png"));
			PImage limg = applet.createImage(bimg.getWidth(), bimg.getHeight(), PConstants.ARGB);
			bimg.getRGB(0, 0, bimg.getWidth(), bimg.getHeight(), limg.pixels, 0, bimg.getWidth());
			loadedPreDefImgs.put(name, limg);
			return limg;
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.err.println("The image \"" + name
				+ ".png\" is missing in the jar or inaccessible. Contact the author of ##library.name## ##library.prettyVersion## if the image should be there.");
		return null;
	}
}
