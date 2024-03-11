package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import jplots.JPlot;
import jplots.axes.JAxis;
import jplots.axes.LogarithmicScale;
import jplots.colour.JColourtable;
import jplots.helper.GeometryTools;
import jplots.maths.AffineBuilder;
import jplots.maths.JDEdge;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import jplots.maths.JDTriangle;
import jplots.maths.JDelaunayTriangulator;
import jplots.maths.JGridTriangulator;
import jplots.maths.JPlotMath;
import jplots.shapes.JDGeometry;
import jplots.shapes.JGroupShape;
import jplots.shapes.JImageShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPolygonShape;
import jplots.transform.JProjection;
import processing.core.PConstants;
import processing.core.PGraphics;

public class JContourLayer extends JPlotsLayer {

	private double EPSILON = Math.pow(2, -52);
	private boolean isFilled, pixelFilling, input2d;
	private double minZ, maxZ, Xin, Xax, Yin, Yax;
	private double[] xarrayx, yarrayy;
	private double[][] xarrayx2, yarrayy2, zarrayz;
	private double[] contourIntervals;
	private int[] startEdge;
	private String[] contourStyle;
	private List<JDPoint> corners, cntCorner;
	private List<JDEdge> edges, contours;
	private List<JDTriangle> triangles;
	private JDGeometry[] fillings;
	
	public JContourLayer(float[] x, float[] y, float[][] z, float zmin, float zmax, int nintervals, JColourtable ct,
			float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		float zin = Float.isNaN(zmin) ? JPlotMath.fmin(z) : zmin;
		float zax = Float.isNaN(zmax) ? JPlotMath.fmax(z) : zmax;
		float[] cntIntervals = new float[nintervals + 1];
		for (int k = 0; k <= nintervals; k++)
			cntIntervals[k] = zin + k * (zax - zin) / nintervals;
		input2d = false;
		JContourLayerFloat(x, y, null, null, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(float[] x, float[] y, float[][] z, float[] intervals, JColourtable ct, float stroke_weight,
			boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = false;
		JContourLayerFloat(x, y, null, null, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(double[] x, double[] y, double[][] z, double zmin, double zmax, int nintervals,
			JColourtable ct, double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		double zin = Double.isNaN(zmin) ? JPlotMath.dmin(z) : zmin;
		double zax = Double.isNaN(zmax) ? JPlotMath.dmax(z) : zmax;
		double[] cntIntervals = new double[nintervals + 1];
		for (int k = 0; k <= nintervals; k++)
			cntIntervals[k] = zin + k * (zax - zin) / nintervals;
		input2d = false;
		JContourLayerDouble(x, y, null, null, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(double[] x, double[] y, double[][] z, double[] intervals, JColourtable ct,
			double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = false;
		JContourLayerDouble(x, y, null, null, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(float[][] x, float[][] y, float[][] z, float zmin, float zmax, int nintervals, JColourtable ct,
			float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		float zin = Float.isNaN(zmin) ? JPlotMath.fmin(z) : zmin;
		float zax = Float.isNaN(zmax) ? JPlotMath.fmax(z) : zmax;
		float[] cntIntervals = new float[nintervals + 1];
		for (int k = 0; k <= nintervals; k++)
			cntIntervals[k] = zin + k * (zax - zin) / nintervals;
		input2d = true;
		JContourLayerFloat(null, null, x, y, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(float[][] x, float[][] y, float[][] z, float[] intervals, JColourtable ct, float stroke_weight,
			boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = true;
		JContourLayerFloat(null, null, x, y, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(double[][] x, double[][] y, double[][] z, double zmin, double zmax, int nintervals,
			JColourtable ct, double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		double zin = Double.isNaN(zmin) ? JPlotMath.dmin(z) : zmin;
		double zax = Double.isNaN(zmax) ? JPlotMath.dmax(z) : zmax;
		double[] cntIntervals = new double[nintervals + 1];
		for (int k = 0; k <= nintervals; k++)
			cntIntervals[k] = zin + k * (zax - zin) / nintervals;
		input2d = true;
		JContourLayerDouble(null, null, x, y, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(double[][] x, double[][] y, double[][] z, double[] intervals, JColourtable ct,
			double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = true;
		JContourLayerDouble(null, null, x, y, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}

	/*
	public JContourLayer(float[] x1, float[] y1, float[][] x2, float[][] y2, float[][] z, float zmin, float zmax, int nintervals, float[] zintervals,
			JColourtable ct, float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		
		
		
	}
	*/
	
	
	private void JContourLayerFloat(float[] x, float[] y, float[][] x2, float[][] y2, float[][] z, float[] intervals,
			JColourtable ct, float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		if (!input2d) {
			xarrayx = new double[x.length];
			for (int i = 0; i < x.length; i++)
				xarrayx[i] = x[i];
			yarrayy = new double[y.length];
			for (int i = 0; i < y.length; i++)
				yarrayy[i] = y[i];
			xarrayx2 = null;
			yarrayy2 = null;
			minX = JPlotMath.dmin(xarrayx);
			maxX = JPlotMath.dmax(xarrayx);
			minY = JPlotMath.dmin(yarrayy);
			maxY = JPlotMath.dmax(yarrayy);
		} else {
			xarrayx = null;
			yarrayy = null;
			xarrayx2 = new double[x2.length][x2[0].length];
			for (int j = 0; j < x2.length; j++)
				for (int i = 0; i < x2[0].length; i++)
					xarrayx2[j][i] = x2[j][i];
			yarrayy2 = new double[y2.length][y2[0].length];
			for (int j = 0; j < y2.length; j++)
				for (int i = 0; i < y2[0].length; i++)
					yarrayy2[j][i] = y2[j][i];
			minX = JPlotMath.dmin(xarrayx2);
			maxX = JPlotMath.dmax(xarrayx2);
			minY = JPlotMath.dmin(yarrayy2);
			maxY = JPlotMath.dmax(yarrayy2);
		}
		zarrayz = new double[z.length][z[0].length];
		for (int j = 0; j < z.length; j++)
			for (int i = 0; i < z[j].length; i++)
				zarrayz[j][i] = z[j][i];
		contourIntervals = new double[intervals.length];
		for (int i = 0; i < intervals.length; i++)
			contourIntervals[i] = intervals[i];
		minZ = contourIntervals[0];
		maxZ = contourIntervals[intervals.length - 1];
		colourtable = ct;
		lw = stroke_weight;
		drawLines = drawContours;
		isFilled = filled;
		pixelFilling = filledAsImage;
		init();
	}

	private void JContourLayerDouble(double[] x, double[] y, double[][] x2, double[][] y2, double[][] z,
			double[] intervals, JColourtable ct, double stroke_weight, boolean drawContours, boolean filled,
			boolean filledAsImage) {
		if (!input2d) {
			xarrayx = x;
			yarrayy = y;
			xarrayx2 = null;
			yarrayy2 = null;
			minX = JPlotMath.dmin(xarrayx);
			maxX = JPlotMath.dmax(xarrayx);
			minY = JPlotMath.dmin(yarrayy);
			maxY = JPlotMath.dmax(yarrayy);
		} else {
			xarrayx = null;
			yarrayy = null;
			xarrayx2 = x2;
			yarrayy2 = y2;
			minX = JPlotMath.dmin(xarrayx2);
			maxX = JPlotMath.dmax(xarrayx2);
			minY = JPlotMath.dmin(yarrayy2);
			maxY = JPlotMath.dmax(yarrayy2);
		}
		zarrayz = z;
		contourIntervals = new double[intervals.length];
		for (int i = 0; i < intervals.length; i++)
			contourIntervals[i] = intervals[i];
		minZ = contourIntervals[0];
		maxZ = contourIntervals[intervals.length - 1];
		colourtable = ct;
		lw = (float) stroke_weight;
		drawLines = drawContours;
		isFilled = filled;
		pixelFilling = filledAsImage;
		init();
	}
	
	private void init() {
		contourStyle = new String[contourIntervals.length];
		for (int cs = 0; cs < contourStyle.length; cs++)
			contourStyle[cs] = "-";

		corners = new ArrayList<>();
		edges = new ArrayList<>();
		triangles = new ArrayList<>();

		cntCorner = new ArrayList<>();
		contours = new ArrayList<>();
		new ArrayList<JDTriangle>();

		fillings = null;
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		//TODO: deal with JMultiAxis
		boolean isXlog = (ax.getScaleX() instanceof LogarithmicScale);
		boolean isYlog = (ax.getScaleY() instanceof LogarithmicScale);
		Xin = isXlog ? Math.log10(minX) : minX;
		Xax = isXlog ? Math.log10(maxX) : maxX;
		Yin = isYlog ? Math.log10(minY) : minY;
		Yax = isYlog ? Math.log10(maxY) : maxY;
		double xs = p[2] / (Xax - Xin), ys = p[3] / (Yax - Yin);
		// double tol = Math.max(Math.abs(maxX-minX), Math.abs(maxY-minY)) * 1.0e-12d;

		// step 1: collect valid corners of grid and project them
		collectValidPoints(ax.getGeoProjection(), isXlog, isYlog, ax.getPlot().isDebug());
		
		// step 2: do delauney-triangulation
		triangulate(ax);
		
		AffineBuilder affine = new AffineBuilder().scale(invertAxisX ? -1d : 1d, invertAxisY ? 1d : -1d)
				.translate(invertAxisX ? maxX : -minX, invertAxisY ? -minY : maxY).scale(xs, ys).translate(p[0], p[1]);
		for(JDTriangle tri: triangles)
			tri.affine(affine.getMatrix());
		
		// step 3: create contours
		if (drawLines || !pixelFilling) {
			createContours(ax.getPlot().isDebug());
		} else {
			if (ax.getPlot().isDebug())
				System.out.println("[DEBUG] JContourLayer: 3] no contours ...");
		}
		
		// step 4: create filling between contours if wished
		if (isFilled) {
			if (pixelFilling) {
				if (ax.getPlot().isDebug())
					System.out.println("[DEBUG] JContourLayer: 4] contour filling pixelwise ...");
				fillPixelByPixel(p, ax, xs, ys, s);
			} else {
				try {
					fillVectorByVector(ax.getGeoProjection(), p, ax, xs, ys, s);
				} catch (Exception e) {
					e.printStackTrace();
					throw new RuntimeException(e);
				}
			}
		} else {
			if (ax.getPlot().isDebug())
				System.out.println("[DEBUG] JContourLayer: 4] no filling ...");
		}

		// step 5: add contours to plot
		if (drawLines) {
			if (ax.getPlot().isDebug())
				System.out.println("[DEBUG] JContourLayer: 5] register " + (contours.size() - startEdge[0])
						+ " contour line segments with " + lw + "px line width...");
			drawContourLines(p, ax, xs, ys, s);
		} else {
			if (ax.getPlot().isDebug())
				System.out.println("[DEBUG] JContourLayer: 5] contours itself will not be drawn ...");
		}
		
		if(ax.getPlot().isDebug()) {
			drawTriangleBorders(p,ax,xs,xs,s);
		}
	}

	//* **************************************** *
	//* ********** GETTER AND SETTER  ********** *
	//* **************************************** *
	
	public double[] getZRange() {
		return new double[] { minZ, maxZ };
	}

	public double[] getLevels() {
		return contourIntervals;
	}
	
	@Override
	public JColourtable getColourtable() {
		if(!isFilled) return null;
		return super.getColourtable();
	}
	
	//* **************************************** *
	//* ********** PUBLIC METHODS     ********** *
	//* **************************************** *
	
	public void setStyles(String[] nstyle) {
		contourStyle = nstyle;
	}
	
	public void drawTriangleBorders(int[] p, JAxis ax, double xs, double ys, JGroupShape s) {
		
	}
	
	

	//* **************************************** *
	//* ********** PRIVATE METHODS    ********** *
	//* **************************************** *
	
	private void collectValidPoints(JProjection outproj, boolean xLog, boolean yLog, boolean debug) {
		int nx = 0, ny = 0;
		if (input2d) {
			ny = xarrayx2.length;
			if (yarrayy2.length != ny)
				throw new IllegalArgumentException("x, y and z are of different shapes!");
			nx = xarrayx2.length;
			if (yarrayy2.length != nx)
				throw new IllegalArgumentException("x, y and z are of different shapes!");
		} else {
			nx = xarrayx.length;
			ny = yarrayy.length;
		}
		corners.clear();
		double[] xy;
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++) {
				if (input2d) {
					xy = inputProj.fromPROJtoLATLON(xLog ? Math.log10(xarrayx2[j][i]) : xarrayx2[j][i],
							yLog ? Math.log10(yarrayy2[j][i]) : yarrayy2[j][i], false, false);
				} else {
					xy = inputProj.fromPROJtoLATLON(xLog ? Math.log10(xarrayx[i]) : xarrayx[i],
							yLog ? Math.log10(yarrayy[j]) : yarrayy[j], false, false);
				}
				xy = outproj.fromLATLONtoPROJ(xy[0], xy[1], false, false);
				if (Double.isFinite(xy[0]) && Double.isFinite(xy[1]))
					corners.add(new JDPoint(xy[0], xy[1], zarrayz[j][i]));
			}
		double xyScale = Math.max(Math.max(-Xin, Xax), Math.max(-Yin, Yax));
		xyScale = Math.max(xyScale, Math.max(Xax - Xin, Yax - Yin));
		EPSILON = xyScale * 1.e-12d;
		if (debug)
			System.out.println("[DEBUG] JContourLayer: 1] has " + corners.size() + " valid sourcepoints");
	}

	private void triangulate(JAxis ax) {
		boolean debug = ax.getPlot().isDebug();
		//TODO: deal with JMultiAxis
		boolean isXlog = (ax.getScaleX() instanceof LogarithmicScale);
		boolean isYlog = (ax.getScaleY() instanceof LogarithmicScale);
		if (debug)
			System.out.println("[DEBUG] JContourLayer: 2] triangulate ...");
		if (input2d || ax.isGeoAxis()) {
			JDelaunayTriangulator delaunay = new JDelaunayTriangulator(corners);
			edges = delaunay.getEdges();
			triangles = delaunay.getTriangles();
		} else {
			JGridTriangulator triangulator = new JGridTriangulator(isXlog ? JPlotMath.log10(xarrayx) : xarrayx,
					isYlog ? JPlotMath.log10(yarrayy) : yarrayy, zarrayz);
			edges = triangulator.getEdges();
			triangles = triangulator.getTriangles();
		}
		if (debug)
			System.out.println("[DEBUG] JContourLayer:   mesh now consists of " + corners.size() + " corners, "
					+ edges.size() + " edges and " + triangles.size() + " triangles");
	}

	private void createContours(boolean debug) {
		if (debug)
			System.out.println("[DEBUG] JContourLayer: 3] create contours from triangle mesh ...");
		startEdge = new int[contourIntervals.length];
		int nanid = -9999;
		contours.clear();
		// first exclude fillvalues
		for (int level = nanid; level < contourIntervals.length; level = Math.max(0, level + 1)) {
			if (debug)
				System.out.println("[DEBUG] JContourLayer:   create contours for level "
						+ (level < 0 ? Double.NaN : contourIntervals[level]));
			for (JDTriangle ct : triangles) {
				JDPoint va = ct.getA();
				JDPoint vb = ct.getB();
				JDPoint vc = ct.getC();
				if (level == nanid) {
					int nancode = 0;
					if (Double.isNaN(va.value))
						nancode |= 1;
					if (Double.isNaN(vb.value))
						nancode |= 2;
					if (Double.isNaN(vc.value))
						nancode |= 4;
					switch (nancode) {
					case 1:
					case 6:
						double x161 = 0.5d * (va.x + vb.x);
						double y161 = 0.5d * (va.y + vb.y);
						int i16s = addCorner(x161, y161, Double.NaN, cntCorner);
						double x162 = 0.5d * (va.x + vc.x);
						double y162 = 0.5d * (va.y + vc.y);
						int i16m = i16s;
						for (int ci = 0; ci < contourIntervals.length; ci++) {
							int ic = vb.value < vc.value ? ci : contourIntervals.length - 1 - ci;
							if ((vb.value < contourIntervals[ic] && contourIntervals[ic] < vc.value)
									|| (vc.value < contourIntervals[ic] && contourIntervals[ic] < vb.value)) {
								double f16 = (contourIntervals[ic] - vb.value) / (vc.value - vb.value);
								i16m = addCorner(x161 + f16 * (x162 - x161), y161 + f16 * (y162 - y161), Double.NaN,
										cntCorner);
								contours.add(new JDEdge(cntCorner.get(i16s), cntCorner.get(i16m)));
								i16s += i16m - i16s;
							}
						}
						int i16e = addCorner(x162, y162, Double.NaN, cntCorner);
						contours.add(new JDEdge(cntCorner.get(i16s), cntCorner.get(i16e)));
						break;
					case 2:
					case 5:
						double x251 = 0.5d * (vb.x + va.x);
						double y251 = 0.5d * (vb.y + va.y);
						int i25s = addCorner(x251, y251, Double.NaN, cntCorner);
						double x252 = 0.5d * (vb.x + vc.x);
						double y252 = 0.5d * (vb.y + vc.y);
						int i25m = i25s;
						for (int ci = 0; ci < contourIntervals.length; ci++) {
							int ic = va.value < vc.value ? ci : contourIntervals.length - 1 - ci;
							if ((va.value < contourIntervals[ic] && contourIntervals[ic] < vc.value)
									|| (vc.value < contourIntervals[ic] && contourIntervals[ic] < va.value)) {
								double f25 = (contourIntervals[ic] - va.value) / (vc.value - va.value);
								i25m = addCorner(x251 + f25 * (x252 - x251), y251 + f25 * (y252 - y251), Double.NaN,
										cntCorner);
								contours.add(new JDEdge(cntCorner.get(i25s), cntCorner.get(i25m)));
								i25s += i25m - i25s;
							}
						}
						int i25e = addCorner(x252, y252, Double.NaN, cntCorner);
						contours.add(new JDEdge(cntCorner.get(i25s), cntCorner.get(i25e)));
						break;
					case 3:
					case 4:
						double x341 = 0.5d * (vc.x + vb.x);
						double y341 = 0.5d * (vc.y + vb.y);
						int i34s = addCorner(x341, y341, Double.NaN, cntCorner);
						double x342 = 0.5d * (vc.x + va.x);
						double y342 = 0.5d * (vc.y + va.y);
						for (int ci = 0; ci < contourIntervals.length; ci++) {
							int ic = vb.value < va.value ? ci : contourIntervals.length - 1 - ci;
							if ((va.value < contourIntervals[ic] && contourIntervals[ic] < vb.value)
									|| (vb.value < contourIntervals[ic] && contourIntervals[ic] < va.value)) {
								double f34 = (contourIntervals[ic] - vb.value) / (va.value - vb.value);
								int i34m = addCorner(x341 + f34 * (x342 - x341), y341 + f34 * (y342 - y341), Double.NaN,
										cntCorner);
								contours.add(new JDEdge(cntCorner.get(i34s), cntCorner.get(i34m)));
								i34s += i34m - i34s;
							}
						}
						int i34e = addCorner(x342, y342, Double.NaN, cntCorner);
						contours.add(new JDEdge(cntCorner.get(i34s), cntCorner.get(i34e)));
						break;
					default:
						break;
					}
				} else {
					double lev = contourIntervals[level];
					double r12 = (lev - va.value) / (vb.value - va.value);
					double r23 = (lev - vb.value) / (vc.value - vb.value);
					double r31 = (lev - vc.value) / (va.value - vc.value);
					int levcode = 0;
					if (Double.isNaN(va.value))
						levcode |= 2;
					if (Double.isNaN(vb.value))
						levcode |= 8;
					if (Double.isNaN(vc.value))
						levcode |= 32;
					if (va.value <= lev)
						levcode |= 1;
					if (vb.value <= lev)
						levcode |= 4;
					if (vc.value <= lev)
						levcode |= 16;
//					if(debug)
//						System.out.println("[DEBUG] JContourplot:   triangle "+t_idx+" -> levcode "+Integer.toBinaryString(levcode));
					switch (levcode) {
					case 1: // x000001
					case 20: // x010100
						double x161 = va.x + r12 * (vb.x - va.x);
						double y161 = va.y + r12 * (vb.y - va.y);
						int i17s = addCorner(x161, y161, lev, cntCorner);
						double x162 = vc.x + r31 * (va.x - vc.x);
						double y162 = vc.y + r31 * (va.y - vc.y);
						int i17e = addCorner(x162, y162, lev, cntCorner);
						if (i17s != i17e)
							contours.add(levcode == 1 ? new JDEdge(cntCorner.get(i17s), cntCorner.get(i17e))
									: new JDEdge(cntCorner.get(i17e), cntCorner.get(i17s)));
						break;
					case 5: // x000101
					case 16: // x010000
						double x341 = vb.x + r23 * (vc.x - vb.x);
						double y341 = vb.y + r23 * (vc.y - vb.y);
						int i20s = addCorner(x341, y341, lev, cntCorner);
						double x342 = vc.x + r31 * (va.x - vc.x);
						double y342 = vc.y + r31 * (va.y - vc.y);
						int i20e = addCorner(x342, y342, lev, cntCorner);
						if (i20s != i20e)
							contours.add(levcode == 5 ? new JDEdge(cntCorner.get(i20s), cntCorner.get(i20e))
									: new JDEdge(cntCorner.get(i20e), cntCorner.get(i20s)));
						break;
					case 4: // x000100
					case 17: // x010001
						double x251 = va.x + r12 * (vb.x - va.x);
						double y251 = va.y + r12 * (vb.y - va.y);
						int i05s = addCorner(x251, y251, lev, cntCorner);
						double x252 = vb.x + r23 * (vc.x - vb.x);
						double y252 = vb.y + r23 * (vc.y - vb.y);
						int i05e = addCorner(x252, y252, lev, cntCorner);
						if (i05s != i05e)
							contours.add(levcode == 17 ? new JDEdge(cntCorner.get(i05s), cntCorner.get(i05e))
									: new JDEdge(cntCorner.get(i05e), cntCorner.get(i05s)));
						break;
					case 6: // x000110
					case 18: // x010010
						double x25c = vb.x + r23 * (vc.x - vb.x);
						double y25c = vb.y + r23 * (vc.y - vb.y);
						int i22s = addCorner(x25c, y25c, lev, cntCorner);
						double x23b = 0.5d * (vb.x + va.x), y23b = 0.5d * (vb.y + va.y);
						double x23c = 0.5d * (vc.x + va.x), y23c = 0.5d * (vc.y + va.y);
						int i22e = addCorner(x23b + r23 * (x23c - x23b), y23b + r23 * (y23c - y23b), Double.NaN,
								cntCorner);
						contours.add(levcode == 6 ? new JDEdge(cntCorner.get(i22s), cntCorner.get(i22e))
								: new JDEdge(cntCorner.get(i22e), cntCorner.get(i22s)));
						break;
					case 9: // x001001
					case 24: // x011000
						double x34c = vc.x + r31 * (va.x - vc.x);
						double y34c = vc.y + r31 * (va.y - vc.y);
						int i26s = addCorner(x34c, y34c, lev, cntCorner);
						double x31c = 0.5d * (vc.x + vb.x), y31c = 0.5d * (vc.y + vb.y);
						double x31a = 0.5d * (va.x + vb.x), y31a = 0.5d * (va.y + vb.y);
						int i26e = addCorner(x31c + r31 * (x31a - x31c), y31c + r31 * (y31a - y31c), Double.NaN,
								cntCorner);
						contours.add(levcode == 24 ? new JDEdge(cntCorner.get(i26s), cntCorner.get(i26e))
								: new JDEdge(cntCorner.get(i26e), cntCorner.get(i26s)));
						break;
					case 33: // x100001
					case 36: // x100100
						double x12c = va.x + r12 * (vb.x - va.x);
						double y12c = va.y + r12 * (vb.y - va.y);
						int i41s = addCorner(x12c, y12c, Double.NaN, cntCorner);
						double x12a = 0.5d * (va.x + vc.x), y12a = 0.5d * (va.y + vc.y);
						double x12b = 0.5d * (vb.x + vc.x), y12b = 0.5d * (vb.y + vc.y);
						int i41e = addCorner(x12a + r12 * (x12b - x12a), y12a + r12 * (y12b - y12a), Double.NaN,
								cntCorner);
						contours.add(levcode == 33 ? new JDEdge(cntCorner.get(i41s), cntCorner.get(i41e))
								: new JDEdge(cntCorner.get(i41e), cntCorner.get(i41s)));
						break;
//						case 21: //x010101
//							int special = 0;
//							if(va.val==lev) special |= 1;
//							if(vb.val==lev) special |= 2;
//							if(vc.val==lev) special |= 4;
//							switch(special) {
//								case 3: //b011
//									int i3a = addCorner(va.x,va.y,lev,-1,cntCorner);
//									int i3b = addCorner(vb.x,vb.y,lev,-1,cntCorner);
//									contours.add(new jedge(i3a, i3b, cntCorner));
//									break;
//								case 5: //b101
//									int i5a = addCorner(va.x,va.y,lev,-1,cntCorner);
//									int i5b = addCorner(vc.x,vc.y,lev,-1,cntCorner);
//									contours.add(new jedge(i5a, i5b, cntCorner));
//									break;
//								case 6: //b110
//									int i6a = addCorner(vb.x,vb.y,lev,-1,cntCorner);
//									int i6b = addCorner(vc.x,vc.y,lev,-1,cntCorner);
//									contours.add(new jedge(i6a, i6b, cntCorner));
//									break;
////								case 7: //b111
////									int i7a = addCorner(va.x,va.y,lev,-1,cntCorner);
////									int i7b = addCorner(vb.x,vb.y,lev,-1,cntCorner);
////									int i7c = addCorner(vc.x,vc.y,lev,-1,cntCorner);
////									int i7d = addCorner((va.x+vb.x+vc.x)/3d,(va.y+vb.y+vc.y)/3d,lev,-1,cntCorner);
////									contours.add(new jedge(i7a, i7d, cntCorner));
////									contours.add(new jedge(i7b, i7d, cntCorner));
////									contours.add(new jedge(i7c, i7d, cntCorner));
////									break;
//								default:
//									break;
//							}
//							break;
					default:
						break;
					}
				}
			}
			if (level + 1 < contourIntervals.length)
				startEdge[Math.max(0, level + 1)] = contours.size();
		}
		// remove duplicates of contour line segments
//		for(int e1=contours.size()-1; e1>0; e1--)
//			for(int e2=e1-1; e2>=0; e2--)
//				if(contours.get(e1).isSame(contours.get(e2))) {
//					for(int se=0; se<startEdge.length; se++)
//						if(startEdge[se]>=e1)
//							startEdge[se]--;
//					contours.remove(e1); break;
//				}
		if (debug) {
			System.out.println("[DEBUG] JContourLayer:    estimated " + contours.size() + " line segments");
			String sts = "";
			for (int element : startEdge)
				sts += ", " + element;
			System.out.println("[DEBUG] JContourLayer:    start indices are: " + sts.substring(2));
		}
//		if(debug) {
//			System.out.println("[DEBUG] JContourLayer:    all corners");
//			for(int c=0; c<cntCorner.size(); c++)
//				System.out.println("    "+PApplet.nf((float)cntCorner.get(c).x,0,4)+
//						"    "+PApplet.nf((float)cntCorner.get(c).y,0,4)+
//						"    "+PApplet.nf((float)cntCorner.get(c).lev,0,4));
//			System.out.println("[DEBUG] JContourLayer:    all edges / line segments");
//			for(int c=0; c<contours.size(); c++)
//				System.out.println("    "+contours.get(c).a+
//						"    "+contours.get(c).b);
//		}
	}

	private int addCorner(double x, double y, double v, List<JDPoint> points) {
		// double tolerance = 1.0e-10d;
		JDPoint np = new JDPoint(x, y, v);
		for (int p = 0; p < points.size(); p++)
			if (points.get(p).equals(np, EPSILON))
				return p;
		points.add(np);
		return points.size() - 1;
	}

	private void fillVectorByVector(JProjection outproj, int[] p, JAxis ax, double xs, double ys, JGroupShape s) {
		if (ax.getPlot().isDebug())
			System.out.println("[DEBUG] JContourLayer: 4] fill contours ...");
		double eps = Math.min(p[2], p[3]) * 1.0e-9d;
		double eps2 = eps * 0.0001d;
		List<JDTriangle> visibleTriangles = new ArrayList<>();
		for (JDTriangle t : triangles) {
//    		if(t.area()<0d)
//    			t.reverse_orientation();
			JDPolygon poly = t.copy().intersectsAABB(p[0], p[1], p[0] + p[2], p[1] + p[2]);
			if (poly == null)
				continue;
//    		if(poly.area()<0d)
//    			poly.reverse_orientation();
			visibleTriangles.addAll(poly.toTriangles());
		}
		fillings = null;
		fillings = new JDGeometry[contourIntervals.length + 1];
		for (int lev = -1; lev < contourIntervals.length; lev++) {
			List<JDPolygon> pl = new ArrayList<>();
			double levmin = -100000000d;
			double levmax = 100000000d;
			if (lev >= 0)
				levmin = contourIntervals[lev];
			if (lev + 1 < contourIntervals.length)
				levmax = contourIntervals[lev + 1];
			if (ax.getPlot().isDebug())
				System.out.println(
						"[DEBUG]                   ... collect geometry for range [" + levmin + " ... " + levmax + "]");
			for (JDTriangle t : visibleTriangles) {
				Object res = cutoutLevelrange(t, levmin, levmax);
				if (res == null)
					continue;
				if (res instanceof JDTriangle) {
					JDTriangle rt = (JDTriangle) res;
					// System.out.println(" ... add "+rt+" with area "+rt.area()+"<>"+eps2);
					if (Math.abs(rt.area()) > eps2)
						addTriangle2polygonList(pl, rt, eps);
				}
				if (res instanceof JDPolygon) {
					JDPolygon rp = (JDPolygon) res;
					// System.out.println(" ... add "+rp+" with area "+rp.area()+"<>"+eps2);
					if (rp.area() > eps2)
						addPolygon2polygonList(pl, rp, eps);
				}
			}
			fillings[lev + 1] = new JDGeometry();
			for (JDPolygon jdp : pl) {
				if (ax.getPlot().isDebug())
					System.out.println("[DEBUG]                       ... use polygon " + jdp);
				fillings[lev + 1].add(jdp);
			}
			if (ax.getPlot().isDebug() && pl.isEmpty())
				System.out.println("[DEBUG]                       ... no polygon created");
		}

		JGroupShape trianglesh = new JGroupShape();
		for (int lev = -1; lev < contourIntervals.length; lev++) {
			double levmin = Double.NEGATIVE_INFINITY;
			double levmax = Double.POSITIVE_INFINITY;
			if (lev >= 0)
				levmin = contourIntervals[lev];
			if (lev + 1 < contourIntervals.length)
				levmax = contourIntervals[lev + 1];
			if (ax.getPlot().isDebug())
				System.out.println(
						"[DEBUG]                   ... draw geometry for range [" + levmin + " ... " + levmax + "]");
			if (fillings[lev + 1] == null)
				continue;
			double pct = minZ - 10d;
			if (lev >= 0 && lev + 1 < contourIntervals.length)
				pct = 0.5d * (contourIntervals[lev] + contourIntervals[lev + 1]);
			if (lev + 1 == contourIntervals.length)
				pct = maxZ + 10d;
			int cct = colourtable.getColour(pct, contourIntervals[0], contourIntervals[contourIntervals.length - 1]);
//			JPlotShape.fill(cct);
//			JPlotShape.noStroke();
//			if(ax.getPlot().isDebug()) {
//				JPlotShape.stroke(0xff999999); JPlotShape.strokeWeight(2f); }
			for (JDPolygon ppp : fillings[lev + 1].getPolygons()) {
				trianglesh.addChild(new JPolygonShape(ppp, cct, cct, 1f, true, true));
			}
		}
		s.addChild(trianglesh);
	}

	private void fillPixelByPixel(int[] p, JAxis ax, double xs, double ys, JGroupShape s) {
		// double us = srcImg.width/(srcExt[2]-srcExt[0]), vs =
		// srcImg.height/(srcExt[3]-srcExt[1]);
		if (img == null) {
			img = ax.getPlot().getApplet().createImage(p[2], p[3], PConstants.ARGB);
		} else if (img.width != p[2] || img.height != p[3]) {
			img = ax.getPlot().getApplet().createImage(p[2], p[3], PConstants.ARGB);
		}
		img.loadPixels();
		double minCI = contourIntervals[0] < 1d ? contourIntervals[0] * 2d : contourIntervals[0] - 10d;
		double maxCI = contourIntervals[contourIntervals.length - 1] > 1d
				? contourIntervals[contourIntervals.length - 1] * 2d
				: contourIntervals[contourIntervals.length - 1] + 10d;
		for (JDTriangle tri : triangles) {
			double	txi = Math.min(tri.x[0], Math.min(tri.x[1], tri.x[2])),
					txa = Math.max(tri.x[0], Math.max(tri.x[1], tri.x[2]));
			double	tyi = Math.min(tri.y[0], Math.min(tri.y[1], tri.y[2])),
					tya = Math.max(tri.y[0], Math.max(tri.y[1], tri.y[2]));
			int ixs =  Math.max((int) txi    - (txi < 0 ? 1 : 0), p[0]),
				ixe = -Math.max((int) (-txa) - (txa > 0 ? 1 : 0), 1 - p[0] - p[2]);
			int iys =  Math.max((int) tyi    - (tyi < 0 ? 1 : 0), p[1]),
				iye = -Math.max((int) (-tya) - (tya > 0 ? 1 : 0), 1 - p[1] - p[3]);
			if ((ixe < ixs) || (iye < iys))
				continue;
			for (int j = iys; j <= iye; j++) {
				for (int i = ixs; i <= ixe; i++) {
					JDPoint ij = new JDPoint(i+0.5d,j+0.5d);
					if(!tri.contains(ij)) continue;
					double val = tri.valueAt(ij);
					if (Double.isNaN(val))
						continue;
					int il = getLevel(val, contourIntervals, -1);
					double dl = il < 1 ? minCI : maxCI;
					if(il>0 && il<contourIntervals.length)
						dl = 0.5d * (contourIntervals[il - 1] + contourIntervals[il]);
					img.pixels[(j-p[1]) * p[2] + i-p[0]] =
							colourtable.getColour(dl,	contourIntervals[0],
														contourIntervals[contourIntervals.length - 1]);
				}
			}
		}
		img.updatePixels();
		s.addChild(new JImageShape(img, p[0], p[1], p[2], p[3]));
	}
	
	private void drawContourLines(int[] p, JAxis ax, double xs, double ys, JGroupShape s) {
		JGroupShape linesh = new JGroupShape();
		for (int c = 0; c < contourIntervals.length; c++) {
			int cs = startEdge[c];
			int ce = contours.size();
			if (c + 1 < contourIntervals.length)
				ce = startEdge[c + 1];
			// TODO edit colour and linestyle info for different contour lines
			int lc = 0xff000000;
			// String ls = "-"; //TODO recreate line-drawing with different linestyles
			String ls = contourStyle[c];
			double lln = 1d, llf = 0d, lpn = 0d, lpf = 0d, loff = 0d;
			if ("-".equals(ls)) {
				lln = 1000 * lw;
				llf = 0;
				lpn = 0;
				lpf = 0;
			}
			if (".".equals(ls)) {
				lln = 0;
				llf = 0;
				lpn = 1 * lw;
				lpf = 3 * lw;
			}
			if (",".equals(ls)) {
				lln = 8 * lw;
				llf = 7 * lw;
				lpn = 0;
				lpf = 0;
			}
			if (";".equals(ls)) {
				lln = 8 * lw;
				llf = 3 * lw;
				lpn = 1 * lw;
				lpf = 3 * lw;
			}
			int li = 0;
			for (int cl = cs; cl < ce; cl++) {
				JDPoint lvs = contours.get(cl).a;
				JDPoint lve = contours.get(cl).b;
				double x1 = lvs.x; //p[0] + xs * (invertAxisX ? Xax - lvs.x : lvs.x - Xin);
				double x2 = lve.x; //p[0] + xs * (invertAxisX ? Xax - lve.x : lve.x - Xin);
				double y1 = lvs.y; //p[1] + ys * (invertAxisY ? lvs.y - Yin : Yax - lvs.y);
				double y2 = lve.y; //p[1] + ys * (invertAxisY ? lve.y - Yin : Yax - lve.y);
				double dx = x2 - x1, dy = y2 - y1;
				double l = Math.sqrt(dx * dx + dy * dy);
				dx /= l;
				dy /= l;
				double lpos = 0d, ldif = 0d;
				while (lpos < l) {
					ldif = 0d;
					switch (li) {
					case 0:
						if (lln == 0d)
							break;
						ldif = Math.min(l - lpos, lln - loff);
						break;
					case 1:
						if (llf == 0d)
							break;
						ldif = Math.min(l - lpos, llf - loff);
						break;
					case 2:
						if (lpn == 0d)
							break;
						ldif = Math.min(l - lpos, lpn - loff);
						break;
					case 3:
						if (lpf == 0d)
							break;
						ldif = Math.min(l - lpos, lpf - loff);
						break;
					}
					float xf1 = (float) (x1 + lpos * dx), yf1 = (float) (y1 + lpos * dy),
							xf2 = (float) (x1 + (lpos + ldif) * dx), yf2 = (float) (y1 + (lpos + ldif) * dy);
					if (xf1 < p[0] && xf2 >= p[0]) {
						yf1 = JPlotMath.map(p[0], xf1, xf2, yf1, yf2);
					}
					if (xf1 > p[0] + p[2] && xf2 <= p[0] + p[2]) {
						yf1 = JPlotMath.map(p[0] + p[2], xf1, xf2, yf1, yf2);
					}
					if (xf2 < p[0] && xf1 >= p[0]) {
						yf2 = JPlotMath.map(p[0], xf1, xf2, yf1, yf2);
					}
					if (xf2 > p[0] + p[2] && xf1 <= p[0] + p[2]) {
						yf2 = JPlotMath.map(p[0] + p[2], xf1, xf2, yf1, yf2);
					}
					if (yf1 < p[1] && yf2 >= p[1]) {
						xf1 = JPlotMath.map(p[1], yf1, yf2, xf1, xf2);
					}
					if (yf1 > p[1] + p[3] && yf2 <= p[1] + p[3]) {
						xf1 = JPlotMath.map(p[1] + p[3], yf1, yf2, xf1, xf2);
					}
					if (yf2 < p[1] && yf1 >= p[1]) {
						xf2 = JPlotMath.map(p[1], yf1, yf2, xf1, xf2);
					}
					if (yf2 > p[1] + p[3] && yf1 <= p[1] + p[3]) {
						xf2 = JPlotMath.map(p[1] + p[3], yf1, yf2, xf1, xf2);
					}
					if (xf1 >= p[0] && xf1 <= p[0] + p[2] && xf2 >= p[0] && xf2 <= p[0] + p[2] && yf1 >= p[1]
							&& yf1 <= p[1] + p[3] && yf2 >= p[1] && yf2 <= p[1] + p[3]) {
						if (li % 2 == 0 && ldif > 0d)
							linesh.addChild(new JLineShape((float)lw, lc, xf1, yf1, xf2, yf2));
						loff += ldif;
						switch (li) {
						case 0:
							if (loff >= lln) {
								loff -= lln;
								li = 1;
							}
							break;
						case 1:
							if (loff >= llf) {
								loff -= llf;
								li = 2;
							}
							break;
						case 2:
							if (loff >= lpn) {
								loff -= lpn;
								li = 3;
							}
							break;
						case 3:
							if (loff >= lpf) {
								loff -= lpf;
								li = 0;
							}
							break;
						}
					}
					lpos += ldif;
				}
			}
		}
		s.addChild(linesh);
	}
	
	//* **************************************** *
	//* ********** STATIC METHODS     ********** *
	//* **************************************** *
	
	public static int getLevel(double value, double[] intervalBorders, int nanLev) {
		if (Double.isNaN(value))
			return nanLev;
		int l = 0;
		for (int cl = 0; cl < intervalBorders.length; cl++)
			if (intervalBorders[cl] < value)
				l = cl + 1;
		return l;
	}

	public Object cutoutLevelrange(JDTriangle tri, double lower, double upper) {
		double x1 = tri.x[0], x2 = tri.x[1], x3 = tri.x[2];
		double y1 = tri.y[0], y2 = tri.y[1], y3 = tri.y[2];
		double v1 = tri.value[0], v2 = tri.value[1], v3 = tri.value[2];
		double vmin = Double.POSITIVE_INFINITY;
		double vmax = Double.NEGATIVE_INFINITY;
		if (Double.isFinite(v1)) {
			if (v1 < vmin)
				vmin = v1;
			if (v1 > vmax)
				vmax = v1;
		}
		if (Double.isFinite(v2)) {
			if (v2 < vmin)
				vmin = v2;
			if (v2 > vmax)
				vmax = v2;
		}
		if (Double.isFinite(v3)) {
			if (v3 < vmin)
				vmin = v3;
			if (v3 > vmax)
				vmax = v3;
		}
		if ((vmax < vmin) || vmax < lower || vmin > upper)
			return null;
		int id = ((Double.isNaN(v1) ? 0 : v1 < lower ? 1 : v1 > upper ? 2 : 3) << 8)
				| ((Double.isNaN(v2) ? 0 : v2 < lower ? 1 : v2 > upper ? 2 : 3) << 4)
				| (Double.isNaN(v3) ? 0 : v3 < lower ? 1 : v3 > upper ? 2 : 3);
//		String s_id = Integer.toHexString(id);
//		while(s_id.length()<3) s_id = "0"+s_id;
//		System.out.println("Triangel "+tri+" -> id="+s_id);
		if (id == 0x000 || id == 0x111 || id == 0x222)
			return null;
		if (id == 0x001 || id == 0x010 || id == 0x100 || id == 0x002 || id == 0x020 || id == 0x200)
			return null;
		if (id == 0x011 || id == 0x101 || id == 0x110 || id == 0x022 || id == 0x202 || id == 0x220)
			return null;
//		JDPoint aaa = tri.getA();
//		JDPoint bbb = tri.getB();
//		JDPoint ccc = tri.getC();
		if (id == 0x333) {
			// System.out.println("[INFO] return full triangle, because it is only 1
			// colour!");
			return new JDTriangle(new JDPoint(x1, y1, v1), new JDPoint(x2, y2, v2), new JDPoint(x3, y3, v3));
		}

//		System.out.println("a="+a+"    b="+b+"    c="+c+
//				"    min="+levmin+"    max="+levmax);
		int nan = ((Double.isNaN(v1) ? 0 : 1) << 8) | ((Double.isNaN(v2) ? 0 : 1) << 4) | (Double.isNaN(v3) ? 0 : 1);
		JDPoint abh = tri.getA().fractionTowards(0.5d, tri.getB());
		JDPoint bch = tri.getB().fractionTowards(0.5d, tri.getC());
		JDPoint cah = tri.getC().fractionTowards(0.5d, tri.getA());
//		if(nan==0x001 || nan==0x010 || nan==0x100) {
//			System.out.println("[INFO] return triangle from NAN-code: "+(nan==0x001?"001":nan==0x010?"010":"100"));
//		}
		if (nan == 0x001)
			return new JDTriangle(tri.getC(), cah, bch);
		if (nan == 0x010)
			return new JDTriangle(tri.getB(), bch, abh);
		if (nan == 0x100)
			return new JDTriangle(tri.getA(), abh, cah);

		// rotate for usage of symmetry
		if (nan == 0x101) {
			double tmp = x1;
			x1 = x2;
			x2 = x3;
			x3 = tmp;
			tmp = y1;
			y1 = y2;
			y2 = y3;
			y3 = tmp;
			tmp = v1;
			v1 = v2;
			v2 = v3;
			v3 = tmp;
			nan = 0x011;
		}
		if (nan == 0x110) {
			double tmp = x1;
			x1 = x3;
			x3 = x2;
			x2 = tmp;
			tmp = y1;
			y1 = y3;
			y3 = y2;
			y2 = tmp;
			tmp = v1;
			v1 = v3;
			v3 = v2;
			v2 = tmp;
			nan = 0x011;
		}

		// TODO
		int count = 3;
		double x4 = Double.NaN, y4 = Double.NaN, v4 = Double.NaN;
		double x5 = Double.NaN, y5 = Double.NaN, v5 = Double.NaN;
		if (nan == 0x011) {
			x4 = JPlotMath.lerp(x1, x3, 0.5d);
			y4 = JPlotMath.lerp(y1, y3, 0.5d);
			v4 = v3;
			x1 = JPlotMath.lerp(x1, x2, 0.5d);
			y1 = JPlotMath.lerp(y1, y2, 0.5d);
			v1 = v2;
			count = 4;
			double fl = (lower - v1) / (v3 - v1);
			if (0.0d < fl && fl < 1.0d) {
				if (v1 < lower) {
					x1 = x1 + fl * (x4 - x1);
					y1 = y1 + fl * (y4 - y1);
					v1 = lower;
					x2 = x2 + fl * (x3 - x2);
					y2 = y2 + fl * (y3 - y2);
					v2 = lower;
				} else {
					x4 = x1 + fl * (x4 - x1);
					y4 = y1 + fl * (y4 - y1);
					v4 = lower;
					x3 = x2 + fl * (x3 - x2);
					y3 = y2 + fl * (y3 - y2);
					v3 = lower;
				}
			}
			double fu = (upper - v1) / (v3 - v1);
			if (0.0d < fu && fu < 1.0d) {
				if (v1 > upper) {
					x1 = x1 + fu * (x4 - x1);
					y1 = y1 + fu * (y4 - y1);
					v1 = upper;
					x2 = x2 + fu * (x3 - x2);
					y2 = y2 + fu * (y3 - y2);
					v2 = upper;
				} else {
					x4 = x1 + fu * (x4 - x1);
					y4 = y1 + fu * (y4 - y1);
					v4 = upper;
					x3 = x2 + fu * (x3 - x2);
					y3 = y2 + fu * (y3 - y2);
					v3 = upper;
				}
			}
		} else {
			int low = (v1 < lower ? 1 : 0) | (v2 < lower ? 2 : 0) | (v3 < lower ? 4 : 0);
			switch (low) {
			case 1:
				count = 4;
				double xf12 = (lower - v1) / (v2 - v1), xf13 = (lower - v1) / (v3 - v1);
				x4 = x1 + xf13 * (x3 - x1);
				y4 = y1 + xf13 * (y3 - y1);
				x1 = x1 + xf12 * (x2 - x1);
				y1 = y1 + xf12 * (y2 - y1);
				v4 = v1 + xf13 * (v3 - v1);
				v1 = v1 + xf12 * (v2 - v1);
				break;
			case 2:
				count = 4;
				double xf21 = (lower - v2) / (v1 - v2), xf23 = (lower - v2) / (v3 - v2);
				x4 = x1;
				y4 = y1;
				x1 = x2 + xf21 * (x1 - x2);
				y1 = y2 + xf21 * (y1 - y2);
				x2 = x2 + xf23 * (x3 - x2);
				y2 = y2 + xf23 * (y3 - y2);
				v4 = v1;
				v1 = v2 + xf21 * (v1 - v2);
				v2 = v2 + xf23 * (v3 - v2);
				break;
			case 4:
				count = 4;
				double xf31 = (lower - v3) / (v1 - v3), xf32 = (lower - v3) / (v2 - v3);
				x4 = x3 + xf31 * (x1 - x3);
				y4 = y3 + xf31 * (y1 - y3);
				x3 = x3 + xf32 * (x2 - x3);
				y3 = y3 + xf32 * (y2 - y3);
				v4 = v3 + xf31 * (v1 - v3);
				v3 = v3 + xf32 * (v2 - v3);
				break;
			case 3:
				count = 3;
				double fx32 = (lower - v3) / (v2 - v3), fx31 = (lower - v3) / (v1 - v3);
				x1 = x3 + fx31 * (x1 - x3);
				y1 = y3 + fx31 * (y1 - y3);
				x2 = x3 + fx32 * (x2 - x3);
				y2 = y3 + fx32 * (y2 - y3);
				v1 = v3 + fx31 * (v1 - v3);
				v2 = v3 + fx32 * (v2 - v3);
				break;
			case 5:
				count = 3;
				double fx21 = (lower - v2) / (v1 - v2), fx23 = (lower - v2) / (v3 - v2);
				x1 = x2 + fx21 * (x1 - x2);
				y1 = y2 + fx21 * (y1 - y2);
				x3 = x2 + fx23 * (x3 - x2);
				y3 = y2 + fx23 * (y3 - y2);
				v1 = v2 + fx21 * (v1 - v2);
				v3 = v2 + fx23 * (v3 - v2);
				break;
			case 6:
				count = 3;
				double fx12 = (lower - v1) / (v2 - v1), fx13 = (lower - v1) / (v3 - v1);
				x2 = x1 + fx12 * (x2 - x1);
				y2 = y1 + fx12 * (y2 - y1);
				x3 = x1 + fx13 * (x3 - x1);
				y3 = y1 + fx13 * (y3 - y1);
				v2 = v1 + fx12 * (v2 - v1);
				v3 = v1 + fx13 * (v3 - v1);
				break;
			default:
				count = 3;
				break;
			}
			int high = (v1 > upper ? 1 : 0) | (v2 > upper ? 2 : 0) | (v3 > upper ? 4 : 0);
			if (count > 3)
				high |= (v4 > upper ? 8 : 0);
			if (count == 3) {
				switch (high) {
				case 1:
					count = 4;
					double xf12 = (upper - v1) / (v2 - v1), xf13 = (upper - v1) / (v3 - v1);
					x4 = x1 + xf13 * (x3 - x1);
					y4 = y1 + xf13 * (y3 - y1);
					x1 = x1 + xf12 * (x2 - x1);
					y1 = y1 + xf12 * (y2 - y1);
					v4 = v1 + xf13 * (v3 - v1);
					v1 = v1 + xf12 * (v2 - v1);
					break;
				case 2:
					count = 4;
					double xf21 = (upper - v2) / (v1 - v2), xf23 = (upper - v2) / (v3 - v2);
					x4 = x1;
					y4 = y1;
					x1 = x2 + xf21 * (x1 - x2);
					y1 = y2 + xf21 * (y1 - y2);
					x2 = x2 + xf23 * (x3 - x2);
					y2 = y2 + xf23 * (y3 - y2);
					v4 = v1;
					v1 = v2 + xf21 * (v1 - v2);
					v2 = v2 + xf23 * (v3 - v2);
					break;
				case 4:
					count = 4;
					double xf31 = (upper - v3) / (v1 - v3), xf32 = (upper - v3) / (v2 - v3);
					x4 = x3 + xf31 * (x1 - x3);
					y4 = y3 + xf31 * (y1 - y3);
					x3 = x3 + xf32 * (x2 - x3);
					y3 = y3 + xf32 * (y2 - y3);
					v4 = v3 + xf31 * (v1 - v3);
					v3 = v3 + xf32 * (v2 - v3);
					break;
				case 3:
					count = 3;
					double fx32 = (upper - v3) / (v2 - v3), fx31 = (upper - v3) / (v1 - v3);
					x1 = x3 + fx31 * (x1 - x3);
					y1 = y3 + fx31 * (y1 - y3);
					x2 = x3 + fx32 * (x2 - x3);
					y2 = y3 + fx32 * (y2 - y3);
					v1 = v3 + fx31 * (v1 - v3);
					v2 = v3 + fx32 * (v2 - v3);
					break;
				case 5:
					count = 3;
					double fx21 = (upper - v2) / (v1 - v2), fx23 = (upper - v2) / (v3 - v2);
					x1 = x2 + fx21 * (x1 - x2);
					y1 = y2 + fx21 * (y1 - y2);
					x3 = x2 + fx23 * (x3 - x2);
					y3 = y2 + fx23 * (y3 - y2);
					v1 = v2 + fx21 * (v1 - v2);
					v3 = v2 + fx23 * (v3 - v2);
					break;
				case 6:
					count = 3;
					double fx12 = (upper - v1) / (v2 - v1), fx13 = (upper - v1) / (v3 - v1);
					x2 = x1 + fx12 * (x2 - x1);
					y2 = y1 + fx12 * (y2 - y1);
					x3 = x1 + fx13 * (x3 - x1);
					y3 = y1 + fx13 * (y3 - y1);
					v2 = v1 + fx12 * (v2 - v1);
					v3 = v1 + fx13 * (v3 - v1);
					break;
				default:
					count = 3;
					break;
				}
			} else if (count == 4) {
				switch (high) {
				case 1:
					count = 5;
					double xf12 = (upper - v1) / (v2 - v1), xf14 = (upper - v1) / (v4 - v1);
					x5 = x1 + xf14 * (x4 - x1);
					y5 = y1 + xf14 * (y4 - y1);
					x1 = x1 + xf12 * (x2 - x1);
					y1 = y1 + xf12 * (y2 - y1);
					v5 = v1 + xf14 * (v4 - v1);
					v1 = v1 + xf12 * (v2 - v1);
					break;
				case 2:
					count = 5;
					double xf21 = (upper - v2) / (v1 - v2), xf23 = (upper - v2) / (v3 - v2);
					x5 = x1;
					y5 = y1;
					x1 = x2 + xf21 * (x1 - x2);
					y1 = y2 + xf21 * (y1 - y2);
					x2 = x2 + xf23 * (x3 - x2);
					y2 = y2 + xf23 * (y3 - y2);
					v5 = v1;
					v1 = x2 + xf21 * (v1 - v2);
					v2 = v2 + xf23 * (v3 - v2);
					break;
				case 4:
					count = 5;
					double xf32 = (upper - v3) / (v2 - v3), xf34 = (upper - v3) / (v4 - v3);
					x5 = x4;
					y5 = y4;
					x4 = x3 + xf34 * (x4 - x3);
					y4 = y3 + xf34 * (y4 - y3);
					x3 = x3 + xf32 * (x2 - x3);
					y3 = y3 + xf32 * (y2 - y3);
					v5 = v4;
					v4 = v3 + xf34 * (v4 - v3);
					v3 = v3 + xf32 * (v2 - v3);
					break;
				case 8:
					count = 5;
					double xf41 = (upper - v4) / (v1 - v4), xf43 = (upper - v4) / (v3 - v4);
					x5 = x4 + xf41 * (x1 - x4);
					y5 = y4 + xf41 * (y1 - y4);
					x4 = x4 + xf43 * (x3 - x4);
					y4 = y4 + xf43 * (y3 - y4);
					v5 = v4 + xf41 * (v1 - v4);
					v4 = v4 + xf43 * (v3 - v4);
					break;
				case 3:
					count = 4;
					double fx14 = (upper - v1) / (v4 - v1), fx23 = (upper - v2) / (v3 - v2);
					x1 = x1 + fx14 * (x4 - x1);
					y1 = y1 + fx14 * (y4 - y1);
					x2 = x2 + fx23 * (x3 - x2);
					y2 = y2 + fx23 * (y3 - y2);
					v1 = v1 + fx14 * (v4 - v1);
					v2 = v2 + fx23 * (v3 - v2);
					break;
				case 6:
					count = 4;
					double fx21 = (upper - v2) / (v1 - v2), fx34 = (upper - v3) / (v4 - v3);
					x2 = x2 + fx21 * (x1 - x2);
					y2 = y2 + fx21 * (y1 - y2);
					x3 = x3 + fx34 * (x4 - x3);
					y3 = y3 + fx34 * (y4 - y3);
					v2 = v2 + fx21 * (v1 - v2);
					v3 = v3 + fx34 * (v4 - v3);
					break;
				case 9:
					count = 4;
					double fx12 = (upper - v1) / (v2 - v1), fx43 = (upper - v4) / (v3 - v4);
					x1 = x1 + fx12 * (x2 - x1);
					y1 = y1 + fx12 * (y2 - y1);
					x4 = x4 + fx43 * (x3 - x4);
					y4 = y4 + fx43 * (y3 - y4);
					v1 = v1 + fx12 * (v2 - v1);
					v4 = v4 + fx43 * (v3 - v4);
					break;
				case 12:
					count = 4;
					double fx41 = (upper - v4) / (v1 - v4), fx32 = (upper - v3) / (v2 - v3);
					x3 = x3 + fx32 * (x2 - x3);
					y3 = y3 + fx32 * (y2 - y3);
					x4 = x4 + fx41 * (x1 - x4);
					y4 = y4 + fx41 * (y1 - y4);
					v3 = v3 + fx32 * (v2 - v3);
					v4 = v4 + fx41 * (v1 - v4);
					break;
				default:
					count = 4;
					break;
				}
			}
		}
		switch (count) {
		case 3:
			return new JDPolygon(new JDPoint(x1, y1, v1), new JDPoint(x2, y2, v2), new JDPoint(x3, y3, v3));
		case 4:
			return new JDPolygon(new JDPoint(x1, y1, v1), new JDPoint(x2, y2, v2), new JDPoint(x3, y3, v3),
					new JDPoint(x4, y4, v4));
		case 5:
			return new JDPolygon(new JDPoint(x1, y1, v1), new JDPoint(x2, y2, v2), new JDPoint(x3, y3, v3),
					new JDPoint(x4, y4, v4), new JDPoint(x5, y5, v5));
		default:
			return null;
		}
	}

	private void addTriangle2polygonList(List<JDPolygon> list, JDTriangle tri, double tol) {
		if (list.isEmpty()) {
			list.add(tri.toPolygon());
		} else {
			boolean failed = true;
			for (int l = 0; l < list.size() && failed; l++) {
				if(list.get(l).c.length>999) continue;
				JDPolygon temp = GeometryTools.union(list.get(l), tri.toPolygon(), tol);
				if(temp==null) continue;
				list.set(l, temp);
				failed = false;
			}
			if (failed)
				list.add(tri.toPolygon());
		}
	}

	private void addPolygon2polygonList(List<JDPolygon> list, JDPolygon poly, double tol) {
		if (list.isEmpty()) {
			list.add(poly.copy());
		} else {
			boolean failed = true;
			for (int l = 0; l < list.size() && failed; l++) {
				if(list.get(l).c.length>999) continue;
				JDPolygon temp = GeometryTools.union(list.get(l), poly, tol);
				if(temp==null) continue;
				list.set(l, temp);
				failed = false;
			}
			if (failed)
				list.add(poly.copy());
		}
	}
}
