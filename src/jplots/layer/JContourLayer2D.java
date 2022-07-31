package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import jplots.JAxis;
import jplots.JPlot;
import jplots.colour.JColourtable;
import jplots.maths.AffineBuilder;
import jplots.maths.JDLine;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import jplots.maths.JDQuad;
import jplots.maths.JDTriangle;
import jplots.maths.JPlotMath;
import jplots.shapes.JDGeometry;
import jplots.shapes.JGroupShape;
import jplots.shapes.JImageShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JPolygonShape;
import jplots.transform.JProjection;
import processing.core.PConstants;
import processing.core.PGraphics;

public class JContourLayer2D extends JPlotsLayer {

//	private double EPSILON = Math.pow(2, -52);
	private boolean isFilled, pixelFilling;
	private double minZ, maxZ, Xin, Xax, Yin, Yax;
	private double[] contourIntervals;
	private String[] contourStyle;
	private JDPoint[][] corners;
	private JDTriangle[] triangles;
	private JDGeometry[] fillings;
	
	public JContourLayer2D(float[] x, float[] y, float[][] z, float zmin, float zmax, int nIntervals, float[] zIntervals) {
		boolean is_valid = (x!=null && y!=null && z!=null);
		if(is_valid) {
			is_valid = y.length == z.length;
			for(int j=0; j<y.length && is_valid; j++)
				is_valid = (x.length == z[j].length);
		}
		if(!is_valid) {
			null_init();
			return;
		}
		corners = new JDPoint[y.length][x.length];
		for(int j=0; j<y.length; j++)
			for(int i=0; i<x.length; i++)
				corners[j][i] = new JDPoint(x[i], y[j], z[j][i]);
		if(zIntervals==null)
			contourIntervals = JPlotMath.linspace((double)zmin, (double)zmax, nIntervals);
		else
			contourIntervals = JPlotMath.toDoubleArray1D(zIntervals);
		init();
	}
	public JContourLayer2D(float[][] x, float[][] y, float[][] z, float zmin, float zmax, int nIntervals, float[] zIntervals) {
		boolean is_valid = (x!=null && y!=null && z!=null);
		if(is_valid) {
			is_valid = (x.length == z.length && y.length == z.length);
			int inner_length = 0; if(is_valid && x.length>0) inner_length = x[0].length;
			for(int j=0; j<x.length && is_valid; j++)
				is_valid = (x[j].length == z[j].length && y[j].length == z[j].length && inner_length == z[j].length);
		}
		if(!is_valid) {
			null_init();
			return;
		}
		corners = new JDPoint[y.length][x.length];
		for(int j=0; j<x.length; j++)
			for(int i=0; i<x[j].length; i++)
				corners[j][i] = new JDPoint(x[j][i], y[j][i], z[j][i]);
		if(zIntervals==null)
			contourIntervals = JPlotMath.linspace((double)zmin, (double)zmax, nIntervals);
		else
			contourIntervals = JPlotMath.toDoubleArray1D(zIntervals);
		init();
	}
	public JContourLayer2D(double[] x, double[] y, double[][] z, double zmin, double zmax, int nIntervals, double[] zIntervals) {
		boolean is_valid = (x!=null && y!=null && z!=null);
		if(is_valid) {
			is_valid = y.length == z.length;
			for(int j=0; j<y.length && is_valid; j++)
				is_valid = (x.length == z[j].length);
		}
		if(!is_valid) {
			null_init();
			return;
		}
		corners = new JDPoint[y.length][x.length];
		for(int j=0; j<y.length; j++)
			for(int i=0; i<x.length; i++)
				corners[j][i] = new JDPoint(x[i], y[j], z[j][i]);
		if(zIntervals==null)
			contourIntervals = JPlotMath.linspace(zmin, zmax, nIntervals);
		else
			contourIntervals = zIntervals;
		init();
	}
	public JContourLayer2D(double[][] x, double[][] y, double[][] z, double zmin, double zmax, int nIntervals, double[] zIntervals) {
		boolean is_valid = (x!=null && y!=null && z!=null);
		if(is_valid) {
			is_valid = (x.length == z.length && y.length == z.length);
			int inner_length = 0; if(is_valid && x.length>0) inner_length = x[0].length;
			for(int j=0; j<x.length && is_valid; j++)
				is_valid = (x[j].length == z[j].length && y[j].length == z[j].length && inner_length == z[j].length);
		}
		if(!is_valid) {
			null_init();
			return;
		}
		corners = new JDPoint[y.length][x.length];
		for(int j=0; j<x.length; j++)
			for(int i=0; i<x[j].length; i++)
				corners[j][i] = new JDPoint(x[j][i], y[j][i], z[j][i]);
		if(zIntervals==null)
			contourIntervals = JPlotMath.linspace(zmin, zmax, nIntervals);
		else
			contourIntervals = zIntervals;
		init();
	}
	
	/*
	public JContourLayer2D(float[] x, float[] y, float[][] z, float zmin, float zmax, int nintervals, JColourtable ct,
			float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		float zin = Float.isNaN(zmin) ? JPlotMath.fmin(z) : zmin;
		float zax = Float.isNaN(zmax) ? JPlotMath.fmax(z) : zmax;
		float[] cntIntervals = new float[nintervals + 1];
		for (int k = 0; k <= nintervals; k++)
			cntIntervals[k] = zin + k * (zax - zin) / nintervals;
		input2d = false;
		JContourLayerFloat(x, y, null, null, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer2D(float[] x, float[] y, float[][] z, float[] intervals, JColourtable ct, float stroke_weight,
			boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = false;
		JContourLayerFloat(x, y, null, null, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer2D(double[] x, double[] y, double[][] z, double zmin, double zmax, int nintervals,
			JColourtable ct, double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		double zin = Double.isNaN(zmin) ? JPlotMath.dmin(z) : zmin;
		double zax = Double.isNaN(zmax) ? JPlotMath.dmax(z) : zmax;
		double[] cntIntervals = new double[nintervals + 1];
		for (int k = 0; k <= nintervals; k++)
			cntIntervals[k] = zin + k * (zax - zin) / nintervals;
		input2d = false;
		JContourLayerDouble(x, y, null, null, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer2D(double[] x, double[] y, double[][] z, double[] intervals, JColourtable ct,
			double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = false;
		JContourLayerDouble(x, y, null, null, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer2D(float[][] x, float[][] y, float[][] z, float zmin, float zmax, int nintervals, JColourtable ct,
			float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		float zin = Float.isNaN(zmin) ? JPlotMath.fmin(z) : zmin;
		float zax = Float.isNaN(zmax) ? JPlotMath.fmax(z) : zmax;
		float[] cntIntervals = new float[nintervals + 1];
		for (int k = 0; k <= nintervals; k++)
			cntIntervals[k] = zin + k * (zax - zin) / nintervals;
		input2d = true;
		JContourLayerFloat(null, null, x, y, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer2D(float[][] x, float[][] y, float[][] z, float[] intervals, JColourtable ct, float stroke_weight,
			boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = true;
		JContourLayerFloat(null, null, x, y, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer2D(double[][] x, double[][] y, double[][] z, double zmin, double zmax, int nintervals,
			JColourtable ct, double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		double zin = Double.isNaN(zmin) ? JPlotMath.dmin(z) : zmin;
		double zax = Double.isNaN(zmax) ? JPlotMath.dmax(z) : zmax;
		double[] cntIntervals = new double[nintervals + 1];
		for (int k = 0; k <= nintervals; k++)
			cntIntervals[k] = zin + k * (zax - zin) / nintervals;
		input2d = true;
		JContourLayerDouble(null, null, x, y, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer2D(double[][] x, double[][] y, double[][] z, double[] intervals, JColourtable ct,
			double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = true;
		JContourLayerDouble(null, null, x, y, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}

	/*
	public JContourLayer(float[] x1, float[] y1, float[][] x2, float[][] y2, float[][] z, float zmin, float zmax, int nintervals, float[] zintervals,
			JColourtable ct, float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		
		
		
	}
	*/
	
	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		if(corners.length==0) return;
		if(corners[0].length==0) return;
		JDPoint[][] cnt2 = new JDPoint[corners.length][corners[0].length];
		for(int j=0; j<corners.length; j++)
			for(int i=0; i<corners[j].length; i++)
				cnt2[j][i] = corners[j][i].copy();
		if(ax.isGeoAxis() && inputProj!=null) {
			for(int j=0; j<corners.length; j++)
				for(int i=0; i<corners[j].length; i++) {
					double[] xy = inputProj.fromPROJtoLATLON(cnt2[j][i].x, cnt2[j][i].y, false);
					xy = ax.getGeoProjection().fromLATLONtoPROJ(xy[0], xy[1], false);
					cnt2[j][i].x = xy[0];
					cnt2[j][i].y = xy[1];
				}
		}
		int[] p = ax.getSize();
		Xin = ax.isXlogAxis() ? Math.log10(minX) : minX;
		Xax = ax.isXlogAxis() ? Math.log10(maxX) : maxX;
		Yin = ax.isYlogAxis() ? Math.log10(minY) : minY;
		Yax = ax.isYlogAxis() ? Math.log10(maxY) : maxY;
		double xs = p[2] / (Xax - Xin), ys = p[3] / (Yax - Yin);
		// double tol = Math.max(Math.abs(maxX-minX), Math.abs(maxY-minY)) * 1.0e-12d;
		AffineBuilder affine = new AffineBuilder().scale(invertAxisX ? -1d : 1d, invertAxisY ? 1d : -1d)
				.translate(invertAxisX ? maxX : -minX, invertAxisY ? -minY : maxY).scale(xs, ys).translate(p[0], p[1]);
		minZ = Double.POSITIVE_INFINITY;
		maxZ = Double.NEGATIVE_INFINITY;
		for(int j=0; j<corners.length; j++)
			for(int i=0; i<corners[j].length; i++) {
				cnt2[j][i].affine(affine.getMatrix());
				if(cnt2[j][i].value<minZ) minZ = cnt2[j][i].value;
				if(cnt2[j][i].value>maxZ) maxZ = cnt2[j][i].value;
			}
		if(contourIntervals[0]<minZ) minZ = contourIntervals[0];
		if(contourIntervals[contourIntervals.length-1]>maxZ) maxZ = contourIntervals[contourIntervals.length-1];
		
		// step 1: create filling between contours if wished
		if (isFilled) {
			if (pixelFilling) {
				if (ax.getPlot().isDebug())
					System.out.println("[DEBUG] JContourLayer2D: 1] contour filling pixelwise ...");
				fillPixelByPixel(p, ax, xs, ys, s, cnt2);
			} else {
				if (ax.getPlot().isDebug())
					System.out.println("[DEBUG] JContourLayer2D: 1] contour filling vectorwise ...");
				try {
					fillVectorByVector(ax.getGeoProjection(), p, ax, xs, ys, s);
				} catch (Exception e) {
					e.printStackTrace();
					throw new RuntimeException(e);
				}
			}
		} else {
			if (ax.getPlot().isDebug())
				System.out.println("[DEBUG] JContourLayer2D: 1] no filling ...");
		}

		// step 2: add contours to plot
		if (drawLines) {
			if (ax.getPlot().isDebug())
				System.out.println("[DEBUG] JContourLayer2D: 2] draw contour line segments with " + lw + "px line width...");
			drawContourLines(p, ax, xs, ys, s, cnt2);
		} else {
			if (ax.getPlot().isDebug())
				System.out.println("[DEBUG] JContourLayer2D: 2] contours itself will not be drawn ...");
		}
		
		if(ax.getPlot().isDebug()) {
			drawTriangleBorders(p,ax,xs,xs,s);
		}
	}

	//* **************************************** *
	//* ********** GETTER AND SETTER  ********** *
	//* **************************************** *
	
	public void setFilled(boolean b) {
		isFilled = b;
	}
	public void setPixelFilling(boolean b) {
		pixelFilling = b;
	}
	
	@Override
	public void setLineColour(int _lc) {
		lc = _lc;
		for(int i=0; i<lcs.length; i++)
			lcs[i] = _lc;
	}
	
	public void setStyles(String[] nstyle) {
		contourStyle = nstyle;
	}
	@Override
	public void setStyle(String lst) {
		for(int i=0; i<contourStyle.length; i++)
			contourStyle[i] = lst;
		ls = lst;
	}
	
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
	
	public void drawTriangleBorders(int[] p, JAxis ax, double xs, double ys, JGroupShape s) {
		
	}
	
	
	
	//* **************************************** *
	//* ********** PRIVATE METHODS    ********** *
	//* **************************************** *
	
	private void null_init() {
		corners = new JDPoint[0][0];
		contourIntervals = new double[] { 0d, 1d };
		//TODO null init
	}
	private void init() {
		drawLines = true;
		isFilled = false;
		pixelFilling = false;
		minX = Double.POSITIVE_INFINITY;
		maxX = Double.NEGATIVE_INFINITY;
		minY = Double.POSITIVE_INFINITY;
		maxY = Double.NEGATIVE_INFINITY;
		for(int j=0; j<corners.length; j++)
			for(int i=0; i<corners[j].length; i++) {
				if(corners[j][i].x<minX) minX = corners[j][i].x;
				if(corners[j][i].x>maxX) maxX = corners[j][i].x;
				if(corners[j][i].y<minY) minY = corners[j][i].y;
				if(corners[j][i].y>maxY) maxY = corners[j][i].y;
			}
		contourStyle = new String[contourIntervals.length];
		lcs = new int[contourIntervals.length];
		for(int i=0; i<contourIntervals.length; i++) {
			contourStyle[i] = "";
			lcs[i] = lc;
		}
	}
	
	private void fillVectorByVector(JProjection outproj, int[] p, JAxis ax, double xs, double ys, JGroupShape s) {
		List<int[]> openset = new ArrayList<>(),
					closedset = new ArrayList<>();
		double[] clev2 = new double[2+contourIntervals.length];
		clev2[0] = Double.NEGATIVE_INFINITY;
		for(int i=0; i<contourIntervals.length; i++) clev2[i+1] = contourIntervals[i];
		clev2[contourIntervals.length+1] = Double.POSITIVE_INFINITY;
		for(int i=1; i<clev2.length; i++) {
			double lower = clev2[i-1];
			double upper = clev2[i];
			
		}
		
		
		
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
	
	private void fillPixelByPixel(int[] p, JAxis ax, double xs, double ys, JGroupShape s, JDPoint[][] points) {
		// double us = srcImg.width/(srcExt[2]-srcExt[0]), vs =
		// srcImg.height/(srcExt[3]-srcExt[1]);
		if(colourtable==null) {
			System.err.println("No colourtable set!");
			return;
		}
		if (img == null) {
			img = ax.getPlot().getApplet().createImage(p[2], p[3], PConstants.ARGB);
		} else if (img.width != p[2] || img.height != p[3]) {
			img = ax.getPlot().getApplet().createImage(p[2], p[3], PConstants.ARGB);
		}
		img.loadPixels();
		for(int j=0; j+1<points.length; j++) {
			for(int i=0; i+1<points[j].length; i++) {
				JDQuad qu = new JDQuad(points[j+1][i], points[j+1][i+1], points[j][i+1], points[j][i]);
				double	txi = Math.min(Math.min(qu.x[0],qu.x[1]), Math.min(qu.x[2], qu.x[3])),
						txa = Math.max(Math.max(qu.x[0],qu.x[1]), Math.max(qu.x[2], qu.x[3]));
				double	tyi = Math.min(Math.min(qu.y[0],qu.y[1]), Math.min(qu.y[2], qu.y[3])),
						tya = Math.max(Math.max(qu.y[0],qu.y[1]), Math.max(qu.y[2], qu.y[3]));
				int ixs =  Math.max((int) txi    - (txi < 0 ? 1 : 0), p[0]),
					ixe = -Math.max((int) (-txa) - (txa > 0 ? 1 : 0), 1 - p[0] - p[2]);
				int iys =  Math.max((int) tyi    - (tyi < 0 ? 1 : 0), p[1]),
					iye = -Math.max((int) (-tya) - (tya > 0 ? 1 : 0), 1 - p[1] - p[3]);
				if ((ixe < ixs) || (iye < iys)) continue;
				for (int v = iys; v <= iye; v++) {
					for (int u = ixs; u <= ixe; u++) {
						JDPoint ij = new JDPoint(u+0.5d,v+0.5d);
						if(!qu.contains(ij)) continue;
						img.pixels[(v-p[1]) * p[2] + u-p[0]] = getColor(qu.valueAt(ij));
					}
				}
			}
		}
		img.updatePixels();
		s.addChild(new JImageShape(img, p[0], p[1], p[2], p[3]));
	}
	
	private void drawContourLines(int[] p, JAxis ax, double xs, double ys, JGroupShape s, JDPoint[][] points) {
		int jlen = points.length-1;
		int ilen = points[0].length-1;
		int[][] visits = new int[jlen][ilen];
		List<JDPoint> path = new ArrayList<>();
		List<JDLine> contours = new ArrayList<>();
		for(int l=0; l<contourIntervals.length; l++) {
			contours.clear();
			for(int j=0; j<jlen; j++)
				for(int i=0; i<ilen; i++)
					visits[j][i] = 0;
			double lev = contourIntervals[l];
			
			for(int j=0; j<jlen; j++) for(int i=0; i<ilen; i++) {
				if(visits[j][i]>1) continue;
				boolean hit_saddle = false;
				int su = i, sv = j;
				int pu = i, pv = j;
				int u = i, v = j;
				int nu = -1, nv = -1;
				path.clear();
				while(su!=nu || sv!=nv) {
					if(visits[v][u]==2) {
						if(path.isEmpty()) break;
						JDLine nl = new JDLine(path.toArray(new JDPoint[0]));
						boolean success = false;
						for(JDLine ol: contours) {
							success = ol.join(nl, 0.01d);
							if(success) break;
						}
						if(success) {
							path.clear();
							break;
						}
					}
					nu = -1; nv = -1;
					JDPoint a = points[ v ][ u ];
					JDPoint b = points[ v ][u+1];
					JDPoint c = points[v+1][u+1];
					JDPoint d = points[v+1][ u ];
					int code =	(Double.isNaN(a.value) ? 0x0007 : a.value<lev ? 0x0000 : 0x0001) |
								(Double.isNaN(b.value) ? 0x0070 : b.value<lev ? 0x0000 : 0x0010) |
								(Double.isNaN(c.value) ? 0x0700 : c.value<lev ? 0x0000 : 0x0100) |
								(Double.isNaN(d.value) ? 0x7000 : d.value<lev ? 0x0000 : 0x1000);
					double ab = (lev-a.value)/(b.value-a.value);
					double ac = (lev-a.value)/(c.value-a.value);
					double ad = (lev-a.value)/(d.value-a.value);
					double ba = (lev-b.value)/(a.value-b.value);
					double bc = (lev-b.value)/(c.value-b.value);
					double bd = (lev-b.value)/(d.value-b.value);
					double ca = (lev-c.value)/(a.value-c.value);
					double cb = (lev-c.value)/(b.value-c.value);
					double cd = (lev-c.value)/(d.value-c.value);
					double da = (lev-d.value)/(a.value-d.value);
					double db = (lev-d.value)/(b.value-d.value);
					double dc = (lev-d.value)/(c.value-d.value);
					JDQuad q = new JDQuad(a,b,c,d);
					switch(code) {
						default: break;
						case 0x0000: break;
						case 0x0001:
							//    1_/   0
							//    /      
							//    0     0
							addCurve(path, q, -1d, 2*ad-1d, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1; break;
						case 0x0007: visits[v][u] = 2; break;
						case 0x0010:
							//    0   \_1
							//          \
							//    0     0
							addCurve(path, q, 1d-2*ba, -1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v; break;
						case 0x0011:
							//    1     1
							//    -------
							//    0     0
							addCurve(path, q, -1d, 2*ad-1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v; break;
						case 0x0017:
							//    N __  1
							//        '--
							//    0     0
							JDPoint p0017 = b.fractionTowards((lev-b.value)/(d.value-b.value), d);
							path.add(p0017.fractionTowards(0.5d, a));
							path.add(b.fractionTowards((lev-b.value)/(c.value-b.value), c)); nu = u+1; nv = v; break;
						case 0x0070: break;
						case 0x0071:
							//    1  __ N
							//    --'    
							//    0     0
							JDPoint p0071 = a.fractionTowards((lev-a.value)/(c.value-a.value), c);
							if(path.isEmpty()) path.add(a.fractionTowards((lev-a.value)/(d.value-a.value), d));
							path.add(p0071.fractionTowards(0.5d, b)); break;
						case 0x0077: break;
						case 0x0100:
							//    0     0
							//         _/
							//    0   / 1
							//show(a,b,c,d);
							addCurve(path, q, 1d, 1d-2*cb, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1; break;
						case 0x0101:
							// du = 0.25d*((value[1]-value[0])*(1-v) + (value[2]-value[3])*(1+v)) = 0
							// dv = 0.25d*((value[3]-value[0])*(1-u) + (value[2]-value[1])*(1+u)) = 0
							/*
							 * (v1-v0)*(1-y)+(v2-v3)*(1+y) = 0
							 * v1-v0+v2-v3 = (v1-v0-v2+v3)*y
							 * 
							 * (v3-v0)*(1-x) + (v2-v1)*(1+x) = 0
							 * v3-v0+v2-v1 = (v3-v0-v2+v1)*x
							 * 
							 */
							hit_saddle = true;
							double u01 = (d.value-a.value+c.value-b.value)/(d.value-a.value-c.value+b.value);
							double v01 = (b.value-a.value+c.value-d.value)/(b.value-a.value-c.value+d.value);
							double m01 = 0.25d*(a.value*(1-u01)*(1-v01)+b.value*(1+u01)*(1-v01)+c.value*(1+u01)*(1+v01)+d.value*(1-u01)*(1+v01));
							boolean takePos01 = visits[v][u]<0;
							if(visits[v][u]==0) takePos01 = (pu>u);
							if(m01<=lev) {
								//    1_/   0
								// -- /  0 _/ +
								//    0   / 1
								if(takePos01) {
									addCurve(path, q, 1d, 1d-2*cb, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1;
									visits[v][u] = visits[v][u]==0 ? 1 : 2;
								} else {
									addCurve(path, q, -1d, 2*ad-1d, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1;
									visits[v][u] = visits[v][u]==0 ? -1 : 2;
								} break;
							} else {
								//    1   \_0
								// -- \_ 1  \ +
								//    0 \   1
								if(takePos01) {
									addCurve(path, q, 1d, 1d-2*cb, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1;
									visits[v][u] = visits[v][u]==0 ? 1 : 2;
								} else {
									addCurve(path, q, -1d, 2*ad-1d, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1;
									visits[v][u] = visits[v][u]==0 ? -1 : 2;
								} break;
							}
						case 0x0107:
							//    N     0
							//         _/
							//    0   / 1
							if(path.isEmpty()) path.add(c.fractionTowards((lev-c.value)/(b.value-c.value), b));
							path.add(c.fractionTowards((lev-c.value)/(d.value-c.value), d)); nu = u; nv = v+1; break;
						case 0x0110:
							//    0  |  1
							//       |   
							//    0  |  1
							addCurve(path, q, 1d-2*ba, -1d, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1; break;
						case 0x0111:
							//    1     1
							//    \_     
							//    0 \   1
							addCurve(path, q, -1d, 2*ad-1d, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1; break;
						case 0x0117:
							//    N     1
							//      \    
							//    0  \  1
							JDPoint p0117 = b.fractionTowards((lev-b.value)/(d.value-b.value), d);
							path.add(p0117.fractionTowards(0.5d,a));
							path.add(c.fractionTowards((lev-c.value)/(d.value-c.value), d)); nu = u; nv = v+1; break;
						case 0x0170:
							//    0     N
							//        /  
							//    0  /  1
							JDPoint p0170 = c.fractionTowards((lev-c.value)/(a.value-c.value), a);
							path.add(p0170.fractionTowards(0.5d,b));
							path.add(c.fractionTowards((lev-c.value)/(d.value-c.value), d)); nu = u; nv = v+1; break;
						case 0x0171:
							//    1     N
							//    \_     
							//    0 \   1
							if(path.isEmpty()) path.add(a.fractionTowards((lev-a.value)/(d.value-a.value), d));
							path.add(c.fractionTowards((lev-c.value)/(d.value-c.value), d)); nu = u; nv = v+1; break;
						case 0x0177:
							//    N     N
							//       |   
							//    0  |  1
							JDPoint p0177 = c.fractionTowards((lev-c.value)/(d.value-c.value), d);
							JDPoint p0177n = b.fractionTowards(0.5d, a);
							path.add(p0177.fractionTowards(0.5d, p0177n)); path.add(p0177); nu = u; nv = v+1; break;
						case 0x0700: break;
						case 0x0701:
							//    1_/   0
							//    /      
							//    0     N
							if(path.isEmpty()) path.add(a.fractionTowards((lev-a.value)/(d.value-a.value), d));
							path.add(a.fractionTowards((lev-a.value)/(b.value-a.value), b)); nu = u; nv = v-1; break;
						case 0x0707: break;
						case 0x0710:
							//    0  \  1
							//        \  
							//    0     N
							JDPoint p0710 = b.fractionTowards((lev-b.value)/(d.value-b.value), d);
							if(path.isEmpty()) path.add(b.fractionTowards((lev-b.value)/(a.value-b.value), a));
							path.add(p0710.fractionTowards(0.5d, c)); break;
						case 0x0711:
							//    1     1
							//    --.__  
							//    0     N
							JDPoint p0711 = b.fractionTowards((lev-b.value)/(d.value-b.value), d);
							if(path.isEmpty()) path.add(a.fractionTowards((lev-a.value)/(d.value-a.value), d));
							path.add(p0711.fractionTowards(0.5d, c)); break;
						case 0x0717:
							//    N _   1
							//       \_  
							//    0     N
							JDPoint p0717 = b.fractionTowards((lev-b.value)/(d.value-b.value), d);
							path.add(p0717.fractionTowards(0.5d, a)); path.add(p0717.fractionTowards(0.5d, c)); break;
						case 0x0770: break;
						case 0x0771:
							//    1     N
							//    ----   
							//    0     N
							JDPoint p0771 = a.fractionTowards((lev-a.value)/(d.value-a.value), d);
							JDPoint p0771n = b.fractionTowards(0.5d, c);
							if(path.isEmpty()) path.add(p0771); path.add(p0771.fractionTowards(0.5d, p0771n)); break;
						case 0x0777: break;
						case 0x1000:
							//    0     0
							//    \_     
							//    1 \   0
							addCurve(path, q, 2*dc-1d, 1d, -1d, 1d-2*da, lev, true); nu = u-1; nv = v; break;
						case 0x1001:
							//    1  |  0
							//       |   
							//    1  |  0
							addCurve(path, q, 2*dc-1d, 1d, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1; break;
						case 0x1007:
							//    N     0
							//      \    
							//    1  \  0
							JDPoint p1007 = d.fractionTowards((lev-d.value)/(b.value-d.value), b);
							if(path.isEmpty()) path.add(d.fractionTowards((lev-d.value)/(c.value-d.value), c));
							path.add(p1007.fractionTowards(0.5d, a)); break;
						case 0x1010:
							hit_saddle = true;
							double u10 = (d.value-a.value+c.value-b.value)/(d.value-a.value-c.value+b.value);
							double v10 = (b.value-a.value+c.value-d.value)/(b.value-a.value-c.value+d.value);
							double m10 = 0.25d*(a.value*(1-u10)*(1-v10)+b.value*(1+u10)*(1-v10)+c.value*(1+u10)*(1+v10)+d.value*(1-u10)*(1+v10));
							boolean takePos10 = visits[v][u]<0;
							if(m10<=lev) {
								//    0   \_1
								// -- \_ 0  \ +
								//    1 \   0
								if(visits[v][u]==0) takePos10 = (pv<v);
								if(takePos10) {
									addCurve(path, q, 1d-2*ba, -1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v;
									visits[v][u] = visits[v][u]==0 ? 1 : 2;
								} else {
									addCurve(path, q, 2*dc-1d, 1d, -1d, 1d-2*da, lev, true); nu = u-1; nv = v;
									visits[v][u] = visits[v][u]==0 ? -1 : 2;
								} break;
							} else {
								//    0_/   1
								// -- /  1 _/ +
								//    1   / 0
								if(visits[v][u]==0) takePos10 = (pv>v);
								if(takePos10) {
									addCurve(path, q, 2*dc-1d, 1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v;
									visits[v][u] = visits[v][u]==0 ? 1 : 2;
								} else {
									addCurve(path, q, 1d-2*ba, -1d, -1d, 1d-2*da, lev, true); nu = u-1; nv = v;
									visits[v][u] = visits[v][u]==0 ? -1 : 2;
								} break;
							}
						case 0x1011:
							//    1     1
							//         _/
							//    1   / 0
							addCurve(path, q, 2*dc-1d, 1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v; break;
						case 0x1017:
							//    N     1
							//         _/
							//    1   / 0
							if(path.isEmpty()) path.add(d.fractionTowards((lev-d.value)/(c.value-d.value), c));
							path.add(b.fractionTowards((lev-b.value)/(c.value-b.value), c)); nu = u+1; nv = v; break;
						case 0x1070:
							//    0     N
							//    \_     
							//    1 \   0
							if(path.isEmpty()) path.add(d.fractionTowards((lev-d.value)/(c.value-d.value), c));
							path.add(d.fractionTowards((lev-d.value)/(a.value-d.value), a)); nu = u-1; nv = v; break;
						case 0x1071:
							//    1     N
							//        /  
							//    1  /  0
							JDPoint p1071 = a.fractionTowards((lev-a.value)/(c.value-a.value), c);
							if(path.isEmpty()) path.add(d.fractionTowards((lev-d.value)/(c.value-d.value), c));
							path.add(p1071.fractionTowards(0.5d, b)); break;
						case 0x1077:
							//    N     N
							//       |   
							//    1  |  0
							JDPoint p1077 = d.fractionTowards((lev-d.value)/(c.value-d.value), c);
							JDPoint p1077n = a.fractionTowards(0.5d, b);
							if(path.isEmpty()) path.add(p1077); path.add(p1077.fractionTowards(0.5d, p1077n)); break;
						case 0x1100:
							//    0     0
							//    -------
							//    1     1
							addCurve(path, q, 1d, 1d-2*cb, -1d, 1d-2*da, lev, true); nu = u-1; nv = v; break;
						case 0x1101:
							//    1   \_0
							//          \
							//    1     1
							addCurve(path, q, 1d, 1d-2*cb, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1; break;
						case 0x1107:
							//    N __  0
							//        '--
							//    1     1
							JDPoint p1107 = d.fractionTowards((lev-d.value)/(b.value-d.value), b);
							if(path.isEmpty()) path.add(c.fractionTowards((lev-c.value)/(b.value-c.value), b));
							path.add(p1107.fractionTowards(0.5d, a)); break;
						case 0x1110:
							//    0_/   1
							//    /      
							//    1     1
							addCurve(path, q, 1d-2*ba, -1d, -1d, 1d-2*da, lev, true); nu = u-1; nv = v; break;
						case 0x1111: break;
						case 0x1117: break;
						case 0x1170:
							//    0  __ N
							//    --'    
							//    1     1
							JDPoint p1170 = c.fractionTowards((lev-c.value)/(a.value-c.value), a);
							path.add(p1170.fractionTowards(0.5d, b));
							path.add(d.fractionTowards((lev-d.value)/(a.value-d.value), a)); nu = u-1; nv = v; break;
						case 0x1171: break;
						case 0x1177: break;
						case 0x1700:
							//    0     0
							//    --.__  
							//    1     N
							JDPoint p1700 = d.fractionTowards((lev-d.value)/(b.value-d.value), b);
							path.add(p1700.fractionTowards(0.5d, c));
							path.add(d.fractionTowards((lev-d.value)/(a.value-d.value), a)); nu = u-1; nv = v; break;
						case 0x1701:
							//    1  \  0
							//        \  
							//    1     N
							JDPoint p1701 = d.fractionTowards((lev-d.value)/(b.value-d.value), b);
							path.add(p1701.fractionTowards(0.5d, c));
							path.add(a.fractionTowards((lev-a.value)/(b.value-a.value), b)); nu = u; nv = v-1; break;
						case 0x1707:
							//    N _   0
							//       \_  
							//    1     N
							JDPoint p1707 = d.fractionTowards((lev-d.value)/(b.value-d.value), b);
							path.add(p1707.fractionTowards(0.5d, c)); path.add(p1707.fractionTowards(0.5d, a)); break;
						case 0x1710:
							//    0_/   1
							//    /
							//    1     N
							if(path.isEmpty()) path.add(b.fractionTowards((lev-b.value)/(a.value-b.value), a));
							path.add(d.fractionTowards((lev-d.value)/(a.value-d.value), a)); nu = u-1; nv = v; break;
						case 0x1711: break;
						case 0x1770:
							//    0     N
							//    ----
							//    1     N
							JDPoint p1770 = d.fractionTowards((lev-d.value)/(a.value-d.value), a);
							JDPoint p1770n = c.fractionTowards(0.5d, b);
							path.add(p1770.fractionTowards(0.5d, p1770n)); path.add(p1770); nu = u-1; nv = v; break;
						case 0x1771: break;
						case 0x1777: break;
						case 0x7000: break;
						case 0x7001:
							//    1  /  0
							//      /    
							//    N     0
							JDPoint p7001 = a.fractionTowards((lev-a.value)/(c.value-a.value), c);
							path.add(p7001.fractionTowards(0.5d, d));
							path.add(a.fractionTowards((lev-a.value)/(b.value-a.value), b)); nu = u; nv = v-1; break;
						case 0x7007: break;
						case 0x7010:
							//    0   \_1
							//          \
							//    N     0
							if(path.isEmpty()) path.add(b.fractionTowards((lev-b.value)/(a.value-b.value), a));
							path.add(b.fractionTowards((lev-b.value)/(c.value-b.value), c)); nu = u+1; nv = v; break;
						case 0x7011:
							//    1     1
							//      __.--
							//    N     0
							JDPoint p7011 = a.fractionTowards((lev-a.value)/(c.value-a.value), c);
							path.add(p7011.fractionTowards(0.5d, d));
							path.add(b.fractionTowards((lev-b.value)/(c.value-b.value), c)); nu = u+1; nv = v; break;
						case 0x7017:
							//    N     1
							//       ----
							//    N     0
							JDPoint p7017 = b.fractionTowards((lev-b.value)/(c.value-b.value), c);
							JDPoint p7017n = a.fractionTowards(0.5d, d);
							path.add(p7017.fractionTowards(0.5d, p7017n)); path.add(p7017); nu = u+1; nv = v; break;
						case 0x7070: break;
						case 0x7071:
							//    1   _ N
							//      _/   
							//    N     0
							JDPoint p7071 = a.fractionTowards((lev-a.value)/(c.value-a.value), c);
							path.add(p7071.fractionTowards(0.5d, d)); path.add(p7071.fractionTowards(0.5d, b)); break;
						case 0x7077: break;
						case 0x7100:
							//    0     0
							//      __.--
							//    N     1
							JDPoint p7100 = c.fractionTowards((lev-c.value)/(a.value-c.value), a);
							if(path.isEmpty()) path.add(c.fractionTowards((lev-c.value)/(b.value-c.value), b));
							path.add(p7100.fractionTowards(0.5d, d)); break;
						case 0x7101:
							//    1   \_0
							//          \
							//    N     1
							if(path.isEmpty()) path.add(c.fractionTowards((lev-c.value)/(b.value-c.value), b));
							path.add(a.fractionTowards((lev-a.value)/(b.value-a.value), b)); nu = u; nv = v-1; break;
						case 0x7107:
							//    N     0
							//       ----
							//    N     1
							JDPoint p7107 = c.fractionTowards((lev-c.value)/(b.value-c.value), b);
							JDPoint p7107n = d.fractionTowards(0.5d, a);
							if(path.isEmpty()) path.add(p7107); path.add(p7107.fractionTowards(0.5d, p7107n)); break;
						case 0x7110:
							//    0  /  1
							//      /
							//    N     1
							JDPoint p7110 = c.fractionTowards((lev-c.value)/(a.value-c.value), a);
							if(path.isEmpty()) path.add(b.fractionTowards((lev-b.value)/(a.value-b.value), a));
							path.add(p7110.fractionTowards(0.5d, d)); break;
						case 0x7111: break;
						case 0x7117: break;
						case 0x7170:
							//    0   _ N
							//      _/   
							//    N     1
							JDPoint p7170 = c.fractionTowards((lev-c.value)/(a.value-c.value), a);
							path.add(p7170.fractionTowards(0.5d, b)); path.add(p7170.fractionTowards(0.5d, d)); break;
						case 0x7171: break;
						case 0x7177: break;
						case 0x7700: break;
						case 0x7701:
							//    1  |  0
							//       |   
							//    N     N
							JDPoint p7701 = a.fractionTowards((lev-a.value)/(b.value-a.value), b);
							JDPoint p7701n = d.fractionTowards(0.5d, c);
							path.add(p7701.fractionTowards(0.5d, p7701n)); path.add(p7701); nu = u; nv = v-1; break;
						case 0x7707: break;
						case 0x7710:
							//    0  |  1
							//       |   
							//    N     N
							JDPoint p7710 = b.fractionTowards((lev-b.value)/(a.value-b.value), a);
							JDPoint p7710n = c.fractionTowards(0.5d, d);
							if(path.isEmpty()) path.add(p7710); path.add(p7710.fractionTowards(0.5d, p7710n)); break;
						case 0x7711: break;
						case 0x7717: break;
						case 0x7770: break;
						case 0x7771: break;
						case 0x7777: break;
					}
					if(visits[v][u]==0) visits[v][u] = 2;
					if(nu<0 || nv<0) break;
					if(nu>=ilen || nv>=jlen) break;
					pu = u; pv = v;
					u = nu; v = nv;
				}
				if(path.size()>1) {
//					System.out.println("  -:- try to add path of size "+path.size());
//					for(JDPoint pt: path) System.out.println("      --> ["+pt.px+", "+pt.py+", "+pt.value+"]");
					JDLine line = new JDLine(path.toArray(new JDPoint[0]));
//					contours.addAll(line.intersectsAABB(p[0],p[0]+p[2], p[1],p[1]+p[3]));
					contours.add(line);
				}
				if(hit_saddle) i--;
			}
			JDLine[] cntArr = contours.toArray(new JDLine[0]);
			contours.clear();
			if(ax.isGeoAxis())
				for(JDLine line: cntArr) contours.add(line);
			else
				for(JDLine line: cntArr) contours.addAll(line.intersectsAABB(p[0], p[0]+p[2], p[1], p[1]+p[3]));
			
			//draw those lines:
//TODO remove
//			if(ax.getPlot().isDebug())
				System.out.println("Found "+contours.size()+" lines for cnt-level "+lev);
			String ls = contourStyle[l];
			double lln = 1d, llf = 0d, lpn = 0d, lpf = 0d;
			if ("-".equals(ls)) { lln = 1000 * lw; llf = 0; lpn = 0; lpf = 0; }
			if (".".equals(ls)) { lln = 0; llf = 0; lpn = 1 * lw; lpf = 3 * lw; }
			if (",".equals(ls)) { lln = 8 * lw; llf = 7 * lw; lpn = 0; lpf = 0; }
			if (";".equals(ls)) { lln = 8 * lw; llf = 3 * lw; lpn = 1 * lw; lpf = 3 * lw; }
			JPlotShape.stroke(lcs[l]);
			JPlotShape.strokeWeight((float)lw);
			for(JDLine line: contours) {
				drawSingleLine(s, line, lln, llf, lpn, lpf);
			}
		}
	}
	private void drawSingleLine(JGroupShape linesh, JDLine line, double lln, double llf, double lpn, double lpf) {
		if(line.getPoints().length<2) return;
		int lc = JPlotShape.strokeColour;
		float lw = JPlotShape.strokeWeight;
		if(lln>999d*lw) { //continuous line
			linesh.addChild(new JLineShape(lw, lc, line.getCoordsAsFloats()));
			return;
		}
		JDPoint[] pnts = line.getPoints();
		int li = 0;
		double lpos = 0d, ldif = 0d, loff = 0d;
		for (int i=0,j=1; j<pnts.length; i=j++) {
			JDPoint lvs = pnts[i];
			JDPoint lve = pnts[j];
			double x1 = lvs.x; //p[0] + xs * (invertAxisX ? Xax - lvs.x : lvs.x - Xin);
			double x2 = lve.x; //p[0] + xs * (invertAxisX ? Xax - lve.x : lve.x - Xin);
			double y1 = lvs.y; //p[1] + ys * (invertAxisY ? lvs.y - Yin : Yax - lvs.y);
			double y2 = lve.y; //p[1] + ys * (invertAxisY ? lve.y - Yin : Yax - lve.y);
			double dx = x2 - x1, dy = y2 - y1;
			double l = Math.sqrt(dx * dx + dy * dy);
			dx /= l;
			dy /= l;
			lpos = 0d;
			while (lpos < l) {
				ldif = 0d;
				switch (li) {
					case 0:
						if (lln == 0d) break;
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
				float	xf1 = (float) (x1 + lpos * dx), yf1 = (float) (y1 + lpos * dy),
						xf2 = (float) (x1 + (lpos + ldif) * dx), yf2 = (float) (y1 + (lpos + ldif) * dy);
				if (li%2 == 0 && ldif > 0d)
					linesh.addChild(new JLineShape(lw, lc, xf1, yf1, xf2, yf2));
				loff += ldif;
				switch (li) {
					case 0: if(loff >= lln) { loff -= lln; li = 1; } break;
					case 1: if(loff >= llf) { loff -= llf; li = 2; } break;
					case 2: if(loff >= lpn) { loff -= lpn; li = 3; } break;
					case 3: if(loff >= lpf) { loff -= lpf; li = 0; } break;
				}
				lpos += ldif;
			}
		}
	}
	private int getColor(double value) {
		if(Double.isNaN(value)) return colourtable.getColour(Double.NaN);
		if(value<contourIntervals[0]) return colourtable.getColour(-1d);
		int cl = contourIntervals.length-1;
		if(value>contourIntervals[cl]) return colourtable.getColour(2d);
		for(int i=1; i<=cl; i++)
			if(value<=contourIntervals[i])
				return colourtable.getColour((i-0.5d)/cl);
		return colourtable.getColour(Double.NaN);
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
	private static void addCurve(List<JDPoint> path, JDQuad q, double u0, double v0, double u4, double v4, double l, boolean check_emptynes) {
		JDPoint p0 = q.pointFromUV(u0, v0);
		double[] m = q.refineUV(0.5d*(u0+u4), 0.5d*(v0+v4), l);
		double[] mn = q.refineUV(0.5d*(u0+m[0]), 0.5d*(v0+m[1]), l);
		double[] mp = q.refineUV(0.5d*(u4+m[0]), 0.5d*(v4+m[1]), l);
		JDPoint p4 = q.pointFromUV(u4, v4);
		if(path.isEmpty() || !check_emptynes) path.add(p0);
		path.add(q.pointFromUV(mn));
		path.add(q.pointFromUV(m));
		path.add(q.pointFromUV(mp));
		path.add(p4);
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
				failed = !list.get(l).union(tri, tol);
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
				failed = !list.get(l).union(poly, tol);
			}
			if (failed)
				list.add(poly.copy());
		}
	}
	
	private static void show(JDPoint A, JDPoint B, JDPoint C, JDPoint D) {
		String o = "//    "+nf(A.value)+" -- "+nf(B.value) + "        ";
		String m = "//    "+  "  |  "  +"    "+  "  |  "   + "   =>   ";
		String u = "//    "+nf(D.value)+" -- "+nf(C.value) + "        ";
		o += "("+nf(A.x)+"/"+nf(A.y)+") -- ("+nf(B.x)+"/"+nf(B.y)+")";
		m += " "+   "     |     "   +"      "+   "     |     "   +" ";
		u += "("+nf(D.x)+"/"+nf(D.y)+") -- ("+nf(C.x)+"/"+nf(C.y)+")";
		System.out.println(o+"\n"+m+"\n"+u);
	}
	private static String nf(double v) {
		double a = Math.abs(v)+0.005d;
		int i = (int)(100d*(a-(int)a));
		return (v<0d?"-":"+")+(int)a+"."+(i<10?"0":"")+i;
	}
}
