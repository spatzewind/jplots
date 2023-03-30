package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import jplots.JPlot;
import jplots.axes.JAxis;
import jplots.colour.JColourtable;
import jplots.maths.AffineBuilder;
import jplots.maths.JDLine;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import jplots.maths.JDQuad;
import jplots.maths.JDTriangle;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JImageShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JPolygonShape;
import processing.core.PConstants;
import processing.core.PGraphics;

public class JPhaseLayer2D extends JPlotsLayer {

//	private double EPSILON = Math.pow(2, -52);
	private boolean isFilled, pixelFilling;
	private double minZ, maxZ, Xin, Xax, Yin, Yax, phScl;
	private double[] contourIntervals;
	private String[] contourStyle;
	private double[][] cosZ,sinZ;
	private JDPoint[][] corners;
	
	public JPhaseLayer2D(float[] x, float[] y, float[][] zc, float[][] zs, float phasemin, float phasemax, int nIntervals, float[] phaseIntervals, float phaseScale) {
		boolean is_valid = (x!=null && y!=null && zc!=null && zs!=null);
		if(is_valid) {
			is_valid = y.length == zc.length && y.length == zs.length;
			for(int j=0; j<y.length && is_valid; j++)
				is_valid = (x.length == zc[j].length && x.length == zs[j].length);
		}
		if(!is_valid) {
			null_init();
			return;
		}
		cosZ = JPlotMath.toDoubleArray2D(zc);
		sinZ = JPlotMath.toDoubleArray2D(zs);
		phScl = phaseScale;
		corners = new JDPoint[y.length][x.length];
		for(int j=0; j<y.length; j++)
			for(int i=0; i<x.length; i++)
				corners[j][i] = new JDPoint(x[i], y[j], 0d);
		if(phaseIntervals==null)
			contourIntervals = JPlotMath.linspace((double)phasemin, (double)phasemax, nIntervals);
		else
			contourIntervals = JPlotMath.toDoubleArray1D(phaseIntervals);
		init();
	}
	public JPhaseLayer2D(float[][] x, float[][] y, float[][] zc, float[][] zs, float phasemin, float phasemax, int nIntervals, float[] phaseIntervals, float phaseScale) {
		boolean is_valid = (x!=null && y!=null && zc!=null && zs!=null);
		if(is_valid) {
			is_valid = (x.length == zc.length && x.length == zs.length && y.length == zc.length && y.length == zs.length);
			int inner_length = 0; if(is_valid && x.length>0) inner_length = x[0].length;
			for(int j=0; j<x.length && is_valid; j++)
				is_valid = (x[j].length == zc[j].length && x[j].length == zs[j].length &&
						    y[j].length == zc[j].length && y[j].length == zs[j].length &&
						    inner_length == zc[j].length && inner_length == zs[j].length);
		}
		if(!is_valid) {
			null_init();
			return;
		}
		cosZ = JPlotMath.toDoubleArray2D(zc);
		sinZ = JPlotMath.toDoubleArray2D(zs);
		phScl = phaseScale;
		corners = new JDPoint[y.length][x.length];
		for(int j=0; j<x.length; j++)
			for(int i=0; i<x[j].length; i++)
				corners[j][i] = new JDPoint(x[j][i], y[j][i], 0d);
		if(phaseIntervals==null)
			contourIntervals = JPlotMath.linspace((double)phasemin, (double)phasemax, nIntervals);
		else
			contourIntervals = JPlotMath.toDoubleArray1D(phaseIntervals);
		init();
	}
	public JPhaseLayer2D(double[] x, double[] y, double[][] zc, double[][] zs, double phasemin, double phasemax, int nIntervals, double[] phaseIntervals, double phaseScale) {
		boolean is_valid = (x!=null && y!=null && zc!=null && zs!=null);
		if(is_valid) {
			is_valid = y.length == zc.length && y.length == zs.length;
			for(int j=0; j<y.length && is_valid; j++)
				is_valid = (x.length == zc[j].length && x.length == zs[j].length);
		}
		if(!is_valid) {
			null_init();
			return;
		}
		cosZ = JPlotMath.toDoubleArray2D(zc);
		sinZ = JPlotMath.toDoubleArray2D(zs);
		phScl = phaseScale;
		corners = new JDPoint[y.length][x.length];
		for(int j=0; j<y.length; j++)
			for(int i=0; i<x.length; i++)
				corners[j][i] = new JDPoint(x[i], y[j], 0d);
		if(phaseIntervals==null)
			contourIntervals = JPlotMath.linspace(phasemin, phasemax, nIntervals);
		else
			contourIntervals = phaseIntervals;
		init();
	}
	public JPhaseLayer2D(double[][] x, double[][] y, double[][] zc, double[][] zs, double phasemin, double phasemax, int nIntervals, double[] phaseIntervals, double phaseScale) {
		boolean is_valid = (x!=null && y!=null && zc!=null && zs!=null);
		if(is_valid) {
			is_valid = (x.length == zc.length && x.length == zs.length && y.length == zc.length && y.length == zs.length);
			int inner_length = 0; if(is_valid && x.length>0) inner_length = x[0].length;
			for(int j=0; j<x.length && is_valid; j++)
				is_valid = (x[j].length == zc[j].length && x[j].length == zs[j].length &&
							y[j].length == zc[j].length && y[j].length == zs[j].length &&
							inner_length == zc[j].length && inner_length == zs[j].length);
		}
		if(!is_valid) {
			null_init();
			return;
		}
		cosZ = JPlotMath.toDoubleArray2D(zc);
		sinZ = JPlotMath.toDoubleArray2D(zs);
		phScl = phaseScale;
		corners = new JDPoint[y.length][x.length];
		for(int j=0; j<x.length; j++)
			for(int i=0; i<x[j].length; i++)
				corners[j][i] = new JDPoint(x[j][i], y[j][i], 0d);
		if(phaseIntervals==null)
			contourIntervals = JPlotMath.linspace(phasemin, phasemax, nIntervals);
		else
			contourIntervals = phaseIntervals;
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
					double[] xy = inputProj.fromPROJtoLATLON(cnt2[j][i].x, cnt2[j][i].y, false, false);
					xy = ax.getGeoProjection().fromLATLONtoPROJ(xy[0], xy[1], false, false);
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
				.translate(invertAxisX ? Xax : -Xin, invertAxisY ? -Yin : Yax).scale(xs, ys).translate(p[0], p[1]);
		minZ = 1d;
		maxZ = 0d;
		if(contourIntervals[0]<minZ) minZ = contourIntervals[0];
		if(contourIntervals[contourIntervals.length-1]>maxZ) maxZ = contourIntervals[contourIntervals.length-1];
		
		// step 1: create filling between contours if wished
		if (isFilled) {
			if (pixelFilling) {
				if (ax.getPlot().isDebug())
					System.out.println("[DEBUG] JContourLayer2D: 1] contour filling pixelwise ...");
				fillPixelByPixel(p, ax, xs, ys, s, cnt2, affine.getMatrix());
			} else {
				if (ax.getPlot().isDebug())
					System.out.println("[DEBUG] JContourLayer2D: 1] contour filling vectorwise ...");
				fillVectorByVector(p, ax, xs, ys, s, cnt2, affine.getMatrix());
			}
		} else {
			if (ax.getPlot().isDebug())
				System.out.println("[DEBUG] JContourLayer2D: 1] no filling ...");
		}

		// step 2: add contours to plot
		if (drawLines) {
			if (ax.getPlot().isDebug())
				System.out.println("[DEBUG] JContourLayer2D: 2] draw contour line segments with " + lw + "px line width...");
			drawContourLines(p, ax, xs, ys, s, cnt2, affine.getMatrix());
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
		minZ = 0d;
		maxZ = 1d;
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
		minZ = 0d;
		maxZ = phScl;
		for(int j=0; j<corners.length; j++)
			for(int i=0; i<corners[j].length; i++) {
				if(corners[j][i].x<minX) minX = corners[j][i].x;
				if(corners[j][i].x>maxX) maxX = corners[j][i].x;
				if(corners[j][i].y<minY) minY = corners[j][i].y;
				if(corners[j][i].y>maxY) maxY = corners[j][i].y;
			}
		contourStyle = new String[contourIntervals.length];
		lw = 2d;
		lc  = 0xff000000;
		lcs = new int[contourIntervals.length];
		ls = "-";
		for(int i=0; i<contourIntervals.length; i++) {
			contourStyle[i] = "-";
			lcs[i] = lc;
		}
	}
	
	private void fillVectorByVector(int[] p, JAxis ax, double xs, double ys, JGroupShape s, JDPoint[][] input, double[][] affine) {
		//TODO implement fillVectorByVector
//		JDPoint[][] points = new JDPoint[input.length][input[0].length];
//		int jlen = input.length;
//		int ilen = input[0].length;
//		for(int j=0; j<jlen; j++)
//			for(int i=0; i<ilen; i++)
//				points[j][i] = input[j][i].copy();
//		double[] clev2 = new double[2+contourIntervals.length];
//		clev2[0] = Double.NEGATIVE_INFINITY;
//		for(int l=0; l<contourIntervals.length; l++) clev2[l+1] = contourIntervals[l];
//		clev2[contourIntervals.length+1] = Double.POSITIVE_INFINITY;
//		for(int l=1; l<clev2.length; l++) {
//			double lower = clev2[l-1];
//			double upper = clev2[l];
//			//find all polygons
//			List<JDPolygon> polys1 = getIntBasedPolygons(input, lower, upper, 0,points.length-1, 0,points[0].length-1);
//			if(polys1==null) continue;
//			List<JDPolygon> polys2 = new ArrayList<JDPolygon>();
//			for(JDPolygon poly: polys1) {
//				JDPoint[] pnts = poly.c;
//				//recalc position for all polygons
//				for(int n=0; n<pnts.length; n++) {
//					double x = pnts[n].x(); int xi = Math.min((int)x, ilen-2);
//					double y = pnts[n].y(); int yj = Math.min((int)y, jlen-2);
//					x -= xi; y -= yj;
//					double xt=points[ yj ][xi].x*(1d-x)+points[ yj ][xi+1].x*x;
//					double yt=points[ yj ][xi].y*(1d-x)+points[ yj ][xi+1].y*x;
//					double xb=points[yj+1][xi].x*(1d-x)+points[yj+1][xi+1].x*x;
//					double yb=points[yj+1][xi].y*(1d-x)+points[yj+1][xi+1].y*x;
//					pnts[n].x = xt*(1d-y)+xb*y;
//					pnts[n].y = yt*(1d-y)+yb*y;
//				}
//				poly.c = pnts;
//				//cutoff what is beyond allowed mapping range
//				polys2.addAll(poly.splitByMapBorder(ax));
//			}
//			int cct = colourtable.getColour((l-1.5d)/(clev2.length-3d));
////			System.out.println("current polygon color: "+Integer.toHexString(cct));
//			JGroupShape cnt_polys = new JGroupShape();
//			for(JDPolygon poly: polys2) {
//				poly.affine(affine);
////				double[] bnds = poly.getBounds();
////				System.out.println("polygon: x={"+bnds[0]+" ... "+bnds[1]+"}  y={"+bnds[2]+" ... "+bnds[3]+"}  [v]="+poly.c.length);
////				for(JDPolygon subpoly: poly.affine(affine).intersectsAABB(p[0], p[0]+p[2], p[1], p[1]+p[3]))
//				cnt_polys.addChild(new JPolygonShape(poly.getCoords(), cct, true, true));
//			}
////			System.out.println("Should be added "+cnt_polys.childCount()+" shapes...");
//			s.addChild(cnt_polys);
//		}
	}
	
	private void fillPixelByPixel(int[] p, JAxis ax, double xs, double ys, JGroupShape s, JDPoint[][] input, double[][] affine) {
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
		JDPoint[][] cpnts = new JDPoint[input.length][input[0].length];
		JDPoint[][] spnts = new JDPoint[input.length][input[0].length];
		for(int j=0; j<input.length; j++)
			for(int i=0; i<input[j].length; i++) {
				cpnts[j][i] = input[j][i].copy().affine(affine);
				spnts[j][i] = cpnts[j][i].copy();
				cpnts[j][i].value = cosZ[j][i];
				spnts[j][i].value = sinZ[j][i];
			}
		img.loadPixels();
		for(int px=0; px<img.pixels.length; px++)
			img.pixels[px] = 0x00999999;
		for(int j=0; j+1<cpnts.length; j++) {
			for(int i=0; i+1<cpnts[j].length; i++) {
				JDQuad quc = new JDQuad(cpnts[j+1][i], cpnts[j+1][i+1], cpnts[j][i+1], cpnts[j][i]);
				JDQuad qus = new JDQuad(spnts[j+1][i], spnts[j+1][i+1], spnts[j][i+1], spnts[j][i]);
				double	txi = Math.min(Math.min(quc.x[0],quc.x[1]), Math.min(quc.x[2], quc.x[3])),
						txa = Math.max(Math.max(quc.x[0],quc.x[1]), Math.max(quc.x[2], quc.x[3]));
				double	tyi = Math.min(Math.min(quc.y[0],quc.y[1]), Math.min(quc.y[2], quc.y[3])),
						tya = Math.max(Math.max(quc.y[0],quc.y[1]), Math.max(quc.y[2], quc.y[3]));
				int ixs =  Math.max((int) txi    - (txi < 0 ? 1 : 0), p[0]),
					ixe = -Math.max((int) (-txa) - (txa > 0 ? 1 : 0), 1 - p[0] - p[2]);
				int iys =  Math.max((int) tyi    - (tyi < 0 ? 1 : 0), p[1]),
					iye = -Math.max((int) (-tya) - (tya > 0 ? 1 : 0), 1 - p[1] - p[3]);
				if ((ixe < ixs) || (iye < iys)) continue;
				for (int v = iys; v <= iye; v++) {
					for (int u = ixs; u <= ixe; u++) {
						JDPoint ij = new JDPoint(u+0.5d,v+0.5d);
						if(!quc.contains(ij)) continue;
						double[] uv = quc.barycentricCoords(ij);
						double cosPart = quc.valueAt(ij);
						double sinPart = qus.valueAt(ij);
						double rad = Math.sqrt(cosPart*cosPart + sinPart*sinPart);
						double phs = 0.5d*phScl*Math.acos(cosPart/(rad+1.0e-10d))/Math.PI;
						if(sinPart<0d) phs = phScl-phs;
						img.pixels[(v-p[1]) * p[2] + u-p[0]] = getColor(phs);
					}
				}
			}
		}
		img.updatePixels();
		s.addChild(new JImageShape(img, p[0], p[1], p[2], p[3]));
	}
	
	private void drawContourLines(int[] p, JAxis ax, double xs, double ys, JGroupShape s, JDPoint[][] points, double[][] affine) {
		//TODO implement drawContourLines
//		int jlen = points.length-1;
//		int ilen = points[0].length-1;
//		int[][] visits = new int[jlen][ilen];
//		List<JDPoint> path = new ArrayList<>();
//		List<JDLine> contours = new ArrayList<>();
//		for(int l=0; l<contourIntervals.length; l++) {
//			contours.clear();
//			for(int j=0; j<jlen; j++)
//				for(int i=0; i<ilen; i++)
//					visits[j][i] = 0;
//			double lev = contourIntervals[l];
//			
//			for(int j=0; j<jlen; j++) for(int i=0; i<ilen; i++) {
//				if(visits[j][i]>1) continue;
//				boolean hit_saddle = false;
//				int su = i, sv = j;
//				int pu = i, pv = j;
//				int u = i, v = j;
//				int nu = -1, nv = -1;
//				path.clear();
//				while(su!=nu || sv!=nv) {
//					if(visits[v][u]==2) {
//						if(path.isEmpty()) break;
//						JDLine nl = new JDLine(path.toArray(new JDPoint[0]));
//						boolean success = false;
//						for(JDLine ol: contours) {
//							success = ol.join(nl, 0.01d);
//							if(success) break;
//						}
//						if(success) {
//							path.clear();
//							break;
//						}
//					}
//					nu = -1; nv = -1;
//					JDPoint a = points[ v ][ u ];
//					JDPoint b = points[ v ][u+1];
//					JDPoint c = points[v+1][u+1];
//					JDPoint d = points[v+1][ u ];
//					int code =	(Double.isNaN(a.value) ? 0x0007 : a.value<lev ? 0x0000 : 0x0001) |
//								(Double.isNaN(b.value) ? 0x0070 : b.value<lev ? 0x0000 : 0x0010) |
//								(Double.isNaN(c.value) ? 0x0700 : c.value<lev ? 0x0000 : 0x0100) |
//								(Double.isNaN(d.value) ? 0x7000 : d.value<lev ? 0x0000 : 0x1000);
//					double ab = (lev-a.value)/(b.value-a.value);
//					double ac = (lev-a.value)/(c.value-a.value);
//					double ad = (lev-a.value)/(d.value-a.value);
//					double ba = (lev-b.value)/(a.value-b.value);
//					double bc = (lev-b.value)/(c.value-b.value);
//					double bd = (lev-b.value)/(d.value-b.value);
//					double ca = (lev-c.value)/(a.value-c.value);
//					double cb = (lev-c.value)/(b.value-c.value);
//					double cd = (lev-c.value)/(d.value-c.value);
//					double da = (lev-d.value)/(a.value-d.value);
//					double db = (lev-d.value)/(b.value-d.value);
//					double dc = (lev-d.value)/(c.value-d.value);
//					JDQuad q = new JDQuad(a,b,c,d);
//					switch(code) {
//						default: break;
//						case 0x0000: break;
//						case 0x0001:
//							//    1_/   0
//							//    /      
//							//    0     0
//							q.addCurve(path, -1d, 2*ad-1d, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1; break;
//						case 0x0007: break;
//						case 0x0010:
//							//    0   \_1
//							//          \
//							//    0     0
//							q.addCurve(path, 1d-2*ba, -1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v; break;
//						case 0x0011:
//							//    1     1
//							//    -------
//							//    0     0
//							q.addCurve(path, -1d, 2*ad-1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v; break;
//						case 0x0017:
//							//    N __  1
//							//        '--
//							//    0     0
//							q.addCurve(path, -bd, bd-1d, 1d, 2*bc-1d, lev, false); nu = u+1; nv = v; break;
//						case 0x0070: break;
//						case 0x0071:
//							//    1  __ N
//							//    --'    
//							//    0     0
//							q.addCurve(path, -1d, 2*ad-1d, ac, ac-1d, lev, true); break;
//						case 0x0077: break;
//						case 0x0100:
//							//    0     0
//							//         _/
//							//    0   / 1
//							//show(a,b,c,d);
//							q.addCurve(path, 1d, 1d-2*cb, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1; break;
//						case 0x0101:
//							// du = 0.25d*((value[1]-value[0])*(1-v) + (value[2]-value[3])*(1+v)) = 0
//							// dv = 0.25d*((value[3]-value[0])*(1-u) + (value[2]-value[1])*(1+u)) = 0
//							/*
//							 * (v1-v0)*(1-y)+(v2-v3)*(1+y) = 0
//							 * v1-v0+v2-v3 = (v1-v0-v2+v3)*y
//							 * 
//							 * (v3-v0)*(1-x) + (v2-v1)*(1+x) = 0
//							 * v3-v0+v2-v1 = (v3-v0-v2+v1)*x
//							 * 
//							 */
//							hit_saddle = true;
//							double u01 = (d.value-a.value+c.value-b.value)/(d.value-a.value-c.value+b.value);
//							double v01 = (b.value-a.value+c.value-d.value)/(b.value-a.value-c.value+d.value);
//							double m01 = 0.25d*(a.value*(1-u01)*(1-v01)+b.value*(1+u01)*(1-v01)+c.value*(1+u01)*(1+v01)+d.value*(1-u01)*(1+v01));
//							boolean takePos01 = visits[v][u]<0;
//							if(visits[v][u]==0) takePos01 = (pu>u);
//							if(m01<=lev) {
//								//    1_/   0
//								// -- /  0 _/ +
//								//    0   / 1
//								if(takePos01) {
//									q.addCurve(path, 1d, 1d-2*cb, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1;
//									visits[v][u] = visits[v][u]==0 ? 1 : 2;
//								} else {
//									q.addCurve(path, -1d, 2*ad-1d, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1;
//									visits[v][u] = visits[v][u]==0 ? -1 : 2;
//								} break;
//							} else {
//								//    1   \_0
//								// -- \_ 1  \ +
//								//    0 \   1
//								if(takePos01) {
//									q.addCurve(path, 1d, 1d-2*cb, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1;
//									visits[v][u] = visits[v][u]==0 ? 1 : 2;
//								} else {
//									q.addCurve(path, -1d, 2*ad-1d, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1;
//									visits[v][u] = visits[v][u]==0 ? -1 : 2;
//								} break;
//							}
//						case 0x0107:
//							//    N     0
//							//         _/
//							//    0   / 1
//							q.addCurve(path, 1d, 1d-2*cb, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1; break;
//						case 0x0110:
//							//    0  |  1
//							//       |   
//							//    0  |  1
//							q.addCurve(path, 1d-2*ba, -1d, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1; break;
//						case 0x0111:
//							//    1     1
//							//    \_     
//							//    0 \   1
//							q.addCurve(path, -1d, 2*ad-1d, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1; break;
//						case 0x0117:
//							//    N     1
//							//      \    
//							//    0  \  1
//							q.addCurve(path, -bd, bd-1d, 1d-2*cd, 1d, lev, false); nu = u; nv = v+1; break;
//						case 0x0170:
//							//    0     N
//							//        /  
//							//    0  /  1
//							q.addCurve(path, 1d-ca, -ca, 1d-2*cd, 1d, lev, false); nu = u; nv = v+1; break;
//						case 0x0171:
//							//    1     N
//							//    \_     
//							//    0 \   1
//							q.addCurve(path, -1d, 2*ad-1d, 1d-2*cd, 1d, lev, true); nu = u; nv = v+1; break;
//						case 0x0177:
//							//    N     N
//							//       |   
//							//    0  |  1
//							path.add(q.pointFromUV(1d-2*cd,0d)); path.add(q.pointFromUV(1d-2*cd,1d)); nu = u; nv = v+1; break;
//						case 0x0700: break;
//						case 0x0701:
//							//    1_/   0
//							//    /      
//							//    0     N
//							q.addCurve(path, -1d, 2*ad-1d, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1; break;
//						case 0x0707: break;
//						case 0x0710:
//							//    0  \  1
//							//        \  
//							//    0     N
//							q.addCurve(path, 1d-2*ba, -1d, 1d-bd, bd, lev, true); break;
//						case 0x0711:
//							//    1     1
//							//    --.__  
//							//    0     N
//							q.addCurve(path, -1d, 2*ad-1d, 1d-bd, bd, lev, true); break;
//						case 0x0717:
//							//    N _   1
//							//       \_  
//							//    0     N
//							q.addCurve(path, -bd, bd-1d, 1d-bd, bd, lev, false); break;
//						case 0x0770: break;
//						case 0x0771:
//							//    1     N
//							//    ----   
//							//    0     N
//							if(path.isEmpty()) path.add(q.pointFromUV(-1d,2*ad-1d)); path.add(q.pointFromUV(0d,2*ad-1d)); break;
//						case 0x0777: break;
//						case 0x1000:
//							//    0     0
//							//    \_     
//							//    1 \   0
//							q.addCurve(path, 2*dc-1d, 1d, -1d, 1d-2*da, lev, true); nu = u-1; nv = v; break;
//						case 0x1001:
//							//    1  |  0
//							//       |   
//							//    1  |  0
//							q.addCurve(path, 2*dc-1d, 1d, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1; break;
//						case 0x1007:
//							//    N     0
//							//      \    
//							//    1  \  0
//							q.addCurve(path, 2*dc-1d, 1d, db-1d, -db, lev, true); break;
//						case 0x1010:
//							hit_saddle = true;
//							double u10 = (d.value-a.value+c.value-b.value)/(d.value-a.value-c.value+b.value);
//							double v10 = (b.value-a.value+c.value-d.value)/(b.value-a.value-c.value+d.value);
//							double m10 = 0.25d*(a.value*(1-u10)*(1-v10)+b.value*(1+u10)*(1-v10)+c.value*(1+u10)*(1+v10)+d.value*(1-u10)*(1+v10));
//							boolean takePos10 = visits[v][u]<0;
//							if(m10<=lev) {
//								//    0   \_1
//								// -- \_ 0  \ +
//								//    1 \   0
//								if(visits[v][u]==0) takePos10 = (pv<v);
//								if(takePos10) {
//									q.addCurve(path, 1d-2*ba, -1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v;
//									visits[v][u] = visits[v][u]==0 ? 1 : 2;
//								} else {
//									q.addCurve(path, 2*dc-1d, 1d, -1d, 1d-2*da, lev, true); nu = u-1; nv = v;
//									visits[v][u] = visits[v][u]==0 ? -1 : 2;
//								} break;
//							} else {
//								//    0_/   1
//								// -- /  1 _/ +
//								//    1   / 0
//								if(visits[v][u]==0) takePos10 = (pv>v);
//								if(takePos10) {
//									q.addCurve(path, 2*dc-1d, 1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v;
//									visits[v][u] = visits[v][u]==0 ? 1 : 2;
//								} else {
//									q.addCurve(path, 1d-2*ba, -1d, -1d, 1d-2*da, lev, true); nu = u-1; nv = v;
//									visits[v][u] = visits[v][u]==0 ? -1 : 2;
//								} break;
//							}
//						case 0x1011:
//							//    1     1
//							//         _/
//							//    1   / 0
//							q.addCurve(path, 2*dc-1d, 1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v; break;
//						case 0x1017:
//							//    N     1
//							//         _/
//							//    1   / 0
//							q.addCurve(path, 2*dc-1d, 1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v; break;
//						case 0x1070:
//							//    0     N
//							//    \_     
//							//    1 \   0
//							q.addCurve(path, 2*dc-1d, 1d, -1d, 1d-2*da, lev, true); nu = u-1; nv = v; break;
//						case 0x1071:
//							//    1     N
//							//        /  
//							//    1  /  0
//							q.addCurve(path, 2*dc-1d, 1d, ac, ac-1d, lev, true); break;
//						case 0x1077:
//							//    N     N
//							//       |   
//							//    1  |  0
//							if(path.isEmpty()) path.add(q.pointFromUV(2*dc-1d,1d)); path.add(q.pointFromUV(2*dc-1d,0d)); break;
//						case 0x1100:
//							//    0     0
//							//    -------
//							//    1     1
//							q.addCurve(path, 1d, 1d-2*cb, -1d, 1d-2*da, lev, true); nu = u-1; nv = v; break;
//						case 0x1101:
//							//    1   \_0
//							//          \
//							//    1     1
//							q.addCurve(path, 1d, 1d-2*cb, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1; break;
//						case 0x1107:
//							//    N __  0
//							//        '--
//							//    1     1
//							q.addCurve(path, 1d, 1d-2*cb, db-1d, -db, lev, true); break;
//						case 0x1110:
//							//    0_/   1
//							//    /      
//							//    1     1
//							q.addCurve(path, 1d-2*ba, -1d, -1d, 1d-2*da, lev, true); nu = u-1; nv = v; break;
//						case 0x1111: break;
//						case 0x1117: break;
//						case 0x1170:
//							//    0  __ N
//							//    --'    
//							//    1     1
//							q.addCurve(path, 1d-ca, -ca, -1d, 1d-2*da, lev, false); nu = u-1; nv = v; break;
//						case 0x1171: break;
//						case 0x1177: break;
//						case 0x1700:
//							//    0     0
//							//    --.__  
//							//    1     N
//							q.addCurve(path, db, 1d-db, -1d, 1d-2*da, lev, false); nu = u-1; nv = v; break;
//						case 0x1701:
//							//    1  \  0
//							//        \  
//							//    1     N
//							q.addCurve(path, db, 1d-db, 2*ab-1d, -1d, lev, false); nu = u; nv = v-1; break;
//						case 0x1707:
//							//    N _   0
//							//       \_  
//							//    1     N
//							q.addCurve(path, db, 1d-db, db-1d, -db, lev, false); break;
//						case 0x1710:
//							//    0_/   1
//							//    /
//							//    1     N
//							q.addCurve(path, 1d-2*ba, -1d, -1d, 1d-2*da, lev, true); nu = u-1; nv = v; break;
//						case 0x1711: break;
//						case 0x1770:
//							//    0     N
//							//    ----
//							//    1     N
//							path.add(q.pointFromUV(0d,1d-2*da)); path.add(q.pointFromUV(-1d,1d-2*da)); nu = u-1; nv = v; break;
//						case 0x1771: break;
//						case 0x1777: break;
//						case 0x7000: break;
//						case 0x7001:
//							//    1  /  0
//							//      /    
//							//    N     0
//							q.addCurve(path, ac-1d, ac, 2*ab-1d, -1d, lev, false); nu = u; nv = v-1; break;
//						case 0x7007: break;
//						case 0x7010:
//							//    0   \_1
//							//          \
//							//    N     0
//							q.addCurve(path, 1d-2*ba, -1d, 1d, 2*bc-1d, lev, true); nu = u+1; nv = v; break;
//						case 0x7011:
//							//    1     1
//							//      __.--
//							//    N     0
//							q.addCurve(path, ac-1d, ac, 1d, 2*bc-1d, lev, false); nu = u+1; nv = v; break;
//						case 0x7017:
//							//    N     1
//							//       ----
//							//    N     0
//							path.add(q.pointFromUV(0d,2*bc-1d)); path.add(q.pointFromUV(1d,2*bc-1d)); nu = u+1; nv = v; break;
//						case 0x7070: break;
//						case 0x7071:
//							//    1   _ N
//							//      _/   
//							//    N     0
//							q.addCurve(path, ac-1d, ac, ac, ac-1d, lev, false); break;
//						case 0x7077: break;
//						case 0x7100:
//							//    0     0
//							//      __.--
//							//    N     1
//							q.addCurve(path, 1d, 1d-2*cb, -ca, 1d-ca, lev, true); break;
//						case 0x7101:
//							//    1   \_0
//							//          \
//							//    N     1
//							q.addCurve(path, 1d, 1d-2*cb, 2*ab-1d, -1d, lev, true); nu = u; nv = v-1; break;
//						case 0x7107:
//							//    N     0
//							//       ----
//							//    N     1
//							if(path.isEmpty()) path.add(q.pointFromUV(1d,1d-2*cb)); path.add(q.pointFromUV(0d,1d-2*cb)); break;
//						case 0x7110:
//							//    0  /  1
//							//      /
//							//    N     1
//							q.addCurve(path, 1d-2*ba, -1d, -ca, 1d-ca, lev, true); break;
//						case 0x7111: break;
//						case 0x7117: break;
//						case 0x7170:
//							//    0   _ N
//							//      _/   
//							//    N     1
//							q.addCurve(path, 1d-ca, -ca, -ca, 1d-ca, lev, false); break;
//						case 0x7171: break;
//						case 0x7177: break;
//						case 0x7700: break;
//						case 0x7701:
//							//    1  |  0
//							//       |   
//							//    N     N
//							path.add(q.pointFromUV(2*ab-1d,0d)); path.add(q.pointFromUV(2*ab-1d,-1d)); nu = u; nv = v-1; break;
//						case 0x7707: break;
//						case 0x7710:
//							//    0  |  1
//							//       |   
//							//    N     N
//							if(path.isEmpty()) path.add(q.pointFromUV(1d-2*ba,-1d)); path.add(q.pointFromUV(1d-2*ba,0d)); break;
//						case 0x7711: break;
//						case 0x7717: break;
//						case 0x7770: break;
//						case 0x7771: break;
//						case 0x7777: break;
//					}
//					if(visits[v][u]==0) visits[v][u] = 2;
//					if(nu<0 || nv<0) break;
//					if(nu>=ilen || nv>=jlen) break;
//					pu = u; pv = v;
//					u = nu; v = nv;
//				}
//				if(path.size()>1) {
//					contours.add(new JDLine(path.toArray(new JDPoint[0])));
//				}
//				if(hit_saddle) i--;
//			}
//			JDLine[] cntArr = contours.toArray(new JDLine[0]);
//			contours.clear();
//			if(ax.isGeoAxis()) {
//				for(JDLine line: cntArr) contours.addAll(ax.getGeoProjection().splitByMapBorder(line));
//				cntArr = contours.toArray(new JDLine[0]);
//				contours.clear();
//			}
//			for(int i=0; i<cntArr.length; i++)
//				cntArr[i].affine(affine);
//			for(JDLine line: cntArr) contours.addAll(line.intersectsAABB(p[0], p[0]+p[2], p[1], p[1]+p[3]));
//			
//			//draw those lines:
//			if(ax.getPlot().isDebug())
//				System.out.println("Found "+contours.size()+" lines for cnt-level "+lev);
//			String ls = contourStyle[l];
//			double lln = 1d, llf = 0d, lpn = 0d, lpf = 0d;
//			if ("-".equals(ls)) { lln = 1000 * lw; llf = 0; lpn = 0; lpf = 0; }
//			if (".".equals(ls)) { lln = 0; llf = 0; lpn = 1 * lw; lpf = 3 * lw; }
//			if (",".equals(ls)) { lln = 8 * lw; llf = 7 * lw; lpn = 0; lpf = 0; }
//			if (";".equals(ls)) { lln = 8 * lw; llf = 3 * lw; lpn = 1 * lw; lpf = 3 * lw; }
//			JPlotShape.stroke(lcs[l]);
//			JPlotShape.strokeWeight((float)lw);
//			for(JDLine line: contours) {
//				drawSingleLine(s, line, lln, llf, lpn, lpf);
//			}
//		}
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
				if (li%2 == 0 && ldif >= 0d)
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
	
	private static List<JDPolygon> getIntBasedPolygons(JDPoint[][] data, double lower, double upper, int i0, int i1, int j0, int j1) {
		List<JDPolygon> res = new ArrayList<JDPolygon>();
		if(i0==i1 || j0==j1)
			return res;
		if(i1-i0 > j1-j0) {
			int im = (i1+1-i0) / 2 + i0;
			for(JDPolygon poly: getIntBasedPolygons(data, lower, upper, i0, im, j0, j1)) {
				boolean isDisjunct = true;
				for(int t=res.size()-1; t>=0 && isDisjunct; t--)
					isDisjunct = !res.get(t).union(poly,0.0001d);
				if(isDisjunct)
					res.add(poly);
			}
			for(JDPolygon poly: getIntBasedPolygons(data, lower, upper, im, i1, j0, j1)) {
				boolean isDisjunct = true;
				for(int t=res.size()-1; t>=0 && isDisjunct; t--)
					isDisjunct = !res.get(t).union(poly,0.0001d);
				if(isDisjunct)
					res.add(poly);
			}
		} else
		if(j1-j0>1) {
			int jm = (j1+1-j0) / 2 + j0;
			for(JDPolygon poly: getIntBasedPolygons(data, lower, upper, i0, i1, j0, jm)) {
				boolean isDisjunct = true;
				for(int t=res.size()-1; t>=0 && isDisjunct; t--)
					isDisjunct = !res.get(t).union(poly,0.0001d);
				if(isDisjunct)
					res.add(poly);
			}
			for(JDPolygon poly: getIntBasedPolygons(data, lower, upper, i0, i1, jm, j1)) {
				boolean isDisjunct = true;
				for(int t=res.size()-1; t>=0 && isDisjunct; t--)
					isDisjunct = !res.get(t).union(poly,0.0001d);
				if(isDisjunct)
					res.add(poly);
			}
		}
		else {
			JDQuad qu = new JDQuad(
					new JDPoint(i0,j0, data[j0][i0].value),
					new JDPoint(i1,j0, data[j0][i1].value),
					new JDPoint(i1,j1, data[j1][i1].value),
					new JDPoint(i0,j1, data[j1][i0].value));
			JDPolygon[] polys = qu.getLevelRangePolygons(lower, upper);
			if(polys!=null) {
//				System.out.println("got "+polys.length+" polygon(s)");
				for(JDPolygon subpoly: polys)
					res.add(subpoly);
			}
		}
		return res;
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

//	private void addTriangle2polygonList(List<JDPolygon> list, JDTriangle tri, double tol) {
//		if (list.isEmpty()) {
//			list.add(tri.toPolygon());
//		} else {
//			boolean failed = true;
//			for (int l = 0; l < list.size() && failed; l++) {
//				if(list.get(l).c.length>999) continue;
//				failed = !list.get(l).union(tri, tol);
//			}
//			if (failed)
//				list.add(tri.toPolygon());
//		}
//	}

//	private void addPolygon2polygonList(List<JDPolygon> list, JDPolygon poly, double tol) {
//		if (list.isEmpty()) {
//			list.add(poly.copy());
//		} else {
//			boolean failed = true;
//			for (int l = 0; l < list.size() && failed; l++) {
//				if(list.get(l).c.length>999) continue;
//				failed = !list.get(l).union(poly, tol);
//			}
//			if (failed)
//				list.add(poly.copy());
//		}
//	}
}
