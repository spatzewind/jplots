package jplots.layer;

import java.util.Arrays;

import jplots.JPlot;
import jplots.axes.JAxis;
import jplots.axes.LogarithmicScale;
import jplots.maths.AffineBuilder;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JPolygonShape;
import processing.core.PGraphics;

public class JHistogramLayer extends JPlotsLayer {
	
	public static final int NORM_COUNT =  1;
	public static final int NORM_MAXONE = 2;
	public static final int NORM_ARRONE = 3;
	public static final int CUM_COUNT =   4;
	public static final int CUM_MAXONE  = 5;

	private double[] yarrayy, classbounds;
	private int normalizationMode;
	
	public JHistogramLayer(float[] y, float[] bounds, int normmode, int colour) {
		yarrayy = JPlotMath.toDoubleArray1D(y);
		classbounds = JPlotMath.toDoubleArray1D(bounds);
		normalizationMode = normmode;
		minX = JPlotMath.dmin(classbounds);
		maxX = JPlotMath.dmax(classbounds);
		pc = colour;
	}
	
	public JHistogramLayer(double[] y, double[] bounds, int normmode, int colour) {
		yarrayy = y;
		classbounds = bounds.clone();
		normalizationMode = normmode;
		minX = JPlotMath.dmin(classbounds);
		maxX = JPlotMath.dmax(classbounds);
		pc = colour;
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		if(classbounds==null || classbounds.length==0) {
			System.err.println("No bounds for histogram bins found/set!");
			return;
		}
		if(classbounds.length==1) {
			double cb = classbounds[0];
			classbounds = new double[] { Math.min(cb-1d, Math.min(0.8d*cb,1.2d*cb)), Math.max(cb+1d, Math.max(0.8d*cb,1.2d*cb)) };
		}
		Arrays.sort(classbounds);
		double[] hist = new double[classbounds.length-1];
		for(int i=0; i<hist.length; i++) hist[i] = 0d;
		for(double y: yarrayy) {
			if(Double.isNaN(y)) continue;
			if(y<classbounds[0]) continue;
			if(y>classbounds[classbounds.length-1]) continue;
			int idx = 0;
			for(idx=0; idx<hist.length; idx++)
				if(y<classbounds[idx+1]) break;
			hist[idx]+=1d;
		}
		if(normalizationMode==CUM_COUNT||normalizationMode==CUM_MAXONE)
			for(int i=1; i<hist.length; i++) hist[i] += hist[i-1];
		double hmax = JPlotMath.dmax(hist);
		if(normalizationMode==NORM_MAXONE||normalizationMode==CUM_MAXONE) {
			for(int i=0; i<hist.length; i++) hist[i] /= hmax;
			hmax = 1d;
		}
		if(normalizationMode==NORM_ARRONE) {
			double hsum = JPlotMath.dsum(hist);
			for(int i=0; i<hist.length; i++) hist[i] /= hsum;
			hmax /= hsum;
		}
		minX=classbounds[0]; maxX=classbounds[classbounds.length-1];
		minY=0d; maxY=JPlotMath.dmax(hist);
		//TODO: deal with JMultiAxis
		boolean isXlog = (ax.getScaleX() instanceof LogarithmicScale);
		boolean isYlog = (ax.getScaleY() instanceof LogarithmicScale);
		double Xin = isXlog ? Math.log10(minX) : minX, Xax = isXlog ? Math.log10(maxX) : maxX;
		double Yin = isYlog ? Math.log10(minY) : minY, Yax = isYlog ? Math.log10(maxY) : maxY;
		double xs = p[2] / (Xax - Xin), ys = p[3] / (Yax - Yin);
		AffineBuilder affine = new AffineBuilder().scale(invertAxisX ? -1d : 1d, invertAxisY ? 1d : -1d)
				.translate(invertAxisX ? Xax : -Xin, invertAxisY ? -Yin : Yax).scale(xs, ys).translate(p[0], p[1]);
		double[][] affmat = affine.getMatrix();
		for(int i=0; i<hist.length; i++) {
			if(hist[i]==0d) continue;
			double cb = classbounds[i+1]-classbounds[i];
			JDPolygon poly = new JDPolygon(
					new JDPoint(  classbounds[i]+0.1d*cb,     0d),
					new JDPoint(  classbounds[i]+0.1d*cb,hist[i]),
					new JDPoint(classbounds[i+1]-0.1d*cb,hist[i]),
					new JDPoint(classbounds[i+1]-0.1d*cb,     0d)
			);
			if(poly.area()<0d) poly.reverse_orientation();
			poly.affine(affmat);
			s.addChild(new JPolygonShape(poly, pc, pc, 0.5f, true, true));
		}
	}
}
