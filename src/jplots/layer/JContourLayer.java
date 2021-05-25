package jplots.layer;

import jplots.JAxis;
import jplots.JPlot;
import jplots.JPlotMath;
import jplots.color.JColourtable;
import jplots.shapes.JGroupShape;
import processing.core.PGraphics;

public class JContourLayer extends JPlotsLayer {

	private boolean hasContours, isFilled;
	private double minZ, maxZ;
	private double[] xarrayx, yarrayy;
	private double[][] zarrayz;
	private JColourtable colourtable;
	private double[] contourIntervals;
	private int[] fillColours;
	
	public JContourLayer(float[] x, float[] y, float[][] z, float zmin, float zmax, int nintervals, JColourtable ct, boolean drawContours, boolean filled) {
		xarrayx = new double[x.length];
		for(int i=0; i<x.length; i++) xarrayx[i] = x[i];
		yarrayy = new double[y.length];
		for(int i=0; i<y.length; i++) yarrayy[i] = y[i];
		zarrayz = new double[z.length][z[0].length];
		for(int j=0; j<z.length; j++)
			for(int i=0; i<z[j].length; i++)
				zarrayz[j][i] = z[j][i];
		minX = JPlotMath.dmin(xarrayx);
		maxX = JPlotMath.dmax(xarrayx);
		minY = JPlotMath.dmin(yarrayy);
		maxY = JPlotMath.dmax(yarrayy);
		minZ = Float.isNaN(zmin) ? JPlotMath.dmin(zarrayz) : zmin;
		maxZ = Float.isNaN(zmax) ? JPlotMath.dmax(zarrayz) : zmax;
		colourtable = ct;
		hasContours = drawContours;
		isFilled = filled;
	}
	

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		//step 1: define contour-intervals -> done in initialisation
		//step 2: detect contours in original array with indexing (level indices)
		int ni=contourIntervals.length, nx=xarrayx.length, ny=yarrayy.length;
		byte[][][] cnt = new byte[ni-1][ny-1][nx-1];
		//step 3: transform via projection
		//step 4: fill contours if with -> needs conversion to triangulare mesh
		//step 5: create (filled) contours from retrieved vectors/triangles
	}

}
