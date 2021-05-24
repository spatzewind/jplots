package pplots.layer;

import pplots.PAxis;
import pplots.PPlot;
import pplots.PPlotMath;
import pplots.color.PColourtable;
import pplots.shapes.PGroupShape;
import processing.core.PGraphics;

public class PContourLayer extends PLayer {

	private boolean isFilled;
	private double minZ, maxZ;
	private double[] xarrayx, yarrayy;
	private double[][] zarrayz;
	private PColourtable colourtable;
	
	public PContourLayer(float[] x, float[] y, float[][] z, float zmin, float zmax, PColourtable ct, boolean filled) {
		xarrayx = new double[x.length];
		for(int i=0; i<x.length; i++) xarrayx[i] = x[i];
		yarrayy = new double[y.length];
		for(int i=0; i<y.length; i++) yarrayy[i] = y[i];
		zarrayz = new double[z.length][z[0].length];
		for(int j=0; j<z.length; j++)
			for(int i=0; i<z[j].length; i++)
				zarrayz[j][i] = z[j][i];
		minX = PPlotMath.dmin(xarrayx);
		maxX = PPlotMath.dmax(xarrayx);
		minY = PPlotMath.dmin(yarrayy);
		maxY = PPlotMath.dmax(yarrayy);
		minZ = zmin;
		maxZ = zmax;
		colourtable = ct;
		isFilled = filled;
	}
	

	@Override
	public void createRasterImg(PPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(PAxis ax, int layernum, PGroupShape s) {
	}

}
