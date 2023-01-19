package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import jplots.JAxis;
import jplots.JPlot;
import jplots.colour.JColourtable;
import jplots.helper.FormatTools;
import jplots.maths.AffineBuilder;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JPolygonShape;
import processing.core.PGraphics;

public class JPColourLayer extends JPlotsLayer {

//	private double EPSILON = Math.pow(2, -52);
	private double minZ, maxZ, Xin, Xax, Yin, Yax;
	private double[][] xCoord,yCoord,zfield;
	
	public JPColourLayer(float[] x, float[] y, float[][] z, float zmin, float zmax) {
		boolean is_valid = (x!=null && y!=null && z!=null);
		if(is_valid) {
			is_valid = (y.length == z.length || y.length-1==z.length);
			for(int j=0; j<z.length && is_valid; j++)
				is_valid = (x.length == z[j].length  || x.length-1==z[j].length);
		}
		if(!is_valid) {
			null_init();
			return;
		}
		xCoord = new double[y.length][x.length];
		yCoord = new double[y.length][x.length];
		for(int j=0; j<y.length; j++)
			for(int i=0; i<x.length; i++) {
				xCoord[j][i] = x[i];
				yCoord[j][i] = y[j];
			}
		zfield = JPlotMath.toDoubleArray2D(z);
		minZ = zmin;
		maxZ = zmax;
		init();
	}
	public JPColourLayer(float[][] x, float[][] y, float[][] z, float zmin, float zmax) {
		boolean is_valid = (x!=null && y!=null && z!=null);
		if(is_valid) {
			is_valid = ((x.length==z.length || x.length-1==z.length) && (y.length==z.length || y.length-1==z.length));
			int inner_length = 0; if(is_valid && x.length>0) inner_length = x[0].length;
			int in_z_len = 0; if(is_valid && z.length>0) in_z_len = z[0].length;
			is_valid = (inner_length==in_z_len || inner_length-1==in_z_len);
			for(int j=0; j<x.length && is_valid; j++) {
				if(j<z.length) is_valid = (z[j].length==in_z_len);
				if(!is_valid) break;
				is_valid = (x[j].length==inner_length && y[j].length==inner_length);
			}
		}
		if(!is_valid) {
			null_init();
			return;
		}
		xCoord = new double[y.length][x[0].length];
		yCoord = new double[y.length][x[0].length];
		for(int j=0; j<y.length; j++)
			for(int i=0; i<x[0].length; i++) {
				xCoord[j][i] = x[j][i];
				yCoord[j][i] = y[j][i];
			}
		zfield = JPlotMath.toDoubleArray2D(z);
		minZ = zmin;
		maxZ = zmax;
		init();
	}
	public JPColourLayer(double[] x, double[] y, double[][] z, double zmin, double zmax) {
		boolean is_valid = (x!=null && y!=null && z!=null);
		if(is_valid) {
			is_valid = (y.length == z.length || y.length-1==z.length);
			for(int j=0; j<z.length && is_valid; j++)
				is_valid = (x.length == z[j].length  || x.length-1==z[j].length);
		}
		if(!is_valid) {
			null_init();
			return;
		}
		xCoord = new double[y.length][x.length];
		yCoord = new double[y.length][x.length];
		for(int j=0; j<y.length; j++)
			for(int i=0; i<x.length; i++) {
				xCoord[j][i] = x[i];
				yCoord[j][i] = y[j];
			}
		zfield = z;
		minZ = zmin;
		maxZ = zmax;
		init();
	}
	public JPColourLayer(double[][] x, double[][] y, double[][] z, double zmin, double zmax) {
		boolean is_valid = (x!=null && y!=null && z!=null);
		if(is_valid) {
			is_valid = ((x.length==z.length || x.length-1==z.length) && (y.length==z.length || y.length-1==z.length));
			int inner_length = 0; if(is_valid && x.length>0) inner_length = x[0].length;
			int in_z_len = 0; if(is_valid && z.length>0) in_z_len = z[0].length;
			is_valid = (inner_length==in_z_len || inner_length-1==in_z_len);
			for(int j=0; j<x.length && is_valid; j++) {
				if(j<z.length) is_valid = (z[j].length==in_z_len);
				if(!is_valid) break;
				is_valid = (x[j].length==inner_length && y[j].length==inner_length);
			}
		}
		if(!is_valid) {
			null_init();
			return;
		}
		xCoord = new double[y.length][x[0].length];
		yCoord = new double[y.length][x[0].length];
		for(int j=0; j<y.length; j++)
			for(int i=0; i<x[0].length; i++) {
				xCoord[j][i] = x[j][i];
				yCoord[j][i] = y[j][i];
			}
		zfield = z;
		minZ = zmin;
		maxZ = zmax;
		init();
	}
	
	
	
	
	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int jlen = xCoord.length;
		int ilen = xCoord[0].length;
		int[] p = ax.getSize();
		boolean debug = ax.getPlot().isDebug();
		
		if(debug) {
			System.out.println("[JPColourLayer] ilen/jlen = "+ilen+"/"+jlen);
			System.out.println("inputX:"); FormatTools.printMat(System.out, xCoord);
			System.out.println("inputY:"); FormatTools.printMat(System.out, yCoord);
			System.out.println("inputZ:"); FormatTools.printMat(System.out, zfield);
		}
		
		JDPoint[][] cnt2 = new JDPoint[jlen][ilen];
		for(int j=0; j<jlen; j++)
			for(int i=0; i<ilen; i++)
				cnt2[j][i] = new JDPoint(ax.isXlogAxis()?Math.log10(xCoord[j][i]):xCoord[j][i],
										 ax.isYlogAxis()?Math.log10(yCoord[j][i]):yCoord[j][i]);
		if(ax.isGeoAxis() && inputProj!=null) {
			for(int j=0; j<jlen; j++)
				for(int i=0; i<ilen; i++) {
					double[] xy = inputProj.fromPROJtoLATLON(cnt2[j][i].x, cnt2[j][i].y, false, false);
					xy = ax.getGeoProjection().fromLATLONtoPROJ(xy[0], xy[1], false, false);
					cnt2[j][i].x = xy[0];
					cnt2[j][i].y = xy[1];
				}
		}
		Xin = ax.isXlogAxis() ? Math.log10(minX) : minX;
		Xax = ax.isXlogAxis() ? Math.log10(maxX) : maxX;
		Yin = ax.isYlogAxis() ? Math.log10(minY) : minY;
		Yax = ax.isYlogAxis() ? Math.log10(maxY) : maxY;
		double xs = p[2] / (Xax - Xin), ys = p[3] / (Yax - Yin);
		// double tol = Math.max(Math.abs(maxX-minX), Math.abs(maxY-minY)) * 1.0e-12d;
		AffineBuilder affine = new AffineBuilder().scale(invertAxisX ? -1d : 1d, invertAxisY ? 1d : -1d)
				.translate(invertAxisX ? Xax : -Xin, invertAxisY ? -Yin : Yax).scale(xs, ys).translate(p[0], p[1]);
		
		// step 1: create polygons
		List<JDPolygon> p2list = new ArrayList<>();
		List<JDPolygon> p3list = new ArrayList<>();
		JGroupShape polys = new JGroupShape();
		for(int j=0; j<jlen-1; j++)
			for(int i=0; i<ilen-1; i++) {
				JDPolygon p1 = new JDPolygon(cnt2[j][i],cnt2[j+1][i],cnt2[j+1][i+1],cnt2[j][i+1]);
//				if(p1.area()<0d) p1.reverse_orientation();
				p1.affine(affine.getMatrix());
				p2list.clear();
				p2list.addAll(ax.getGeoProjection().splitByMapBorder(p1));
				p3list.clear();
				for(JDPolygon p2: p2list)
					p3list.addAll(p2.intersectsAABB(p[0], p[0]+p[2], p[1], p[1]+p[3]));
				int colour = colourtable.getColour(zfield[j][i], minZ, maxZ);
				for(JDPolygon p3: p3list)
					polys.addChild(new JPolygonShape(p3, colour, colour, 1f, true, false));
			}
		if(debug) System.out.println("[JPColourLayer] draw "+polys.childCount()+" polygons");
		s.addChild(polys);
	}

	//* **************************************** *
	//* ********** GETTER AND SETTER  ********** *
	//* **************************************** *
	
	public double[] getZRange() {
		return new double[] { minZ, maxZ };
	}
	
	@Override
	public JColourtable getColourtable() {
		return super.getColourtable();
	}
	
	//* **************************************** *
	//* ********** PUBLIC METHODS     ********** *
	//* **************************************** *
	
	
	
	
	//* **************************************** *
	//* ********** PRIVATE METHODS    ********** *
	//* **************************************** *
	
	private void null_init() {
		xCoord = new double[][] { { 0d } };
		yCoord = new double[][] { { 0d } };
		zfield = new double[0][0];
		//TODO null init
	}
	private void init() {
		minX = JPlotMath.dmin(xCoord);
		maxX = JPlotMath.dmax(xCoord);
		minY = JPlotMath.dmin(yCoord);
		maxY = JPlotMath.dmax(yCoord);
	}
}
