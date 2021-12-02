package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.AffineTransformation;
import jplots.JAxis;
import jplots.JPlot;
import jplots.colour.JColourtable;
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
import jplots.shapes.JPlotShape;
import jplots.shapes.JPolygonShape;
import jplots.transform.JProjection;
import processing.core.PApplet;
import processing.core.PGraphics;

public class JContourLayer extends JPlotsLayer {

	private double EPSILON  = Math.pow(2, -52);
	private boolean isFilled, pixelFilling, input2d;
	private double minZ, maxZ;
	private double[] xarrayx, yarrayy;
	private double[][] xarrayx2,yarrayy2,zarrayz;
	private double[] contourIntervals;
	private int[] startEdge;
	private String[] contourStyle;
	private List<JDPoint> corners, cntCorner;
	private List<JDEdge> edges, contours;
	private List<JDTriangle> triangles;
	private JDGeometry[] fillings;
	
	public JContourLayer(float[] x, float[] y, float[][] z, float zmin, float zmax, int nintervals, JColourtable ct, float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		float zin = Float.isNaN(zmin) ? JPlotMath.fmin(z) : zmin;
		float zax = Float.isNaN(zmax) ? JPlotMath.fmax(z) : zmax;
		float[] cntIntervals = new float[nintervals+1];
		for(int k=0; k<=nintervals; k++)
			cntIntervals[k] = zin + k*(zax-zin)/nintervals;
		input2d = false;
		JContourLayerFloat(x, y, null, null, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(float[] x, float[] y, float[][] z, float[] intervals, JColourtable ct, float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = false;
		JContourLayerFloat(x, y, null, null, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(double[] x, double[] y, double[][] z, double zmin, double zmax, int nintervals, JColourtable ct, double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		double zin = Double.isNaN(zmin) ? JPlotMath.dmin(z) : zmin;
		double zax = Double.isNaN(zmax) ? JPlotMath.dmax(z) : zmax;
		double[] cntIntervals = new double[nintervals+1];
		for(int k=0; k<=nintervals; k++)
			cntIntervals[k] = zin + k*(zax-zin)/nintervals;
		input2d = false;
		JContourLayerDouble(x, y, null, null, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(double[] x, double[] y, double[][] z, double[] intervals, JColourtable ct, double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = false;
		JContourLayerDouble(x, y, null, null, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}

	public JContourLayer(float[][] x, float[][] y, float[][] z, float zmin, float zmax, int nintervals, JColourtable ct, float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		float zin = Float.isNaN(zmin) ? JPlotMath.fmin(z) : zmin;
		float zax = Float.isNaN(zmax) ? JPlotMath.fmax(z) : zmax;
		float[] cntIntervals = new float[nintervals+1];
		for(int k=0; k<=nintervals; k++)
			cntIntervals[k] = zin + k*(zax-zin)/nintervals;
		input2d = true;
		JContourLayerFloat(null, null, x, y, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(float[][] x, float[][] y, float[][] z, float[] intervals, JColourtable ct, float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = true;
		JContourLayerFloat(null, null, x, y, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(double[][] x, double[][] y, double[][] z, double zmin, double zmax, int nintervals, JColourtable ct, double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		double zin = Double.isNaN(zmin) ? JPlotMath.dmin(z) : zmin;
		double zax = Double.isNaN(zmax) ? JPlotMath.dmax(z) : zmax;
		double[] cntIntervals = new double[nintervals+1];
		for(int k=0; k<=nintervals; k++)
			cntIntervals[k] = zin + k*(zax-zin)/nintervals;
		input2d = true;
		JContourLayerDouble(null, null, x, y, z, cntIntervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}
	public JContourLayer(double[][] x, double[][] y, double[][] z, double[] intervals, JColourtable ct, double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		input2d = true;
		JContourLayerDouble(null, null, x, y, z, intervals, ct, stroke_weight, drawContours, filled, filledAsImage);
	}

	private void JContourLayerFloat(float[] x, float[] y, float[][] x2, float[][] y2, float[][] z, float[] intervals, JColourtable ct, float stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		if(!input2d) {
			xarrayx = new double[x.length];
			for(int i=0; i<x.length; i++) xarrayx[i] = x[i];
			yarrayy = new double[y.length];
			for(int i=0; i<y.length; i++) yarrayy[i] = y[i];
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
			for(int j=0; j<x2.length; j++) for(int i=0; i<x2[0].length; i++) xarrayx2[j][i] = x2[j][i];
			yarrayy2 = new double[y2.length][y2[0].length];
			for(int j=0; j<y2.length; j++) for(int i=0; i<y2[0].length; i++) yarrayy2[j][i] = y2[j][i];
			minX = JPlotMath.dmin(xarrayx2);
			maxX = JPlotMath.dmax(xarrayx2);
			minY = JPlotMath.dmin(yarrayy2);
			maxY = JPlotMath.dmax(yarrayy2);
		}
		zarrayz = new double[z.length][z[0].length];
		for(int j=0; j<z.length; j++)
			for(int i=0; i<z[j].length; i++)
				zarrayz[j][i] = z[j][i];
		contourIntervals = new double[intervals.length];
		for(int i=0; i<intervals.length; i++)
			contourIntervals[i] = intervals[i];
		minZ = contourIntervals[0];
		maxZ = contourIntervals[intervals.length-1];
		colourtable = ct;
		lw = stroke_weight;
		drawLines = drawContours;
		isFilled = filled;
		pixelFilling = filledAsImage;
		init();
	}
	private void JContourLayerDouble(double[] x, double[] y, double[][] x2, double[][] y2, double[][] z, double[] intervals, JColourtable ct, double stroke_weight, boolean drawContours, boolean filled, boolean filledAsImage) {
		if(!input2d) {
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
		for(int i=0; i<intervals.length; i++)
			contourIntervals[i] = intervals[i];
		minZ = contourIntervals[0];
		maxZ = contourIntervals[intervals.length-1];
		colourtable = ct;
		lw = (float) stroke_weight;
		drawLines = drawContours;
		isFilled = filled;
		pixelFilling = filledAsImage;
		init();
	}
	private void init() {
		contourStyle = new String[contourIntervals.length];
		for(int cs=0; cs<contourStyle.length; cs++)
			contourStyle[cs] = "-";
		
		corners   = new ArrayList<JDPoint>();
		edges     = new ArrayList<JDEdge>();
		triangles = new ArrayList<JDTriangle>();
		
		cntCorner    = new ArrayList<JDPoint>();
		contours     = new ArrayList<JDEdge>();
		new ArrayList<JDTriangle>();
		
		fillings = null;
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		double xs = p[2]/(maxX-minX), ys = p[3]/(maxY-minY);
		double tol = Math.max(Math.abs(maxX-minX), Math.abs(maxY-minY)) * 1.0e-12d;
		
		//step 1: collect valid corners of grid and project them
		collectValidPoints(ax.getGeoProjection(), ax.getPlot().isDebug());
        
		//step 2: do delauney-triangulation
		triangulate(ax.getPlot().isDebug(), ax);
		
		//step 3: create contours
		if(drawLines || !pixelFilling) {
			createContours(ax.getPlot().isDebug());
		} else {
			if(ax.getPlot().isDebug())
				System.out.println("[DEBUG] JContourLayer: 3] no contours ...");
		}
		
		//step 4: create filling between contours if wished
		if(isFilled) {
			if(pixelFilling) {
				if(ax.getPlot().isDebug())
					System.out.println("[DEBUG] JContourLayer: 4] contour filling pixelwise ...");
				fillPixelByPixel(p, ax, xs, ys, s);
			} else {
//				fillContours(ax.getGeoProjection(), ax.getPlot().isDebug());
//				fillByPolygons(p, ax, xs, ys, s);
				try {
					fillContours(ax.getGeoProjection(), p, ax, xs, ys, s, tol);
				} catch(Exception e) {
					e.printStackTrace();
					throw new RuntimeException(e);
				}
			}
		} else {
            if(ax.getPlot().isDebug())
                System.out.println("[DEBUG] JContourLayer: 4] no filling ...");
        }
		
		//step 5: add contours to plot
		if(drawLines) {
            if(ax.getPlot().isDebug())
                System.out.println("[DEBUG] JContourLayer: 5] register "+(contours.size()-startEdge[0])+" contour line segments with "+lw+"px line width...");
            drawContourLines(p, ax, xs, ys, s);
		} else {
            if(ax.getPlot().isDebug())
                System.out.println("[DEBUG] JContourLayer: 5] contours itself will not be drawn ...");
        }
	}
	
	public double[] getZRange() {
		return new double[] {minZ, maxZ}; }
	public double[] getLevels() {
		return contourIntervals; }
	
	private void collectValidPoints(JProjection outproj, boolean debug) {
		int nx=0, ny=0;
		if(input2d) {
			ny = xarrayx2.length;
			if(yarrayy2.length!=ny)
				throw new IllegalArgumentException("x, y and z are of different shapes!");
			nx = xarrayx2.length;
			if(yarrayy2.length!=nx)
				throw new IllegalArgumentException("x, y and z are of different shapes!");
		} else {
			nx=xarrayx.length;
			ny=yarrayy.length;
		}
		corners.clear();
		double[] xy;
		for(int j=0; j<ny; j++)
			for(int i=0; i<nx; i++) {
				if(input2d) {
					xy = inputProj.fromPROJtoLATLON(xarrayx2[j][i], yarrayy2[j][i], false);
				} else {
					xy = inputProj.fromPROJtoLATLON(xarrayx[i], yarrayy[j], false);
				}
				xy = outproj.fromLATLONtoPROJ(xy[0], xy[1], false);
				if(Double.isFinite(xy[0]) && Double.isFinite(xy[1]))
					corners.add(new JDPoint(xy[0],xy[1],zarrayz[j][i]));
			}
		double xyScale = Math.max(Math.max(-minX,maxX), Math.max(-minY, maxY));
		xyScale = Math.max(xyScale, Math.max(maxX-minX, maxY-minY));
		EPSILON  =     xyScale     * 1.e-12d;
		if(debug)
            System.out.println("[DEBUG] JContourLayer: 1] has "+corners.size()+" valid sourcepoints");
	}
	private void triangulate(boolean debug, JAxis ax) {
        if(debug)
            System.out.println("[DEBUG] JContourLayer: 2] triangulate ...");
        if(input2d || ax.isGeoAxis()) {
	        JDelaunayTriangulator delaunay = new JDelaunayTriangulator(corners);
	        edges = delaunay.getEdges();
	        triangles = delaunay.getTriangles();
        } else {
        	JGridTriangulator triangulator = new JGridTriangulator(xarrayx, yarrayy, zarrayz);
        	edges = triangulator.getEdges();
        	triangles = triangulator.getTriangles();
        }
        if(debug)
            System.out.println("[DEBUG] JContourLayer:   mesh now consists of "+corners.size()+" corners, "+edges.size()+" edges and "+triangles.size()+" triangles");
	}
	private void createContours(boolean debug) {
        if(debug)
            System.out.println("[DEBUG] JContourLayer: 3] create contours from triangle mesh ...");
		startEdge = new int[contourIntervals.length];
		int nanid = -9999;
		contours.clear();
		//first exclude fillvalues
		for(int level=nanid; level<contourIntervals.length; level = Math.max(0, level+1)) {
            if(debug)
                System.out.println("[DEBUG] JContourLayer:   create contours for level "+(level<0 ? Double.NaN : contourIntervals[level]));
			for(JDTriangle ct: triangles) {
				JDPoint va = ct.a;
				JDPoint vb = ct.b;
				JDPoint vc = ct.c;
				if(level==nanid) {
					int nancode = 0;
					if(Double.isNaN(va.value)) nancode |= 1;
					if(Double.isNaN(vb.value)) nancode |= 2;
					if(Double.isNaN(vc.value)) nancode |= 4;
					switch(nancode) {
	                    case 1:
	                    case 6:
	                        double x161 = 0.5d*(va.x+vb.x);
	                        double y161 = 0.5d*(va.y+vb.y);
	                        int i16s = addCorner(x161,y161,Double.NaN,cntCorner);
	                        double x162 = 0.5d*(va.x+vc.x);
	                        double y162 = 0.5d*(va.y+vc.y);
	                        int i16m = i16s;
	                        for(int ci=0; ci<contourIntervals.length; ci++) {
	                            int ic = vb.value<vc.value ? ci : contourIntervals.length-1-ci;
	                            if((vb.value<contourIntervals[ic] && contourIntervals[ic]<vc.value) ||
	                                    (vc.value<contourIntervals[ic] && contourIntervals[ic]<vb.value)) {
	                                double f16 = (contourIntervals[ic]-vb.value) / (vc.value-vb.value);
	                                i16m = addCorner(x161+f16*(x162-x161),y161+f16*(y162-y161),Double.NaN,cntCorner);
	                                contours.add(new JDEdge(cntCorner.get(i16s), cntCorner.get(i16m)));
	                                i16s += i16m-i16s;
	                            }
	                        }
	                        int i16e = addCorner(x162,y162,Double.NaN,cntCorner);
	                        contours.add(new JDEdge(cntCorner.get(i16s), cntCorner.get(i16e)));
	                        break;
	                    case 2:
	                    case 5:
	                        double x251 = 0.5d*(vb.x+va.x);
	                        double y251 = 0.5d*(vb.y+va.y);
	                        int i25s = addCorner(x251,y251,Double.NaN,cntCorner);
	                        double x252 = 0.5d*(vb.x+vc.x);
	                        double y252 = 0.5d*(vb.y+vc.y);
	                        int i25m = i25s;
	                        for(int ci=0; ci<contourIntervals.length; ci++) {
	                            int ic = va.value<vc.value ? ci : contourIntervals.length-1-ci;
	                            if((va.value<contourIntervals[ic] && contourIntervals[ic]<vc.value) ||
	                                    (vc.value<contourIntervals[ic] && contourIntervals[ic]<va.value)) {
	                                double f25 = (contourIntervals[ic]-va.value) / (vc.value-va.value);
	                                i25m = addCorner(x251+f25*(x252-x251),y251+f25*(y252-y251),Double.NaN,cntCorner);
	                                contours.add(new JDEdge(cntCorner.get(i25s), cntCorner.get(i25m)));
	                                i25s += i25m - i25s;
	                            }
	                        }
	                        int i25e = addCorner(x252,y252,Double.NaN,cntCorner);
	                        contours.add(new JDEdge(cntCorner.get(i25s), cntCorner.get(i25e)));
	                        break;
	                    case 3:
	                    case 4:
	                        double x341 = 0.5d*(vc.x+vb.x);
	                        double y341 = 0.5d*(vc.y+vb.y);
	                        int i34s = addCorner(x341,y341,Double.NaN,cntCorner);
	                        double x342 = 0.5d*(vc.x+va.x);
	                        double y342 = 0.5d*(vc.y+va.y);
	                        for(int ci=0; ci<contourIntervals.length; ci++) {
	                            int ic = vb.value<va.value ? ci : contourIntervals.length-1-ci;
	                            if((va.value<contourIntervals[ic] && contourIntervals[ic]<vb.value) ||
	                                    (vb.value<contourIntervals[ic] && contourIntervals[ic]<va.value)) {
	                                double f34 = (contourIntervals[ic]-vb.value) / (va.value-vb.value);
	                                int i34m = addCorner(x341+f34*(x342-x341),y341+f34*(y342-y341),Double.NaN,cntCorner);
	                                contours.add(new JDEdge(cntCorner.get(i34s), cntCorner.get(i34m)));
	                                i34s += i34m-i34s;
	                            }
	                        }
	                        int i34e = addCorner(x342,y342,Double.NaN,cntCorner);
	                        contours.add(new JDEdge(cntCorner.get(i34s), cntCorner.get(i34e)));
	                        break;
						default:
							break;
					}
				} else {
					double lev = contourIntervals[level];
					double r12 = (lev-va.value)/(vb.value-va.value);
					double r23 = (lev-vb.value)/(vc.value-vb.value);
					double r31 = (lev-vc.value)/(va.value-vc.value);
					int levcode = 0;
					if(Double.isNaN(va.value)) levcode |=  2;
					if(Double.isNaN(vb.value)) levcode |=  8;
					if(Double.isNaN(vc.value)) levcode |= 32;
					if(va.value <= lev) levcode |=  1;
					if(vb.value <= lev) levcode |=  4;
					if(vc.value <= lev) levcode |= 16;
//					if(debug)
//						System.out.println("[DEBUG] JContourplot:   triangle "+t_idx+" -> levcode "+Integer.toBinaryString(levcode));
					switch(levcode) {
	                    case  1: //x000001
	                    case 20: //x010100
	                        double x161 = va.x+r12*(vb.x-va.x);
	                        double y161 = va.y+r12*(vb.y-va.y);
	                        int i17s = addCorner(x161,y161,lev,cntCorner);
	                        double x162 = vc.x+r31*(va.x-vc.x);
	                        double y162 = vc.y+r31*(va.y-vc.y);
	                        int i17e = addCorner(x162,y162,lev,cntCorner);
	                        if(i17s!=i17e)
	                            contours.add(levcode==1 ? new JDEdge(cntCorner.get(i17s), cntCorner.get(i17e)) :
	                                                      new JDEdge(cntCorner.get(i17e), cntCorner.get(i17s)));
	                        break;
	                    case  5: //x000101
	                    case 16: //x010000
	                        double x341 = vb.x+r23*(vc.x-vb.x);
	                        double y341 = vb.y+r23*(vc.y-vb.y);
	                        int i20s = addCorner(x341,y341,lev,cntCorner);
	                        double x342 = vc.x+r31*(va.x-vc.x);
	                        double y342 = vc.y+r31*(va.y-vc.y);
	                        int i20e = addCorner(x342,y342,lev,cntCorner);
	                        if(i20s!=i20e)
	                            contours.add(levcode==5 ? new JDEdge(cntCorner.get(i20s), cntCorner.get(i20e)) :
	                                                      new JDEdge(cntCorner.get(i20e), cntCorner.get(i20s)));
	                        break;
	                    case  4: //x000100
	                    case 17: //x010001
	                        double x251 = va.x+r12*(vb.x-va.x);
	                        double y251 = va.y+r12*(vb.y-va.y);
	                        int i05s = addCorner(x251,y251,lev,cntCorner);
	                        double x252 = vb.x+r23*(vc.x-vb.x);
	                        double y252 = vb.y+r23*(vc.y-vb.y);
	                        int i05e = addCorner(x252,y252,lev,cntCorner);
	                        if(i05s!=i05e)
	                            contours.add(levcode==17 ? new JDEdge(cntCorner.get(i05s), cntCorner.get(i05e)) :
	                                                       new JDEdge(cntCorner.get(i05e), cntCorner.get(i05s)));
	                        break;
	                    case  6: //x000110
	                    case 18: //x010010
	                        double x25c = vb.x+r23*(vc.x-vb.x);
	                        double y25c = vb.y+r23*(vc.y-vb.y);
	                        int i22s = addCorner(x25c,y25c,lev,cntCorner);
	                        double x23b = 0.5d*(vb.x+va.x), y23b = 0.5d*(vb.y+va.y);
	                        double x23c = 0.5d*(vc.x+va.x), y23c = 0.5d*(vc.y+va.y);
	                        int i22e = addCorner(x23b+r23*(x23c-x23b),y23b+r23*(y23c-y23b),Double.NaN,cntCorner);
	                        contours.add(levcode==6 ? new JDEdge(cntCorner.get(i22s), cntCorner.get(i22e)) :
	                                                  new JDEdge(cntCorner.get(i22e), cntCorner.get(i22s)));
	                        break;
	                    case  9: //x001001
	                    case 24: //x011000
	                        double x34c = vc.x+r31*(va.x-vc.x);
	                        double y34c = vc.y+r31*(va.y-vc.y);
	                        int i26s = addCorner(x34c,y34c,lev,cntCorner);
	                        double x31c = 0.5d*(vc.x+vb.x), y31c = 0.5d*(vc.y+vb.y);
	                        double x31a = 0.5d*(va.x+vb.x), y31a = 0.5d*(va.y+vb.y);
	                        int i26e = addCorner(x31c+r31*(x31a-x31c),y31c+r31*(y31a-y31c),Double.NaN,cntCorner);
	                        contours.add(levcode==24 ? new JDEdge(cntCorner.get(i26s), cntCorner.get(i26e)) :
	                                                   new JDEdge(cntCorner.get(i26e), cntCorner.get(i26s)));
	                        break;
	                    case 33: //x100001
	                    case 36: //x100100
	                        double x12c = va.x+r12*(vb.x-va.x);
	                        double y12c = va.y+r12*(vb.y-va.y);
	                        int i41s = addCorner(x12c,y12c,Double.NaN,cntCorner);
	                        double x12a = 0.5d*(va.x+vc.x), y12a = 0.5d*(va.y+vc.y);
	                        double x12b = 0.5d*(vb.x+vc.x), y12b = 0.5d*(vb.y+vc.y);
	                        int i41e = addCorner(x12a+r12*(x12b-x12a),y12a+r12*(y12b-y12a),Double.NaN,cntCorner);
	                        contours.add(levcode==33 ? new JDEdge(cntCorner.get(i41s), cntCorner.get(i41e)) :
	                                                   new JDEdge(cntCorner.get(i41e), cntCorner.get(i41s)));
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
			if(level+1<contourIntervals.length)
				startEdge[Math.max(0,level+1)] = contours.size();
		}
		//remove duplicates of contour line segments
//		for(int e1=contours.size()-1; e1>0; e1--)
//			for(int e2=e1-1; e2>=0; e2--)
//				if(contours.get(e1).isSame(contours.get(e2))) {
//					for(int se=0; se<startEdge.length; se++)
//						if(startEdge[se]>=e1)
//							startEdge[se]--;
//					contours.remove(e1); break;
//				}
		if(debug) {
			System.out.println("[DEBUG] JContourLayer:    estimated "+contours.size()+" line segments");
			String sts = "";
			for(int i=0; i<startEdge.length; i++)
				sts += ", "+startEdge[i];
			System.out.println("[DEBUG] JContourLayer:    start indices are: "+sts.substring(2));
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
        //double tolerance = 1.0e-10d;
        JDPoint np = new JDPoint(x,y,v);
        for(int p=0; p<points.size(); p++)
            if(points.get(p).equals(np, EPSILON))
                return p;
        points.add(np);
        return points.size()-1;
	}
	private void fillContours(JProjection outproj, int[] p, JAxis ax, double xs, double ys, JGroupShape s, double eps) {
        if(ax.getPlot().isDebug())
            System.out.println("[DEBUG] JContourLayer: 4] fill contours ...");
        double eps2 = eps*eps;
        GeometryFactory gf = new GeometryFactory();
        Geometry rahmen = gf.createPolygon(new Coordinate[] {
        		new Coordinate(p[0],p[1]),
        		new Coordinate(p[0]+p[2],p[1]),
        		new Coordinate(p[0]+p[2],p[1]+p[3]),
        		new Coordinate(p[0],p[1]+p[3]),
        		new Coordinate(p[0],p[1])
        });
        AffineTransformation affine = new AffineTransformation()
        		.setToIdentity()
        		.scale(invertAxisX?-1d:1d, invertAxisY?1d:-1d)
        		.translate(invertAxisX?maxX:-minX, invertAxisY?-minY:maxY)
        		.scale(xs, ys)
        		.translate(p[0], p[1])
        		;
        fillings = null;
        fillings = new JDGeometry[contourIntervals.length+1];
        for(int lev=-1; lev<contourIntervals.length; lev++) {
        	List<JDPolygon> pl = new ArrayList<JDPolygon>();
        	double levmin = -100000000d;
        	double levmax =  100000000d;
        	if(lev>=0) levmin = contourIntervals[lev];
        	if(lev+1<contourIntervals.length) levmax = contourIntervals[lev+1];
            if(ax.getPlot().isDebug())
                System.out.println("[DEBUG]                   ... collect geometry for range ["+levmin+" ... "+levmax+"]");
        	for(JDTriangle t: triangles) {
        		double vmin = Double.POSITIVE_INFINITY;
        		double vmax = Double.NEGATIVE_INFINITY;
        		double a = t.a.value();
        		if(Double.isFinite(a) && a<vmin) vmin = a; if(Double.isFinite(a) && a>vmax) vmax = a;
        		double b = t.b.value();
        		if(Double.isFinite(b) && b<vmin) vmin = b; if(Double.isFinite(b) && b>vmax) vmax = b;
        		double c = t.c.value();
        		if(Double.isFinite(c) && c<vmin) vmin = c; if(Double.isFinite(c) && c>vmax) vmax = c;
        		if(vmax<vmin) continue;
        		if(vmin>levmax || vmax<levmin) continue;
        		int id = ((Double.isNaN(a) ? 0 : a<levmin ? 1 : a>levmax ? 2 : 3)<<8) |
        				 ((Double.isNaN(b) ? 0 : b<levmin ? 1 : b>levmax ? 2 : 3)<<4) |
        				  (Double.isNaN(c) ? 0 : c<levmin ? 1 : c>levmax ? 2 : 3);
        		if(id==0x000 || id==0x111 || id==0x222)
        			continue;
        		if(id==0x001 || id==0x010 || id==0x100 || id==0x002 || id==0x020 || id==0x200)
        			continue;
        		if(id==0x011 || id==0x101 || id==0x110 || id==0x022 || id==0x202 || id==0x220)
        			continue;
        		Coordinate aaa = new Coordinate(t.a.x, t.a.y);
        		Coordinate bbb = new Coordinate(t.b.x, t.b.y);
        		Coordinate ccc = new Coordinate(t.c.x, t.c.y);
        		if(id==0x333) {
        			if(ax.getPlot().isDebug())
        				System.out.println("[DEBUG] triangle "+t+" seems to be fully one color (lev=["+levmin+"..."+levmax+"])");
        			Geometry res = gf.createPolygon(new Coordinate[] {aaa,bbb,ccc,aaa});
            		res = affine.transform(res).intersection(rahmen);
            		switch(res.getCoordinates().length) {
            			case 4:  addTriangle2polygonList(pl, new JDTriangle(res.getCoordinates())); break;
            			default: addPolygon2polygonList(pl, new JDPolygon(res.getCoordinates())); break;
            		}
        			continue;
        		}
//        		System.out.println("a="+a+"    b="+b+"    c="+c+
//						"    min="+levmin+"    max="+levmax);
        		int nan = ((Double.isNaN(a) ? 0 : 1)<<8) |
      				  ((Double.isNaN(b) ? 0 : 1)<<4) |
      				   (Double.isNaN(c) ? 0 : 1);
        		Coordinate abh = new Coordinate( t.a.x+0.5d*(t.b.x-t.a.x), t.a.y+0.5d*(t.b.y-t.a.y) );
        		Coordinate bch = new Coordinate( t.b.x+0.5d*(t.c.x-t.b.x), t.b.y+0.5d*(t.c.y-t.b.y) );
        		Coordinate cah = new Coordinate( t.c.x+0.5d*(t.a.x-t.c.x), t.c.y+0.5d*(t.a.y-t.c.y) );
        		Polygon nanGeom = null;
        		switch(nan) {
        			case 0x001: nanGeom = gf.createPolygon(new Coordinate[] {ccc, cah, bch, ccc}); break;
        			case 0x010: nanGeom = gf.createPolygon(new Coordinate[] {bbb, bch, abh, bbb}); break;
        			case 0x011: nanGeom = gf.createPolygon(new Coordinate[] {ccc, cah, abh, bbb, ccc}); break;
        			case 0x100: nanGeom = gf.createPolygon(new Coordinate[] {aaa, abh, cah, aaa}); break;
        			case 0x101: nanGeom = gf.createPolygon(new Coordinate[] {aaa, abh, bch, ccc, aaa}); break;
        			case 0x110: nanGeom = gf.createPolygon(new Coordinate[] {bbb, bch, cah, aaa, bbb}); break;
        			case 0x111: nanGeom = gf.createPolygon(new Coordinate[] {aaa, bbb, ccc, aaa}); break;
        			default: break;
        		}
        		if(nanGeom==null) {
        			System.err.println("[ERROR] FilledContours(Vec): In triangle "+t+" no lowGeometry could be created! (id="+
        					Integer.toHexString(id)+",nan="+Integer.toHexString(nan)+")");
        			continue;
        		}
        		if(nanGeom.getArea()<eps2) continue;
        		if((id&0xf00)==0) a = levmin;
        		if((id&0x0f0)==0) b = levmin;
        		if((id&0x00f)==0) c = levmin;
        		int imin = ((a<levmin ? 0 : 1)<<8) | ((b<levmin ? 0 : 1)<<4) | (c<levmin ? 0 : 1);
        		double  abfi = (levmin-a) / (b-a),
        				bcfi = (levmin-b) / (c-b),
        				cafi = (levmin-c) / (a-c);
        		Coordinate abi = new Coordinate( t.a.x+abfi*(t.b.x-t.a.x), t.a.y+abfi*(t.b.y-t.a.y) );
        		Coordinate bci = new Coordinate( t.b.x+bcfi*(t.c.x-t.b.x), t.b.y+bcfi*(t.c.y-t.b.y) );
        		Coordinate cai = new Coordinate( t.c.x+cafi*(t.a.x-t.c.x), t.c.y+cafi*(t.a.y-t.c.y) );
        		Polygon lowGeom = null;
        		switch(imin) {
	    			case 0x001: lowGeom = gf.createPolygon(new Coordinate[] {ccc, cai, bci, ccc}); break;
	    			case 0x010: lowGeom = gf.createPolygon(new Coordinate[] {bbb, bci, abi, bbb}); break;
	    			case 0x011: lowGeom = gf.createPolygon(new Coordinate[] {ccc, cai, abi, bbb, ccc}); break;
	    			case 0x100: lowGeom = gf.createPolygon(new Coordinate[] {aaa, abi, cai, aaa}); break;
	    			case 0x101: lowGeom = gf.createPolygon(new Coordinate[] {aaa, abi, bci, ccc, aaa}); break;
	    			case 0x110: lowGeom = gf.createPolygon(new Coordinate[] {bbb, bci, cai, aaa, bbb}); break;
        			case 0x111: lowGeom = gf.createPolygon(new Coordinate[] {aaa, bbb, ccc, aaa}); break;
	    			default: break;
        		}
        		if(lowGeom==null) {
        			System.err.println("[ERROR] FilledContours(Vec): In triangle "+t+" no lowGeometry could be created! (id="+
        					Integer.toHexString(id)+",nan="+Integer.toHexString(nan)+",min="+Integer.toHexString(imin)+")");
        			continue;
        		}
        		if(lowGeom.getArea()<eps2) continue;
        		if((id&0xf00)==0) a = levmax;
        		if((id&0x0f0)==0) b = levmax;
        		if((id&0x00f)==0) c = levmax;
        		int imax = ((a>levmax ? 0 : 1)<<8) | ((b>levmax ? 0 : 1)<<4) | (c>levmax ? 0 : 1);
        		double  abfa = (levmax-a) / (b-a),
        				bcfa = (levmax-b) / (c-b),
        				cafa = (levmax-c) / (a-c);
        		Coordinate aba = new Coordinate( t.a.x+abfa*(t.b.x-t.a.x), t.a.y+abfa*(t.b.y-t.a.y) );
        		Coordinate bca = new Coordinate( t.b.x+bcfa*(t.c.x-t.b.x), t.b.y+bcfa*(t.c.y-t.b.y) );
        		Coordinate caa = new Coordinate( t.c.x+cafa*(t.a.x-t.c.x), t.c.y+cafa*(t.a.y-t.c.y) );
        		Polygon higGeom = null;
        		switch(imax) {
	    			case 0x001: higGeom = gf.createPolygon(new Coordinate[] {ccc, caa, bca, ccc}); break;
	    			case 0x010: higGeom = gf.createPolygon(new Coordinate[] {bbb, bca, aba, bbb}); break;
	    			case 0x011: higGeom = gf.createPolygon(new Coordinate[] {ccc, caa, aba, bbb, ccc}); break;
	    			case 0x100: higGeom = gf.createPolygon(new Coordinate[] {aaa, aba, caa, aaa}); break;
	    			case 0x101: higGeom = gf.createPolygon(new Coordinate[] {aaa, aba, bca, ccc, aaa}); break;
	    			case 0x110: higGeom = gf.createPolygon(new Coordinate[] {bbb, bca, caa, aaa, bbb}); break;
        			case 0x111: higGeom = gf.createPolygon(new Coordinate[] {aaa, bbb, ccc, aaa}); break;
	    			default: break;
        		}
        		if(higGeom==null) {
        			System.err.println("[ERROR] FilledContours(Vec): In triangle "+t+" no highGeometry could be created! (id="+
        					Integer.toHexString(id)+",nan="+Integer.toHexString(nan)+",min="+Integer.toHexString(imin)+",max="+
        					Integer.toHexString(imax)+")");
        			continue;
        		}
        		if(higGeom.getArea()<eps2) continue;
//        		System.out.println("minfraction{ab="+abfi+",bc="+bcfi+",ca="+cafi+"}  maxfraction{ab="+abfa+",bc="+bcfa+",ca="+cafa+"}");
        		Geometry res = nanGeom.intersection(lowGeom.intersection(higGeom));
        		res = affine.transform(res).intersection(rahmen);
        		if(res.getArea()<eps2)
        			continue;
//        		System.out.println("found "+printGeom(res)+" from\n"+
//									"        "+printGeom(nanGeom)+"  (nan="+Integer.toHexString(nan)+")"+
//									"\n      /\\"+printGeom(lowGeom)+"  (low="+Integer.toHexString(imin)+")"+
//									"\n      /\\"+printGeom(higGeom)+"  (hig="+Integer.toHexString(imax)+")");
        		switch(res.getCoordinates().length) {
        			case 4: addTriangle2polygonList(pl, new JDTriangle(res.getCoordinates())); break;
        			default: addPolygon2polygonList(pl, new JDPolygon(res.getCoordinates())); break;
        		}
        	}
        	fillings[lev+1] = new JDGeometry();
        	for(JDPolygon jdp: pl) {
        		if(ax.getPlot().isDebug())
        			System.out.println("[DEBUG]                       ... use polygon "+jdp);
        		fillings[lev+1].add(jdp);
        	}
    		if(ax.getPlot().isDebug() && pl.isEmpty())
    			System.out.println("[DEBUG]                       ... no polygon created");
        }
        
		//add corner points to fill last bits in normal plot // ussually not necessary for complex geographical projections
//		int nx=xarrayx.length, ny=yarrayy.length;
//		double[] xy00 = inputProj.fromPROJtoLATLON(xarrayx[0], yarrayy[0], false);
//		xy00 = outproj.fromLATLONtoPROJ(xy00[0], xy00[1], false);
//		if(Double.isFinite(xy00[0]) && Double.isFinite(xy00[1]))
//			addCorner(xy00[0], xy00[1], zarrayz[0][0], cntCorner);
//		double[] xy01 = inputProj.fromPROJtoLATLON(xarrayx[nx-1], yarrayy[0], false);
//		xy01 = outproj.fromLATLONtoPROJ(xy01[0], xy01[1], false);
//		if(Double.isFinite(xy01[0]) && Double.isFinite(xy01[1]))
//			addCorner(xy01[0], xy01[1], zarrayz[0][nx-1], cntCorner);
//		double[] xy10 = inputProj.fromPROJtoLATLON(xarrayx[0], yarrayy[ny-1], false);
//		xy10 = outproj.fromLATLONtoPROJ(xy10[0], xy10[1], false);
//		if(Double.isFinite(xy10[0]) && Double.isFinite(xy10[1]))
//			addCorner(xy10[0], xy10[1], zarrayz[ny-1][0], cntCorner);
//		double[] xy11 = inputProj.fromPROJtoLATLON(xarrayx[nx-1], yarrayy[ny-1], false);
//		xy11 = outproj.fromLATLONtoPROJ(xy11[0], xy11[1], false);
//		if(Double.isFinite(xy11[0]) && Double.isFinite(xy11[1]))
//			addCorner(xy11[0], xy11[1], zarrayz[ny-1][nx-1], cntCorner);
//		if(debug)
//			System.out.println("[DEBUG] JContourLayer:    added corner points.");
//		
//		//create delaunay triangulation
//		JConstrainedDelaunayTriangulator delaunayFill = new JConstrainedDelaunayTriangulator(contours, cntCorner);
//		//constrain edges afterwards to respect contours
//		//delaunayFill.constrain(contours);
//		if(debug)
//			System.out.println("[DEBUG] JContourLayer:    CDT successful.");
//		
//		//order triangles with respect to contour level (or by x-axis amount when level is equal)
//		cntTriangles = delaunayFill.getTriangles();
//		List<JDTriangle> equalLevelTriangle = new ArrayList<JDTriangle>();
//		double[] ci = new double[contourIntervals.length];
//		ci[0] = contourIntervals[0]*(contourIntervals[0]<0d ? 1.00000001d : 0.99999999d);
//		for(int c=1; c<ci.length; c++)
//				ci[c] = 0.00000001d*contourIntervals[c-1]+0.99999999d*contourIntervals[c];
//		for(JDTriangle cntTri: cntTriangles) {
//			int[] l =  {getLevel(cntTri.a.value, contourIntervals, -9999),
//						getLevel(cntTri.b.value, contourIntervals, -9999),
//						getLevel(cntTri.c.value, contourIntervals, -9999)};
//			if(l[0]<0 || l[1]<0 || l[2]<0) {
//				cntTri.lev = -1;
//				equalLevelTriangle.add(cntTri);
//			} else
//			if(l[0]==l[1] && l[1]==l[2]) {
//				cntTri.lev = l[0];
//				equalLevelTriangle.add(cntTri);
//			} else {
//				cntTri.lev = (l[0]+l[1]+l[2])/3;
//			}
//		}
//		if(debug)
//			System.out.println("[DEBUG] JContourLayer:    Simple colouring done.");
//		int iter=0;
//		while(equalLevelTriangle.size()>0 && iter<1000) {
//			iter++;
//			for(int t=equalLevelTriangle.size()-1; t>=0; t--) {
//				JDTriangle cntTri = equalLevelTriangle.get(t);
//				for(JDEdge e: new JDEdge[] {cntTri.ab, cntTri.bc, cntTri.ca}) {
//					if(e.has2triangles()) {
//						JDTriangle[] w = e.getWing();
//						JDTriangle next = w[0].equals(cntTri) ? w[1] : w[0];
//						if(next==null)
//							continue;
//						if(equalLevelTriangle.contains(next))
//							continue;
//						if(contours.contains(e)) {
//							double ev = 0.5d*(e.a.value+e.b.value);
//							if(Double.isNaN(e.a.value)) ev = e.b.value;
//							if(Double.isNaN(e.b.value)) ev = e.a.value;
//							int elev = getLevel(ev, contourIntervals, -1);
//							if(elev<0)
//								continue;
//							if(next.lev>=0) {
//								cntTri.lev = 2*elev - 1 - next.lev;
//							} else {
//								cntTri.lev = -1;
//							}
//						} else {
//							cntTri.lev = next.lev;
//						}
//						equalLevelTriangle.remove(t);
//						break;
//					}
//				}
//			}
//		}
//		if(debug)
//			System.out.println("[DEBUG] JContourLayer:    colors to triangle assigned.");
//		cntTriangles.sort(new Comparator<JDTriangle>() {
//			@Override
//			public int compare(JDTriangle o1, JDTriangle o2) {
//				if(o1==null)
//					return -1;
//				if(o2==null)
//					return  1;
//				if(o1.lev!=o2.lev)
//					return o1.lev - o2.lev;
//				return Double.compare(o1.a.x+o1.b.x+o1.c.x, o2.a.x+o2.b.x+o2.c.x);
//			}
//		});
//		startTriangle = new int[contourIntervals.length+1];
//		int lastLev = -1;
//		for(int c=1; c<cntTriangles.size(); c++) {
//			if(cntTriangles.get(c-1).lev<cntTriangles.get(c).lev) {
//				lastLev = Math.max(lastLev, cntTriangles.get(c).lev);
//				startTriangle[cntTriangles.get(c).lev] = c;
//			}
//		}
//		for(int c=1; c<startTriangle.length; c++)
//			if(startTriangle[c]<startTriangle[c-1]) {
//				int bigger = c+1;
//				for(bigger=c+1; bigger<startTriangle.length; bigger++)
//					if(startTriangle[bigger]>startTriangle[c-1])
//						break;
//				int bc = cntTriangles.size();
//				if(bigger<startTriangle.length) {
//					bc = startTriangle[bigger];
//				} else {
//					bigger = startTriangle.length-1;
//				}
//				for(bigger=bigger+0; bigger>=c; bigger--)
//					startTriangle[bigger] = bc;
//			}
////		for(int c=lastLev+1; c<startTriangle.length; c++)
////			startTriangle[c] = cntTriangles.size();
//		if(debug) {
//			System.out.println("[DEBUG] JContourLayer:    estimated "+cntTriangles.size()+" triangles between contourlines");
//			String sts = "";
//			for(int i=0; i<startTriangle.length; i++)
//				sts += ", "+startTriangle[i];
//			System.out.println("[DEBUG] JContourLayer:    start indices are: "+sts.substring(2));
//		}

		JGroupShape trianglesh = new JGroupShape();
        for(int lev=-1; lev<contourIntervals.length; lev++) {
        	double levmin = Double.NEGATIVE_INFINITY;
        	double levmax = Double.POSITIVE_INFINITY;
        	if(lev>=0) levmin = contourIntervals[lev];
        	if(lev+1<contourIntervals.length) levmax = contourIntervals[lev+1];
            if(ax.getPlot().isDebug())
                System.out.println("[DEBUG]                   ... draw geometry for range ["+levmin+" ... "+levmax+"]");
        	if(fillings[lev+1]==null)
        		continue;
			double pct = minZ-10d;
			if(lev>=0 && lev+1<contourIntervals.length)
				pct = 0.5d*(contourIntervals[lev]+contourIntervals[lev+1]);
			if(lev+1==contourIntervals.length)
				pct = maxZ + 10d;
			int cct = colourtable.getColour(pct, contourIntervals[0], contourIntervals[contourIntervals.length-1]);
//			JPlotShape.fill(cct);
//			JPlotShape.noStroke();
//			if(ax.getPlot().isDebug()) {
//				JPlotShape.stroke(0xff999999); JPlotShape.strokeWeight(2f); }
			for(JDPolygon ppp: fillings[lev+1].getPolygons()) {
				trianglesh.addChild(new JPolygonShape(ppp, cct, 0xff999999, 2f, true, ax.getPlot().isDebug()));
//				Geometry triGeom = polygons.getGeometryN(gi);
//				if(triGeom.getCoordinates().length<3)
//					continue;
//				
//				Coordinate tv1 = triGeom.getCoordinates()[0];
//				Coordinate tv2 = triGeom.getCoordinates()[1];
//				Coordinate tv3 = triGeom.getCoordinates()[2];
//				System.out.println(tv1);
//				System.out.println(tv2);
//				System.out.println(tv3);
//				double x1 = p[0]+xs*(invertAxisX ? maxX-tv1.x : tv1.x-minX);
//				double x2 = p[0]+xs*(invertAxisX ? maxX-tv2.x : tv2.x-minX);
//				double x3 = p[0]+xs*(invertAxisX ? maxX-tv3.x : tv3.x-minX);
//				double y1 = p[1]+ys*(invertAxisY ? tv1.y-minY : maxY-tv1.y);
//				double y2 = p[1]+ys*(invertAxisY ? tv2.y-minY : maxY-tv2.y);
//				double y3 = p[1]+ys*(invertAxisY ? tv3.y-minY : maxY-tv3.y);
//				for(JDTriangle ttt: cutoff(x1,y1,x2,y2,x3,y3, p[0],p[1],p[0]+p[2],p[1]+p[3])) {
//					trianglesh.addChild(
//							new JTriangleShape((float)ttt.a.x, (float)ttt.a.y, (float)ttt.b.x, (float)ttt.b.y, (float)ttt.c.x, (float)ttt.c.y)
//					);
//				}
			}
        }
		s.addChild(trianglesh);
	}
	private void fillPixelByPixel(int[] p, JAxis ax, double xs, double ys, JGroupShape s) {
		//double us = srcImg.width/(srcExt[2]-srcExt[0]), vs = srcImg.height/(srcExt[3]-srcExt[1]);
		if(img==null) {
			img = ax.getPlot().getApplet().createImage(p[2], p[3], PApplet.ARGB);
		} else if(img.width!=p[2] || img.height!=p[3]) {
			img = ax.getPlot().getApplet().createImage(p[2], p[3], PApplet.ARGB);
		}
		img.loadPixels();
		for(JDTriangle tri: triangles) {
			JDPoint va = tri.a;
			JDPoint vb = tri.b;
			JDPoint vc = tri.c;
			double x1 = xs*(invertAxisX ? maxX-va.x : va.x-minX);
			double x2 = xs*(invertAxisX ? maxX-vb.x : vb.x-minX);
			double x3 = xs*(invertAxisX ? maxX-vc.x : vc.x-minX);
			double y1 = ys*(invertAxisY ? va.y-minY : maxY-va.y);
			double y2 = ys*(invertAxisY ? vb.y-minY : maxY-vb.y);
			double y3 = ys*(invertAxisY ? vc.y-minY : maxY-vc.y);
			double txi = Math.min(x1,Math.min(x2,x3)), txa = Math.max(x1,Math.max(x2,x3));
			double tyi = Math.min(y1,Math.min(y2,y3)), tya = Math.max(y1,Math.max(y2,y3));
			int ixs = Math.max((int) txi - (txi<0 ? 1 : 0), 0),
				ixe = -Math.max((int) (-txa) - (txa>0 ? 1 : 0), 1-p[2]);
			int iys = Math.max((int) tyi - (tyi<0 ? 1 : 0), 0),
				iye = -Math.max((int) (-tya) - (tya>0 ? 1 : 0), 1-p[3]);
			if(ixe<ixs) continue;
			if(iye<iys) continue;
			double minCI = contourIntervals[0]<1d ? contourIntervals[0]*2d : contourIntervals[0]-10d;
			double maxCI = contourIntervals[contourIntervals.length-1]>1d ? contourIntervals[contourIntervals.length-1]*2d : contourIntervals[contourIntervals.length-1]+10d;
			for(int j=iys; j<=iye; j++)
				for(int i=ixs; i<=ixe; i++) {
					double vx = i+0.5d, vy = j+0.5d;
					double det     = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3);
				    double lambda1 = ((y2-y3)*(vx-x3) + (x3-x2)*(vy-y3)) / det;
				    if(lambda1<0d || lambda1>1d) continue;
					double lambda2 = ((y3-y1)*(vx-x3) + (x1-x3)*(vy-y3)) / det;
				    if(lambda2<0d || lambda2>1d) continue;
					double lambda3 = 1d - lambda1 - lambda2;
				    if(lambda3<0d || lambda3>1d) continue;
					int nancode = 0;
					if(Double.isNaN(va.value)) nancode |= 1;
					if(Double.isNaN(vb.value)) nancode |= 2;
					if(Double.isNaN(vc.value)) nancode |= 4;
					double val = 0d;
					switch(nancode) {
						case 0:
							val = lambda1*va.value + lambda2*vb.value + lambda3*vc.value;
							break;
						case 1:
						case 2:
						case 4:
							if(nancode==1) val = lambda1<=0.5d ? (lambda2*vb.value + lambda3*vc.value) / (lambda2+lambda3) : Double.NaN;
							if(nancode==2) val = lambda2<=0.5d ? (lambda1*va.value + lambda3*vc.value) / (lambda1+lambda3) : Double.NaN;
							if(nancode==4) val = lambda3<=0.5d ? (lambda1*va.value + lambda2*vb.value) / (lambda1+lambda2) : Double.NaN;
							break;
						case 3:
						case 5:
						case 6:
							if(nancode==3) val = lambda3>=0.5d ? vc.value : Double.NaN;
							if(nancode==5) val = lambda2>=0.5d ? vb.value : Double.NaN;
							if(nancode==6) val = lambda1>=0.5d ? va.value : Double.NaN;
							break;
						default:
							val = Double.NaN;
							break;
					}
					if(Double.isNaN(val))
						continue;
					int il = getLevel(val, contourIntervals, -1);
					double dl = il<1 ? minCI : il>contourIntervals.length-1 ? maxCI : 0.5d*(contourIntervals[il-1]+contourIntervals[il]);
					img.pixels[j*p[2]+i] = colourtable.getColour(dl, contourIntervals[0], contourIntervals[contourIntervals.length-1]);
				}
		}
		img.updatePixels();
		s.addChild(new JImageShape(img, p[0], p[1], p[2], p[3]));
	}
	private void drawContourLines(int[] p, JAxis ax, double xs, double ys, JGroupShape s) {
        JGroupShape linesh = new JGroupShape();
		for(int c=0; c<contourIntervals.length; c++) {
			int cs = startEdge[c];
			int ce = contours.size();
			if(c+1<contourIntervals.length)
				ce = startEdge[c+1];
			//TODO edit colour and linestyle info for different contour lines
			int lc = 0xff000000;
			//String ls = "-"; //TODO recreate line-drawing with different linestyles
			JPlotShape.stroke(lc);
			JPlotShape.strokeWeight((float)lw);
			String ls = contourStyle[c];
			double lln=1d, llf=0d, lpn=0d, lpf=0d, loff = 0d;
			if("-".equals(ls)) { lln=1000*lw; llf=0; lpn=0; lpf=0; }
			if(".".equals(ls)) { lln=0; llf=0; lpn=1*lw; lpf=3*lw; }
			if(",".equals(ls)) { lln=8*lw; llf=7*lw; lpn=0; lpf=0; }
			if(";".equals(ls)) { lln=8*lw; llf=3*lw; lpn=1*lw; lpf=3*lw; }
			int li = 0;
			for(int cl=cs; cl<ce; cl++) {
				JDPoint lvs = contours.get(cl).a;
				JDPoint lve = contours.get(cl).b;
				double x1 = p[0]+xs*(invertAxisX ? maxX-lvs.x : lvs.x-minX);
				double x2 = p[0]+xs*(invertAxisX ? maxX-lve.x : lve.x-minX);
				double y1 = p[1]+ys*(invertAxisY ? lvs.y-minY : maxY-lvs.y);
				double y2 = p[1]+ys*(invertAxisY ? lve.y-minY : maxY-lve.y);
				double dx = x2-x1, dy = y2-y1;
				double l = Math.sqrt(dx*dx+dy*dy);
				dx /= l; dy /= l;
				double lpos = 0d, ldif = 0d;
				while(lpos<l) {
					ldif = 0d;
					switch(li) {
						case 0: if(lln==0d) break; ldif = Math.min(l-lpos, lln-loff); break;
						case 1: if(llf==0d) break; ldif = Math.min(l-lpos, llf-loff); break;
						case 2: if(lpn==0d) break; ldif = Math.min(l-lpos, lpn-loff); break;
						case 3: if(lpf==0d) break; ldif = Math.min(l-lpos, lpf-loff); break;
					}
					float xf1 = (float)(x1+lpos*dx), yf1 = (float)(y1+lpos*dy),
							  xf2 = (float)(x1+(lpos+ldif)*dx), yf2 = (float)(y1+(lpos+ldif)*dy);
					if(xf1<p[0]      && xf2>=p[0])      { yf1 = JPlotMath.flerp(p[0],      xf1, xf2, yf1, yf2); }
					if(xf1>p[0]+p[2] && xf2<=p[0]+p[2]) { yf1 = JPlotMath.flerp(p[0]+p[2], xf1, xf2, yf1, yf2); }
					if(xf2<p[0]      && xf1>=p[0])      { yf2 = JPlotMath.flerp(p[0],      xf1, xf2, yf1, yf2); }
					if(xf2>p[0]+p[2] && xf1<=p[0]+p[2]) { yf2 = JPlotMath.flerp(p[0]+p[2], xf1, xf2, yf1, yf2); }
					
					if(yf1<p[1]      && yf2>=p[1])      { xf1 = JPlotMath.flerp(p[1],      yf1, yf2, xf1, xf2); }
					if(yf1>p[1]+p[3] && yf2<=p[1]+p[3]) { xf1 = JPlotMath.flerp(p[1]+p[3], yf1, yf2, xf1, xf2); }
					if(yf2<p[1]      && yf1>=p[1])      { xf2 = JPlotMath.flerp(p[1],      yf1, yf2, xf1, xf2); }
					if(yf2>p[1]+p[3] && yf1<=p[1]+p[3]) { xf2 = JPlotMath.flerp(p[1]+p[3], yf1, yf2, xf1, xf2); }
					if(xf1>=p[0] && xf1<=p[0]+p[2] && xf2>=p[0] && xf2<=p[0]+p[2] &&
							yf1>=p[1] && yf1<=p[1]+p[3] && yf2>=p[1] && yf2<=p[1]+p[3]) {
						if(li%2==0 && ldif>0d)
							linesh.addChild(new JLineShape(xf1,yf1,xf2,yf2));
						loff += ldif;
						switch(li) {
							case 0: if(loff>=lln) {loff -= lln; li=1; } break;
							case 1: if(loff>=llf) {loff -= llf; li=2; } break;
							case 2: if(loff>=lpn) {loff -= lpn; li=3; } break;
							case 3: if(loff>=lpf) {loff -= lpf; li=0; } break;
						}
					}
					lpos += ldif;
				}
			}
		}
		s.addChild(linesh);
	}

	public static int getLevel(double value, double[] intervalBorders, int nanLev) {
		if(Double.isNaN(value))
			return nanLev;
		int l = 0;
		for(int cl=0; cl<intervalBorders.length; cl++)
			if(intervalBorders[cl]<value)
				l = cl+1;
		return l;
	}
	private void addTriangle2polygonList(List<JDPolygon> list, JDTriangle tri) {
		if(list.isEmpty()) {
			list.add(new JDPolygon(tri.a, tri.b, tri.c));
		} else {
			boolean failed = true;
			for(int l=0; l<list.size() && failed; l++)
				failed = !list.get(l).union(tri, 1.0e-9d);
			if(failed)
				list.add(new JDPolygon(tri.a, tri.b, tri.c));
		}
	}
	private void addPolygon2polygonList(List<JDPolygon> list, JDPolygon poly) {
		if(list.isEmpty()) {
			list.add(poly.copy());
		} else {
			boolean failed = true;
			for(int l=0; l<list.size() && failed; l++)
				failed = !list.get(l).union(poly, 1.0e-9d);
			if(failed)
				list.add(poly.copy());
		}
	}
}
