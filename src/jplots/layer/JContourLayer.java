package jplots.layer;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import jplots.JAxis;
import jplots.JPlot;
import jplots.colour.JColourtable;
import jplots.maths.JConstrainedDelaunayTriangulator;
import jplots.maths.JDEdge;
import jplots.maths.JDPoint;
import jplots.maths.JDTriangle;
import jplots.maths.JDelaunayTriangulator;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JImageShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JTriangleShape;
import jplots.transform.JProjection;
import processing.core.PApplet;
import processing.core.PGraphics;

public class JContourLayer extends JPlotsLayer {

	private double EPSILON  = Math.pow(2, -52);
	private double EPSILON2 = Math.pow(2, -52);
	private boolean isFilled, pixelFilling, input2d;
	private double minZ, maxZ;
	private double[] xarrayx, yarrayy;
	private double[][] xarrayx2,yarrayy2,zarrayz;
	private double[] contourIntervals;
	private int[] startEdge, startTriangle;
	private String[] contourStyle;
	private List<JDPoint> corners, cntCorner;
	private List<JDEdge> edges, contours;
	private List<JDTriangle> triangles, cntTriangles;
	
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
		cntTriangles = new ArrayList<JDTriangle>();
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		double xs = p[2]/(maxX-minX), ys = p[3]/(maxY-minY);
		
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
		
//		if(ax.getPlot().isDebug()) {
//			JColourtable dct = JColourtable.pctables.get("default");
//			int ti = 0; double tl = triangles.size()-1;
//			for(JDTriangle ttt: triangles) {
//				JDPoint[] abc = {ttt.a, ttt.b, ttt.c};
//				double[][] uvw = new double[3][2];
//				for(int k=0; k<3; k++) {
//					uvw[k][0] = p[0]+xs*(abc[k].x-minX);
//					uvw[k][1] = p[1]+ys*(maxY-abc[k].y);
//				}
//				JPlotShape.fill(dct.getColour(ti/tl));
//				JPlotShape.noStroke();
//				s.addChild(new JTriangleShape((float)uvw[0][0], (float)uvw[0][1], (float)uvw[1][0], (float)uvw[1][1], (float)uvw[2][0], (float)uvw[2][1]));
//				ti++;
//			}
//		}
		
		//step 4: create filling between contours if wished
		if(isFilled) {
			if(pixelFilling) {
				if(ax.getPlot().isDebug())
					System.out.println("[DEBUG] JContourLayer: 4] contour filling pixelwise ...");
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
					//TODO simple pixel colouring!
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
			} else {
				fillContours(ax.getGeoProjection(), ax.getPlot().isDebug());
				JGroupShape trianglesh = new JGroupShape();
				for(int c=0; c<contourIntervals.length; c++) {
					int cs = startTriangle[c];
					int ce = cntTriangles.size();
					if(c+1<contourIntervals.length)
						ce = startTriangle[c+1];
					//TODO edit colour and linestyle info for different contour lines
					double pct = minZ-10d;
					if(c>0 && c<contourIntervals.length)
						pct = 0.5d*(contourIntervals[c-1]+contourIntervals[c]);
					if(c==contourIntervals.length)
						pct = maxZ + 10d;
					int cct = colourtable.getColour(pct, contourIntervals[0], contourIntervals[contourIntervals.length-1]);
					if(ax.getPlot().isDebug())
						System.out.println("[DEBUG] JContourLayer:    color for level "+pct+" (between "+contourIntervals[0]+" and "+
										   contourIntervals[contourIntervals.length-1]+") is "+Integer.toHexString(cct));
					JPlotShape.fill(cct);
					JPlotShape.noStroke();
					if(ax.getPlot().isDebug()) {
						JPlotShape.stroke(0xff999999); JPlotShape.strokeWeight(2f); }
					for(int cl=cs; cl<ce; cl++) {
						JDPoint tv1 = cntTriangles.get(cl).a;
						JDPoint tv2 = cntTriangles.get(cl).b;
						JDPoint tv3 = cntTriangles.get(cl).c;
						double x1 = p[0]+xs*(invertAxisX ? maxX-tv1.x : tv1.x-minX);
						double x2 = p[0]+xs*(invertAxisX ? maxX-tv2.x : tv2.x-minX);
						double x3 = p[0]+xs*(invertAxisX ? maxX-tv3.x : tv3.x-minX);
						double y1 = p[1]+ys*(invertAxisY ? tv1.y-minY : maxY-tv1.y);
						double y2 = p[1]+ys*(invertAxisY ? tv2.y-minY : maxY-tv2.y);
						double y3 = p[1]+ys*(invertAxisY ? tv3.y-minY : maxY-tv3.y);
						for(JDTriangle ttt: cutoff(x1,y1,x2,y2,x3,y3, p[0],p[1],p[0]+p[2],p[1]+p[3])) {
							trianglesh.addChild(
									new JTriangleShape((float)ttt.a.x, (float)ttt.a.y, (float)ttt.b.x, (float)ttt.b.y, (float)ttt.c.x, (float)ttt.c.y)
							);
						}
					}
				}
				s.addChild(trianglesh);
			}
		} else {
            if(ax.getPlot().isDebug())
                System.out.println("[DEBUG] JContourLayer: 4] no filling ...");
        }
		
		//step 5: add contours to plot
		if(drawLines) {
            if(ax.getPlot().isDebug())
                System.out.println("[DEBUG] JContourLayer: 5] register "+(contours.size()-startEdge[0])+" contour line segments with "+lw+"px line width...");
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
		EPSILON2 = xyScale*xyScale * 1.e-18d;
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
        	JDPoint[][] points = new JDPoint[yarrayy.length][xarrayx.length];
        	edges.clear();
        	triangles.clear();
        	for(int j=0; j<yarrayy.length; j++)
        		for(int i=0; i<xarrayx.length; i++)
        			points[j][i] = new JDPoint(xarrayx[i], yarrayy[j], zarrayz[j][i]);
        	for(int j=1; j<yarrayy.length; j++)
        		for(int i=1; i<xarrayx.length; i++)
        			if((2*i>=xarrayx.length?1:0)==(2*j>=yarrayy.length?1:0)) {
        				if(i==1) edges.add(new JDEdge(points[j-1][0], points[j][0]));
        				if(j==1) edges.add(new JDEdge(points[0][i-1], points[0][i]));
        				edges.add(new JDEdge(points[j-1][i], points[j][i]));
        				edges.add(new JDEdge(points[j][i-1], points[j][i]));
        				edges.add(new JDEdge(points[j][i-1], points[j-1][i]));
        				triangles.add(new JDTriangle(points[j-1][i-1], points[j-1][i], points[j][i-1]));
        				triangles.add(new JDTriangle(points[j][i-1], points[j-1][i], points[j][i]));
        			} else {
        				if(i==1) edges.add(new JDEdge(points[j-1][0], points[j][0]));
        				if(j==1) edges.add(new JDEdge(points[0][i-1], points[0][i]));
        				edges.add(new JDEdge(points[j-1][i], points[j][i]));
        				edges.add(new JDEdge(points[j][i-1], points[j][i]));
        				edges.add(new JDEdge(points[j-1][i-1], points[j][i]));
        				triangles.add(new JDTriangle(points[j-1][i-1], points[j-1][i], points[j][i]));
        				triangles.add(new JDTriangle(points[j][i-1], points[j-1][i-1], points[j][i]));
        			}
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
	private void fillContours(JProjection outproj, boolean debug) {
        if(debug)
            System.out.println("[DEBUG] JContourLayer: 4] fill contours ...");
        
		//add corner points to fill last bits in normal plot // ussually not necessary for complex geographical projections
		int nx=xarrayx.length, ny=yarrayy.length;
		double[] xy00 = inputProj.fromPROJtoLATLON(xarrayx[0], yarrayy[0], false);
		xy00 = outproj.fromLATLONtoPROJ(xy00[0], xy00[1], false);
		if(Double.isFinite(xy00[0]) && Double.isFinite(xy00[1]))
			addCorner(xy00[0], xy00[1], zarrayz[0][0], cntCorner);
		double[] xy01 = inputProj.fromPROJtoLATLON(xarrayx[nx-1], yarrayy[0], false);
		xy01 = outproj.fromLATLONtoPROJ(xy01[0], xy01[1], false);
		if(Double.isFinite(xy01[0]) && Double.isFinite(xy01[1]))
			addCorner(xy01[0], xy01[1], zarrayz[0][nx-1], cntCorner);
		double[] xy10 = inputProj.fromPROJtoLATLON(xarrayx[0], yarrayy[ny-1], false);
		xy10 = outproj.fromLATLONtoPROJ(xy10[0], xy10[1], false);
		if(Double.isFinite(xy10[0]) && Double.isFinite(xy10[1]))
			addCorner(xy10[0], xy10[1], zarrayz[ny-1][0], cntCorner);
		double[] xy11 = inputProj.fromPROJtoLATLON(xarrayx[nx-1], yarrayy[ny-1], false);
		xy11 = outproj.fromLATLONtoPROJ(xy11[0], xy11[1], false);
		if(Double.isFinite(xy11[0]) && Double.isFinite(xy11[1]))
			addCorner(xy11[0], xy11[1], zarrayz[ny-1][nx-1], cntCorner);
		if(debug)
			System.out.println("[DEBUG] JContourLayer:    added corner points.");
		
		//create delaunay triangulation
		JConstrainedDelaunayTriangulator delaunayFill = new JConstrainedDelaunayTriangulator(contours, cntCorner);
		//constrain edges afterwards to respect contours
		//delaunayFill.constrain(contours);
		if(debug)
			System.out.println("[DEBUG] JContourLayer:    CDT successful.");
		
		//order triangles with respect to contour level (or by x-axis amount when level is equal)
		cntTriangles = delaunayFill.getTriangles();
		List<JDTriangle> equalLevelTriangle = new ArrayList<JDTriangle>();
		double[] ci = new double[contourIntervals.length];
		ci[0] = contourIntervals[0]*(contourIntervals[0]<0d ? 1.00000001d : 0.99999999d);
		for(int c=1; c<ci.length; c++)
				ci[c] = 0.00000001d*contourIntervals[c-1]+0.99999999d*contourIntervals[c];
		for(JDTriangle cntTri: cntTriangles) {
			int[] l =  {getLevel(cntTri.a.value, contourIntervals, -9999),
						getLevel(cntTri.b.value, contourIntervals, -9999),
						getLevel(cntTri.c.value, contourIntervals, -9999)};
			if(l[0]<0 || l[1]<0 || l[2]<0) {
				cntTri.lev = -1;
				equalLevelTriangle.add(cntTri);
			} else
			if(l[0]==l[1] && l[1]==l[2]) {
				cntTri.lev = l[0];
				equalLevelTriangle.add(cntTri);
			} else {
				cntTri.lev = (l[0]+l[1]+l[2])/3;
			}
		}
		if(debug)
			System.out.println("[DEBUG] JContourLayer:    Simple colouring done.");
		int iter=0;
		while(equalLevelTriangle.size()>0 && iter<1000) {
			iter++;
			for(int t=equalLevelTriangle.size()-1; t>=0; t--) {
				JDTriangle cntTri = equalLevelTriangle.get(t);
				for(JDEdge e: new JDEdge[] {cntTri.ab, cntTri.bc, cntTri.ca}) {
					if(e.has2triangles()) {
						JDTriangle[] w = e.getWing();
						JDTriangle next = w[0].equals(cntTri) ? w[1] : w[0];
						if(next==null)
							continue;
						if(equalLevelTriangle.contains(next))
							continue;
						if(contours.contains(e)) {
							double ev = 0.5d*(e.a.value+e.b.value);
							if(Double.isNaN(e.a.value)) ev = e.b.value;
							if(Double.isNaN(e.b.value)) ev = e.a.value;
							int elev = getLevel(ev, contourIntervals, -1);
							if(elev<0)
								continue;
							if(next.lev>=0) {
								cntTri.lev = 2*elev - 1 - next.lev;
							} else {
								cntTri.lev = -1;
							}
						} else {
							cntTri.lev = next.lev;
						}
						equalLevelTriangle.remove(t);
						break;
					}
				}
			}
		}
		if(debug)
			System.out.println("[DEBUG] JContourLayer:    colors to triangle assigned.");
		cntTriangles.sort(new Comparator<JDTriangle>() {
			@Override
			public int compare(JDTriangle o1, JDTriangle o2) {
				if(o1==null)
					return -1;
				if(o2==null)
					return  1;
				if(o1.lev!=o2.lev)
					return o1.lev - o2.lev;
				return Double.compare(o1.a.x+o1.b.x+o1.c.x, o2.a.x+o2.b.x+o2.c.x);
			}
		});
		startTriangle = new int[contourIntervals.length+1];
		int lastLev = -1;
		for(int c=1; c<cntTriangles.size(); c++) {
			if(cntTriangles.get(c-1).lev<cntTriangles.get(c).lev) {
				lastLev = Math.max(lastLev, cntTriangles.get(c).lev);
				startTriangle[cntTriangles.get(c).lev] = c;
			}
		}
		for(int c=1; c<startTriangle.length; c++)
			if(startTriangle[c]<startTriangle[c-1]) {
				int bigger = c+1;
				for(bigger=c+1; bigger<startTriangle.length; bigger++)
					if(startTriangle[bigger]>startTriangle[c-1])
						break;
				int bc = cntTriangles.size();
				if(bigger<startTriangle.length) {
					bc = startTriangle[bigger];
				} else {
					bigger = startTriangle.length-1;
				}
				for(bigger=bigger+0; bigger>=c; bigger--)
					startTriangle[bigger] = bc;
			}
//		for(int c=lastLev+1; c<startTriangle.length; c++)
//			startTriangle[c] = cntTriangles.size();
		if(debug) {
			System.out.println("[DEBUG] JContourLayer:    estimated "+cntTriangles.size()+" triangles between contourlines");
			String sts = "";
			for(int i=0; i<startTriangle.length; i++)
				sts += ", "+startTriangle[i];
			System.out.println("[DEBUG] JContourLayer:    start indices are: "+sts.substring(2));
		}
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
	private JDTriangle[] cutoff(double x1, double y1, double x2, double y2, double x3, double y3, double le,double to, double ri,double bt) {
		double xi = Math.min(x1,Math.min(x2,x3));
		double xa = Math.max(x1,Math.max(x2,x3));
		if(xa<=le || xi>=ri)
			return new JDTriangle[0];
		double yi = Math.min(y1,Math.min(y2,y3));
		double ya = Math.max(y1,Math.max(y2,y3));
		if(ya<=to || yi>=bt)
			return new JDTriangle[0];
		
		//recalc triangle, if it intersects with vertical borders
		double x4 = Double.NaN, y4 = Double.NaN;
		double x5 = Double.NaN, y5 = Double.NaN;
		int count = 3;
		if(xi<le) {
			int xcode = (x1<le?1:0) | (x2<le?2:0) | (x3<le?4:0);
			switch(xcode) {
				case 1: count=4; double xf12=JPlotMath.dlerp(le,x1,x2,0d,1d), xf13=JPlotMath.dlerp(le,x1,x3,0d,1d);
					x4=x1+xf12*(x2-x1); y4=y1+xf12*(y2-y1); x1=x1+xf13*(x3-x1); y1=y1+xf13*(y3-y1); break;
				case 2: count=4; double xf21=JPlotMath.dlerp(le,x2,x1,0d,1d), xf23=JPlotMath.dlerp(le,x2,x3,0d,1d);
					x4=x1; y4=y1; x1=x2+xf21*(x1-x2); y1=y2+xf21*(y1-y2); x1=x2+xf23*(x3-x2); y2=y1+xf23*(y3-y2); break;
				case 4: count=4; double xf31=JPlotMath.dlerp(le,x3,x1,0d,1d), xf32=JPlotMath.dlerp(le,x3,x2,0d,1d);
					x4=x3+xf31*(x1-x3); y4=y3+xf31*(y1-y3); x3=x3+xf32*(x2-x3); y3=y3+xf32*(y2-y3); break;
				case 3: count=3; double fx32=JPlotMath.dlerp(le,x3,x2,0d,1d), fx31=JPlotMath.dlerp(le,x3,x1,0d,1d);
					x1=x3+fx31*(x1-x3); y1=y3+fx31*(y1-y3); x2=x3+fx32*(x2-x3); y2=y3+fx31*(y2-y3); break;
				case 5: count=3; double fx21=JPlotMath.dlerp(le,x2,x1,0d,1d), fx23=JPlotMath.dlerp(le,x2,x3,0d,1d);
					x1=x2+fx21*(x1-x2); y1=y2+fx21*(y1-y2); x3=x2+fx23*(x3-x2); y3=y2+fx23*(y3-y2); break;
				case 6: count=3; double fx12=JPlotMath.dlerp(le,x1,x2,0d,1d), fx13=JPlotMath.dlerp(le,x1,x3,0d,1d);
					x2=x1+fx12*(x2-x1); y2=y1+fx12*(y2-y1); x3=x1+fx13*(x3-x1); y3=y1+fx13*(y3-y1); break;
				default: count = 3; break;
			}
		}
		if(xa>ri) {
			int xcode = (x1>ri?1:0) | (x2>ri?2:0) | (x3>ri?4:0) | (x4>ri?8:0);
			if(count==3) switch(xcode) {
				case 1: count=4; double xf12=JPlotMath.dlerp(ri,x1,x2,0d,1d), xf13=JPlotMath.dlerp(ri,x1,x3,0d,1d);
					x4=x1+xf12*(x2-x1); y4=y1+xf12*(y2-y1); x1=x1+xf13*(x3-x1); y1=y1+xf13*(y3-y1); break;
				case 2: count=4; double xf21=JPlotMath.dlerp(ri,x2,x1,0d,1d), xf23=JPlotMath.dlerp(ri,x2,x3,0d,1d);
					x4=x1; y4=y1; x1=x2+xf21*(x1-x2); y1=y2+xf21*(y1-y2); x1=x2+xf23*(x3-x2); y2=y1+xf23*(y3-y2); break;
				case 4: count=4; double xf31=JPlotMath.dlerp(ri,x3,x1,0d,1d), xf32=JPlotMath.dlerp(ri,x3,x2,0d,1d);
					x4=x3+xf31*(x1-x3); y4=y3+xf31*(y1-y3); x3=x3+xf32*(x2-x3); y3=y3+xf32*(y2-y3); break;
				case 3: count=3; double fx32=JPlotMath.dlerp(ri,x3,x2,0d,1d), fx31=JPlotMath.dlerp(ri,x3,x1,0d,1d);
					x1=x3+fx31*(x1-x3); y1=y3+fx31*(y1-y3); x2=x3+fx32*(x2-x3); y2=y3+fx31*(y2-y3); break;
				case 5: count=3; double fx21=JPlotMath.dlerp(ri,x2,x1,0d,1d), fx23=JPlotMath.dlerp(ri,x2,x3,0d,1d);
					x1=x2+fx21*(x1-x2); y1=y2+fx21*(y1-y2); x3=x2+fx23*(x3-x2); y3=y2+fx23*(y3-y2); break;
				case 6: count=3; double fx12=JPlotMath.dlerp(ri,x1,x2,0d,1d), fx13=JPlotMath.dlerp(ri,x1,x3,0d,1d);
					x2=x1+fx12*(x2-x1); y2=y1+fx12*(y2-y1); x3=x1+fx13*(x3-x1); y3=y1+fx13*(y3-y1); break;
				default: count=3; break;
			}
			if(count==4) switch(xcode) {
				case 1: count=5; double xf12=JPlotMath.dlerp(ri,x1,x2,0d,1d), xf14=JPlotMath.dlerp(ri,x1,x4,0d,1d);
					x5=x1+xf14*(x4-x1); y5=y1+xf14*(y4-y1); x1=x1+xf12*(x2-x1); y1=y1+xf12*(y2-y1); break;
				case 2: count=5; double xf21=JPlotMath.dlerp(ri,x2,x1,0d,1d), xf23=JPlotMath.dlerp(ri,x2,x3,0d,1d);
					x5=x1; y5=y1; x1=x2+xf21*(x1-x2); y1=y2+xf21*(y1-y2); x2=x2+xf23*(x3-x2); y2=y2+xf23*(y3-y2); break;
				case 4: count=5; double xf32=JPlotMath.dlerp(ri,x3,x2,0d,1d), xf34=JPlotMath.dlerp(ri,x3,x4,0d,1d);
					x5=x4; y5=y4; x4=x3+xf34*(x4-x3); y4=y3+xf34*(y4-y3); x3=x3+xf32*(x2-x3); y3=y3+xf32*(y2-y3); break;
				case 8: count=5; double xf41=JPlotMath.dlerp(ri,x4,x1,0d,1d), xf43=JPlotMath.dlerp(ri,x4,x3,0d,1d);
					x5=x4+xf41*(x1-x4); y5=y4+xf41*(y1-y4); x4=x4+xf43*(x3-x4); y4=y4+xf43*(y3-y4); break;
				case 3: count=4; double fx14=JPlotMath.dlerp(ri,x1,x4,0d,1d), fx23=JPlotMath.dlerp(ri,x2,x3,0d,1d);
					x1=x1+fx14*(x4-x1); y1=y1+fx14*(y4-y1); x2=x2+fx23*(x3-x2); y2=y2+fx23*(y3-y2); break;
				case 6: count=4; double fx21=JPlotMath.dlerp(ri,x2,x1,0d,1d), fx34=JPlotMath.dlerp(ri,x3,x4,0d,1d);
					x2=x2+fx21*(x1-x2); y2=y2+fx21*(y1-y2); x3=x3+fx34*(x4-x3); y3=y3+fx34*(y4-y3); break;
				case 9: count=4; double fx12=JPlotMath.dlerp(ri,x1,x2,0d,1d), fx43=JPlotMath.dlerp(ri,x4,x3,0d,1d);
					x1=x1+fx12*(x2-x1); y1=y1+fx12*(y2-y1); x4=x4+fx43*(x3-x4); y4=y4+fx43*(y3-y4); break;
				case 12: count=4; double fx41=JPlotMath.dlerp(ri,x4,x1,0d,1d), fx32=JPlotMath.dlerp(ri,x3,x2,0d,1d);
					x3=x3+fx32*(x2-x3); y3=y3+fx32*(y2-y3); x4=x4+fx41*(x1-x4); y4=y4+fx41*(y1-y4); break;
				default: count=4; break;
			}
		}
		yi = Math.min(y1,Math.min(y2,y3)); if(count>3 && yi>y4) yi=y4; if(count>4 && yi>y5) yi=y5;
		ya = Math.max(y1,Math.max(y2,y3)); if(count>3 && ya<y4) ya=y4; if(count>4 && ya<y5) ya=y5;
		if(ya<=to || yi>=bt)
			return new JDTriangle[0];

		//recalc triangle, if it intersects with horizontal borders
		double x6 = Double.NaN, y6 = Double.NaN;
		double x7 = Double.NaN, y7 = Double.NaN;
		if(yi<to) {
			int ycode = (y1<to?1:0) | (y2<to?2:0) | (y3<to?4:0) | (y4<to?8:0) | (y5<to?16:0);
			if(count==3) switch(ycode) {
				case 1: count=4; double yy12=JPlotMath.dlerp(to,y1,y2,0d,1d), yy13=JPlotMath.dlerp(to,y1,y3,0d,1d);
					x4=x1+yy12*(x2-x1); y4=y1+yy12*(y2-y1); x1=x1+yy13*(x3-x1); y1=y1+yy13*(y3-y1); break;
				case 2: count=4; double yy21=JPlotMath.dlerp(to,y2,y1,0d,1d), yy23=JPlotMath.dlerp(to,y2,y3,0d,1d);
					x4=x1; y4=y1; x1=x2+yy21*(x1-x2); y1=y2+yy21*(y1-y2); x1=x2+yy23*(x3-x2); y2=y1+yy23*(y3-y2); break;
				case 4: count=4; double yy31=JPlotMath.dlerp(to,y3,y1,0d,1d), yy32=JPlotMath.dlerp(to,y3,y2,0d,1d);
					x4=x3+yy31*(x1-x3); y4=y3+yy31*(y1-y3); x3=x3+yy32*(x2-x3); y3=y3+yy32*(y2-y3); break;
				case 3: count=3; double ff32=JPlotMath.dlerp(to,y3,y2,0d,1d), ff31=JPlotMath.dlerp(to,y3,y1,0d,1d);
					x1=x3+ff31*(x1-x3); y1=y3+ff31*(y1-y3); x2=x3+ff32*(x2-x3); y2=y3+ff31*(y2-y3); break;
				case 5: count=3; double ff21=JPlotMath.dlerp(to,y2,y1,0d,1d), ff23=JPlotMath.dlerp(to,y2,y3,0d,1d);
					x1=x2+ff21*(x1-x2); y1=y2+ff21*(y1-y2); x3=x2+ff23*(x3-x2); y3=y2+ff23*(y3-y2); break;
				case 6: count=3; double ff12=JPlotMath.dlerp(to,y1,y2,0d,1d), ff13=JPlotMath.dlerp(to,y1,y3,0d,1d);
					x2=x1+ff12*(x2-x1); y2=y1+ff12*(y2-y1); x3=x1+ff13*(x3-x1); y3=y1+ff13*(y3-y1); break;
				default: count=3; break;
			}
			if(count==4) switch(ycode) {
				case 1: count=5; double yf12=JPlotMath.dlerp(to,y1,y2,0d,1d), yf14=JPlotMath.dlerp(to,y1,y4,0d,1d);
					x5=x1+yf14*(x4-x1); y5=y1+yf14*(y4-y1); x1=x1+yf12*(x2-x1); y1=y1+yf12*(y2-y1); break;
				case 2: count=5; double yf21=JPlotMath.dlerp(to,y2,y1,0d,1d), yf23=JPlotMath.dlerp(to,y2,y3,0d,1d);
					x5=x1; y5=y1; x1=x2+yf21*(x1-x2); y1=y2+yf21*(y1-y2); x2=x2+yf23*(x3-x2); y2=y2+yf23*(y3-y2); break;
				case 4: count=5; double yf32=JPlotMath.dlerp(to,y3,y2,0d,1d), yf34=JPlotMath.dlerp(to,y3,y4,0d,1d);
					x5=x4; y5=y4; x4=x3+yf34*(x4-x3); y4=y3+yf34*(y4-y3); x3=x3+yf32*(x2-x3); y3=y3+yf32*(y2-y3); break;
				case 8: count=5; double yf41=JPlotMath.dlerp(to,y4,y1,0d,1d), yf43=JPlotMath.dlerp(to,y4,y3,0d,1d);
					x5=x4+yf41*(x1-x4); y5=y4+yf41*(y1-y4); x4=x4+yf43*(x3-x4); y4=y4+yf43*(y3-y4); break;
				case 3: count=4; double fy14=JPlotMath.dlerp(to,y1,y4,0d,1d), fy23=JPlotMath.dlerp(to,y2,y3,0d,1d);
					x1=x1+fy14*(x4-x1); y1=y1+fy14*(y4-y1); x2=x2+fy23*(x3-x2); y2=y2+fy23*(y3-y2); break;
				case 6: count=4; double fy21=JPlotMath.dlerp(to,y2,y1,0d,1d), fy34=JPlotMath.dlerp(to,y3,y4,0d,1d);
					x2=x2+fy21*(x1-x2); y2=y2+fy21*(y1-y2); x3=x3+fy34*(x4-x3); y3=y3+fy34*(y4-y3); break;
				case 9: count=4; double fy12=JPlotMath.dlerp(to,y1,y2,0d,1d), fy43=JPlotMath.dlerp(to,y4,y3,0d,1d);
					x1=x1+fy12*(x2-x1); y1=y1+fy12*(y2-y1); x4=x4+fy43*(x3-x4); y4=y4+fy43*(y3-y4); break;
				case 12: count=4; double fy41=JPlotMath.dlerp(to,y4,y1,0d,1d), fy32=JPlotMath.dlerp(to,y3,y2,0d,1d);
					x3=x3+fy32*(x2-x3); y3=y3+fy32*(y2-y3); x4=x4+fy41*(x1-x4); y4=y4+fy41*(y1-y4); break;
				case 7: count=3; double yy41=JPlotMath.dlerp(to,y4,y1,0d,1d), yy43=JPlotMath.dlerp(to,y4,y3,0d,1d);
					x2=x4+yy41*(x1-x4); y2=y4+yy41*(y1-y4); x3=x4+yy43*(x3-x4); y3=y4+yy43*(y3-y4); x1=x4; y1=y4; x4=Double.NaN; y4=Double.NaN; break;
				case 11: count=3; double yy32=JPlotMath.dlerp(to,y3,y2,0d,1d), yy34=JPlotMath.dlerp(to,y3,y4,0d,1d);
					x1=x3+yy34*(x4-x3); y1=y3+yy34*(y4-y3); x2=x3+yy32*(x2-x3); y2=y3+yy32*(y2-y3); x4=Double.NaN; y4=Double.NaN; break;
				case 13: count=3; double yy21=JPlotMath.dlerp(to,y2,y1,0d,1d), yy23=JPlotMath.dlerp(to,y2,y3,0d,1d);
					x1=x2+yy21*(x1-x2); y1=y2+yy21*(y1-y2); x3=x2+yy23*(x3-x2); y3=y2+yy23*(y3-y2); x4=Double.NaN; y4=Double.NaN; break;
				case 14: count=3; double yy12=JPlotMath.dlerp(to,y1,y2,0d,1d), yy14=JPlotMath.dlerp(to,y1,y4,0d,1d);
					x2=x1+yy12*(x2-x1); y2=y1+yy12*(y2-y1); x3=x1+yy14*(x4-x1); y3=y1+yy14*(y4-y1); x4=Double.NaN; y4=Double.NaN; break;
				default: count=4; break;
			}
			if(count==5) switch(ycode) {
				case 1: count=6; double yf12=JPlotMath.dlerp(to,y1,y2,0d,1d), yf15=JPlotMath.dlerp(to,y1,y5,0d,1d);
					x6=x1+yf15*(x5-x1); y6=y1+yf15*(y5-y1); x1=x1+yf12*(x2-x1); y1=y1+yf12*(y2-y1); break;
				case 2: count=6; double yf21=JPlotMath.dlerp(to,y2,y1,0d,1d), yf23=JPlotMath.dlerp(to,y2,y3,0d,1d);
					x6=x1; y6=y1; x1=x2+yf21*(x1-x2); y1=y2+yf21*(y1-y2); x2=x2+yf23*(x3-x2); y2=y2+yf23*(y3-y2); break;
				case 4: count=6; double yf32=JPlotMath.dlerp(to,y3,y2,0d,1d), yf34=JPlotMath.dlerp(to,y3,y4,0d,1d);
					x6=x1; y6=y1; x1=x2; y1=y2; x2=x3+yf32*(x2-x3); y2=y3+yf32*(y2-y3); x3=x3+yf34*(x4-x3); y3=y3+yf34*(y4-y3); break;
				case 8: count=6; double yf43=JPlotMath.dlerp(to,y4,y3,0d,1d), yf45=JPlotMath.dlerp(to,y4,y5,0d,1d);
					x6=x5; y6=y5; x5=x4+yf45*(x5-x4); y5=y4+yf45*(y5-y4); x4=x4+yf43*(x3-x4); y4=y4+yf43*(y3-y4); break;
				case 16: count=6; double yf51=JPlotMath.dlerp(to,y5,y1,0d,1d), yf54=JPlotMath.dlerp(to,y5,y4,0d,1d);
					x6=x5+yf51*(x1-x5); y6=y5+yf51*(y1-y5); x5=x5+yf54*(x4-x5); y5=y5+yf54*(y4-y5); break;
				case 3: count=5; double fy15=JPlotMath.dlerp(to,y1,y5,0d,1d), fy23=JPlotMath.dlerp(to,y2,y3,0d,1d);
					x1=x1+fy15*(x5-x1); y1=y1+fy15*(y5-y1); x2=x2+fy23*(x3-x2); y2=y2+fy23*(y3-y2); break;
				case 6: count=5; double fy21=JPlotMath.dlerp(to,y2,y1,0d,1d), fy34=JPlotMath.dlerp(to,y3,y4,0d,1d);
					x2=x2+fy21*(x1-x2); y2=y2+fy21*(y1-y2); x3=x3+fy34*(x4-x3); y3=y3+fy34*(y4-y3); break;
				case 12: count=5; double fy32=JPlotMath.dlerp(to,y3,y2,0d,1d), fy45=JPlotMath.dlerp(to,y4,y5,0d,1d);
					x3=x3+fy32*(x2-x3); y3=y3+fy32*(y2-y3); x4=x4+fy45*(x5-x4); y4=y4+fy45*(y5-y4); break;
				case 17: count=5; double fy12=JPlotMath.dlerp(to,y1,y2,0d,1d), fy54=JPlotMath.dlerp(to,y5,y4,0d,1d);
					x1=x1+fy12*(x2-x1); y1=y1+fy12*(y2-y1); x5=x5+fy54*(x4-x5); y5=y5+fy54*(y4-y5); break;
				case 24: count=5; double fy43=JPlotMath.dlerp(to,y4,y3,0d,1d), fy51=JPlotMath.dlerp(to,y5,y1,0d,1d);
					x4=x4+fy43*(x3-x4); y4=y4+fy43*(y3-y4); x5=x5+fy51*(x1-x5); y5=y5+fy51*(y1-y5); break;
				case 7: count=4; double ff15=JPlotMath.dlerp(to,y1,y5,0d,1d), ff34=JPlotMath.dlerp(to,y3,y4,0d,1d);
					x1=x1+ff15*(x5-x1); y1=y1+ff15*(y5-y1); x2=x3+ff34*(x4-x3); y2=y3+ff34*(y4-y3); x3=x4; y3=y4; x4=x5; y4=y5;
					x5=Double.NaN; y5=Double.NaN; break;
				case 14: count=4; double ff21=JPlotMath.dlerp(to,y2,y1,0d,1d), ff45=JPlotMath.dlerp(to,y4,y5,0d,1d);
					x2=x2+ff21*(x1-x2); y2=y2+ff21*(y1-y2); x3=x4+ff45*(x5-x4); y3=y4+ff45*(y5-y4); x4=x5; y4=y5;
					x5=Double.NaN; y5=Double.NaN; break;
				case 19: count=4; double ff23=JPlotMath.dlerp(to,y2,y3,0d,1d), ff54=JPlotMath.dlerp(to,y5,y4,0d,1d);
					x1=x5+ff54*(x4-x5); y1=y5+ff54*(y4-y5); x2=x2+ff23*(x3-x2); y2=y2+ff23*(y3-y2); x5=Double.NaN; y5=Double.NaN; break;
				case 25: count=4; double ff12=JPlotMath.dlerp(to,y1,y2,0d,1d), ff43=JPlotMath.dlerp(to,y4,y3,0d,1d);
					x1=x1+ff12*(x2-x1); y1=y1+ff12*(y2-y1); x4=x4+ff43*(x3-x4); y4=y4+ff43*(y3-y4); x5=Double.NaN; y5=Double.NaN; break;
				case 28: count=4; double ff32=JPlotMath.dlerp(to,y3,y2,0d,1d), ff51=JPlotMath.dlerp(to,y5,y1,0d,1d);
					x3=x3+ff32*(x2-x3); y3=y3+ff32*(y2-y3); x4=x5+ff51*(x1-x5); y4=y5+ff51*(y1-y5); x5=Double.NaN; y5=Double.NaN; break;
				case 15: count=3; double yy15=JPlotMath.dlerp(to,y1,y5,0d,1d), yy45=JPlotMath.dlerp(to,y4,y5,0d,1d);
					x1=x1+yy15*(x5-x1); y1=y1+yy15*(y5-y1); x2=x4+yy45*(x5-x4); y2=y4+yy45*(y5-y4); x3=x5; y3=y5;
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; break;
				case 23: count=3; double yy34=JPlotMath.dlerp(to,y3,y4,0d,1d), yy54=JPlotMath.dlerp(to,y5,y4,0d,1d);
					x1=x4; y1=y4; x2=x5+yy54*(x4-x5); y2=y5+yy54*(y4-y5); x3=x3+yy34*(x4-x3); y3=y3+yy34*(y4-y3);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; break;
				case 27: count=3; double yy23=JPlotMath.dlerp(to,y2,y3,0d,1d), yy43=JPlotMath.dlerp(to,y4,y3,0d,1d);
					x1=x4+yy43*(x3-x4); y1=y4+yy43*(y3-y4); x2=x2+yy23*(x3-x2); y2=y2+yy23*(y3-y2);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; break;
				case 29: count=3; double yy12=JPlotMath.dlerp(to,y1,y2,0d,1d), yy32=JPlotMath.dlerp(to,y3,y2,0d,1d);
					x1=x1+yy12*(x2-x1); y1=y1+yy12*(y2-y1); x3=x3+yy32*(x2-x3); y3=y3+yy32*(y2-y3);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; break;
				case 30: count=3; double yy21=JPlotMath.dlerp(to,y2,y1,0d,1d), yy51=JPlotMath.dlerp(to,y5,y1,0d,1d);
					x2=x2+yy21*(x1-x2); y2=y2+yy21*(y1-y2); x3=x5+yy51*(x1-x5); y3=y5+yy51*(y1-y5);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; break;
				default: count=5; break;
			}
		}
		if(ya>bt) {
			int ycode = (y1>bt?1:0) | (y2>bt?2:0) | (y3>bt?4:0) | (y4>bt?8:0) | (y5>bt?16:0) | (y6>bt?32:0);
			if(count==3) switch(ycode) {
				case 1: count=4; double yy12=JPlotMath.dlerp(bt,y1,y2,0d,1d), yy13=JPlotMath.dlerp(bt,y1,y3,0d,1d);
					x4=x1+yy12*(x2-x1); y4=y1+yy12*(y2-y1); x1=x1+yy13*(x3-x1); y1=y1+yy13*(y3-y1); break;
				case 2: count=4; double yy21=JPlotMath.dlerp(bt,y2,y1,0d,1d), yy23=JPlotMath.dlerp(bt,y2,y3,0d,1d);
					x4=x1; y4=y1; x1=x2+yy21*(x1-x2); y1=y2+yy21*(y1-y2); x1=x2+yy23*(x3-x2); y2=y1+yy23*(y3-y2); break;
				case 4: count=4; double yy31=JPlotMath.dlerp(bt,y3,y1,0d,1d), yy32=JPlotMath.dlerp(bt,y3,y2,0d,1d);
					x4=x3+yy31*(x1-x3); y4=y3+yy31*(y1-y3); x3=x3+yy32*(x2-x3); y3=y3+yy32*(y2-y3); break;
				case 3: count=3; double ff32=JPlotMath.dlerp(bt,y3,y2,0d,1d), ff31=JPlotMath.dlerp(bt,y3,y1,0d,1d);
					x1=x3+ff31*(x1-x3); y1=y3+ff31*(y1-y3); x2=x3+ff32*(x2-x3); y2=y3+ff31*(y2-y3); break;
				case 5: count=3; double ff21=JPlotMath.dlerp(bt,y2,y1,0d,1d), ff23=JPlotMath.dlerp(bt,y2,y3,0d,1d);
					x1=x2+ff21*(x1-x2); y1=y2+ff21*(y1-y2); x3=x2+ff23*(x3-x2); y3=y2+ff23*(y3-y2); break;
				case 6: count=3; double ff12=JPlotMath.dlerp(bt,y1,y2,0d,1d), ff13=JPlotMath.dlerp(bt,y1,y3,0d,1d);
					x2=x1+ff12*(x2-x1); y2=y1+ff12*(y2-y1); x3=x1+ff13*(x3-x1); y3=y1+ff13*(y3-y1); break;
				default: count=3; break;
			}
			if(count==4) switch(ycode) {
				case 1: count=5; double yf12=JPlotMath.dlerp(bt,y1,y2,0d,1d), yf14=JPlotMath.dlerp(bt,y1,y4,0d,1d);
					x5=x1+yf14*(x4-x1); y5=y1+yf14*(y4-y1); x1=x1+yf12*(x2-x1); y1=y1+yf12*(y2-y1); break;
				case 2: count=5; double yf21=JPlotMath.dlerp(bt,y2,y1,0d,1d), yf23=JPlotMath.dlerp(bt,y2,y3,0d,1d);
					x5=x1; y5=y1; x1=x2+yf21*(x1-x2); y1=y2+yf21*(y1-y2); x2=x2+yf23*(x3-x2); y2=y2+yf23*(y3-y2); break;
				case 4: count=5; double yf32=JPlotMath.dlerp(bt,y3,y2,0d,1d), yf34=JPlotMath.dlerp(bt,y3,y4,0d,1d);
					x5=x4; y5=y4; x4=x3+yf34*(x4-x3); y4=y3+yf34*(y4-y3); x3=x3+yf32*(x2-x3); y3=y3+yf32*(y2-y3); break;
				case 8: count=5; double yf41=JPlotMath.dlerp(bt,y4,y1,0d,1d), yf43=JPlotMath.dlerp(bt,y4,y3,0d,1d);
					x5=x4+yf41*(x1-x4); y5=y4+yf41*(y1-y4); x4=x4+yf43*(x3-x4); y4=y4+yf43*(y3-y4); break;
				case 3: count=4; double fy14=JPlotMath.dlerp(bt,y1,y4,0d,1d), fy23=JPlotMath.dlerp(bt,y2,y3,0d,1d);
					x1=x1+fy14*(x4-x1); y1=y1+fy14*(y4-y1); x2=x2+fy23*(x3-x2); y2=y2+fy23*(y3-y2); break;
				case 6: count=4; double fy21=JPlotMath.dlerp(bt,y2,y1,0d,1d), fy34=JPlotMath.dlerp(bt,y3,y4,0d,1d);
					x2=x2+fy21*(x1-x2); y2=y2+fy21*(y1-y2); x3=x3+fy34*(x4-x3); y3=y3+fy34*(y4-y3); break;
				case 9: count=4; double fy12=JPlotMath.dlerp(bt,y1,y2,0d,1d), fy43=JPlotMath.dlerp(bt,y4,y3,0d,1d);
					x1=x1+fy12*(x2-x1); y1=y1+fy12*(y2-y1); x4=x4+fy43*(x3-x4); y4=y4+fy43*(y3-y4); break;
				case 12: count=4; double fy41=JPlotMath.dlerp(bt,y4,y1,0d,1d), fy32=JPlotMath.dlerp(bt,y3,y2,0d,1d);
					x3=x3+fy32*(x2-x3); y3=y3+fy32*(y2-y3); x4=x4+fy41*(x1-x4); y4=y4+fy41*(y1-y4); break;
				case 7: count=3; double yy41=JPlotMath.dlerp(bt,y4,y1,0d,1d), yy43=JPlotMath.dlerp(bt,y4,y3,0d,1d);
					x2=x4+yy41*(x1-x4); y2=y4+yy41*(y1-y4); x3=x4+yy43*(x3-x4); y3=y4+yy43*(y3-y4); x1=x4; y1=y4; x4=Double.NaN; y4=Double.NaN; break;
				case 11: count=3; double yy32=JPlotMath.dlerp(bt,y3,y2,0d,1d), yy34=JPlotMath.dlerp(bt,y3,y4,0d,1d);
					x1=x3+yy34*(x4-x3); y1=y3+yy34*(y4-y3); x2=x3+yy32*(x2-x3); y2=y3+yy32*(y2-y3); x4=Double.NaN; y4=Double.NaN; break;
				case 13: count=3; double yy21=JPlotMath.dlerp(bt,y2,y1,0d,1d), yy23=JPlotMath.dlerp(bt,y2,y3,0d,1d);
					x1=x2+yy21*(x1-x2); y1=y2+yy21*(y1-y2); x3=x2+yy23*(x3-x2); y3=y2+yy23*(y3-y2); x4=Double.NaN; y4=Double.NaN; break;
				case 14: count=3; double yy12=JPlotMath.dlerp(bt,y1,y2,0d,1d), yy14=JPlotMath.dlerp(bt,y1,y4,0d,1d);
					x2=x1+yy12*(x2-x1); y2=y1+yy12*(y2-y1); x3=x1+yy14*(x4-x1); y3=y1+yy14*(y4-y1); x4=Double.NaN; y4=Double.NaN; break;
				default: count=4; break;
			}
			if(count==5) switch(ycode) {
				case 1: count=6; double yf12=JPlotMath.dlerp(bt,y1,y2,0d,1d), yf15=JPlotMath.dlerp(bt,y1,y5,0d,1d);
					x6=x1+yf15*(x5-x1); y6=y1+yf15*(y5-y1); x1=x1+yf12*(x2-x1); y1=y1+yf12*(y2-y1); break;
				case 2: count=6; double yf21=JPlotMath.dlerp(bt,y2,y1,0d,1d), yf23=JPlotMath.dlerp(bt,y2,y3,0d,1d);
					x6=x1; y6=y1; x1=x2+yf21*(x1-x2); y1=y2+yf21*(y1-y2); x2=x2+yf23*(x3-x2); y2=y2+yf23*(y3-y2); break;
				case 4: count=6; double yf32=JPlotMath.dlerp(bt,y3,y2,0d,1d), yf34=JPlotMath.dlerp(bt,y3,y4,0d,1d);
					x6=x1; y6=y1; x1=x2; y1=y2; x2=x3+yf32*(x2-x3); y2=y3+yf32*(y2-y3); x3=x3+yf34*(x4-x3); y3=y3+yf34*(y4-y3); break;
				case 8: count=6; double yf43=JPlotMath.dlerp(bt,y4,y3,0d,1d), yf45=JPlotMath.dlerp(bt,y4,y5,0d,1d);
					x6=x5; y6=y5; x5=x4+yf45*(x5-x4); y5=y4+yf45*(y5-y4); x4=x4+yf43*(x3-x4); y4=y4+yf43*(y3-y4); break;
				case 16: count=6; double yf51=JPlotMath.dlerp(bt,y5,y1,0d,1d), yf54=JPlotMath.dlerp(bt,y5,y4,0d,1d);
					x6=x5+yf51*(x1-x5); y6=y5+yf51*(y1-y5); x5=x5+yf54*(x4-x5); y5=y5+yf54*(y4-y5); break;
				case 3: count=5; double fy15=JPlotMath.dlerp(bt,y1,y5,0d,1d), fy23=JPlotMath.dlerp(bt,y2,y3,0d,1d);
					x1=x1+fy15*(x5-x1); y1=y1+fy15*(y5-y1); x2=x2+fy23*(x3-x2); y2=y2+fy23*(y3-y2); break;
				case 6: count=5; double fy21=JPlotMath.dlerp(bt,y2,y1,0d,1d), fy34=JPlotMath.dlerp(bt,y3,y4,0d,1d);
					x2=x2+fy21*(x1-x2); y2=y2+fy21*(y1-y2); x3=x3+fy34*(x4-x3); y3=y3+fy34*(y4-y3); break;
				case 12: count=5; double fy32=JPlotMath.dlerp(bt,y3,y2,0d,1d), fy45=JPlotMath.dlerp(bt,y4,y5,0d,1d);
					x3=x3+fy32*(x2-x3); y3=y3+fy32*(y2-y3); x4=x4+fy45*(x5-x4); y4=y4+fy45*(y5-y4); break;
				case 17: count=5; double fy12=JPlotMath.dlerp(bt,y1,y2,0d,1d), fy54=JPlotMath.dlerp(bt,y5,y4,0d,1d);
					x1=x1+fy12*(x2-x1); y1=y1+fy12*(y2-y1); x5=x5+fy54*(x4-x5); y5=y5+fy54*(y4-y5); break;
				case 24: count=5; double fy43=JPlotMath.dlerp(bt,y4,y3,0d,1d), fy51=JPlotMath.dlerp(bt,y5,y1,0d,1d);
					x4=x4+fy43*(x3-x4); y4=y4+fy43*(y3-y4); x5=x5+fy51*(x1-x5); y5=y5+fy51*(y1-y5); break;
				case 7: count=4; double ff15=JPlotMath.dlerp(bt,y1,y5,0d,1d), ff34=JPlotMath.dlerp(bt,y3,y4,0d,1d);
					x1=x1+ff15*(x5-x1); y1=y1+ff15*(y5-y1); x2=x3+ff34*(x4-x3); y2=y3+ff34*(y4-y3); x3=x4; y3=y4; x4=x5; y4=y5;
					x5=Double.NaN; y5=Double.NaN; break;
				case 14: count=4; double ff21=JPlotMath.dlerp(bt,y2,y1,0d,1d), ff45=JPlotMath.dlerp(bt,y4,y5,0d,1d);
					x2=x2+ff21*(x1-x2); y2=y2+ff21*(y1-y2); x3=x4+ff45*(x5-x4); y3=y4+ff45*(y5-y4); x4=x5; y4=y5;
					x5=Double.NaN; y5=Double.NaN; break;
				case 19: count=4; double ff23=JPlotMath.dlerp(bt,y2,y3,0d,1d), ff54=JPlotMath.dlerp(bt,y5,y4,0d,1d);
					x1=x5+ff54*(x4-x5); y1=y5+ff54*(y4-y5); x2=x2+ff23*(x3-x2); y2=y2+ff23*(y3-y2); x5=Double.NaN; y5=Double.NaN; break;
				case 25: count=4; double ff12=JPlotMath.dlerp(bt,y1,y2,0d,1d), ff43=JPlotMath.dlerp(bt,y4,y3,0d,1d);
					x1=x1+ff12*(x2-x1); y1=y1+ff12*(y2-y1); x4=x4+ff43*(x3-x4); y4=y4+ff43*(y3-y4); x5=Double.NaN; y5=Double.NaN; break;
				case 28: count=4; double ff32=JPlotMath.dlerp(bt,y3,y2,0d,1d), ff51=JPlotMath.dlerp(bt,y5,y1,0d,1d);
					x3=x3+ff32*(x2-x3); y3=y3+ff32*(y2-y3); x4=x5+ff51*(x1-x5); y4=y5+ff51*(y1-y5); x5=Double.NaN; y5=Double.NaN; break;
				case 15: count=3; double yy15=JPlotMath.dlerp(bt,y1,y5,0d,1d), yy45=JPlotMath.dlerp(bt,y4,y5,0d,1d);
					x1=x1+yy15*(x5-x1); y1=y1+yy15*(y5-y1); x2=x4+yy45*(x5-x4); y2=y4+yy45*(y5-y4); x3=x5; y3=y5;
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; break;
				case 23: count=3; double yy34=JPlotMath.dlerp(bt,y3,y4,0d,1d), yy54=JPlotMath.dlerp(bt,y5,y4,0d,1d);
					x1=x4; y1=y4; x2=x5+yy54*(x4-x5); y2=y5+yy54*(y4-y5); x3=x3+yy34*(x4-x3); y3=y3+yy34*(y4-y3);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; break;
				case 27: count=3; double yy23=JPlotMath.dlerp(bt,y2,y3,0d,1d), yy43=JPlotMath.dlerp(bt,y4,y3,0d,1d);
					x1=x4+yy43*(x3-x4); y1=y4+yy43*(y3-y4); x2=x2+yy23*(x3-x2); y2=y2+yy23*(y3-y2);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; break;
				case 29: count=3; double yy12=JPlotMath.dlerp(bt,y1,y2,0d,1d), yy32=JPlotMath.dlerp(bt,y3,y2,0d,1d);
					x1=x1+yy12*(x2-x1); y1=y1+yy12*(y2-y1); x3=x3+yy32*(x2-x3); y3=y3+yy32*(y2-y3);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; break;
				case 30: count=3; double yy21=JPlotMath.dlerp(bt,y2,y1,0d,1d), yy51=JPlotMath.dlerp(bt,y5,y1,0d,1d);
					x2=x2+yy21*(x1-x2); y2=y2+yy21*(y1-y2); x3=x5+yy51*(x1-x5); y3=y5+yy51*(y1-y5);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; break;
				default: count=5; break;
			}
			if(count==6) switch(ycode) {
				case 1: count=7; double yf12=JPlotMath.dlerp(bt,y1,y2,0d,1d), yf16=JPlotMath.dlerp(bt,y1,y6,0d,1d);
					x7=x1+yf16*(x6-x1); y7=y1+yf16*(y6-y1); x1=x1+yf12*(x2-x1); y1=y1+yf12*(y2-y1); break;
				case 2: count=7; double yf21=JPlotMath.dlerp(bt,y2,y1,0d,1d), yf23=JPlotMath.dlerp(bt,y2,y3,0d,1d);
					x7=x1; y7=y1; x1=x2+yf21*(x1-x2); y1=y2+yf21*(y1-y2); x2=x2+yf23*(x3-x2); y2=y2+yf23*(y3-y2); break;
				case 4: count=7; double yf32=JPlotMath.dlerp(bt,y3,y2,0d,1d), yf34=JPlotMath.dlerp(bt,y3,y4,0d,1d);
					x7=x1; y7=y1; x1=x2; y1=y2; x2=x3+yf32*(x2-x3); y2=y3+yf32*(y2-y3); x3=x3+yf34*(x4-x3); y3=y3+yf34*(y4-y3); break;
				case 8: count=7; double yf43=JPlotMath.dlerp(bt,y4,y3,0d,1d), yf45=JPlotMath.dlerp(bt,y4,y5,0d,1d);
					x7=x6; y7=y6; x6=x5; y6=y5; x5=x4+yf45*(x5-x4); y5=y4+yf45*(y5-y4); x4=x4+yf43*(x3-x4); y4=y4+yf43*(y3-y4); break;
				case 16: count=7; double yf54=JPlotMath.dlerp(bt,y5,y4,0d,1d), yf56=JPlotMath.dlerp(bt,y5,y6,0d,1d);
					x7=x6; y7=y6; x6=x5+yf56*(x6-x5); y6=y5+yf56*(y6-y5); x5=x5+yf54*(x4-x5); y5=y5+yf54*(y4-y5); break;
				case 32: count=7; double yf61=JPlotMath.dlerp(bt,y6,y1,0d,1d), yf65=JPlotMath.dlerp(bt,y6,y5,0d,1d);
					x7=x6+yf61*(x1-x6); y7=y6+yf61*(y1-y6); x6=x6+yf65*(x5-x6); y6=y6+yf65*(y5-y6); break;
				case 3: count=6; double fy16=JPlotMath.dlerp(bt,y1,y6,0d,1d), fy23=JPlotMath.dlerp(bt,y2,y3,0d,1d);
					x1=x1+fy16*(x6-x1); y1=y1+fy16*(y6-y1); x2=x2+fy23*(x3-x2); y2=y2+fy23*(y3-y2); break;
				case 6: count=6; double fy21=JPlotMath.dlerp(bt,y2,y1,0d,1d), fy34=JPlotMath.dlerp(bt,y3,y4,0d,1d);
					x2=x2+fy21*(x1-x2); y2=y2+fy21*(y1-y2); x3=x3+fy34*(x4-x3); y3=y3+fy34*(y4-y3); break;
				case 12: count=6; double fy32=JPlotMath.dlerp(bt,y3,y2,0d,1d), fy45=JPlotMath.dlerp(bt,y4,y5,0d,1d);
					x3=x3+fy32*(x2-x3); y3=y3+fy32*(y2-y3); x4=x4+fy45*(x5-x4); y4=y4+fy45*(y5-y4); break;
				case 24: count=6; double fy43=JPlotMath.dlerp(bt,y4,y3,0d,1d), fy56=JPlotMath.dlerp(bt,y5,y6,0d,1d);
					x4=x4+fy43*(x3-x4); y4=y4+fy43*(y3-y4); x5=x5+fy56*(x6-x5); y5=y5+fy56*(y6-y5); break;
				case 48: count=6; double fy54=JPlotMath.dlerp(bt,y5,y4,0d,1d), fy61=JPlotMath.dlerp(bt,y6,y1,0d,1d);
					x5=x5+fy54*(x4-x5); y5=y5+fy54*(y4-y5); x6=x6+fy61*(x1-x6); y6=y6+fy61*(y1-y6); break;
				case 33: count=6; double fy12=JPlotMath.dlerp(bt,y1,y2,0d,1d), fy65=JPlotMath.dlerp(bt,y6,y5,0d,1d);
					x1=x1+fy12*(x2-x1); y1=y1+fy12*(y2-y1); x6=x6+fy65*(x5-x6); y6=y6+fy65*(y5-y6); break;
				case 7: count=5; double yy16=JPlotMath.dlerp(bt,y1,y6,0d,1d), yy34=JPlotMath.dlerp(bt,y3,y4,0d,1d);
					x2=x1+yy16*(x6-x1); y2=y1+yy16*(y6-y1); x3=x3+yy34*(x4-x3); y3=y3+yy34*(y4-y3); x1=x6; y1=y6;
					x6=Double.NaN; y6=Double.NaN; break;
				case 14: count=5; double yy21=JPlotMath.dlerp(bt,y2,y1,0d,1d), yy45=JPlotMath.dlerp(bt,y4,y5,0d,1d);
					x2=x2+yy21*(x1-x2); y2=y2+yy21*(y1-y2); x3=x4+yy45*(x5-x4); y3=y4+yy45*(y5-y4); x4=x5; y4=y5; x5=x6; y5=y6;
					x6=Double.NaN; y6=Double.NaN; break;
				case 28: count=5; double yy32=JPlotMath.dlerp(bt,y3,y2,0d,1d), yy56=JPlotMath.dlerp(bt,y5,y6,0d,1d);
					x3=x3+yy32*(x2-x3); y3=y3+yy32*(y2-y3); x4=x5+yy56*(x6-x5); y4=y5+yy56*(y6-y5); x5=x6; y5=y6;
					x6=Double.NaN; y6=Double.NaN; break;
				case 36: count=5; double yy23=JPlotMath.dlerp(bt,y2,y3,0d,1d), yy65=JPlotMath.dlerp(bt,y6,y5,0d,1d);
					x2=x2+yy23*(x3-x2); y2=y2+yy23*(y3-y2); x1=x6+yy65*(x5-x6); y1=y6+yy65*(y5-y6); x6=Double.NaN; y6=Double.NaN; break;
				case 49: count=5; double yy12=JPlotMath.dlerp(bt,y1,y2,0d,1d), yy54=JPlotMath.dlerp(bt,y5,y4,0d,1d);
					x1=x1+yy12*(x2-x1); y1=y1+yy12*(y2-y1); x5=x5+yy54*(x4-x5); y5=y5+yy54*(y4-y5); x6=Double.NaN; y6=Double.NaN; break;
				case 56: count=5; double yy43=JPlotMath.dlerp(bt,y4,y3,0d,1d), yy61=JPlotMath.dlerp(bt,y6,y1,0d,1d);
					x4=x4+yy43*(x3-x4); y4=y4+yy43*(y3-y4); x5=x6+yy61*(x1-x6); y5=y6+yy61*(y1-y6); x6=Double.NaN; y6=Double.NaN; break;
				case 30: count=4; double ff21=JPlotMath.dlerp(bt,y2,y1,0d,1d), ff56=JPlotMath.dlerp(bt,y5,y6,0d,1d);
					x2=x2+ff21*(x1-x2); y2=y2+ff21*(y2-y1); x3=x5+ff56*(x6-x5); y3=y5+ff56*(y6-y5); x4=x6; y4=y6;
					x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 15: count=4; double ff16=JPlotMath.dlerp(bt,y1,y6,0d,1d), ff45=JPlotMath.dlerp(bt,y4,y5,0d,1d);
					x1=x1+ff16*(x6-x1); y1=y1+ff16*(y6-y1); x2=x4+ff45*(x5-x4); y2=y4+ff45*(y5-y4); x3=x5; y3=y5; x4=x6; y4=y6;
					x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 39: count=4; double ff34=JPlotMath.dlerp(bt,y3,y4,0d,1d), ff65=JPlotMath.dlerp(bt,y6,y5,0d,1d);
					x1=x5; y1=y5; x2=x6+ff65*(x5-x6); y2=y6+ff65*(y5-y6); x3=x3+ff34*(x4-x3); y3=y3+ff34*(y4-y3);
					x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 51: count=4; double ff23=JPlotMath.dlerp(bt,y2,y3,0d,1d), ff54=JPlotMath.dlerp(bt,y5,y4,0d,1d);
					x1=x5+ff54*(x4-x5); y1=y5+ff54*(y4-y5); x2=x2+ff23*(x3-x2); y2=y2+ff23*(y3-y2);
					x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 57: count=4; double ff12=JPlotMath.dlerp(bt,y1,y2,0d,1d), ff43=JPlotMath.dlerp(bt,y4,y3,0d,1d);
					x1=x1+ff12*(x2-x1); y1=y1+ff12*(y2-y1); x4=x4+ff43*(x3-x4); y4=y4+ff43*(y3-y4);
					x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 60: count=4; double ff61=JPlotMath.dlerp(bt,y6,y1,0d,1d), ff32=JPlotMath.dlerp(bt,y3,y2,0d,1d);
					x3=x3+ff32*(x2-x3); y3=y3+ff32*(y2-y3); x4=x6+ff61*(x1-x6); y4=y6+ff61*(y1-y6);
					x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 31: count=3; double jj16=JPlotMath.dlerp(bt,y1,y6,0d,1d), jj56=JPlotMath.dlerp(bt,y5,y6,0d,1d);
					x1=x1+jj16*(x6-x1); y1=y1+jj16*(y6-y1); x2=x5+jj56*(x6-x5); y2=y5+jj56*(y6-y5); x3=x6; y3=y6;
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 47: count=3; double jj45=JPlotMath.dlerp(bt,y4,y5,0d,1d), jj65=JPlotMath.dlerp(bt,y6,y5,0d,1d);
					x1=x4+jj45*(x5-x4); y1=y4+jj45*(y5-y4); x2=x5; y2=y5; x3=x6+jj65*(x5-x6); y3=y6+jj65*(y5-y6);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 55: count=3; double jj34=JPlotMath.dlerp(bt,y3,y4,0d,1d), jj54=JPlotMath.dlerp(bt,y5,y4,0d,1d);
					x1=x4; y1=y4; x2=x5+jj54*(x4-x5); y2=y5+jj54*(y4-y5); x3=x3+jj34*(x4-x3); y3=y3+jj34*(y4-y3);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 59: count=3; double jj23=JPlotMath.dlerp(bt,y2,y3,0d,1d), jj43=JPlotMath.dlerp(bt,y4,y3,0d,1d);
					x1=x4+jj43*(x3-x4); y1=y4+jj43*(y3-y4); x2=x2+jj23*(x3-x2); y2=y2+jj23*(y3-y2);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 61: count=3; double jj12=JPlotMath.dlerp(bt,y1,y2,0d,1d), jj32=JPlotMath.dlerp(bt,y3,y2,0d,1d);
					x1=x1+jj12*(x2-x1); y1=y1+jj12*(y2-y1); x3=x3+jj32*(x2-x3); y3=y3+jj32*(y2-y3);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				case 62: count=3; double jj21=JPlotMath.dlerp(bt,y2,y1,0d,1d), jj61=JPlotMath.dlerp(bt,y6,y1,0d,1d);
					x2=x2+jj21*(x1-x2); y2=y2+jj21*(y1-y2); x3=x6+jj61*(x1-x6); y3=y6+jj61*(y1-y6);
					x4=Double.NaN; y4=Double.NaN; x5=Double.NaN; y5=Double.NaN; x6=Double.NaN; y6=Double.NaN; break;
				default: count=6; break;
			}
		}

		JDPoint p1 = new JDPoint(x1, y1, 0d);
		JDPoint p2 = new JDPoint(x2, y2, 0d);
		JDPoint p3 = new JDPoint(x3, y3, 0d);
		JDPoint p4 = new JDPoint(x4, y4, 0d);
		JDPoint p5 = new JDPoint(x5, y5, 0d);
		JDPoint p6 = new JDPoint(x6, y6, 0d);
		JDPoint p7 = new JDPoint(x7, y7, 0d);
		switch(count) {
			default:
			case 3: return new JDTriangle[] {
					new JDTriangle(p1, p2, p3)};
			case 4: return new JDTriangle[] {
					new JDTriangle(p1, p2, p3),
					new JDTriangle(p1, p2, p4)};
			case 5: return new JDTriangle[] {
					new JDTriangle(p1, p2, p3),
					new JDTriangle(p1, p2, p4),
					new JDTriangle(p1, p2, p5)};
			case 6: return new JDTriangle[] {
					new JDTriangle(p1, p2, p3),
					new JDTriangle(p1, p2, p4),
					new JDTriangle(p1, p2, p5),
					new JDTriangle(p1, p2, p6)};
			case 7: return new JDTriangle[] {
					new JDTriangle(p1, p2, p3),
					new JDTriangle(p1, p2, p4),
					new JDTriangle(p1, p2, p5),
					new JDTriangle(p1, p2, p6),
					new JDTriangle(p1, p2, p7)};
		}
	}
}
