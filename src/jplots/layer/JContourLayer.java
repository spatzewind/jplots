package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import jplots.JAxis;
import jplots.JPlot;
import jplots.color.JColourtable;
import jplots.maths.JDEdge;
import jplots.maths.JDPoint;
import jplots.maths.JDTriangle;
import jplots.maths.JDelaunayTriangulator;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JTriangleShape;
import jplots.transform.JProjection;
import processing.core.PApplet;
import processing.core.PGraphics;

public class JContourLayer extends JPlotsLayer {

	private boolean hasContours, isFilled;
	private double minZ, maxZ;
	private double[] xarrayx, yarrayy;
	private double[][] zarrayz;
	private double[] contourIntervals;
	private int[] startEdge, startTriangle;
	private String[] contourStyle;
	private List<JDPoint> corners, cntCorner;
	private List<JDEdge> edges, contours;
	private List<JDTriangle> triangles, cntTriangles;
	
	public JContourLayer(float[] x, float[] y, float[][] z, float zmin, float zmax, int nintervals, JColourtable ct, float stroke_weight, boolean drawContours, boolean filled) {
		float zin = Float.isNaN(zmin) ? JPlotMath.fmin(z) : zmin;
		float zax = Float.isNaN(zmax) ? JPlotMath.fmax(z) : zmax;
		float[] cntIntervals = new float[nintervals+1];
		cntIntervals[nintervals] = zax;
		for(int k=0; k<nintervals; k++)
			cntIntervals[k] = zin + k*(zax-zin)/nintervals;
		JContourLayerFloat(x, y, z, cntIntervals, ct, stroke_weight, drawContours, filled);
	}
	public JContourLayer(float[] x, float[] y, float[][] z, float[] intervals, JColourtable ct, float stroke_weight, boolean drawContours, boolean filled) {
		JContourLayerFloat(x, y, z, intervals, ct, stroke_weight, drawContours, filled);
	}
	public JContourLayer(double[] x, double[] y, double[][] z, double zmin, double zmax, int nintervals, JColourtable ct, double stroke_weight, boolean drawContours, boolean filled) {
		double zin = Double.isNaN(zmin) ? JPlotMath.dmin(z) : zmin;
		double zax = Double.isNaN(zmax) ? JPlotMath.dmax(z) : zmax;
		double[] cntIntervals = new double[nintervals+1];
		cntIntervals[nintervals] = zax;
		for(int k=0; k<nintervals; k++)
			cntIntervals[k] = zin + k*(zax-zin)/nintervals;
		JContourLayerDouble(x, y, z, cntIntervals, ct, stroke_weight, drawContours, filled);
	}
	public JContourLayer(double[] x, double[] y, double[][] z, double[] intervals, JColourtable ct, double stroke_weight, boolean drawContours, boolean filled) {
		JContourLayerDouble(x, y, z, intervals, ct, stroke_weight, drawContours, filled);
	}

	private void JContourLayerFloat(float[] x, float[] y, float[][] z, float[] intervals, JColourtable ct, float stroke_weight, boolean drawContours, boolean filled) {
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
		contourIntervals = new double[intervals.length];
		for(int i=0; i<intervals.length; i++)
			contourIntervals[i] = intervals[i];
		minZ = contourIntervals[0];
		maxZ = contourIntervals[intervals.length-1];
		colourtable = ct;
		lw = stroke_weight;
		hasContours = drawContours;
		isFilled = filled;
		init();
	}
	private void JContourLayerDouble(double[] x, double[] y, double[][] z, double[] intervals, JColourtable ct, double stroke_weight, boolean drawContours, boolean filled) {
		xarrayx = x;
		yarrayy = y;
		zarrayz = z;
		minX = JPlotMath.dmin(xarrayx);
		maxX = JPlotMath.dmax(xarrayx);
		minY = JPlotMath.dmin(yarrayy);
		maxY = JPlotMath.dmax(yarrayy);
		contourIntervals = new double[intervals.length];
		for(int i=0; i<intervals.length; i++)
			contourIntervals[i] = intervals[i];
		minZ = contourIntervals[0];
		maxZ = contourIntervals[intervals.length-1];
		colourtable = ct;
		lw = (float) stroke_weight;
		hasContours = drawContours;
		isFilled = filled;
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
		triangulate(ax.getPlot().isDebug());
		
		//step 3: create contours
        if(ax.getPlot().isDebug())
            System.out.println("[DEBUG] JContourLayer: 3] create contours from triangle mesh ...");
		createContours(ax.getPlot().isDebug());
		if(ax.getPlot().isDebug()) {
			System.out.println("[DEBUG] JContourLayer:    all corners");
			for(int c=0; c<cntCorner.size(); c++)
				System.out.println("    "+PApplet.nf((float)cntCorner.get(c).x,0,4)+
						"    "+PApplet.nf((float)cntCorner.get(c).y,0,4)+
						"    "+PApplet.nf((float)cntCorner.get(c).lev,0,4));
			System.out.println("[DEBUG] JContourLayer:    all edges / line segments");
			for(int c=0; c<contours.size(); c++)
				System.out.println("    "+contours.get(c).a+
						"    "+contours.get(c).b);
		}
		
		if(ax.getPlot().isDebug()) {
			JColourtable dct = JColourtable.pctables.get("default");
			int ti = 0; double tl = triangles.size()-1;
			for(JDTriangle ttt: triangles) {
				JDPoint[] abc = {ttt.a, ttt.b, ttt.c};
				double[][] uvw = new double[3][2];
				for(int k=0; k<3; k++) {
					double[] xy = inputProj.fromPROJtoLATLON(abc[k].x, abc[k].y, false);
					uvw[k] = ax.getGeoProjection().fromLATLONtoPROJ(xy[0], xy[1], false);
					uvw[k][0] = p[0]+xs*(uvw[k][0]-minX);
					uvw[k][1] = p[1]+ys*(maxY-uvw[k][1]);
				}
				JPlotShape.fill(dct.getColour(ti/tl));
				JPlotShape.noStroke();
				s.addChild(new JTriangleShape((float)uvw[0][0], (float)uvw[0][1], (float)uvw[1][0], (float)uvw[1][1], (float)uvw[2][0], (float)uvw[2][1]));
				ti++;
			}
		}
		
		//step 4: create filling between contours if wished
		if(isFilled) {
            if(ax.getPlot().isDebug())
                System.out.println("[DEBUG] JContourLayer: 4] fill contours ...");
			fillContours();
			for(int c=-1; c<contourIntervals.length; c++) {
				int cs = 0;
				if(c>=0)
					cs = startTriangle[c];
				int ce = cntTriangles.size();
				if(c+1<contourIntervals.length)
					ce = startTriangle[c+1];
				//TODO edit colour and linestyle info for different contour lines
				double pct = minZ-10d;
				if(c>=0 && c+1<contourIntervals.length)
					pct = 0.5d*(contourIntervals[c]+contourIntervals[c+1]);
				if(c==contourIntervals.length-1)
					pct = maxZ + 10d;
				JPlotShape.fill(colourtable.getColour(pct, minZ, maxZ));
				JPlotShape.noStroke();
				for(int cl=cs; cl<ce; cl++) {
					JDPoint lc1 = cntTriangles.get(cl).a;
					JDPoint lc2 = cntTriangles.get(cl).b;
					JDPoint lc3 = cntTriangles.get(cl).c;
					s.addChild(new JTriangleShape((float)lc1.x, (float)lc1.y,
							(float)lc2.x, (float)lc2.y, (float)lc3.x, (float)lc3.y));
				}
			}
		} else {
            if(ax.getPlot().isDebug())
                System.out.println("[DEBUG] JContourLayer: 4] no filling ...");
        }
		
		//step 5: add contours to plot
		if(hasContours) {
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
					double[] xy0 = inputProj.fromPROJtoLATLON(lvs.x, lvs.y, false);
					double[] xy1 = inputProj.fromPROJtoLATLON(lve.x, lve.y, false);
					if(ax.isGeoAxis()) {
						xy0 = ax.getGeoProjection().fromLATLONtoPROJ(xy0[0], xy0[1], false);
						xy1 = ax.getGeoProjection().fromLATLONtoPROJ(xy1[0], xy1[1], false);
					}
					double x1 = p[0]+xs*(xy0[0]-minX);
					double x2 = p[0]+xs*(xy1[0]-minX);
					double y1 = p[1]+ys*(maxY-xy0[1]);
					double y2 = p[1]+ys*(maxY-xy1[1]);
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
	
	private void collectValidPoints(JProjection outproj, boolean debug) {
		int nx=xarrayx.length, ny=yarrayy.length;
		corners.clear();
		for(int j=0; j<ny; j++)
			for(int i=0; i<nx; i++) {
				double[] xy = inputProj.fromPROJtoLATLON(xarrayx[i], yarrayy[j], false);
				xy = outproj.fromLATLONtoPROJ(xy[0], xy[1], false);
				if(Double.isFinite(xy[0]) && Double.isFinite(xy[1]))
					corners.add(new JDPoint(xy[0],xy[1],zarrayz[j][i]));
			}
        if(debug)
            System.out.println("[DEBUG] JContourLayer: 1] has "+corners.size()+" valid sourcepoints");
	}
	private void triangulate(boolean debug) {
        if(debug)
            System.out.println("[DEBUG] JContourLayer: 2] triangulate ...");
        JDelaunayTriangulator delaunay = new JDelaunayTriangulator(corners);
        edges = delaunay.getEdges();
        triangles = delaunay.getTriangles();
        if(debug)
            System.out.println("[DEBUG] JContourLayer:   mesh now consists of "+corners.size()+" corners, "+edges.size()+" edges and "+triangles.size()+" triangles");
	}
	private void createContours(boolean debug) {
		startEdge = new int[contourIntervals.length];
		int nanid = -9999;
		contours.clear();
		//first exclude fillvalues
		for(int level=nanid; level<contourIntervals.length; level = Math.max(0, level+1)) {
            if(debug)
                System.out.println("[DEBUG] JContourLayer:   create contours for level "+(level<0 ? Double.NaN : contourIntervals[level]));
			int t_idx = -1;
			for(JDTriangle ct: triangles) {
				t_idx++;
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
							JDPoint i16s = addCorner(x161,y161,Double.NaN,cntCorner);
							double x162 = 0.5d*(va.x+vc.x);
							double y162 = 0.5d*(va.y+vc.y);
							JDPoint i16e = addCorner(x162,y162,Double.NaN,cntCorner);
							if(!i16s.equals(i16e))
								contours.add(new JDEdge(i16s, i16e));
							break;
						case 2:
						case 5:
							double x251 = 0.5d*(vb.x+va.x);
							double y251 = 0.5d*(vb.y+va.y);
							JDPoint i25s = addCorner(x251,y251,Double.NaN,cntCorner);
							double x252 = 0.5d*(vb.x+vc.x);
							double y252 = 0.5d*(vb.y+vc.y);
							JDPoint i25e = addCorner(x252,y252,Double.NaN,cntCorner);
							if(!i25s.equals(i25e))
								contours.add(new JDEdge(i25s, i25e));
							break;
						case 3:
						case 4:
							double x341 = 0.5d*(vc.x+vb.x);
							double y341 = 0.5d*(vc.y+vb.y);
							JDPoint i34s = addCorner(x341,y341,Double.NaN,cntCorner);
							double x342 = 0.5d*(vc.x+va.x);
							double y342 = 0.5d*(vc.y+va.y);
							JDPoint i34e = addCorner(x342,y342,Double.NaN,cntCorner);
							if(!i34s.equals(i34e))
								contours.add(new JDEdge(i34s, i34e));
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
					if(debug)
						System.out.println("[DEBUG] JContourplot:   triangle "+t_idx+" -> levcode "+Integer.toBinaryString(levcode));
					switch(levcode) {
						case  1: //x000001
						case 20: //x010100
							double x161 = va.x+r12*(vb.x-va.x);
							double y161 = va.y+r12*(vb.y-va.y);
							JDPoint i17s = addCorner(x161,y161,lev,cntCorner);
							double x162 = vc.x+r31*(va.x-vc.x);
							double y162 = vc.y+r31*(va.y-vc.y);
							JDPoint i17e = addCorner(x162,y162,lev,cntCorner);
							if(!i17s.equals(i17e))
								contours.add(new JDEdge(i17s, i17e));
							break;
						case  5: //x000101
						case 16: //x010000
							double x341 = vb.x+r23*(vc.x-vb.x);
							double y341 = vb.y+r23*(vc.y-vb.y);
							JDPoint i20s = addCorner(x341,y341,lev,cntCorner);
							double x342 = vc.x+r31*(va.x-vc.x);
							double y342 = vc.y+r31*(va.y-vc.y);
							JDPoint i20e = addCorner(x342,y342,lev,cntCorner);
							if(!i20s.equals(i20e))
								contours.add(new JDEdge(i20s, i20e));
							break;
						case  4: //x000100
						case 17: //x010001
							double x251 = va.x+r12*(vb.x-va.x);
							double y251 = va.y+r12*(vb.y-va.y);
							JDPoint i05s = addCorner(x251,y251,lev,cntCorner);
							double x252 = vb.x+r23*(vc.x-vb.x);
							double y252 = vb.y+r23*(vc.y-vb.y);
							JDPoint i05e = addCorner(x252,y252,lev,cntCorner);
							if(!i05s.equals(i05e))
								contours.add(new JDEdge(i05s, i05e));
							break;
						case 33: //x100001
						case 36: //x100100
							double x12c = va.x+r12*(vb.x-va.x);
							double y12c = va.y+r12*(vb.y-va.y);
							JDPoint i41s = addCorner(x12c,y12c,Double.NaN,cntCorner);
							double x12a = 0.5d*(va.x+vc.x), y12a = 0.5d*(va.y+vc.y);
							double x12b = 0.5d*(vb.x+vc.x), y12b = 0.5d*(vb.y+vc.y);
							JDPoint i41e = addCorner(x12a+r12*(x12b-x12a),y12a+r12*(y12b-y12a),Double.NaN,cntCorner);
							contours.add(new JDEdge(i41s, i41e));
							break;
						case  6: //x000110
						case 18: //x010010
							double x25c = vb.x+r23*(vc.x-vb.x);
							double y25c = vb.y+r23*(vc.y-vb.y);
							JDPoint i22s = addCorner(x25c,y25c,lev,cntCorner);
							double x23b = 0.5d*(vb.x+va.x), y23b = 0.5d*(vb.y+va.y);
							double x23c = 0.5d*(vc.x+va.x), y23c = 0.5d*(vc.y+va.y);
							JDPoint i22e = addCorner(x23b+r23*(x23c-x23b),y23b+r23*(y23c-y23b),Double.NaN,cntCorner);
							contours.add(new JDEdge(i22s, i22e));
							break;
						case  7: //x001001
						case 24: //x011000
							double x34c = vc.x+r31*(va.x-vc.x);
							double y34c = vc.y+r31*(va.y-vc.y);
							JDPoint i26s = addCorner(x34c,y34c,lev,cntCorner);
							double x31c = 0.5d*(vc.x+vb.x), y31c = 0.5d*(vc.y+vb.y);
							double x31a = 0.5d*(va.x+vb.x), y31a = 0.5d*(va.y+vb.y);
							JDPoint i26e = addCorner(x31c+r31*(x31a-x31c),y31c+r31*(y31a-y31c),Double.NaN,cntCorner);
							contours.add(new JDEdge(i26s, i26e));
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
	}
	private JDPoint addCorner(double x, double y, double v, List<JDPoint> points) {
		//double tolerance = 1.0e-10d;
		JDPoint np = new JDPoint(x,y,v);
		for(JDPoint p: points)
			if(p.equals(np))
				return p;
		points.add(np);
		return np;
	}
	private void fillContours() {
		startTriangle = new int[contourIntervals.length];
	}

}
