package jplots.layer;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import jplots.JAxis;
import jplots.JPlot;
import jplots.maths.AffineBuilder;
import jplots.maths.JDEdge;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import jplots.maths.JDTriangle;
import jplots.maths.JDelaunayTriangulator;
import jplots.maths.JGridTriangulator;
import jplots.maths.JPlotMath;
import jplots.shapes.JEllipseShape;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JPolygonShape;
import jplots.transform.JProjection;
import processing.core.PGraphics;

public class JHatchLayer extends JPlotsLayer {

	private boolean input2d;
	private double minZ, maxZ, Xin, Xax, Yin, Yax;
	private double[] xarrayx, yarrayy;
	private double[][] xarrayx2, yarrayy2, zarrayz;
	private String pattern;
	private List<JDPoint> corners;
	private List<JDEdge> edges;
	private List<JDTriangle> triangles;
	private Comparator<JDEdge> comparator;

	public JHatchLayer(float[] x, float[] y, float[][] x2, float[][] y2, float[][] z, float stroke_weight, float lowerBound, float upperBound, String pattern) {
		input2d = true;
		if(x!=null || y!=null) input2d = false;
		else if(x2==null || y2==null)
			throw new IllegalArgumentException("Either 1D or 2D coordinate arrays x,y have to be given, but all are null");
		if (!input2d) {
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
			for(int j=0; j<x2.length; j++)
				for(int i=0; i<x2[j].length; i++)
					xarrayx2[j][i] = x2[j][i];
			yarrayy2 = new double[y2.length][y2[0].length];
			for(int j=0; j<y2.length; j++)
				for(int i=0; i<y2[j].length; i++)
					yarrayy2[j][i] = y2[j][i];
			minX = JPlotMath.dmin(xarrayx2);
			maxX = JPlotMath.dmax(xarrayx2);
			minY = JPlotMath.dmin(yarrayy2);
			maxY = JPlotMath.dmax(yarrayy2);
		}
		zarrayz = new double[z.length][z[0].length];
		for(int j=0; j<z.length; j++)
			for(int i=0; i<z[j].length; i++)
				zarrayz[j][i] = z[j][i];
		minZ = lowerBound;
		maxZ = upperBound;
		if(Double.isNaN(upperBound)) maxZ = 1d + JPlotMath.fmax(z);
		lw = (float) stroke_weight;
		this.pattern = pattern;
		init();
	}

	public JHatchLayer(double[] x, double[] y, double[][] x2, double[][] y2, double[][] z, double stroke_weight, double lowerBound, double upperBound, String pattern) {
		input2d = true;
		if(x!=null || y!=null) input2d = false;
		else if(x2==null || y2==null)
			throw new IllegalArgumentException("Either 1D or 2D coordinate arrays x,y have to be given, but all are null");
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
		minZ = lowerBound;
		maxZ = upperBound;
		if(Double.isNaN(upperBound)) maxZ = 1d + JPlotMath.dmax(z);
		lw = (float) stroke_weight;
		this.pattern = pattern;
		init();
	}
	
	private void init() {
		corners = new ArrayList<>();
		edges = new ArrayList<>();
		triangles = new ArrayList<>();
		
		comparator = new Comparator<JDEdge>() {
			@Override
			public int compare(JDEdge e1, JDEdge e2) {
				return e1.a.compareTo(e2.a);
			}
		};
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
//		System.out.println("Do hatching!");
		int[] p = ax.getSize();
		Xin = ax.isXlogAxis() ? Math.log10(minX) : minX;
		Xax = ax.isXlogAxis() ? Math.log10(maxX) : maxX;
		Yin = ax.isYlogAxis() ? Math.log10(minY) : minY;
		Yax = ax.isYlogAxis() ? Math.log10(maxY) : maxY;
		double xs = p[2] / (Xax - Xin), ys = p[3] / (Yax - Yin);
		// double tol = Math.max(Math.abs(maxX-minX), Math.abs(maxY-minY)) * 1.0e-12d;

		// step 1: collect valid corners of grid and project them
		collectValidPoints(ax.getGeoProjection(), ax.isXlogAxis(), ax.isYlogAxis(), ax.getPlot().isDebug());

		// step 2: do delauney-triangulation
		triangulate(ax);

		// step 4: create filling between contours if wished
		int lineCode = -1, pointCode = -1;
		if(pattern.equals("|"))    lineCode = 0x01;
		if(pattern.equals("||"))   lineCode = 0x11;
		if(pattern.equals("-"))    lineCode = 0x02;
		if(pattern.equals("--"))   lineCode = 0x12;
		if(pattern.equals("+"))    lineCode = 0x03;
		if(pattern.equals("++"))   lineCode = 0x13;
		if(pattern.equals("/"))    lineCode = 0x04;
		if(pattern.equals("//"))   lineCode = 0x14;
		if(pattern.equals("\\"))   lineCode = 0x08;
		if(pattern.equals("\\\\")) lineCode = 0x18;
		if(pattern.equals("x"))    lineCode = 0x0c;
		if(pattern.equals("xx"))   lineCode = 0x1c;
		if(pattern.equals("X"))    lineCode = 0x0c;
		if(pattern.equals("XX"))   lineCode = 0x1c;
		if(pattern.equals("o"))    pointCode = 0x01;
		if(pattern.equals("O"))    pointCode = 0x02;
		if(pattern.equals("Oo"))   pointCode = 0x03;
		if(pattern.equals("oO"))   pointCode = 0x03;
		if(pattern.equals("."))    pointCode = 0x04;
		if(pattern.equals(".."))   pointCode = 0x14;
		if(pattern.equals("*"))    pointCode = 0x05;
		if(pattern.equals("**"))   pointCode = 0x15;
		
		if(lineCode<0 && pointCode<0) {
			System.err.println("Cannot recognize patter '"+pattern+"'");
			return;
		}
		

		AffineBuilder affine = new AffineBuilder().scale(invertAxisX ? -1d : 1d, invertAxisY ? 1d : -1d)
				.translate(invertAxisX ? maxX : -minX, invertAxisY ? -minY : maxY).scale(xs, ys).translate(p[0], p[1]);
		for(JDTriangle t: triangles) t.affine(affine.getMatrix());
		JPlot parent = ax.getPlot();
		int[] psize = parent.getSize();
		double pdist = Math.min((double)p[0]/(double)parent.getNumColumns(), (double)p[1]/(double)parent.getNumRows()) / 10d;
		psize = ax.getSize();
		
		if(pointCode>0) {
			boolean dblDens = (pointCode & 0x10) != 0;
			float pdx = (float) (dblDens ? 0.5d*pdist : pdist);
			float pdy = pdx*(float)Math.sqrt(0.75d);
			JPlotShape.fill(lc); JPlotShape.stroke(lc); JPlotShape.strokeWeight(Math.min(0.15f*pdx, 2f*(float)lw));
//			if(parent.isDebug())
//			System.out.println("Parent size:        "+parent.getSize()[0]+"x"+parent.getSize()[1]);
//			System.out.println("Axis-pos&size:      "+psize[0]+"/"+psize[1]+" ... "+(psize[0]+psize[2])+"/"+(psize[1]+psize[3]));
//			System.out.println("Point-stepdistance: "+pdx+"/"+pdy);
			int si=(int)(psize[0]/pdx), ei=1+(int)((psize[0]+psize[2])/pdx);
			int sj=(int)(psize[1]/pdy), ej=1+(int)((psize[1]+psize[3])/pdy);
			for(int j=sj; j<=ej; j++) {
				float y = j*(float)pdy;
				if(y<psize[1]) continue;
				if(y>psize[1]+psize[3]) continue;
				for(int i=si; i<=ei; i++) {
					float x = (i + (j%2==1?0.5f:0f))*pdx;
					if(x<psize[0]) continue;
					if(x>psize[0]+psize[2]) continue;
					JDPoint pt = new JDPoint(x,y);
					for(JDTriangle t: triangles) {
						if(!t.contains(pt, 0.05d))
							continue;
						double v = t.valueAt(pt);
						if(Double.isNaN(v)) continue;
						if(v<minZ) continue;
						if(v>maxZ) continue;
						switch(pointCode&0x0f) {
							case 1: s.addChild(new JEllipseShape(x, y, 0.2f*pdx, 0.2f*pdx, false)); break;
							case 2: s.addChild(new JEllipseShape(x, y, 0.4f*pdx, 0.4f*pdx, false)); break;
							case 3: s.addChild(new JEllipseShape(x, y, 0.2f*pdx, 0.2f*pdx, false));
									s.addChild(new JEllipseShape(x, y, 0.4f*pdx, 0.4f*pdx, false)); break;
							case 4: s.addChild(new JEllipseShape(x, y, 0.2f*pdx, 0.2f*pdx, true)); break;
							case 5: s.addChild(new JPolygonShape(new float[][] {
								{x+0.000f*pdx,y-0.600f*pdx}, {x-0.134f*pdx,y-0.185f*pdx}, {x-0.571f*pdx,y-0.185f*pdx},
								{x-0.217f*pdx,y+0.071f*pdx}, {x-0.353f*pdx,y+0.485f*pdx}, {x+0.000f*pdx,y+0.228f*pdx},
								{x+0.353f*pdx,y+0.485f*pdx}, {x+0.217f*pdx,y+0.071f*pdx}, {x+0.571f*pdx,y-0.185f*pdx},
								{x+0.134f*pdx,y-0.185f*pdx}}, true, false)); break;
							default: break;
						}
					}
				}
			}
		}
		if(lineCode>0) {
			boolean dblDens = (lineCode & 0x10) != 0;
			float pd = (float) (dblDens ? 0.5d*pdist : pdist);
			JPlotShape.stroke(lc); JPlotShape.strokeWeight(Math.min(0.15f*pd, 2f*(float)lw));
//			if(parent.isDebug())
//			System.out.println("Parent size:        "+parent.getSize()[0]+"x"+parent.getSize()[1]);
//			System.out.println("Axis-pos&size:      "+psize[0]+"/"+psize[1]+" ... "+(psize[0]+psize[2])+"/"+(psize[1]+psize[3]));
//			System.out.println("Point-stepdistance: "+pd);
			List<JDEdge> lines = new ArrayList<>();
			
			JPlotShape.stroke(lc); JPlotShape.strokeWeight(Math.min(0.3f*pd, (float)lw));
			if((lineCode&0x01)==1) {
				lines.clear();
				int st = (int) (p[0]/pd), en = (int) ((p[0]+p[2])/pd+1);
				for(int i=st; i<=en; i++) {
					if(i*pd<p[0]) continue;
					if(i*pd>p[0]+p[2]) continue;
					JDPoint top = new JDPoint(i*pd, p[1]);
					JDPoint bottom = new JDPoint(i*pd, p[1]+p[3]);
					for(JDTriangle t: triangles)
						cutoutLevelrange(t, top, bottom, lines, false);
				}
//				System.out.println("befor: "+lines.size()+"x |");
				lineCombination(lines);
//				System.out.println("after: "+lines.size()+"x |");
				lineCombination(lines);
//				System.out.println("after: "+lines.size()+"x |");
				for(JDEdge l: lines) {
//					JPlotShape.stroke(l.b.y>l.a.y?0xff00ff00:0xffff0000);
					s.addChild(new JLineShape((float)l.a.x, (float)l.a.y, (float)l.b.x, (float)l.b.y));
				}
			}

			if((lineCode&0x02)==2) {
				lines.clear();
				int st = (int) (p[1]/pd), en = (int) ((p[1]+p[3])/pd+1);
				for(int i=st; i<=en; i++) {
					if(i*pd<p[1]) continue;
					if(i*pd>p[1]+p[3]) continue;
					JDPoint left  = new JDPoint(p[0], i*pd);
					JDPoint right = new JDPoint(p[0]+p[2], i*pd);
					for(JDTriangle t: triangles)
						cutoutLevelrange(t, left, right, lines, false);
				}
//				System.out.println("befor: "+lines.size()+"x --");
				lineCombination(lines);
//				System.out.println("after: "+lines.size()+"x --");
				lineCombination(lines);
//				System.out.println("after: "+lines.size()+"x --");
				for(JDEdge l: lines) {
//					JPlotShape.stroke(l.b.x>l.a.x?0xff00ff00:0xffff0000);
					s.addChild(new JLineShape((float)l.a.x, (float)l.a.y, (float)l.b.x, (float)l.b.y));
				}
			}
			
			pd *= (float) Math.sqrt(2d);
			
			if((lineCode&0x04)==4) {
				lines.clear();
				int o = p[0]+p[1];
				int st = (int) (o/pd), en = (int) ((o+p[2]+p[3])/pd+1);
				for(int i=st; i<=en; i++) {
					float m = i*pd - o;
					if(m<0) continue;
					if(m>p[2]+p[3]) continue;
					JDPoint topRight = new JDPoint(p[0]+Math.min(m, p[2]), p[1]+Math.max(m-p[2], 0));
					JDPoint leftBott = new JDPoint(p[0]+Math.max(m-p[3], 0), p[1]+Math.min(m, p[3]));
					for(JDTriangle t: triangles)
						cutoutLevelrange(t, topRight, leftBott, lines, true);
				}
//				System.out.println("befor: "+lines.size()+"x --");
				lineCombination(lines);
//				System.out.println("after: "+lines.size()+"x --");
				lineCombination(lines);
//				System.out.println("after: "+lines.size()+"x --");
				for(JDEdge l: lines) {
//					JPlotShape.stroke(l.b.x>l.a.x?0xff00ff00:0xffff0000);
					s.addChild(new JLineShape((float)l.a.x, -(float)l.a.y, (float)l.b.x, -(float)l.b.y));
				}
			}
			

			if((lineCode&0x08)==8) {
				lines.clear();
				int o = p[0]-p[1]-p[3];
				int st = (int) (o/pd)-(o<0?1:0), en = (int) ((o+p[2]+p[3])/pd)-(o+p[2]+p[3]<0?0:1);
				for(int i=st; i<=en; i++) {
					float fwd = i*pd - o;
					float bwd = p[2]+p[3]-fwd;
					if(fwd<0 || bwd<0) continue;
					JDPoint LeftTop  = new JDPoint(p[0]+Math.max(fwd-p[3], 0), p[1]+Math.max(p[3]-fwd, 0));
					JDPoint botRight = new JDPoint(p[0]+Math.min(fwd, p[2]),   p[1]+Math.min(bwd, p[3]));
					for(JDTriangle t: triangles)
						cutoutLevelrange(t, LeftTop, botRight, lines, false);
				}
//				System.out.println("befor: "+lines.size()+"x --");
				lineCombination(lines);
//				System.out.println("after: "+lines.size()+"x --");
				lineCombination(lines);
//				System.out.println("after: "+lines.size()+"x --");
				for(JDEdge l: lines) {
//					JPlotShape.stroke(l.b.x>l.a.x?0xff00ff00:0xffff0000);
					s.addChild(new JLineShape((float)l.a.x, (float)l.a.y, (float)l.b.x, (float)l.b.y));
				}
			}
		}
	}

	public double[] getZRange() {
		return new double[] { minZ, maxZ };
	}

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
							yLog ? Math.log10(yarrayy2[j][i]) : yarrayy2[j][i], false);
				} else {
					xy = inputProj.fromPROJtoLATLON(xLog ? Math.log10(xarrayx[i]) : xarrayx[i],
							yLog ? Math.log10(yarrayy[j]) : yarrayy[j], false);
				}
				xy = outproj.fromLATLONtoPROJ(xy[0], xy[1], false);
				if (Double.isFinite(xy[0]) && Double.isFinite(xy[1]))
					corners.add(new JDPoint(xy[0], xy[1], zarrayz[j][i]));
			}
		double xyScale = Math.max(Math.max(-Xin, Xax), Math.max(-Yin, Yax));
		xyScale = Math.max(xyScale, Math.max(Xax - Xin, Yax - Yin));
		if (debug)
			System.out.println("[DEBUG] JContourLayer: 1] has " + corners.size() + " valid sourcepoints");
	}

	private void triangulate(JAxis ax) {
		boolean debug = ax.getPlot().isDebug();
		if (debug)
			System.out.println("[DEBUG] JContourLayer: 2] triangulate ...");
		if (input2d || ax.isGeoAxis()) {
			JDelaunayTriangulator delaunay = new JDelaunayTriangulator(corners);
			edges = delaunay.getEdges();
			triangles = delaunay.getTriangles();
		} else {
			JGridTriangulator triangulator = new JGridTriangulator(ax.isXlogAxis() ? JPlotMath.log10(xarrayx) : xarrayx,
					ax.isYlogAxis() ? JPlotMath.log10(yarrayy) : yarrayy, zarrayz);
			edges = triangulator.getEdges();
			triangles = triangulator.getTriangles();
		}
		if (debug)
			System.out.println("[DEBUG] JContourLayer:   mesh now consists of " + corners.size() + " corners, "
					+ edges.size() + " edges and " + triangles.size() + " triangles");
	}
	public List<JDTriangle> getTriangles() {
		return triangles;
	}
	
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

	private void cutoutLevelrange(JDTriangle t, JDPoint p1, JDPoint p2, List<JDEdge> edges, boolean reverseY) {
		// p.x = w1*t.x[0] + w2*t.x[1] + w3*t.x[2];
		// p.y = w1*t.y[0] + w2*t.y[1] + w3*t.y[2];
		// 1.0 = w1        + w2        + w3;
		JDPoint[] line = t.intersectsLine(p1,p2,0.0001d);
		if(line==null) return;
		double ry = reverseY ? -1 : 1;
		
		JDTriangle sub1 = null;
		JDTriangle sub2 = null;
		int nancode = (Double.isNaN(t.value[0])?1:0) | (Double.isNaN(t.value[1])?2:0) | (Double.isNaN(t.value[2])?4:0);
		JDPoint ab = t.getA().fractionTowards(0.5d, t.getB());
		JDPoint bc = t.getB().fractionTowards(0.5d, t.getC());
		JDPoint ca = t.getC().fractionTowards(0.5d, t.getA());
		switch(nancode) {
			case 0: break;
			case 1: ab.value = t.value[1]; ca.value = t.value[2];
					sub1 = new JDTriangle(t.getB(),t.getC(),ca); sub2 = new JDTriangle(t.getB(),ca,ab);
					break;
			case 2: ab.value = t.value[0]; bc.value = t.value[2];
					sub1 = new JDTriangle(t.getC(),t.getA(),ab); sub2 = new JDTriangle(t.getC(),ab,bc);
					break;
			case 3: bc.value = t.value[2]; ca.value = t.value[2];
					sub1 = new JDTriangle(t.getC(),ca,bc);
					break;
			case 4: bc.value = t.value[1]; ca.value = t.value[0];
					sub1 = new JDTriangle(t.getA(),t.getB(),bc); sub2 = new JDTriangle(t.getA(),bc,ca);
					break;
			case 5: ab.value = t.value[1]; bc.value = t.value[1];
					sub1 = new JDTriangle(t.getB(),bc,ab);
					break;
			case 6: ab.value = t.value[0]; ca.value = t.value[0];
					sub1 = new JDTriangle(t.getA(),ab,ca);
					break;
			default:
			case 7: return;
		}
		
		if(sub1==null) {
			double v0 = t.valueAt(line[0]);
			if(Double.isNaN(v0)) return;
			double v1 = t.valueAt(line[1]);
			if(Double.isNaN(v1)) return;
			if(v0<minZ && v1<minZ) return;
			if(v0>maxZ && v1>maxZ) return;
			if(v0<minZ) line[0] = line[0].fractionTowards((minZ-v0)/(v1-v0), line[1]);
			if(v0>maxZ) line[0] = line[0].fractionTowards((maxZ-v0)/(v1-v0), line[1]);
		    v0 = t.valueAt(line[0]);
			if(v1<minZ) line[1] = line[1].fractionTowards((minZ-v1)/(v0-v1), line[0]);
			if(v1>maxZ) line[1] = line[1].fractionTowards((maxZ-v1)/(v0-v1), line[0]);
			edges.add(new JDEdge(new JDPoint(line[0].x, ry*line[0].y), new JDPoint(line[1].x, ry*line[1].y)));
			return;
		}
		
		line = sub1.intersectsLine(p1, p2);
		if(line!=null) {
			double v0 = sub1.valueAt(line[0]);
			if(Double.isNaN(v0)) return;
			double v1 = sub1.valueAt(line[1]);
			if(Double.isNaN(v1)) return;
			if(v0<minZ && v1<minZ) return;
			if(v0>maxZ && v1>maxZ) return;
			if(v0<minZ) line[0] = line[0].fractionTowards((minZ-v0)/(v1-v0), line[1]);
			if(v0>maxZ) line[0] = line[0].fractionTowards((maxZ-v0)/(v1-v0), line[1]);
		    v0 = sub1.valueAt(line[0]);
			if(v1<minZ) line[1] = line[1].fractionTowards((minZ-v1)/(v0-v1), line[0]);
			if(v1>maxZ) line[1] = line[1].fractionTowards((maxZ-v1)/(v0-v1), line[0]);
			edges.add(new JDEdge(new JDPoint(line[0].x, ry*line[0].y), new JDPoint(line[1].x, ry*line[1].y)));
		}
		
		if(sub2==null) return;
		
		line = sub2.intersectsLine(p1, p2);
		if(line!=null) {
			double v0 = sub2.valueAt(line[0]);
			if(Double.isNaN(v0)) return;
			double v1 = sub2.valueAt(line[1]);
			if(Double.isNaN(v1)) return;
			if(v0<minZ && v1<minZ) return;
			if(v0>maxZ && v1>maxZ) return;
			if(v0<minZ) line[0] = line[0].fractionTowards((minZ-v0)/(v1-v0), line[1]);
			if(v0>maxZ) line[0] = line[0].fractionTowards((maxZ-v0)/(v1-v0), line[1]);
		    v0 = sub2.valueAt(line[0]);
			if(v1<minZ) line[1] = line[1].fractionTowards((minZ-v1)/(v0-v1), line[0]);
			if(v1>maxZ) line[1] = line[1].fractionTowards((maxZ-v1)/(v0-v1), line[0]);
			edges.add(new JDEdge(new JDPoint(line[0].x, ry*line[0].y), new JDPoint(line[1].x, ry*line[1].y)));
		}
	}
	private void lineCombination(List<JDEdge> lines) {
		lines.sort(comparator);
		for(int j=lines.size()-1; j>0; j--) {
			JDEdge ej = lines.get(j);
			double tj = Math.abs(ej.a.x-ej.b.x) + Math.abs(ej.a.y-ej.b.y);
			boolean[] rem = new boolean[j];
			for(int i=0; i<j; i++) rem[i] = false;
			for(int i=j-1; i>=0; i--) {
				JDEdge ei = lines.get(i);
				double ti = Math.abs(ei.a.x-ei.b.x) + Math.abs(ei.a.y-ei.b.y); 
				double tol = 0.000001d*Math.min(tj, ti);
//				if(ej.a.equals(ei.a, tol)) {
//					lines.set(j,new JDEdge(ei.b, ej.b));
//					ej = lines.get(j); rem[i] = true; continue;
//				}
				if(ej.a.equals(ei.b, tol)) {
					lines.set(j,new JDEdge(ei.a, ej.b));
					ej = lines.get(j); rem[i] = true; continue;
				}
				if(ej.b.equals(ei.a, tol)) {
					lines.add(new JDEdge(ei.b, ej.a));
					ej = lines.get(j); rem[i] = true; continue;
				}
//				if(ej.b.equals(ei.b, tol)) {
//					lines.add(new JDEdge(ei.a, ej.a));
//					ej = lines.get(j); rem[i] = true; continue;
//				}
			}
			for(int i=j-1; i>=0; i--) if(rem[i]) {
				lines.remove(i);
				j--;
			}
		}
	}
}
