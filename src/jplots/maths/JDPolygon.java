package jplots.maths;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;

public class JDPolygon {

	public int idx;
	public int lev;
	public JDPoint[] c;
	public JDEdge[] e;

	private Integer hash = null;

	public JDPolygon(JDPoint... abc) {
		this(abc, 0);
	}
	public JDPolygon(JDPoint[] abc, int index) {
		if(abc==null || abc.length<3) {
			System.err.println("Cannot create polygon of less than 3 points!");
			c = null;
			e = null;
			idx = 0;
			return;
		}
		boolean clockwise = (area(abc)<0d);
		c = new JDPoint[abc.length];
		for(int i=0; i<abc.length; i++)
			c[clockwise?abc.length-1-i:i] = new JDPoint(abc[i].x, abc[i].y, abc[i].value);
		e = new JDEdge[abc.length];
		for(int i=0,j=e.length-1; i<e.length; j=i++)
			e[i] = new JDEdge(c[j],c[i]);
		idx = index;
	}
	
	public JDPolygon(Coordinate[] coords) {
		this(coords, 0);
	}
	public JDPolygon(Coordinate[] coords, int index) {
		if(coords==null || coords.length<4) {
			System.err.println("Cannot create polygon of less than 3 points!");
			c = null;
			e = null;
			idx = 0;
			return;
		}
		c = new JDPoint[coords.length-1];
		for(int i=0; i<c.length; i++)
			c[i] = new JDPoint(coords[i]);
		double a = area(c);
		if(a<0d)
			for(int i=0; i<c.length; i++)
				c[i] = new JDPoint(coords[c.length-1-i]);
		e = new JDEdge[c.length];
		for(int i=0,j=c.length-1; i<c.length; j=i++)
			e[i] = new JDEdge(c[j],c[i]);
		idx = index;
	}

//	public void edges(JDEdge ab, JDEdge bc, JDEdge ca) {
//		this.ab = ab.equals(a, b) ? ab : bc.equals(a, b) ? bc : ca;
//		this.bc = ab.equals(b, c) ? ab : bc.equals(b, c) ? bc : ca;
//		this.ca = ab.equals(c, a) ? ab : bc.equals(c, a) ? bc : ca;
//	}
	
	public JDPolygon copy() {
		return new JDPolygon(c==null?new JDPoint[0]:c);
	}
	public Geometry toGeometry(GeometryFactory factory) {
		if(c==null)
			return null;
		Coordinate[] ccc = new Coordinate[c.length+1];
		for(int i=0; i<c.length; i++)
			ccc[i] = new Coordinate(c[i].x, c[i].y, c[i].value);
		ccc[c.length] = new Coordinate(c[0].x, c[0].y, c[0].value);
		return factory.createPolygon(ccc);
	}
	
	public List<JDTriangle> toTriangles() {
		List<JDTriangle> triangles = new ArrayList<JDTriangle>();
		List<JDPoint> points = new LinkedList<JDPoint>(Arrays.asList(c));
		while(points.size()>3) {
			int smallestIdx = -1;
			double smallestArea = Double.POSITIVE_INFINITY;
			int slen = points.size();
			for(int s0=0; s0<slen-1; s0++) {
				int s1 = (s0+1)%slen;
				int s2 = (s0+2)%slen;
				double a = area(points.get(s0), points.get(s1), points.get(s2));
				if(a<0d) continue;
				if(a<smallestArea) {
					smallestIdx = s0;
					smallestArea = a;
				}
			}
			if(smallestIdx==-1)
				throw new RuntimeException("Triangulation of polygon failed!");
			int next0 = (smallestIdx+1)%slen;
			int next1 = (smallestIdx+2)%slen;
			triangles.add(new JDTriangle(points.get(smallestIdx), points.get(next0), points.get(next1)));
			points.remove(next0);
		}
		triangles.add(new JDTriangle(points.get(0), points.get(1), points.get(2)));
		return triangles;
	}
	
	public boolean union(JDTriangle tri, double tolerance) {
		double ar = area(tri.a,tri.b,tri.c);
		boolean triClockwise = (ar<0d);
		//System.out.println("area of triangle to add is "+ar);
		int ai = containsCorner(tri.a, tolerance);
		int bi = containsCorner(triClockwise?tri.c:tri.b, tolerance);
		int ci = containsCorner(triClockwise?tri.b:tri.c, tolerance);
		//System.out.println("Corner ids: a="+ai+" b="+bi+" c="+ci);
		if(ai==-1) {
			if(bi==-1 || ci==-1) return false;
			return union(bi,ci, tri.a);
		}
		if(bi==-1) {
			if(ci==-1) return false;
			return union(ci,ai,triClockwise?tri.c:tri.b);
		}
		if(ci==-1) {
			return union(ai,bi,triClockwise?tri.b:tri.c);
		}
		//System.out.println("Triangle seams to be connected via 3 of 3 corners.");
		int ac = 2*ai<c.length?ai+c.length:ai;
		int bc = 2*bi<c.length?bi+c.length:bi;
		int cc = 2*ci<c.length?ci+c.length:ci;
		if(ac-bc==1) {
			if(cc-ac==1) return remove(ac);
			if(bc-cc==1) return remove(bc);
		}
		if(bc-cc==1 && cc-ac==1) return remove(cc);
		return false;
	}
	public boolean union(JDPolygon other, double tolerance) {
		List<JDTriangle> triangles = other.toTriangles();
//		System.out.println("polygon "+other+"\nis separated in:");
//		for(JDTriangle tri: triangles)
//			System.out.println("    "+tri);
		while(!triangles.isEmpty()) {
			boolean noInsert = true;
			for(int t=triangles.size()-1; t>=0; t--) {
				boolean success = union(triangles.get(t), tolerance);
				if(success) {
					triangles.remove(t);
					noInsert = false;
				}
			}
			if(noInsert) break;
		}
		return triangles.isEmpty();
	}

	@Override
	public String toString() {
		String s = "g[";
		for(int i=0; i<c.length; i++)
			s += (i==0?"(":" - (")+c[i].x+","+c[i].y+")";
		return s + "]";
	}

	@Override
	public int hashCode() {
		if (hash != null) {
			return hash;
		}
		return 1;//hash = hash(a, b, c);
	}

	public static int hash(JDPoint a, JDPoint b, JDPoint c) {
		final int prime = 31;
		int hash = 1;
		hash = prime * hash + a.hashCode();
		hash = prime * hash + b.hashCode();
		hash = prime * hash + c.hashCode();
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		JDPoint[] A = this.c;
		JDPoint[] B = ((JDPolygon) obj).c;
		if(A.length!=B.length)
			return false;
		int len = A.length;
		for(int o=0; o<len; o++) {
			boolean is_equal = true;
			for(int i=0; i<len && is_equal; i++)
				is_equal = A[i].equals(B[(i+o)%len]);
			if(is_equal)
				return true;
		}
		return false;
	}
	public boolean equals(Object obj, double tolerance) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		JDPoint[] A = this.c;
		JDPoint[] B = ((JDPolygon) obj).c;
		if(A.length!=B.length)
			return false;
		int len = A.length;
		for(int o=0; o<len; o++) {
			boolean is_equal = true;
			for(int i=0; i<len && is_equal; i++)
				is_equal = A[i].equals(B[(i+o)%len], tolerance);
			if(is_equal)
				return true;
		}
		return false;
	}
	
	private double area(JDPoint... p) {
		double a=0d;
		for(int i=0,j=p.length-1; i<p.length; j=i++) {
			JDPoint u = p[i];
			JDPoint v = p[j];
			a += (v.x-u.x) * (v.y+u.y);
		}
		return 0.5d * a;
	}
	@SuppressWarnings("unused")
	private int containsCorner(JDPoint p) {
		for(int i=0; i<c.length; i++)
			if(c[i].equals(p))
				return i;
		return -1;
	}
	private int containsCorner(JDPoint p, double tolerance) {
		for(int i=0; i<c.length; i++)
			if(c[i].equals(p, tolerance))
				return i;
		return -1;
	}
	private boolean union(int a, int b, JDPoint p) {
		//System.out.println(" ... {UNION} try to add config: a="+a+" b="+b+" p="+p);
		if(a==b) return false;
		int q = -1;
		if(a>b) {
			if(a-b!=1) return false;
			q = a;
		} else {
			if(a!=0 || b!=c.length-1) return false;
			q = c.length;
		}
		//System.out.println(" ... {UNION} current polygon is "+this.toString());
		//System.out.println(" ... {UNION} found q="+q);
		JDPoint[] temp = new JDPoint[c.length+1];
		for(int i=0; i<temp.length; i++) {
			int j = (i>=q ? i-1 : i);
			temp[i] = new JDPoint(c[j].x, c[j].y, c[j].value);
		}
		temp[q].x     = p.x;
		temp[q].y     = p.y;
		temp[q].value = p.value;
		//System.out.println(" ... {UNION} try to override old polygon coords.");
		c = null;
		c = temp;
		//System.out.println(" ... {UNION} new polygon is "+this.toString());
		return true;
	}
	private boolean remove(int q) {
		if(c.length<4)
			return false;
		int r = q%c.length;
		JDPoint[] temp = new JDPoint[c.length-1];
		for(int i=0; i<temp.length; i++) {
			int j = (i>=r ? i+1 : i);
			temp[i] = new JDPoint(c[j].x, c[j].y, c[j].value);
		}
		c = null;
		c = temp;
		return true;
	}
}
