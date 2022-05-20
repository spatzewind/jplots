package jplots.maths;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;

import jplots.JAxis;
import jplots.helper.GeometryTools;

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
		if (abc == null || abc.length < 3) {
			// System.err.println("Cannot create polygon of less than 3 points!");
			c = new JDPoint[0];
			e = new JDEdge[0];
			idx = 0;
			return;
		}
		c = new JDPoint[abc.length];
		boolean clockwise = false; // (area(abc)<0d);
		for (int i = 0; i < abc.length; i++)
			c[clockwise ? abc.length - 1 - i : i] = new JDPoint(abc[i].x, abc[i].y, abc[i].value);
		e = new JDEdge[abc.length];
		for (int i = 0, j = e.length - 1; i < e.length; j = i++)
			e[i] = new JDEdge(c[j], c[i]);
		idx = index;
	}

	public JDPolygon(Coordinate[] coords) {
		this(coords, 0);
	}

	public JDPolygon(Coordinate[] coords, int index) {
		if (coords == null || coords.length < 4) {
			// System.err.println("Cannot create polygon of less than 3 points!");
			c = new JDPoint[0];
			e = new JDEdge[0];
			idx = 0;
			return;
		}
		c = new JDPoint[coords.length - 1];
		for (int i = 0; i < c.length; i++)
			c[i] = new JDPoint(coords[i]);
//		double a = area(c);
//		if(a<0d)
//			for(int i=0; i<c.length; i++)
//				c[i] = new JDPoint(coords[c.length-1-i]);
		e = new JDEdge[c.length];
		for (int i = 0, j = c.length - 1; i < c.length; j = i++)
			e[i] = new JDEdge(c[j], c[i]);
		idx = index;
	}

	public JDPolygon(double[]... coords) {
		this(coords, 0);
	}

	public JDPolygon(double[][] coords, int index) {
		if (coords == null || coords.length < 3) {
			// System.err.println("Cannot create polygon of less than 3 points!");
			c = new JDPoint[0];
			e = new JDEdge[0];
			idx = 0;
			return;
		}
		c = new JDPoint[coords.length];
		for (int i = 0; i < c.length; i++)
			c[i] = new JDPoint(coords[i][0], coords[i][1], 0d);
//		boolean clockwise = (area(c)<0d);
//		if(clockwise)
//			for(int i=0; i<c.length; i++)
//				c[c.length-1-i] = new JDPoint(coords[i][0], coords[i][1], 0d);
		e = new JDEdge[c.length];
		for (int i = 0, j = e.length - 1; i < e.length; j = i++)
			e[i] = new JDEdge(c[j], c[i]);
		idx = index;
	}
//	public void edges(JDEdge ab, JDEdge bc, JDEdge ca) {
//		this.ab = ab.equals(a, b) ? ab : bc.equals(a, b) ? bc : ca;
//		this.bc = ab.equals(b, c) ? ab : bc.equals(b, c) ? bc : ca;
//		this.ca = ab.equals(c, a) ? ab : bc.equals(c, a) ? bc : ca;
//	}

	public JDPolygon copy() {
		return new JDPolygon(c == null ? new JDPoint[0] : c);
	}

	public Geometry toGeometry(GeometryFactory factory) {
		if (c == null)
			return null;
		Coordinate[] ccc = new Coordinate[c.length + 1];
		for (int i = 0; i < c.length; i++)
			ccc[i] = new Coordinate(c[i].x, c[i].y, c[i].value);
		ccc[c.length] = new Coordinate(c[0].x, c[0].y, c[0].value);
		return factory.createPolygon(ccc);
	}

	public List<JDTriangle> toTriangles() {
		List<JDTriangle> triangles = new ArrayList<>();
		List<JDPoint> points = new LinkedList<>(Arrays.asList(c));
		int failCount = 0, maxFails = points.size();
		double sign = area() < 0d ? -1d : 1d;
		while (points.size() > 3) {
			int smallestIdx = -1;
			double smallestArea = Double.POSITIVE_INFINITY;
			int slen = points.size();
			for (int s0 = 0; s0 < slen - 1; s0++) {
				int s1 = (s0 + 1) % slen;
				int s2 = (s0 + 2) % slen;
				double a = sign * GeometryTools.area(points.get(s0), points.get(s1), points.get(s2));
				if (a < 0d)
					continue;
				if (a < smallestArea) {
					smallestIdx = s0;
					smallestArea = a;
				}
			}
			if (smallestIdx == -1) {
//				failCount++;
//				if(failCount>maxFails)
//					break;
//				points.add(points.get(0));
//				points.remove(0);
//				continue;
				System.out.println("[POLY] failed to triangulate: [" + c.length + " points, " + area() + " area]");
				break;
//				throw new RuntimeException("Triangulation of polygon failed!");
			}
			int next0 = (smallestIdx + 1) % slen;
			int next1 = (smallestIdx + 2) % slen;
			triangles.add(new JDTriangle(points.get(smallestIdx), points.get(next0), points.get(next1)));
			points.remove(next0);
		}
		triangles.add(new JDTriangle(points.get(0), points.get(1), points.get(2)));
		return triangles;
	}

	public boolean union(JDTriangle tri, double tolerance) {
		double ar = tri.area();
		boolean triClockwise = (ar < 0d);
		// System.out.println("area of triangle to add is "+ar);
		int ai = containsCorner(tri.getA(), tolerance);
		int bi = containsCorner(triClockwise ? tri.getC() : tri.getB(), tolerance);
		int ci = containsCorner(triClockwise ? tri.getB() : tri.getC(), tolerance);
		// System.out.println("Corner ids: a="+ai+" b="+bi+" c="+ci);
		if (ai == -1) {
			if (bi == -1 || ci == -1)
				return false;
			return union(bi, ci, tri.getA());
		}
		if (bi == -1) {
			if (ci == -1)
				return false;
			return union(ci, ai, triClockwise ? tri.getC() : tri.getB());
		}
		if (ci == -1) {
			return union(ai, bi, triClockwise ? tri.getB() : tri.getC());
		}
		// System.out.println("Triangle seams to be connected via 3 of 3 corners.");
		int ac = 2 * ai < c.length ? ai + c.length : ai;
		int bc = 2 * bi < c.length ? bi + c.length : bi;
		int cc = 2 * ci < c.length ? ci + c.length : ci;
		if (ac - bc == 1) {
			if (cc - ac == 1)
				return remove(ac);
			if (bc - cc == 1)
				return remove(bc);
		}
		if (bc - cc == 1 && cc - ac == 1)
			return remove(cc);
		return false;
	}

	public boolean union(JDPolygon other, double tolerance) {
		List<JDTriangle> triangles = other.toTriangles();
//		System.out.println("polygon "+other+"\nis separated in:");
//		for(JDTriangle tri: triangles)
//			System.out.println("    "+tri);
		while (!triangles.isEmpty()) {
			boolean noInsert = true;
			for (int t = triangles.size() - 1; t >= 0; t--) {
				boolean success = union(triangles.get(t), tolerance);
				if (success) {
					triangles.remove(t);
					noInsert = false;
				}
			}
			if (noInsert)
				break;
		}
		return triangles.isEmpty();
	}

	@Override
	public String toString() {
		String s = "g[";
		for (int i = 0; i < c.length; i++)
			s += (i == 0 ? "(" : " - (") + c[i].x + "," + c[i].y + ")";
		return s + "]";
	}

	@Override
	public int hashCode() {
		if (hash != null) {
			return hash;
		}
		return 1;// hash = hash(a, b, c);
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
		if ((obj == null) || (getClass() != obj.getClass())) {
			return false;
		}
		JDPoint[] A = this.c;
		JDPoint[] B = ((JDPolygon) obj).c;
		if (A.length != B.length)
			return false;
		int len = A.length;
		for (int o = 0; o < len; o++) {
			boolean is_equal = true;
			for (int i = 0; i < len && is_equal; i++)
				is_equal = A[i].equals(B[(i + o) % len]);
			if (is_equal)
				return true;
		}
		return false;
	}

	public boolean equals(Object obj, double tolerance) {
		if (this == obj) {
			return true;
		}
		if ((obj == null) || (getClass() != obj.getClass())) {
			return false;
		}
		JDPoint[] A = this.c;
		JDPoint[] B = ((JDPolygon) obj).c;
		if (A.length != B.length)
			return false;
		int len = A.length;
		for (int o = 0; o < len; o++) {
			boolean is_equal = true;
			for (int i = 0; i < len && is_equal; i++)
				is_equal = A[i].equals(B[(i + o) % len], tolerance);
			if (is_equal)
				return true;
		}
		return false;
	}

	public float[][] getCoords() {
		float[][] coords = new float[c.length][2];
		for (int i = 0; i < c.length; i++) {
			coords[i][0] = (float) c[i].x;
			coords[i][1] = (float) c[i].y;
		}
		return coords;
	}

	public double getDefaultTolerance() {
		double minX = Double.POSITIVE_INFINITY;
		double maxX = Double.NEGATIVE_INFINITY;
		double minY = Double.POSITIVE_INFINITY;
		double maxY = Double.NEGATIVE_INFINITY;
		for (JDPoint p : c) {
			if (p.x < minX)
				minX = p.x;
			if (p.x > maxX)
				maxX = p.x;
			if (p.y < minX)
				minY = p.y;
			if (p.y > maxX)
				maxY = p.y;
		}
		return Math.max(maxX - minX, maxY - minY) * 1.0e-10d;
	}

	/**
	 * apply affine transformation to polygon </b> (the orientation stays the same,
	 * even if the determinant of the transformation is negative)
	 * 
	 * @param tm transformation matrix (2x3)
	 * @return itself after transformation
	 */
	public JDPolygon affine(double[][] tm) {
		for (JDPoint element : c)
			element.affine(tm);
		if (tm[0][0] * tm[1][1] - tm[0][1] * tm[1][0] < 0d)
			reverse_orientation();
		return this;
	}

	public JDPolygon reverse_orientation() {
		double[] xx = new double[c.length];
		double[] yy = new double[c.length];
		double[] vv = new double[c.length];
		for (int i = 0; i < c.length; i++) {
			xx[i] = c[i].x;
			yy[i] = c[i].y;
			vv[i] = c[i].value;
		}
		for (int i = 0; i < c.length; i++) {
			c[c.length - 1 - i] = new JDPoint(xx[i], yy[i], vv[i]);
		}
		for (int i = 0, j = e.length - 1; i < e.length; j = i++)
			e[i] = new JDEdge(c[j], c[i]);
		return this;
	}

	public List<JDPolygon> splitByMapBorder(JAxis ax) {
		List<JDPolygon> res = new ArrayList<>();
		if (!ax.isGeoAxis())
			return res;
		return ax.getGeoProjection().splitByMapBorder(this);
	}

	public List<JDPolygon> intersectsAABB(double le, double to, double ri, double bt) {
		double minx = Math.min(le, ri), maxx = Math.max(le, ri);
		double miny = Math.min(to, bt), maxy = Math.max(to, bt);
		double eps = Math.max(Math.abs(le - ri), Math.abs(to - bt)) * 1.0e-10d;
//		int[] counts = { 0, 0, 0, 0 };
		double[] normal = { -1d, 0d };
		List<JDPoint[]> first = GeometryTools.SutherlandHodgmanAlgorithm(c, normal, -minx, eps);
//		counts[0] = first.size();
		List<JDPoint[]> second = new ArrayList<>();
		normal[0] = 1d;
		normal[1] = 0d;
		for (JDPoint[] p : first)
			second.addAll(GeometryTools.SutherlandHodgmanAlgorithm(p, normal, maxx, eps));
//		counts[1] = second.size();

		first.clear();
		normal[0] = 0d;
		normal[1] = -1d;
		for (JDPoint[] p : second)
			first.addAll(GeometryTools.SutherlandHodgmanAlgorithm(p, normal, -miny, eps));
//		counts[2] = first.size();
		second.clear();
		normal[0] = 0d;
		normal[1] = 1d;
		for (JDPoint[] p : first)
			second.addAll(GeometryTools.SutherlandHodgmanAlgorithm(p, normal, maxy, eps));
//		counts[3] = second.size();
//		System.out.println("[JDPOLYGON] intersectsAABB: created "+counts[0]+" -> "+counts[1]+" -> "+counts[2]+" -> "+counts[3]+" polygon(s)");

		List<JDPolygon> res = new ArrayList<>();
		for (JDPoint[] p : second)
			res.add(new JDPolygon(p));
		return res;
	}
//	public List<JDPolygon> intersectsAABB(double le, double to,  double ri, double bt) {
//		if(ri<le) return intersectsAABB(ri, to, le, bt);
//		if(bt<to) return intersectsAABB(le, bt, ri, to);
//		JDPoint center = new JDPoint(0.5d*(le+ri), 0.5d*(to+bt));
//		double ve = Math.min(ri-le, bt-to) * 1.0e-10d;
//		List<JDPolygon> res = new ArrayList<>();
//		boolean foundPointOutside = false;
//		for(int i=0; i<c.length && !foundPointOutside; i++) {
//			if(c[i].x < le) foundPointOutside = true;
//			if(c[i].x > ri) foundPointOutside = true;
//			if(c[i].y < to) foundPointOutside = true;
//			if(c[i].y > bt) foundPointOutside = true;
//		}
//		if(area()<0d && !foundPointOutside) {
//			int nearest_to_edge = -1, nearest_side = -1;
//			double nearest_x = Double.NaN, nearest_y = Double.NaN;
//			double nearest_dist = Double.POSITIVE_INFINITY, off = 100d;
//			for(int i=0; i<c.length; i++) {
//				if(c[i].x >= le && c[i].x <= ri &&
//					c[i].y >= to && c[i].y <= bt) {
//					if(c[i].x-le < nearest_dist) { nearest_to_edge = i;
//						nearest_dist = c[i].x-le; nearest_side = 1; nearest_x = c[i].x; nearest_y = c[i].y; }
//					if(ri-c[i].x < nearest_dist) { nearest_to_edge = i;
//						nearest_dist = ri-c[i].x; nearest_side = 2; nearest_x = c[i].x; nearest_y = c[i].y; }
//					if(c[i].y-to < nearest_dist) { nearest_to_edge = i;
//						nearest_dist = c[i].y-to; nearest_side = 3; nearest_x = c[i].x; nearest_y = c[i].y; }
//					if(bt-c[i].y < nearest_dist) { nearest_to_edge = i;
//						nearest_dist = bt-c[i].y; nearest_side = 4; nearest_x = c[i].x; nearest_y = c[i].y; }
//				}
//				double d = 10d+Math.max(Math.max(le-c[i].x, c[i].x-ri), Math.max(to-c[i].y, c[i].y-bt));
//				if(off<d) off = d;
//			}
//			List<JDPoint> temp = new ArrayList<>();
//			switch(nearest_side) {
//				case 1: temp.add(new JDPoint(le-off,nearest_y)); temp.add(new JDPoint(le-off,bt+off));
//						temp.add(new JDPoint(ri+off,bt+off));    temp.add(new JDPoint(ri+off,to-off));
//						temp.add(new JDPoint(le-off,to-off));    temp.add(new JDPoint(le-off,nearest_y));
//						break;
//				case 2: temp.add(new JDPoint(ri+off,nearest_y)); temp.add(new JDPoint(ri+off,to-off));
//						temp.add(new JDPoint(le-off,to-off));    temp.add(new JDPoint(le-off,bt+off));
//						temp.add(new JDPoint(ri+off,bt+off));    temp.add(new JDPoint(ri+off,nearest_y));
//						break;
//				case 3: temp.add(new JDPoint(nearest_x,to-off)); temp.add(new JDPoint(le-off,to-off));
//						temp.add(new JDPoint(le-off,bt+off));    temp.add(new JDPoint(ri+off,bt+off));
//						temp.add(new JDPoint(ri+off,to-off));    temp.add(new JDPoint(nearest_x,to-off));
//						break;
//				case 4: temp.add(new JDPoint(nearest_x,bt+off)); temp.add(new JDPoint(ri+off,bt+off));
//						temp.add(new JDPoint(ri+off,to-off));    temp.add(new JDPoint(le-off,to-off));
//						temp.add(new JDPoint(le-off,bt+off));    temp.add(new JDPoint(nearest_x,bt+off));
//						break;
//				default: return res;
//			}
//			for(int i=0; i<c.length; i++)
//				temp.add(c[(i+nearest_to_edge)%c.length]);
//			res.add(new JDPolygon(temp.toArray(new JDPoint[0]), this.idx));
//		} else {
//			res.add(this);
//			if(!foundPointOutside)
//				return res;
//		}
//		List<JDPoint> points = new ArrayList<>();
//		List<Double> cutids = new ArrayList<>();
//		for(int[] side: new int[][] {{0,-1},{0,1},{1,-1},{1,1}}) {
//			int polycount = res.size();
//			double vc = side[0]==0?(side[1]<0?-le:ri):(side[1]<0?-to:bt);
//			double[] normal = { side[0]==0?side[1]:0d, side[0]==0?0d:side[1] };
//			for(int p_i=0; p_i<polycount; p_i++) {
//				points.clear();
//				cutids.clear();
//				if(res.get(p_i)==null) continue;
//				JDPoint[] pnts = res.get(p_i).c;
//				if(pnts==null) continue;
//				for(int i=pnts.length-1,j=0; j<pnts.length; i=j++) {
//					double v1 = (side[0]==0?pnts[i].x:pnts[i].y)*side[1];
//					double v2 = (side[0]==0?pnts[j].x:pnts[j].y)*side[1];
//					if(v1<=vc) {
//						if(v2<=vc) {
//							points.add(pnts[j]);
//						} else {
//							double vf = (vc-v1) / (v2-v1);
//							points.add(pnts[i].fractionTowards(vf, pnts[j]));
//							GeometryTools.checkCut(points, cutids, normal, vc, ve, center);
//						}
//					} else {
//						if(v2<=vc) {
//							double vf = (vc-v1) / (v2-v1);
//							points.add(pnts[i].fractionTowards(vf, pnts[j]));
//							GeometryTools.checkCut(points, cutids, normal, vc, ve, center);
//							points.add(pnts[j]);
//						} else {
//							//nothing to add
//						}
//					}
//					if(points.size()>1) {
//						int k = points.size()-2;
//						int l = points.size()-1;
//						double uk = (side[0]==0 ? points.get(k).x : points.get(k).y) * side[1];
//						double ul = (side[0]==0 ? points.get(l).x : points.get(l).y) * side[1];
//						if( Math.abs(uk-vc)<ve && Math.abs(ul-vc)<ve ) {
//							boolean isNegative = GeometryTools.area(center,points.get(k),points.get(l))<0d;
//							if(isNegative) cutids.add( -1d-points.size() );
//							else {
//								boolean isLastPositive = false;
//								if(!cutids.isEmpty()) isLastPositive = cutids.get(cutids.size()-1)>0d;
//								if(!isLastPositive) cutids.add(1d+cutids.size());
//							}
//						}
//					}
//				}
//				if(points.size()>0) {
//					points.add(points.get(0));
//					GeometryTools.checkCut(points, cutids, normal, vc, ve, center);
//					points.remove(points.size()-1);
//				}
//				if(cutids.size()>0) {
//					boolean isLastPositive = false;
//					if(!cutids.isEmpty()) isLastPositive = cutids.get(cutids.size()-1)>0d;
//					if(isLastPositive) {
//						if(cutids.get(0)>0d) cutids.remove(0);
//					}
//				}
//				if(cutids.size()>0) {
//					cutids.add(cutids.get(0));
//					int pl = points.size();
//					for(int n=1; n<cutids.size(); n++) {
//						int ci1 = (int) (Math.abs(cutids.get(n-1))-0.5d);
//						int ci2 = (int) (Math.abs(cutids.get( n ))-0.5d);
//						if(ci2<ci1) ci2 += pl;
//						JDPoint[] temp = new JDPoint[ci2-ci1];
//						for(int t=0; t<temp.length; t++)
//							temp[t] = points.get((ci1+t)%pl);
//						res.add(new JDPolygon(temp));
//					}
//				} else {
//					res.add(new JDPolygon(points.toArray(new JDPoint[0])));
//				}
//			}
//			for(int p_i=polycount-1; p_i>=0; p_i--)
//				res.remove(p_i);
//		}
//		//TODO create AABB intersection
//		return res;
//	}

	public double area() {
		return GeometryTools.area(c);
	}

	@SuppressWarnings("unused")
	private int containsCorner(JDPoint p) {
		for (int i = 0; i < c.length; i++)
			if (c[i].equals(p))
				return i;
		return -1;
	}

	private int containsCorner(JDPoint p, double tolerance) {
		for (int i = 0; i < c.length; i++)
			if (c[i].equals(p, tolerance))
				return i;
		return -1;
	}

	private boolean union(int a, int b, JDPoint p) {
		// System.out.println(" ... {UNION} try to add config: a="+a+" b="+b+" p="+p);
		if (a == b)
			return false;
		int q = -1;
		if (a > b) {
			if (a - b != 1)
				return false;
			q = a;
		} else {
			if (a != 0 || b != c.length - 1)
				return false;
			q = c.length;
		}
		// System.out.println(" ... {UNION} current polygon is "+this.toString());
		// System.out.println(" ... {UNION} found q="+q);
		JDPoint[] temp = new JDPoint[c.length + 1];
		for (int i = 0; i < temp.length; i++) {
			int j = (i >= q ? i - 1 : i);
			temp[i] = new JDPoint(c[j].x, c[j].y, c[j].value);
		}
		temp[q].x = p.x;
		temp[q].y = p.y;
		temp[q].value = p.value;
		// System.out.println(" ... {UNION} try to override old polygon coords.");
		c = null;
		c = temp;
		// System.out.println(" ... {UNION} new polygon is "+this.toString());
		return true;
	}

	private boolean remove(int q) {
		if (c.length < 4)
			return false;
		int r = q % c.length;
		JDPoint[] temp = new JDPoint[c.length - 1];
		for (int i = 0; i < temp.length; i++) {
			int j = (i >= r ? i + 1 : i);
			temp[i] = new JDPoint(c[j].x, c[j].y, c[j].value);
		}
		c = null;
		c = temp;
		return true;
	}
}
