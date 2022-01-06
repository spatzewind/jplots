package jplots.maths;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class JGridTriangulator {

	public int[] triangles;
	public int[] halfedges;
	public JDPoint[] points;

	private double[] coords;
	public int[] hull;
	
	private int xlen;
	private int ylen;

	public JGridTriangulator(double[] x, double[] y, double[][] z) {
		if (x.length < 2 || y.length < 2) {
			System.err.println("Need at least 2x2 points");
			return;
		}

		xlen = x.length;
		ylen = y.length;
		points = new JDPoint[xlen*ylen];
		coords = new double[xlen*ylen*2];
		for(int j=0; j<ylen; j++) {
			for(int i=0; i<xlen; i++) {
				int p = j*xlen + i;
				points[p] = new JDPoint(x[i], y[j], z[j][i]);
				coords[2*p+0] = x[i];
				coords[2*p+1] = y[j];
			}
		}

		//create convex hull
		findConvexHull();

		//generate delaunay triangulation
		generate();
	}


	private List<JDTriangle> trias = null;
	private List<JDEdge> edges = null;
	private List<JDPoint> poinz = null;
	private List<JDEdge> hulls = null;
	private List<JDEdge> voron = null;
	private List<JDEdge> vhull = null;

	public List<JDTriangle> getTriangles() {
		return trias;
	}

	public List<JDEdge> getEdges() {
		return edges;
	}

	public List<JDPoint> getPoints() {
		return poinz;
	}

	public List<JDEdge> getHullEdges() {
		return hulls;
	}

	public List<JDEdge> getVoronoiEdges() {
		return voron;
	}

	public List<JDEdge> getVoronoiHullEdges() {
		return vhull;
	}

	private void generate() {
		Map<JDEdge, JDEdge> unq_edges = new HashMap<>();
		Map<JDTriangle, JDTriangle> unq_trias = new HashMap<>();

		for (int n = 0; triangles != null && n < triangles.length / 3; n++) {
			int i = triangles[3 * n + 0];
			int j = triangles[3 * n + 1];
			int k = triangles[3 * n + 2];
			JDPoint a = points[i];
			JDPoint b = points[j];
			JDPoint c = points[k];

			JDPoint[] tmp = { a, b, c };
			Arrays.sort(tmp);
			a = tmp[0];
			b = tmp[1];
			c = tmp[2];

			JDEdge ab = new JDEdge(a, b);
			JDEdge bc = new JDEdge(b, c);
			JDEdge ca = new JDEdge(c, a);

			if (unq_edges.containsKey(ab)) {
				ab = unq_edges.get(ab);
			} else {
				unq_edges.put(ab, ab);
			}

			if (unq_edges.containsKey(bc)) {
				bc = unq_edges.get(bc);
			} else {
				unq_edges.put(bc, bc);
			}

			if (unq_edges.containsKey(ca)) {
				ca = unq_edges.get(ca);
			} else {
				unq_edges.put(ca, ca);
			}

			JDTriangle t = new JDTriangle(a, b, c);
			if (unq_trias.containsKey(t)) {
				System.err.println("[ERR] duplicated triangle " + t + " vs " + unq_trias.get(t));
				t = unq_trias.get(t);
			} else {
				t.idx = unq_trias.size();
				unq_trias.put(t, t);
			}

			t.edges(ab, bc, ca);
			ab.wing(t);
			bc.wing(t);
			ca.wing(t);
		}

		////

		List<JDEdge> hulledges = new ArrayList<>();
		for (int i = 0; i < hull.length; i++) {
			JDPoint a = points[hull[i]];
			JDPoint b = points[hull[(i + 1) % hull.length]];
			JDEdge e = new JDEdge(a, b);

			if (unq_edges.containsKey(e)) {
				e = unq_edges.get(e);
			} else {
				System.err.println("[ERR] cannot found valid edge " + e);
			}
			hulledges.add(e);
		}

		////

		List<JDEdge> voronoiedges = new ArrayList<>();
		List<JDEdge> voronoihulledges = new ArrayList<>();
//		Set<Integer> tmp = new LinkedHashSet<Integer>();
//		for (int i = 0; i < triangles.length; i++) {
//			int id = triangles[nextHalfEdge(i)];
//			if (!tmp.contains(id)) {
//				tmp.add(id);
//				List<Integer> edges = edgesAroundPoint(i);
//				List<JDPoint> point = new ArrayList<>();
//
//				for (int j = 0; j < edges.size(); j++) {
//					JDPoint[] pnt = getTrianglePoints(triangleOfEdge(edges.get(j)));
//					JDPoint cen = getCentroid(pnt);
//					if (cen == null) {
//						continue;
//					}
//					point.add(cen);
//				}
//
//				for (int j = 0; j < point.size(); j++) {
//					JDPoint a = point.get(j);
//					JDPoint b = point.get((j + 1) % point.size());
//					JDEdge e = new JDEdge(a, b);
//
//					voronoiedges.add(e);
//				}
//			}
//		}

		this.trias = new ArrayList<JDTriangle>(unq_trias.values());
		this.edges = new ArrayList<JDEdge>(unq_edges.values());
		this.poinz = Arrays.asList(points);
		this.hulls = hulledges;
		this.voron = voronoiedges;
		this.vhull = hulledges; // voronoihulledges; //TODO it must be implement !!!
	}

	///////////////////////////////////////////////////////////////////////////

	private void findConvexHull() {
		//create triangle set
		int count = 2 * (xlen-1) * (ylen-1);
		triangles = new int[3 * count];
		for(int j=0; j<ylen-1; j++)
			for(int i=0; i<xlen-1; i++) {
				int n = 2 * ( j * (xlen-1) +  i );
				int p = j * xlen + i;
				boolean mirror = Double.isNaN(points[p+1].value) || Double.isNaN(points[p+xlen].value);
				if(mirror && (Double.isNaN(points[p].value) || Double.isNaN(points[p+xlen+1].value)))
					mirror = false;
				if(Math.abs(points[p].value-points[p+xlen+1].value)<Math.abs(points[p+1].value-points[p+xlen].value))
					mirror = true;
				if(mirror) {
					triangles[3 * n + 0] = p;
					triangles[3 * n + 1] = p        + 1;
					triangles[3 * n + 2] = p + xlen + 1;
					n += 1;
					triangles[3 * n + 0] = p;
					triangles[3 * n + 1] = p + xlen + 1;
					triangles[3 * n + 2] = p + xlen;
				} else {
					triangles[3 * n + 0] = p;
					triangles[3 * n + 1] = p        + 1;
					triangles[3 * n + 2] = p + xlen;
					n += 1;
					triangles[3 * n + 0] = p + xlen;
					triangles[3 * n + 1] = p        + 1;
					triangles[3 * n + 2] = p + xlen + 1;
				}
			}
		
		//create hull set
		count = 2 * ( ( xlen - 1 ) + ( ylen - 1) );
		hull = new int[count];
		int halfcount = count / 2;
		for(int i=0; i<xlen-1; i++) {
			int p = i;
			int q = xlen*ylen - 1 - i;
			hull[i]             = p;
			hull[i + halfcount] = q;
		}
		for(int j=0; j<ylen-1; j++) {
			int p = j * xlen + xlen - 1;
			int q = (ylen-1 - j) * xlen;
			hull[xlen-1 + j]             = p;
			hull[xlen-1 + j + halfcount] = q;
		}
	}

	///////////////////////////////////////////////////////////////////////////

	public int[] edgesOfTriangle(int t) {
		return new int[] { 3 * t, 3 * t + 1, 3 * t + 2 };
	}

	public int triangleOfEdge(int e) {
		return (int) Math.floor(e / 3);
	}

	public int nextHalfEdge(int e) {
		return (e % 3 == 2) ? e - 2 : e + 1;
	}

	public int prevHalfEdge(int e) {
		return (e % 3 == 0) ? e + 2 : e - 1;
	}

	public int[] pointsOfTriangle(int t) {
		int[] e = edgesOfTriangle(t);
		int a = triangles[e[0]];
		int b = triangles[e[1]];
		int c = triangles[e[2]];
		return new int[] { a, b, c };
	}

	public JDPoint[] getTrianglePoints(int t) {
		int[] p = pointsOfTriangle(t);
		JDPoint a = points[p[0]];
		JDPoint b = points[p[1]];
		JDPoint c = points[p[2]];

		return new JDPoint[] { a, b, c };
	}

	public JDPoint[] getEdgesOfTriangle(int t) {
		int[] e = edgesOfTriangle(t);
		JDPoint a = points[e[0]];
		JDPoint b = points[e[1]];
		JDPoint c = points[e[2]];
		return new JDPoint[] { a, b, c };
	}

	public JDPoint[] GetHullPoints() {
		JDPoint[] hullpoint = new JDPoint[hull.length];
		for (int i = 0; i < hull.length; i++) {
			hullpoint[i] = points[hull[i]];
		}
		return hullpoint;
	}

	public JDPoint getCentroid(int t) {
		JDPoint[] vertices = getTrianglePoints(t);
		return getCentroid(vertices);
	}

	public JDPoint getCentroid(JDPoint[] points) {
		double accumulatedArea = 0.0f;
		double centerX = 0.0d;
		double centerY = 0.0d;
		double centerV = 0.0d;

		for (int i = 0, j = points.length - 1; i < points.length; j = i++) {
			double temp = points[i].x * points[j].y - points[j].x * points[i].y;
			accumulatedArea += temp;
			centerX += (points[i].x + points[j].x) * temp;
			centerY += (points[i].y + points[j].y) * temp;
			centerV += points[i].value;
		}

		if (Math.abs(accumulatedArea) < 1E-7f) {
			return null; // new JDPoint(0, 0, Double.NaN); // ???
		}

		accumulatedArea *= 3d;
		return new JDPoint(centerX / accumulatedArea, centerY / accumulatedArea, centerV / points.length);
	}

	public List<Integer> edgesAroundPoint(int start) {
		List<Integer> lst = new ArrayList<>();
		int incoming = start;
		do {
			lst.add(incoming);
			int outgoing = nextHalfEdge(incoming);
			incoming = halfedges[outgoing];
		} while (incoming != -1 && incoming != start);

		return lst;
	}

	public List<Integer> trianglesAdjacentToTriangle(int t) {
		List<Integer> adjacentTriangles = new ArrayList<>();
		int[] triangleEdges = edgesOfTriangle(t);

		for (int i = 0; i < triangleEdges.length; i++) {
			int e = triangleEdges[i];
			int opposite = halfedges[e];
			if (opposite >= 0) {
				adjacentTriangles.add(triangleOfEdge(opposite));
			}
		}
		return adjacentTriangles;
	}
}
