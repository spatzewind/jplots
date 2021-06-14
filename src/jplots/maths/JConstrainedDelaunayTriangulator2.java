package jplots.maths;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import processing.data.IntList;

public class JConstrainedDelaunayTriangulator2 {

	private static double EPSILON  = Math.pow(2, -52);
	private static double EPSILON2 = Math.pow(2, -52);
	private static final int    SUCCESS        = 0;
	private static final int    ERROR_CENTER   = 1;
	private int[] EDGE_STACK = new int[512];

	public int[] triangles;
	public int[] halfedges;
	public JDPoint[] ep;
	public JDPoint[] points;
	public JDEdge[] constrain;

	private int hashSize;
	private int[] hullPrev;
	private int[] hullNext;
	private int[] hullTria;
	private int[] hullHash;

	private double minX=-1d, maxX=1d, minY=-1d, maxY=1d;
	private double i0x,i0y, i1x,i1y, i2x,i2y;
	private double cx, cy;
	private int   i0, i1, i2;
	private int[] ids;
	private JDPoint center;

	private int trianglesLen;
	private double[] coords;
	private int hullStart;
	private int hullSize;
	public int[] hull;

	public JConstrainedDelaunayTriangulator2(List<JDEdge> constrains, List<JDPoint> extra_points) {
		this(uniquee(constrains), uniquep(extra_points));
	}

	public JConstrainedDelaunayTriangulator2(JDEdge[] constrains, JDPoint[] extra_points) {
		double range = Math.max(Math.max(-minX,maxX), Math.max(-minY, maxY));
		range = Math.max(range, Math.max(maxX-minX, maxY-minY));
		EPSILON  =    range    * 1.d-10d;
		EPSILON2 = range*range * 1.d-15d;
		ep = extra_points;
		merge(constrains, null, extra_points, null);
		
		//initialise and find seed-point;
		initialise();
		int n = this.points.length;

		//find seed as point closest to mass-center of points
		int err = SUCCESS;
		if((err=findCenter())!=SUCCESS) {
			System.err.print(err);
			return;
		}

		// sort the points by distance from the seed triangle circumcenter
		double[] dists = new double[n];
		for (int i = 0; i < n; i++) {
			dists[i] = dist(coords[2 * i], coords[2 * i + 1], center.x, center.y);
		}
		quicksort(ids, dists, 0, n - 1);

		//create convex hull
		//check delaunay condition and flip edges automatically
		findConvexHull();

		//remember the constrains
		//flip edges accordingly to satisfy the constrains
		System.out.println("all edges:");
		for(JDEdge e: constrain)
			System.out.println("    "+e);
		int iter = 0;
		boolean checkConstrains = true;
		while(checkConstrains && iter<2) {
			checkConstrains = false;
			iter++;
			if(checkConstrainCrossings()) {
				initialise();
				err = findCenter();
				if(err!=SUCCESS) {
					System.err.print(err);
					return;
				}
				n = this.points.length;
				dists = new double[n];
				for (int i = 0; i < n; i++) {
					dists[i] = dist(coords[2 * i], coords[2 * i + 1], center.x, center.y);
				}
				quicksort(ids, dists, 0, n - 1);
				findConvexHull();
				checkConstrains = true;
			}
		}
		System.out.println("all edges:");
		for(JDEdge e: constrain)
			System.out.println("    "+e);

		// generates Lists of
		// -- all triangles as JDTriangle with JDPoint's and JDEdge's
		// -- all edges of the convex hull as JDEdge with JDPoint's
		// -- voronoi edges as JDEdge with JDPoint's
		// -- voronoi hull edges (must be implement) as JDEdge with JDPoint's
		generate();

//		initialise();
//		
//		generateConvexHull();
//		
//		addConstrains();
//		
//		orderEdgesViaHashtable();
//		
//		makeDelaunayOfNonConstrainedEdges();
//		
//		generate();
	}

	///////////////////////////////////////////////////////////////////////////

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

	public static List<JDPoint> convert(List<Point2D> points) {
		List<JDPoint> lst = new ArrayList<>();
		points.forEach(p -> lst.add(new JDPoint(p.getX(), p.getY(), Double.NaN)));
		return lst;
	}

	/**
	 * make every point unique with offsetting by as small amount when overlap
	 * 
	 * @param points list of JDPoint objects to make unique
	 * @return
	 */
	public static JDPoint[] uniquep(List<JDPoint> points) {
		Set<JDPoint> unq = new LinkedHashSet<>();
		for (int i = 0; i < points.size(); i++) {
			JDPoint v = points.get(i);
			if (unq.contains(v)) {
				while (unq.contains(v) == false) {
					System.err.printf("[INF] found duplicated point (%f, %f), fix it will be plus +1e-6... \n", (float) v.x, (float) v.y);
					v = new JDPoint(v.x() + 1e-6, v.y() + 1e-6, v.value);
				}
			}
			unq.add(v);
		}

		return unq.toArray(new JDPoint[] {});
	}

	/**
	 * make every edge unique with offsetting points by as small amount when overlap
	 * 
	 * @param edges list of JDEdge objects to make unique
	 * @return
	 */
	public static JDEdge[] uniquee(List<JDEdge> edges) {
		Set<JDEdge> unq = new LinkedHashSet<>();
		for (int i = 0; i < edges.size(); i++) {
			JDEdge e = edges.get(i);
			if (unq.contains(e)) {
				while (unq.contains(e) == false) {
					System.err.printf("[INF] found duplicated edge (%f/%f, %f/%f), fix it will be plus +1e-6... \n", (float) e.a.x, (float) e.a.y, (float) e.b.x, (float) e.b.y);
					e = new JDEdge(new JDPoint(e.a.x() + 1e-6, e.a.y() + 1e-6, e.a.value), new JDPoint(e.b.x() + 1e-6, e.b.y() + 1e-6, e.b.value));
				}
			}
			unq.add(e);
		}

		return unq.toArray(new JDEdge[] {});
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
		Set<Integer> tmp = new LinkedHashSet<Integer>();
		for (int i = 0; i < triangles.length; i++) {
			int id = triangles[nextHalfEdge(i)];
			if (!tmp.contains(id)) {
				tmp.add(id);
				List<Integer> edges = edgesAroundPoint(i);
				List<JDPoint> point = new ArrayList<>();

				for (int j = 0; j < edges.size(); j++) {
					JDPoint[] pnt = getTrianglePoints(triangleOfEdge(edges.get(j)));
					JDPoint cen = getCentroid(pnt);
					if (cen == null) {
						continue;
					}
					point.add(cen);
				}

				for (int j = 0; j < point.size(); j++) {
					JDPoint a = point.get(j);
					JDPoint b = point.get((j + 1) % point.size());
					JDEdge e = new JDEdge(a, b);

					voronoiedges.add(e);
				}
			}
		}

		this.trias = new ArrayList<JDTriangle>(unq_trias.values());
		this.edges = new ArrayList<JDEdge>(unq_edges.values());
		this.poinz = Arrays.asList(points);
		this.hulls = hulledges;
		this.voron = voronoiedges;
		this.vhull = hulledges; // voronoihulledges; //it must be implement !!!
	}

	///////////////////////////////////////////////////////////////////////////

	private int legalize(int a) {
		int i = 0;
		int ar;

		// recursion eliminated with a fixed-size stack
		while (true) {
			int b = halfedges[a];

			/*
			 * if the pair of triangles doesn't satisfy the Delaunay condition (p1 is inside
			 * the circumcircle of [p0, pl, pr]), flip them, then do the same check/flip
			 * recursively for the new pair of triangles
			 *
			 *       pl                  pl
			 *      /||\                /  \
			 *   al/ || \bl          al/    \a
			 *    /  ||  \            /      \
			 *   /  a||b  \   flip   /___ar___\
			 * p0\   ||   /p1  =>  p0\---bl---/p1
			 *    \  ||  /            \      /
			 *   ar\ || /br           b\    /br
			 *      \||/                \  /
			 *       pr                  pr
			 */
			int a0 = a - a % 3;
			ar = a0 + (a + 2) % 3;
			// ar = a - (a%3) + (a+2)%3 = a + { a%3==0 ? 2 : -1 }

			if (b == -1) { // convex hull edge
				if (i == 0) {
					break;
				}
				a = EDGE_STACK[--i];
				continue;
			}

			int b0 = b - b % 3;
			int al = a0 + (a + 1) % 3;
			int bl = b0 + (b + 2) % 3;

			int p0 = triangles[ar];
			int pr = triangles[a];
			int pl = triangles[al];
			int p1 = triangles[bl];

			boolean illegal = inCircle(coords[2 * p0], coords[2 * p0 + 1], coords[2 * pr], coords[2 * pr + 1], coords[2 * pl],
						coords[2 * pl + 1], coords[2 * p1], coords[2 * p1 + 1]);

			if (illegal) {
				triangles[a] = p1;
				triangles[b] = p0;

				int hbl = halfedges[bl];

				// edge swapped on the other side of the hull (rare); fix the halfedge reference
				if (hbl == -1) {
					int e = hullStart;
					do {
						if (hullTria[e] == bl) {
							hullTria[e] = a;
							break;
						}
						e = hullPrev[e];
					} while (e != hullStart);
				}
				link(a, hbl);
				link(b, halfedges[ar]);
				link(ar, bl);

				int br = b0 + (b + 1) % 3;

				// don't worry about hitting the cap: it can only happen on extremely degenerate
				// input
				if (i < EDGE_STACK.length) {
					EDGE_STACK[i++] = br;
				}
			} else {
				if (i == 0) {
					break;
				}
				a = EDGE_STACK[--i];
			}
		}

		return ar;
	}
	private void legalize(double ax, double ay, double bx, double by) {
		/*
		 * if the pair of triangles doesn't satisfy the Delaunay condition (p1 is inside
		 * the circumcircle of [p0, pl, pr]), flip them, then do the same check/flip
		 * recursively for the new pair of triangles
		 *
		 *       pl                  pl
		 *      /||\                /  \
		 *   al/ || \bl          al/    \a
		 *    /  ||  \            /      \
		 *   /  a||b  \   flip   /___ar___\
		 * p0\   ||   /p1  =>  p0\---bl---/p1
		 *    \  ||  /            \      /
		 *   ar\ || /br           b\    /br
		 *      \||/                \  /
		 *       pr                  pr
		 */
		boolean foundCrossing = true;
		int iter = 0;
		while(foundCrossing && iter<200) {
			iter++;
			foundCrossing = false;
			int[] halfs = iter<200 ? halfedges : shuffle(halfedges);
			for(int a: halfs) {
				if(a<0)
					continue;
				int b = halfedges[a];
				if(b<a)
					continue;

				int a0 = a - a % 3;
				int b0 = b - b % 3;
				int al = a0 + (a + 1) % 3;
				int bl = b0 + (b + 2) % 3;
				int ar = a0 + (a + 2) % 3;

				int p0 = triangles[ar];
				int pr = triangles[a];
				int pl = triangles[al];
				int p1 = triangles[bl];

				boolean illegal = cross(ax, ay, bx, by, coords[2 * pl], coords[2 * pl + 1], coords[2 * pr], coords[2 * pr + 1]);
				
				if(illegal) {
					foundCrossing = true;
					if(inLine(coords[2 * p0], coords[2 * p0 + 1], coords[2 * pl], coords[2 * pl + 1], coords[2 * p1], coords[2 * p1 + 1]) ||
							inLine(coords[2 * p0], coords[2 * p0 + 1], coords[2 * pr], coords[2 * pr + 1], coords[2 * p1], coords[2 * p1 + 1]))
						illegal = false;
					if(illegal && (!orient(coords[2 * p0], coords[2 * p0 + 1], coords[2 * pl], coords[2 * pl + 1], coords[2 * p1], coords[2 * p1 + 1]) &&
							!orient(coords[2 * p0], coords[2 * p0 + 1], coords[2 * p1], coords[2 * p1 + 1], coords[2 * pr], coords[2 * pr + 1])))
						illegal = false;
				}

				if (illegal) {
					triangles[a] = p1;
					triangles[b] = p0;

					int hbl = halfedges[bl];

					// edge swapped on the other side of the hull (rare); fix the halfedge reference
					if (hbl == -1) {
						int e = hullStart;
						do {
							if (hullTria[e] == bl) {
								hullTria[e] = a;
								break;
							}
							e = hullPrev[e];
						} while (e != hullStart);
					}
					link(a, hbl);
					link(b, halfedges[ar]);
					link(ar, bl);
				}
			}
		}
		if(iter>1)
			System.out.println("iter="+iter+"  for  constrain ["+ax+", "+ay+", "+bx+", "+by+"]");
	}
	private int addTriangle(int i0, int i1, int i2, int a, int b, int c) {
		int t = trianglesLen;

		triangles[t] = i0;
		triangles[t + 1] = i1;
		triangles[t + 2] = i2;

		link(t, a);
		link(t + 1, b);
		link(t + 2, c);

		trianglesLen += 3;
		return t;
	}

	private void link(int a, int b) {
		halfedges[a] = b;
		if (b != -1) {
			halfedges[b] = a;
		}
	}

	private int hashKey(double x, double y) {
		return (int) (Math.floor(pseudoAngle(x - cx, y - cy) * hashSize) % hashSize);
	}

	private static boolean inCircle(double ax, double ay, double bx, double by, double cx, double cy, double px, double py) {
		double dx = ax - px;
		double dy = ay - py;
		double ex = bx - px;
		double ey = by - py;
		double fx = cx - px;
		double fy = cy - py;

		double ap = dx * dx + dy * dy;
		double bp = ex * ex + ey * ey;
		double cp = fx * fx + fy * fy;

		return dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx) < 0;
	}

	private static boolean cross(double asx, double asy, double aex, double aey, double bsx, double bsy, double bex, double bey) {
		double adx = aex - asx;
		double ady = aey - asy;
		double bdx = bex - bsx;
		double bdy = bey - bsy;
		double ppx = bsx - asx;
		double ppy = bsy - asy;
		
		if( Math.abs(bsx-asx)<EPSILON && Math.abs(bsy-asy)<EPSILON)
			return false;
		if( Math.abs(bsx-aex)<EPSILON && Math.abs(bsy-aey)<EPSILON)
			return false;
		if( Math.abs(bex-asx)<EPSILON && Math.abs(bey-asy)<EPSILON)
			return false;
		if( Math.abs(bex-aex)<EPSILON && Math.abs(bey-aey)<EPSILON)
			return false;
		
		double u = ( ppx*ady - ppy*adx ) / ( adx*bdy - ady*bdx );
		double v = ( ppx*bdy - ppy*bdx ) / ( adx*bdy - ady*bdx );
		
		if( 0d<u && u<1d && 0d<v && v<1d )
			return true;
		
		return false;
	}
	private static boolean inLine(double ax, double ay, double bx, double by, double cx, double cy) {
		double abx = ax - bx;
		double aby = ay - by;
		double bcx = cx - bx;
		double bcy = cy - by;
		
		double s = Math.abs(bcx*aby - abx*bcy);
		
		if(Math.abs(s)<EPSILON2)
			return true;
		
		return false;
	}

	private static double pseudoAngle(double dx, double dy) {
		double p = dx / (Math.abs(dx) + Math.abs(dy));
		return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
	}

	private static void quicksort(int[] ids, double[] dists, int left, int right) {
		if (right - left <= 20) {
			for (int i = left + 1; i <= right; i++) {
				int temp = ids[i];
				double tempDist = dists[temp];
				int j = i - 1;
				while (j >= left && dists[ids[j]] > tempDist) {
					ids[j + 1] = ids[j--];
				}
				ids[j + 1] = temp;
			}
		} else {
			int median = (left + right) >> 1;
			int i = left + 1;
			int j = right;
			swap(ids, median, i);
			if (dists[ids[left]] > dists[ids[right]]) {
				swap(ids, left, right);
			}
			if (dists[ids[i]] > dists[ids[right]]) {
				swap(ids, i, right);
			}
			if (dists[ids[left]] > dists[ids[i]]) {
				swap(ids, left, i);
			}

			int temp = ids[i];
			double tempDist = dists[temp];
			while (true) {
				do {
					i++;
				} while (dists[ids[i]] < tempDist);
				do {
					j--;
				} while (dists[ids[j]] > tempDist);
				if (j < i) {
					break;
				}
				swap(ids, i, j);
			}
			ids[left + 1] = ids[j];
			ids[j] = temp;

			if (right - i + 1 >= j - left) {
				quicksort(ids, dists, i, right);
				quicksort(ids, dists, left, j - 1);
			} else {
				quicksort(ids, dists, left, j - 1);
				quicksort(ids, dists, i, right);
			}
		}
	}

	private static int[] shuffle(int[] arr) {
		int[] shuf = new int[arr.length];
		Random rand = new Random();
		for(int s=0; s<arr.length; s++)
			shuf[s] = arr[s];
		int a = (int) (0.5d*(arr.length + Math.sqrt(arr.length)));
		for(int t=0; t*(t-1)<a; t++)
			if(rand.nextDouble()>0.5d) {
				int s = rand.nextInt(arr.length);
				int e = rand.nextInt(arr.length);
				swap(shuf, s, (s+e)%arr.length);
			}
		return shuf;
	}
	private static void swap(int[] arr, int i, int j) {
		int tmp = arr[i];
		arr[i] = arr[j];
		arr[j] = tmp;
	}

	private static boolean orient(double px, double py, double qx, double qy, double rx, double ry) {
		return (qy - py) * (rx - qx) - (qx - px) * (ry - qy) < 0;
	}

	private static double circumradius(double ax, double ay, double bx, double by, double cx, double cy) {
		double dx = bx - ax;
		double dy = by - ay;
		double ex = cx - ax;
		double ey = cy - ay;
		double bl = dx * dx + dy * dy;
		double cl = ex * ex + ey * ey;
		double d = 0.5 / (dx * ey - dy * ex);
		double x = (ey * bl - dy * cl) * d;
		double y = (dx * cl - ex * bl) * d;
		return x * x + y * y;
	}

	protected static JDPoint circumcenter(double ax, double ay, double bx, double by, double cx, double cy) {
		double dx = bx - ax;
		double dy = by - ay;
		double ex = cx - ax;
		double ey = cy - ay;
		double bl = dx * dx + dy * dy;
		double cl = ex * ex + ey * ey;
		double d = 0.5 / (dx * ey - dy * ex);
		double x = ax + (ey * bl - dy * cl) * d;
		double y = ay + (dx * cl - ex * bl) * d;

		return new JDPoint(x, y, Double.NaN);
	}

	private static double dist(double ax, double ay, double bx, double by) {
		double dx = ax - bx;
		double dy = ay - by;
		return dx * dx + dy * dy;
	}

	private int findID(JDPoint p) {
		for(int i=0; i<points.length; i++)
			if(p.x==points[i].x && p.y==points[i].y)
				return i;
		return -1;
	}

	///////////////////////////////////////////////////////////////////////////


	private void merge(JDEdge[] edges1, List<JDEdge> edges2, JDPoint[] points1, List<JDPoint> points2) {
		List<JDPoint> pointsN = new ArrayList<>();
		List<JDEdge> edgesN   = new ArrayList<>();
		if(edges1 != null)
			for(JDEdge e: edges1) {
				pointsN.add(e.a); pointsN.add(e.b); edgesN.add(e); }
		if(edges2 != null)
			for(JDEdge e: edges2) {
				pointsN.add(e.a); pointsN.add(e.b); edgesN.add(e); }
		if(points1 != null)
			for(JDPoint p: points1)
				pointsN.add(p);
		if(points2 != null)
			for(JDPoint p: points2)
				pointsN.add(p);
		
		//remove duplicates
		for(int p1=pointsN.size()-1; p1>0; p1--)
			for(int p2=p1-1; p2>=0; p2--)
				if(pointsN.get(p1).equals(pointsN.get(p2))) {
					pointsN.remove(p1); break; }
		
		if (pointsN.size() < 3) {
			System.err.println("Need at least 3 points");
			points = new JDPoint[0];
			constrain = new JDEdge[0];
			return;
		}
		
		points = uniquep(pointsN);
		constrain = uniquee(edgesN);
	}
	private void initialise() {
		this.coords = new double[points.length * 2];

		for (int i = 0; i < points.length; i++) {
			JDPoint p = points[i];
			coords[2 * i] = p.x;
			coords[2 * i + 1] = p.y;
		}

		int n = coords.length >> 1;
		int maxTriangles = 2 * n - 5;

		triangles = new int[maxTriangles * 3];

		halfedges = new int[maxTriangles * 3];
		hashSize = (int) Math.ceil(Math.sqrt(n));

		hullPrev = new int[n];
		hullNext = new int[n];
		hullTria = new int[n];
		hullHash = new int[hashSize];

		ids = new int[n];

		minX = Double.MAX_VALUE;
		minY = Double.MAX_VALUE;
		maxX = -Double.MAX_VALUE;
		maxY = -Double.MAX_VALUE;

		for (int i = 0; i < n; i++) {
			double x = coords[2 * i];
			double y = coords[2 * i + 1];
			if (x < minX) {
				minX = x;
			}
			if (y < minY) {
				minY = y;
			}
			if (x > maxX) {
				maxX = x;
			}
			if (y > maxY) {
				maxY = y;
			}
			ids[i] = i;
		}
	}

	private int findCenter() {
		int n = points.length;

		double cx = (minX + maxX) / 2;
		double cy = (minY + maxY) / 2;

		double minDist = Double.MAX_VALUE;
		i0 = 0;
		i1 = 0;
		i2 = 0;

		for (int i = 0; i < n; i++) {// pick a seed point close to the center
			double d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);
			if (d < minDist) {
				i0 = i;
				minDist = d;
			}
		}
		i0x = coords[2 * i0];
		i0y = coords[2 * i0 + 1];

		minDist = Double.MAX_VALUE;

		for (int i = 0; i < n; i++) {// find the point closest to the seed
			if (i == i0) {
				continue;
			}
			double d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
			if (d < minDist && d > 0) {
				i1 = i;
				minDist = d;
			}
		}

		i1x = coords[2 * i1];
		i1y = coords[2 * i1 + 1];
		double minRadius = Double.MAX_VALUE;

		for (int i = 0; i < n; i++) {// find the third point which forms the smallest circumcircle with the first two
			if (i == i0 || i == i1) {
				continue;
			}
			double r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);
			if (r < minRadius) {
				i2 = i;
				minRadius = r;
			}
		}
		i2x = coords[2 * i2];
		i2y = coords[2 * i2 + 1];

		if (minRadius == Double.MAX_VALUE) {
			System.err.println("No Delaunay triangulation exists for this input.");
			return ERROR_CENTER;
		}

		if (orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
			int i = i1;
			double x = i1x;
			double y = i1y;
			i1 = i2;
			i1x = i2x;
			i1y = i2y;
			i2 = i;
			i2x = x;
			i2y = y;
		}

		center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
		this.cx = center.x;
		this.cy = center.y;
		
		return SUCCESS;
	}
	private void findConvexHull() {
		// set up the seed triangle as the starting hull
		hullStart = i0;
		hullSize = 3;

		hullNext[i0] = hullPrev[i2] = i1;
		hullNext[i1] = hullPrev[i0] = i2;
		hullNext[i2] = hullPrev[i1] = i0;

		hullTria[i0] = 0;
		hullTria[i1] = 1;
		hullTria[i2] = 2;

		hullHash[hashKey(i0x, i0y)] = i0;
		hullHash[hashKey(i1x, i1y)] = i1;
		hullHash[hashKey(i2x, i2y)] = i2;

		trianglesLen = 0;
		addTriangle(i0, i1, i2, -1, -1, -1);

		double xp = 0;
		double yp = 0;

		for (int k = 0; k < ids.length; k++) {
			int i = ids[k];
			double x = coords[2 * i];
			double y = coords[2 * i + 1];

			// skip near-duplicate points
			if (k > 0 && Math.abs(x - xp) <= EPSILON && Math.abs(y - yp) <= EPSILON) {
				continue;
			}
			xp = x;
			yp = y;

			// skip seed triangle points
			if (i == i0 || i == i1 || i == i2) {
				continue;
			}

			// find a visible edge on the convex hull using edge hash
			int start = 0;
			for (int j = 0; j < hashSize; j++) {
				int key = hashKey(x, y);
				start = hullHash[(key + j) % hashSize];
				if (start != -1 && start != hullNext[start]) {
					break;
				}
			}

			start = hullPrev[start];
			int e = start;
			int q = hullNext[e];

			while (!orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) {
				e = q;
				if (e == start) {
					e = Integer.MAX_VALUE;
					break;
				}

				q = hullNext[e];
			}

			if (e == Integer.MAX_VALUE) {
				continue; // likely a near-duplicate point; skip it
			}

			// add the first triangle from the point
			int t = addTriangle(e, i, hullNext[e], -1, -1, hullTria[e]);

			// recursively flip triangles from the point until they satisfy the Delaunay
			// condition
			hullTria[i] = legalize(t + 2);
			hullTria[e] = t; // keep track of boundary triangles on the hull
			hullSize++;

			// walk forward through the hull, adding more triangles and flipping recursively
			int next = hullNext[e];
			q = hullNext[next];

			while (orient(x, y, coords[2 * next], coords[2 * next + 1], coords[2 * q], coords[2 * q + 1])) {
				t = addTriangle(next, i, q, hullTria[i], -1, hullTria[next]);
				hullTria[i] = legalize(t + 2);
				hullNext[next] = next; // mark as removed
				hullSize--;
				next = q;

				q = hullNext[next];
			}

			if (e == start) {// walk backward from the other side, adding more triangles and flipping
				q = hullPrev[e];

				while (orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
					t = addTriangle(q, i, e, -1, hullTria[e], hullTria[q]);
					legalize(t + 2);
					hullTria[q] = t;
					hullNext[e] = e; // mark as removed
					hullSize--;
					e = q;
					q = hullPrev[e];
				}
			}

			// update the hull indices
			hullStart = hullPrev[i] = e;
			hullNext[e] = hullPrev[next] = i;
			hullNext[i] = next;

			// save the two new edges in the hash table
			hullHash[hashKey(x, y)] = i;
			hullHash[hashKey(coords[2 * e], coords[2 * e + 1])] = e;
		}

		hull = new int[hullSize];
		int s = hullStart;
		for (int i = 0; i < hullSize; i++) {
			hull[i] = s;
			s = hullNext[s];
		}
		
		//consider constrains!
//		for(int k=0; k<10; k++) {
//			for(int c=0; c<constrain.length; c++) {
//				JDPoint u = constrain[c].a;
//				JDPoint v = constrain[c].b;
//				legalize(u.x, u.y, v.x, v.y);
//				//System.out.println("Legalize constrain ["+coords[2*u]+", "+coords[2*u+1]+", "+coords[2*v]+", "+coords[2*v+1]+"]");
//			}
//		}

		hullPrev = hullNext = hullTria = null; // get rid of temporary arrays

		// trim typed triangle mesh arrays
		int[] tempTriangles = new int[trianglesLen];
		System.arraycopy(triangles, 0, tempTriangles, 0, trianglesLen);
		triangles = tempTriangles;

		int[] tempHalfedges = new int[trianglesLen];
		System.arraycopy(halfedges, 0, tempHalfedges, 0, trianglesLen);
		halfedges = tempHalfedges;
	}
	
	private boolean checkConstrainCrossings() {
		List<JDEdge> constrainCopy = new ArrayList<JDEdge>();
		for(JDEdge c: constrain)
			constrainCopy.add(c);
		System.out.println("Check edges:");
		for(JDEdge e: constrainCopy)
			System.out.println(e);
		boolean atLeastOneCrossing = false;
		List<JDPoint> newPoints = new ArrayList<JDPoint>();
		for(int cc=constrainCopy.size()-1; cc>=0; cc--) {
			List<JDEdge> newEdges = splitConstrain(constrainCopy.get(cc), newPoints);
			if(newEdges.size()>0) {
				System.out.println("List \"newEdges\" is not empty:");
				System.out.println("    remove edge: "+constrainCopy.get(cc));
				constrainCopy.remove(cc);
				System.out.println("    add new edges:");
				for(JDEdge e: newEdges)
					System.out.println("        "+e);
				constrainCopy.addAll(newEdges);
				atLeastOneCrossing = true;
			}
		}

		if(atLeastOneCrossing)
			merge(null, constrainCopy, ep, newPoints);

		return atLeastOneCrossing;
	}
	private List<JDEdge> splitConstrain(JDEdge ce, List<JDPoint> newpoints) {
		double ccx = ce.b.x - ce.a.x;
		double ccy = ce.b.y - ce.a.y;
		List<JDEdge> splitEdges = new ArrayList<JDEdge>();
		JDPoint csp = new JDPoint(ce.a.x,ce.a.y,ce.a.value);
		for(int a=0; a<halfedges.length; a++) {
			int b = halfedges[a];
			if(halfedges[a]<a)
				continue;

			int a0 = a - a % 3;
			int b0 = b - b % 3;
			int al = a0 + (a + 1) % 3;
			int bl = b0 + (b + 2) % 3;
			int ar = a0 + (a + 2) % 3;

			int p0 = triangles[ar];
			int pr = triangles[a];
			int pl = triangles[al];
			int p1 = triangles[bl];

			double plx = coords[2 * pl];
			double ply = coords[2 * pl + 1];
			double prx = coords[2 * pr];
			double pry = coords[2 * pr + 1];
			
			//check crossing:
			double lrx = plx - prx;
			double lry = ply - pry;
			double adx = prx - ce.a.x;
			double ady = pry - ce.a.y;
			
			if( Math.abs(prx-ce.a.x)<EPSILON && Math.abs(pry-ce.a.y)<EPSILON)
				continue;
			if( Math.abs(prx-ce.b.x)<EPSILON && Math.abs(pry-ce.b.y)<EPSILON)
				continue;
			if( Math.abs(plx-ce.a.x)<EPSILON && Math.abs(ply-ce.a.y)<EPSILON)
				continue;
			if( Math.abs(plx-ce.b.x)<EPSILON && Math.abs(ply-ce.b.y)<EPSILON)
				continue;
			
			//pr  +  u * (pl-pr)  =  ca  +  v * (cb-ca)
			//
			//(pr-ca)  +  u * (pl-pr)  =  v * (cb-ca)
			//
			// x(prca) * y(cbca)  +  u * x(plpr) * y(cbca)
			//-y(prca) * x(cbca)  -  u * y(plpr) * x(cbca)  =  0
			double u = ( adx*ccy - ady*ccx ) / ( ccx*lry - ccy*lrx ); //crossing on halfedge
			double v = ( adx*lry - ady*lrx ) / ( ccx*lry - ccy*lrx ); //crossing on constrain
			
			if(0d<u && u<1d && 0d<v && v<1d) {
				System.out.println("Found intersection on constrain:  e["+ce.a+" - "+ce.b+"]");
				double mx = ce.a.x  +  v * ccx;
				double my = ce.a.y  +  v * ccy;
				double mv = Double.isNaN(ce.a.value) ? ce.b.value : Double.isNaN(ce.b.value) ? ce.a.value : 0.5d*(ce.a.value+ce.b.value);
				System.out.println("   new point inserted:  p["+mx+", "+my+", "+mv+"]");
				JDPoint newPoint = new JDPoint(mx,my,mv);
				newpoints.add(newPoint);
				splitEdges.add(new JDEdge(new JDPoint(csp.x,csp.y,csp.value), newPoint));
				System.out.println("   new edge on stack:  "+splitEdges.get(splitEdges.size()-1));
				csp = new JDPoint(newPoint.x, newPoint.y, newPoint.value);
			}
		}
		if(splitEdges.size()>0) {
			splitEdges.add(new JDEdge(csp, new JDPoint(ce.b.x,ce.b.y,ce.b.value)));
			System.out.println("   new edge on stack:  "+splitEdges.get(splitEdges.size()-1));
		}
		return splitEdges;
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

	public JDPoint getTriangleCircumcenter(int t) {
		JDPoint[] vertices = getTrianglePoints(t);
		return getCircumcenter(vertices[0], vertices[1], vertices[2]);
	}

	public JDPoint getCentroid(int t) {
		JDPoint[] vertices = getTrianglePoints(t);
		return getCentroid(vertices);
	}

	public JDPoint getCircumcenter(JDPoint a, JDPoint b, JDPoint c) {
		return circumcenter(a.x, a.y, b.x, b.y, c.x, c.y);
	}

	public JDPoint getCentroid(JDPoint[] points) {
		double accumulatedArea = 0.0f;
		double centerX = 0.0f;
		double centerY = 0.0f;
		double centerV = 0.0f;

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

		accumulatedArea *= 3f;
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

	///////////////////////////////////////////////////////////////////////////

	private void generateConvexHull() {
		int n = points.length;

		//get left most point:
		double xref = Double.MAX_VALUE;
		int lmIDX = -9999;
		for(int c=0; c<points.length; c++) {
			if ( points[c].x < xref ) {
				lmIDX = c;
				xref = points[c].x;
			}
		}

		int currIDX = lmIDX, newIDX = -9999;
		double[] dir = { 0, -1 };
		int hullcount = 0;
		int[] hullids = new int[n];
		while(newIDX != lmIDX) {
			double minDist = (maxX-minX)*(maxX-minX) + (maxY-minY)*(maxY-minY);
			double maxDirEqual = -minDist;
			JDPoint current = points[currIDX];
			for(int c=0; c<points.length; c++) {
				if ( c==currIDX ) continue;
				JDPoint cc = points[c];
				double dx = cc.x-current.x, dy = cc.y - current.y;
				double distToNewPoint = Math.sqrt(dx*dx+dy*dy);
				double[] dirToNewPoint = { dx/distToNewPoint, dy/distToNewPoint };
				double dt = dir[0]*dirToNewPoint[0] + dir[1]*dirToNewPoint[1];
				if ( dt > maxDirEqual - 0.00000001d ) {
					boolean newFound = true;
					if ( Math.abs(dt-maxDirEqual) < 0.000000011d ) newFound = ( distToNewPoint < minDist );
					if ( newFound ) {
						minDist = distToNewPoint;
						maxDirEqual = dt;
						newIDX = c;
						hullids[hullcount] = newIDX;
						hullHash[hullcount] = hashKey(cc.x,cc.y);
						hullcount++;
					}
				}
			}
		}
		hullids = null;
	}

	private void addConstrains() {
		
	}

	private void orderEdgesViaHashtable() {
		
	}

	private void makeDelaunayOfNonConstrainedEdges() {
		
	}
}
