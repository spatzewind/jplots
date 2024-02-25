package jplots.helper;

import java.util.ArrayList;
import java.util.List;

import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;

public class GeometryTools {
	
	/**
	 * calculates the signed area from a polygon defined by the series of x,y coordinates
	 * @param x array of all x-coordinates
	 * @param y array of all y-coordinates
	 * @return negative area when corners are counterclockwise in right-handed coordinate system
	 */
	public static float area(float[] x, float[] y) {
		if ((x == null) || (y == null) || (x.length < 2))
			return Float.NaN;
		if (x.length == 2)
			return 0f;
		float a = 0f;
		for (int s = x.length - 1, e = 0; e < x.length; s = e++)
			a += (x[e] - x[s]) * (y[e] + y[s]);
		return 0.5f * a;
	}
	/**
	 * calculates the signed area from a polygon defined by the series of x,y coordinates
	 * @param x array of all x-coordinates
	 * @param y array of all y-coordinates
	 * @return negative area when corners are counterclockwise in right-handed coordinate system
	 */
	public static double area(double[] x, double[] y) {
		if ((x == null) || (y == null) || (x.length < 2))
			return Double.NaN;
		if (x.length == 2)
			return 0d;
		double a = 0d;
		for (int s = x.length - 1, e = 0; e < x.length; s = e++)
			a += (x[e] - x[s]) * (y[e] + y[s]);
		return 0.5d * a;
	}
	/**
	 * calculates the signed area of a polygon
	 * 
	 * @param points sequence of corner points of the polygon
	 * @return negative area when corners are counterclockwise in right-handed coordinate system
	 */
	public static double area(JDPoint... points) {
		if ((points == null) || (points.length < 2))
			return Double.NaN;
		if (points.length == 2)
			return 0d;
		double a = 0d;
		for (int s = points.length - 1, e = 0; e < points.length; s = e++)
			a += (points[e].x - points[s].x) * (points[e].y + points[s].y);
		return 0.5d * a;
	}
	/**
	 * calculates the signed area of a polygon
	 * 
	 * @param  points   array of corner points of the polygon
	 * @param  edges    array of pairs of indices of points, which defines edges of the polygon
	 * @return negative area when corners are counterclockwise in right-handed coordinate system
	 */
	public static double area(JDPoint[] points, int[][] edges) {
		if (points==null || edges==null)
			return Double.NaN;
		if (edges.length == 2)
			return 0d;
		double a = 0d;
		for(int[] e: edges)
			a += (points[e[1]].x - points[e[0]].x) * (points[e[1]].y + points[e[0]].y);
		return 0.5d * a;
	}
	
	
	public static double atan(JDPoint p1, JDPoint p2) {
		return atan(p1, p2, new double[] { 1d, 0d });
	}
	public static double atan(JDPoint p1, JDPoint p2, double[] null_dir) {
		double u = p2.x - p1.x, v = p2.y - p1.y;
		double x = u * null_dir[0] + v * null_dir[1];
		double y = u * null_dir[1] - v * null_dir[0];
		double r = Math.sqrt(x * x + y * y) + 1.0e-10d;
		double c = Math.acos(x / r);
		if (y < 0d)
			c = 2 * Math.PI - c;
		return c;
	}
	
	public static boolean segmentContainsPoint(JDPoint es, JDPoint ee, JDPoint p, double tol) {
		double dx = ee.x-es.x, dy = ee.y-es.y;
		double dr = Math.sqrt(dx*dx+dy*dy), dt = tol*dr;
//		double or = Math.abs((p.x-es.x)*dy-(p.y-es.y)*dx);
//		double ta = (p.x-es.x)*dx+(p.y-es.y)*dy;
//		System.out.println("SCP:   ta/or="+ta+"/"+or+"  with  dr/dt="+dr+"/"+dt);
		if(Math.abs((p.x-es.x)*dy-(p.y-es.y)*dx)>=dt) return false;
		double dl = (p.x-es.x)*dx+(p.y-es.y)*dy;
		if(dl<=-dt) return false;
		if(dl>=dr*dr+dt) return false;
		return true;
	}
	public static JDPoint segmentIntersection(JDPoint as, JDPoint ae, JDPoint bs, JDPoint be, double tol) {
		// as + u*(ae-as)  =  bs + v*(be-bs)
		// 
		// as.x + u*(ae.x-as.x) = bs.x + v*(be.x-bs.x)
		// as.y + u*(ae.y-as.y) = bs.y + v*(be.y-bs.y)
		// 
		// u*(ae.x-as.x) + v*(bs.x-be.x) = bs.x-as.x
		// u*(ae.y-as.y) + v*(bs.y-be.y) = bs.y-as.y
		double den = (as.x-ae.x)*(bs.y-be.y) - (as.y-ae.y)*(bs.x-be.x);
		double u = ((as.x-bs.x)*(bs.y-be.y) - (as.y-bs.y)*(bs.x-be.x)) / den;
		if( u < -tol || u > 1d+tol ) return null;
		double v = ((as.x-bs.x)*(as.y-ae.y) - (as.y-bs.y)*(as.x-ae.x)) / den;
		v = Math.max(0d, Math.min(1d, v));
		return bs.fractionTowards(v, be);
	}
	public static JDPoint intersects(JDPoint as, JDPoint ae, JDPoint bs, JDPoint be, double tolerance) {
		if(as.equals(bs, tolerance)) return null;
		if(ae.equals(be, tolerance)) return null;
		
		if(GeometryTools.segmentContainsPoint(as, ae, bs, tolerance))
			return bs;
		if(GeometryTools.segmentContainsPoint(bs, be, as, tolerance))
			return as;
		if(GeometryTools.segmentContainsPoint(as, ae, be, tolerance))
			return null;
		if(GeometryTools.segmentContainsPoint(bs, be, ae, tolerance))
			return null;
	    // As + u*(Ae-As) = Bs + v*(Be-Bs)
	    // u*(Ae-As) + v*(Bs-Be) = Bs - As
	    double abx = bs.x-as.x, aby = bs.y-as.y;
	    double aax = ae.x-as.x, aay = ae.y-as.y;
	    double bbx = be.x-bs.x, bby = be.y-bs.y;
	    //println("A: "+aax+"/"+aay+"\nB: "+bbx+"/"+bby+"\nAB: "+abx+"/"+aby);
	    // u*aax - v*bbx = abx
	    // u*aay - v*bby = aby
	    double u = (abx*bby - aby*bbx) / (aax*bby - aay*bbx);
	    double v = (abx*aay - aby*aax) / (aax*bby - aay*bbx);
	    //println("u = "+u+"\nv = "+v);
	    if( u<0d || u>1d ) return null;
	    if( v<0d || v>1d ) return null;
	    return new JDPoint( 0.5d*(as.x+u*aax + bs.x+v*bbx), 0.5d*(as.y+u*aay + bs.y+v*bby) );
	}

	public static JDPolygon union(JDPolygon poly1, JDPolygon poly2) {
		return union(poly1, poly2, Double.NaN);
	}
	public static JDPolygon union(JDPolygon poly1, JDPolygon poly2, double tolerance) {
	    if(poly1==null) return poly2;
	    if(poly2==null) return poly1;
	    
	    double[] pb = poly1.getBounds(), tb = poly2.getBounds();
	    if(pb[1]<tb[0] || pb[0]>tb[1] || pb[3]<tb[2] || pb[2]>tb[3]) return null;
	    
	    if(Double.isNaN(tolerance)) {
	    	double ix = Math.min(pb[0],tb[0]), iy = Math.min(pb[2],tb[2]);
	    	double ax = Math.max(pb[1],tb[1]), ay = Math.max(pb[3],tb[3]);
	    	tolerance = Math.max(ax-ix,ay-iy) * 1e-10d;
	    }
	    
	    List<Intersection> intersections = new ArrayList<>();
	    for(int i=0; i<poly1.edges.length; i++) {
	    	int[] ee = poly1.edges[i];
	    	for(int j=0; j<poly2.edges.length; j++) {
	    		int[] ff = poly2.edges[j];
	            JDPoint X = intersects(poly1.c[ee[0]],poly1.c[ee[1]], poly2.c[ff[0]],poly2.c[ff[1]], tolerance);
	            if(X!=null)
	                intersections.add(new Intersection(ee[0],ee[1],ff[0],ff[1], X, i,j));
	        }
	    }
	    
	    if(intersections.isEmpty()) {
	        if(poly2.contains(poly1.c[0], tolerance)) {
	            return poly2;
	        }
	        if(poly1.contains(poly2.c[0], tolerance)) {
	            return poly1;
	        }
	        return null;
	    }
	    
        boolean[] bp = new boolean[poly1.c.length];
        boolean[] bq = new boolean[poly2.c.length];
        for(int i=Math.max(poly1.c.length,poly2.c.length)-1; i>=0; i--) {
            if(i<poly1.c.length)
                bp[i] = !poly2.contains(poly1.c[i], tolerance);
            if(i<poly2.c.length)
                bq[i] = !poly1.contains(poly2.c[i], tolerance);
        }
        
        List<JDPoint> new_points = new ArrayList<>();
        
        //find left most point.
        double leftmost = Double.POSITIVE_INFINITY;
        int id=-1; int currPoly = 0;
        for(int i=0; i<poly1.edges.length; i++) {
            JDPoint pp = poly1.c[poly1.edges[i][0]];
            if(pp.x<leftmost) { leftmost = pp.x;
                id = i; currPoly = 1; }
        }
        for(int i=0; i<poly2.edges.length; i++) {
            JDPoint qq = poly2.c[poly2.edges[i][0]];
            if(qq.x<leftmost) { leftmost = qq.x;
                id = i; currPoly = 2; }
        }
        
        JDPoint lastpoint = poly1.c[0];
        double vx=0d, vy=0d;
        if(currPoly==1) { lastpoint = poly1.c[poly1.edges[id][0]]; JDPoint n = poly1.c[poly1.edges[id][1]];
            vx = n.x-lastpoint.x; vy = n.y-lastpoint.y; bp[poly1.edges[id][0]] = false; }
        if(currPoly==2) { lastpoint = poly2.c[poly2.edges[id][0]]; JDPoint n = poly2.c[poly2.edges[id][1]];
            vx = n.x-lastpoint.x; vy = n.y-lastpoint.y; bq[poly2.edges[id][0]] = false; }
//        println((currPoly==1?"P1["+poly1.edges[id][0]:"P2["+poly2.edges[id][0])+"]: "+lastpoint);
        JDPoint startPoint = lastpoint.copy();
        JDPolygon res = null;
        int hi = 0;
        while(lastpoint!=null) {
        	//trace boundary
	        while(true) {
	            new_points.add(lastpoint.copy());
	            int isct=-1, isct_cnt=0;
	            double dd = Double.POSITIVE_INFINITY;
	            int lpi = -1;
	            if(currPoly==1) lpi = poly1.edges[id][0];
	            if(currPoly==2) lpi = poly2.edges[id][0];
	            for(int k=0; k<intersections.size(); k++) {
	                Intersection ii = intersections.get(k);
	                boolean find = (ii.p1s==lpi && currPoly==1) || (ii.p2s==lpi && currPoly==2);
	                if(find) {
	                    JDPoint X = ii.p;
	                    double d = (X.x-lastpoint.x)*vx + (X.y-lastpoint.y)*vy;
	                    if(d>=0d) {
	                        if(d<dd) {
	                            dd = d;
	                            isct = k;
	                        }
	                        isct_cnt++;
	                    }
	                }
	            }
	            if(isct_cnt==0) {
	                if(currPoly==1) {
	                	lastpoint = poly1.c[poly1.edges[id][1]];
	                    bp[poly1.edges[id][1]] = false;
	                	id = poly1.edges[id][2];
	                	JDPoint n = poly1.c[poly1.edges[id][1]];
	                    vx = n.x-lastpoint.x; vy = n.y-lastpoint.y;
	                }
	                if(currPoly==2) {
	                	lastpoint = poly2.c[poly2.edges[id][1]];
	                    bq[poly2.edges[id][1]] = false;
	                	id = poly2.edges[id][2];
	                	JDPoint n = poly2.c[poly2.edges[id][1]];
	                    vx = n.x-lastpoint.x; vy = n.y-lastpoint.y;
	                }
//	                println((currPoly==1?"P1["+poly1.edges[id][0]:"P2["+poly2.edges[id][0])+"]: "+lastpoint);
	            } else {
	            	Intersection ii = intersections.get(isct);
	                lastpoint = ii.p;
	                intersections.remove(isct);
	                int nextPoly = 0;
	                if(currPoly==1) { nextPoly = 2; id = ii.n2; JDPoint n = poly2.c[poly2.edges[id][1]];
	                    vx = n.x-lastpoint.x; vy = n.y-lastpoint.y; }
	                if(currPoly==2) { nextPoly = 1; id = ii.n1; JDPoint n = poly1.c[poly1.edges[id][1]];
	                    vx = n.x-lastpoint.x; vy = n.y-lastpoint.y; }
	                currPoly = nextPoly;
//	                println("isect: "+lastpoint+"  "+ii.p1s+"/"+ii.p1e+"/"+ii.p2s+"/"+ii.p2e+"  (-> "+(currPoly==1?"P1":"P2")+")");
	            }
	            //TODO rem break;
//	            if(new_points.size()>100) break;
	            if(lastpoint.equals(startPoint)) break;
	        }
	        if(res==null) {
//	        	System.out.println("Polygon created...\n\n");
	        	res = new JDPolygon(new_points.toArray(new JDPoint[0]));
//	        	fill(0xffff9900); noStroke();
//	            for(Intersection i: intersections)
//	            	circle((float)i.p.x, height-(float)i.p.y, diam);
	        } else {
//	        	System.out.println("Hole added to polygon...");
	        	res.addHole(new_points.toArray(new JDPoint[0]));
//	        	//TODO rem break;
//	        	break;
	        }
	        lastpoint = null;
	        id = -1;
	        currPoly = -1;
	        new_points.clear();
	        for(int i=0; i<poly1.edges.length && lastpoint==null; i++)
	        	if(bp[poly1.edges[i][0]]) {
	        		lastpoint = poly1.c[poly1.edges[i][0]]; JDPoint n = poly1.c[poly1.edges[i][1]];
	        		vx = n.x-lastpoint.x; vy = n.y-lastpoint.y;
	        		id = i; currPoly = 1; }
		    for(int i=0; i<poly2.edges.length && lastpoint==null; i++)
		    	if(bq[poly2.edges[i][0]]) {
			        lastpoint = poly2.c[poly2.edges[i][0]]; JDPoint n = poly2.c[poly2.edges[i][1]];
	        		vx = n.x-lastpoint.x; vy = n.y-lastpoint.y;
			        id = i; currPoly = 2; }
		    if(lastpoint==null && !intersections.isEmpty()) {
		    	Intersection ii = intersections.get(0);
		    	vx = poly1.c[ii.p1e].x-poly1.c[ii.p1s].x;
		    	vy = poly1.c[ii.p1e].y-poly1.c[ii.p1s].y;
		    	double os = vx*(poly2.c[ii.p2s].y-poly1.c[ii.p1s].y) - vy*(poly2.c[ii.p2s].x-poly1.c[ii.p1s].x);
		    	double oe = vx*(poly2.c[ii.p2e].y-poly1.c[ii.p1s].y) - vy*(poly2.c[ii.p2e].x-poly1.c[ii.p1s].x);
		    	boolean n2_i2o = oe>os;
		    	vx = poly2.c[ii.p2e].x-poly2.c[ii.p2s].x;
		    	vy = poly2.c[ii.p2e].y-poly2.c[ii.p2s].y;
		    	os = vx*(poly1.c[ii.p1s].y-poly2.c[ii.p2s].y) - vy*(poly1.c[ii.p1s].x-poly2.c[ii.p2s].x);
		    	oe = vx*(poly1.c[ii.p1e].y-poly2.c[ii.p2s].y) - vy*(poly1.c[ii.p1e].x-poly2.c[ii.p2s].x);
		    	boolean n1_i2o = oe>os;
		    	if(n1_i2o && !n2_i2o) {
		    		lastpoint = ii.p.copy(); JDPoint n = poly1.c[ii.p1e];
		    		vx = n.x-lastpoint.x; vy = n.y-lastpoint.y;
		    		id = ii.n1; currPoly = 1;
		    	} else
		    	if(n2_i2o && !n1_i2o) {
		    		lastpoint = ii.p.copy(); JDPoint n = poly2.c[ii.p2e];
		    		vx = n.x-lastpoint.x; vy = n.y-lastpoint.y;
		    		id = ii.n2; currPoly = 2;
		    	}
		    }
		    if(lastpoint!=null) {
		    	startPoint = lastpoint.copy();
//		    	println((currPoly==1?"P1["+poly1.edges[id][0]:"P2["+poly2.edges[id][0])+"]: "+lastpoint);
		    }
        }
        
        return res;
	}
	
	
	
	
	
	public static void checkCut(List<JDPoint> p, List<Double> c, double[] normal, double crit, double eps) {
		if (p.size() < 2)
			return;
		int k = p.size() - 2;
		int l = p.size() - 1;
		double uk = p.get(k).x * normal[0] + p.get(k).y * normal[1];// (xy=='x' ? p.get(k).x : p.get(k).y) * sign;
		double ul = p.get(l).x * normal[0] + p.get(l).y * normal[1];// (xy=='x' ? p.get(l).x : p.get(l).y) * sign;
		if (Math.abs(uk - crit) < eps && Math.abs(ul - crit) < eps) {
			double vk = p.get(k).x * normal[1] - p.get(k).y * normal[0];
			double vl = p.get(l).x * normal[1] - p.get(l).y * normal[0];
			if (vl < vk)
				c.add(-1d - l);
			else
				c.add(1d + l);
		}
	}

	public static void checkCutC(List<JDPoint> p, List<double[]> c, JDPoint center, double radius,
			double[] minus_x_dir) {
		if (p.size() < 2)
			return;
		int k = p.size() - 2;
		int l = p.size() - 1;
		double uk = center.dist2(p.get(k));
		double ul = center.dist2(p.get(l));
		double crit = radius * radius;
		double eps = crit * 1.0e-10d;
		if (Math.abs(uk - crit) < eps && Math.abs(ul - crit) < eps) {
			double vk = 0d;
			if (minus_x_dir[0] == 0d && minus_x_dir[1] == 0d) {
				minus_x_dir[0] = p.get(k).x - center.x;
				minus_x_dir[1] = p.get(k).y - center.y;
				double mxdR = 1d / Math.sqrt(minus_x_dir[0] * minus_x_dir[0] + minus_x_dir[1] * minus_x_dir[1]);
				minus_x_dir[0] *= mxdR;
				minus_x_dir[1] *= mxdR;
			} else {
				vk = atan(center, p.get(k), minus_x_dir);
			}
			double vl = atan(center, p.get(l), minus_x_dir);
			boolean isNegative = vl < vk;
			if (radius < 0d)
				isNegative = !isNegative;
			c.add(new double[] { isNegative ? -1d - l : 1d + l, vk, vl });
		}
	}

	private static void insert_edgepoints(List<JDPoint> sub_polygon, double[] normal, double crit, double eps,
			double step) {
		double delta = step * Math.sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
		for (int s = sub_polygon.size() - 1, e = 0; s >= 0; e = s--) {
			JDPoint ps = sub_polygon.get(s), pe = sub_polygon.get(e);
			double us = ps.x * normal[0] + ps.y * normal[1], ue = pe.x * normal[0] + pe.y * normal[1];
			if (Math.abs(us - crit) < eps && Math.abs(ue - crit) < eps) {
//				int idx = 0;
				double vs = ps.x * normal[1] - ps.y * normal[0];
				double ve = pe.x * normal[1] - pe.y * normal[0];
				int cnt_insertes = (int) (Math.abs(ve-vs)/delta+0.0001d);
				if(cnt_insertes>0) {
					for(int i=0; i<cnt_insertes; i++)
						sub_polygon.add(s + 1, pe.fractionTowards((i+1d)/(cnt_insertes+1d), ps));
				}
//				while (idx >= 0) {
//					idx = -1;
//					double maxf = 0d;
//					double vc = Double.NaN;
//					for (int i = inserts.size() - 1; i >= 0; i--) {
//						JDPoint pi = inserts.get(i);
//						double vi = pi.x * normal[1] - pi.y * normal[0];
//						double f = (vi - vs) / (ve - vs);
//						if (maxf < f && f < 1d) {
//							maxf = f;
//							idx = i;
//							vc = vi;
//						}
//					}
//					if (idx >= 0) {
////						System.out.println("[GEOTOOLS] insert point "+inserts.get(idx)+" with s="+s+"/"+sub_polygon.size()+
////								" v="+vs+"|"+ve+" f="+maxf);
//						sub_polygon.add(s + 1, inserts.get(idx));
//						ve = vc;// +(vs<ve?-1.0e-8d:-1.0e-8d);
//					}
//				}
			}
		}
	}

	/*
	private static JDPoint circle_line_intersection(JDPoint start_point, JDPoint end_point, JDPoint circle_center,
			double circle_radius) {
		double dx = end_point.x - start_point.x;
		double dy = end_point.y - start_point.x;
		double dr = dx * dx + dy * dy;
		double di = (start_point.x - circle_center.x) * (end_point.y - circle_center.y)
				- (end_point.x - circle_center.x) * (start_point.y - circle_center.y);
		double w = Math.sqrt(circle_radius * circle_radius * dr - di * di);
		double u = (dy < 0d ? -1d : 1d) * dx * w;
		double v = Math.abs(dy) * w;
		double rx = (di * dy + u) / dr + circle_center.x, ry = (-di * dx + v) / dr + circle_center.y;
		double f = (rx - start_point.x) * dx + (ry - start_point.y) * dy;
		if (f < 0d || f > 1d) {
			rx = (di * dy - u) / dr + circle_center.x;
			ry = (-di * dx - v) / dr + circle_center.y;
			f = (rx - start_point.x) * dx + (ry - start_point.y) * dy;
		}
		f /= dr;
		double rv = Double.isNaN(start_point.value) ? end_point.value
				: Double.isNaN(end_point.value) ? start_point.value
						: start_point.value + f * (end_point.value - start_point.value);
		return new JDPoint(rx, ry, rv);
	}
	*/

	private static void insert_circlepoints(List<JDPoint> sub_polygon, JDPoint center, double radius,
			double angle_resolution, double[] x_dir) {
		double ar = Math.max(0.1d, Math.min(15d, angle_resolution)) * Math.PI / 180d;
		double crit = radius * radius;
		double eps = crit * 1.0e-10d;
		for (int s = sub_polygon.size() - 1, e = 0; s >= 0; e = s--) {
			JDPoint ps = sub_polygon.get(s), pe = sub_polygon.get(e);
			double us = ps.dist2(center), ue = pe.dist2(center);
			if (Math.abs(us - crit) < eps && Math.abs(ue - crit) < eps) {
				double vs = atan(center, ps, x_dir);
				double ve = atan(center, ps, x_dir);
				int idxS = (int) (Math.ceil(vs / ar) + 0.5d);
				int idxE = (int) (Math.floor(ve / ar) + 0.5d);
				if (vs > ve) {
					idxS = (int) (Math.ceil(ve / ar) + 0.5d);
					idxE = (int) (Math.floor(vs / ar) + 0.5d);
				}
				for (int i = idxS; i <= idxE; i++) {
					double u = Math.cos(ar * i);
					double v = -Math.sin(ar * i);
					double x = Math.abs(radius) * (u * x_dir[0] + v * x_dir[1]);
					double y = Math.abs(radius) * (u * x_dir[1] - v * x_dir[0]);
					sub_polygon.add(s + 1, new JDPoint(center.x + x, center.y + y));
				}
			}
		}
	}
	
	/*
	private static void add_if_negative(List<JDPoint> points, List<Double> ids, double[] norm, double value,
			double epsilon) {
		JDPoint p = points.get(points.size() - 2);
		JDPoint q = points.get(points.size() - 1);
		boolean p_on_line = Math.abs(p.x * norm[0] + p.y * norm[1] - value) < epsilon;
		boolean q_on_line = Math.abs(q.x * norm[0] + q.y * norm[1] - value) < epsilon;
		if (p_on_line && q_on_line) {
			double u = p.x * norm[1] - p.y * norm[0];
			double v = q.x * norm[1] - q.y * norm[0];
			if (v < u)
				ids.add(0d + points.size());
		}
	}
	*/

	public static List<JDPoint[]> SutherlandHodgmanAlgorithm(JDPoint[] points, double[] normal, double crit,
			double eps) {
		return SutherlandHodgmanAlgorithm(points, normal, crit, eps, -1d);
	}

	public static List<JDPoint[]> SutherlandHodgmanAlgorithm(JDPoint[] points, double[] normal, double crit, double eps,
			double fillstep) {
		List<JDPoint[]> res = new ArrayList<>();
		if (points == null) return res;
		List<JDPoint> new_points = new ArrayList<>();
		List<Double> cutids = new ArrayList<>();
		for (int s = points.length - 1, e = 0; e < points.length; s = e++) {
			double vs = points[s].x * normal[0] + points[s].y * normal[1];
			double ve = points[e].x * normal[0] + points[e].y * normal[1];
			if (ve - crit < eps) {
				if (vs - crit >= eps) {
					double vf = (crit - vs) / (ve - vs);
					new_points.add(points[s].fractionTowards(vf, points[e]));
					if (new_points.size() > 1)
						checkCut(new_points, cutids, normal, crit, eps);
				}
				new_points.add(points[e]);
			} else if (vs - crit < eps) {
				double vf = (crit - vs) / (ve - vs);
				new_points.add(points[s].fractionTowards(vf, points[e]));
			}
		}
		if (new_points.size() > 1) {
			new_points.add(new_points.get(0));
			checkCut(new_points, cutids, normal, crit, eps);
			new_points.remove(new_points.size() - 1);
		}

		int pn = new_points.size();
		for (int c = cutids.size() - 1; c >= 0; c--) {
			Double ci = cutids.get(c);
			if (ci < 0d) continue;
			int e = (int) (ci - 0.5d) % pn;
			int s = (e + pn - 1) % pn;
			double vs = new_points.get(s).x * normal[1] - new_points.get(s).y * normal[0];
			double ve = new_points.get(e).x * normal[1] - new_points.get(e).y * normal[0];
			boolean overlaps = false;
			for (Double cj : cutids) {
				if (cj.equals(ci))
					continue;
				int n = (int) (Math.abs(cj) - 0.5d) % pn;
				int m = (n + pn - 1) % pn;
				double vm = new_points.get(m).x * normal[1] - new_points.get(m).y * normal[0];
				double fm = (vm - vs) / (ve - vs);
				if (0d < fm && fm < 1d) {
					overlaps = true;
					break;
				}
				double vn = new_points.get(n).x * normal[1] - new_points.get(n).y * normal[0];
				double fn = (vn - vs) / (ve - vs);
				if (0d < fn && fn < 1d) {
					overlaps = true;
					break;
				}
			}
			if (!overlaps)
				cutids.remove(c);
		}

		if (cutids.size() == 0) {
			if (fillstep > 0d) insert_edgepoints(new_points, normal, crit, eps, fillstep);
			res.add(new_points.toArray(new JDPoint[0]));
			return res;
		}
		Double ci0 = cutids.get(0);
		Double ciS = ci0;
		List<JDPoint> subpoly = new ArrayList<>();
		for (int c = 1; c < cutids.size(); c++) {
			subpoly.clear();
			Double ci1 = cutids.get(c);
			int s = (int) (Math.abs(ci0) - 0.5d);
			int e = (int) (Math.abs(ci1) - 0.5d);
			int l = (e - s + 2 * pn) % pn;
			for (int n = 0; n < l; n++)
				subpoly.add(new_points.get((n + s) % pn));
			if (fillstep > 0d)
				insert_edgepoints(subpoly, normal, crit, eps, fillstep);
			res.add(subpoly.toArray(new JDPoint[0]));
			ci0 = ci1;
		}
		subpoly.clear();
		int s = (int) (Math.abs(ci0) - 0.5d);
		int e = (int) (Math.abs(ciS) - 0.5d);
		int l = (e - s + 2 * pn) % pn;
		for (int n = 0; n < l; n++)
			subpoly.add(new_points.get((n + s) % pn));
		if (fillstep > 0d)
			insert_edgepoints(subpoly, normal, crit, eps, fillstep);
		res.add(subpoly.toArray(new JDPoint[0]));

		return res;
	}

//	public static List<JDPoint[]> WeilerAthertonAlgorithm(JDPoint[] points, JDPoint[] clipper) {
//		List<JDPoint[]> res = new ArrayList<>();
//		if(points==null) return res;
//		if(clipper==null) return res;
//		boolean[] used = new boolean[points.length];
//		for(int i=0; i<points.length; i++) used[i] = false;
//		//try for each vertex of subject polygon as start vertex a new polygon creation
//		for(int s=0; s<points.length; s++) {
//			boolean isBeginning = true;
//			if()
//		}
//
//		return res;
//	}
	/*
	private static double isInside(JDPoint test_point, JDPoint[] polygon) {
		double r = -1d;
		for (int i = polygon.length - 1, j = 0; j < polygon.length; i = j++) {
			boolean checkY = (polygon[i].y > test_point.y) != (polygon[j].y > test_point.y);
			double f = (polygon[j].x - polygon[i].x) / (polygon[j].y - polygon[i].y);
			boolean checkX = test_point.x < polygon[i].x + (test_point.y - polygon[i].y) * f;
			if (checkY && checkX)
				r = -r;
		}
		return r;
	}
	*/

	public static List<JDPoint[]> SutherlandHodgmanAlgorithmC(JDPoint[] points, double circle_radius,
			JDPoint circle_center) {
		return SutherlandHodgmanAlgorithmC(points, circle_radius, circle_center, 0.5d);
	}

	public static List<JDPoint[]> SutherlandHodgmanAlgorithmC(JDPoint[] points, double radius, JDPoint center,
			double angle_resolution_degree) {
		List<JDPoint[]> res = new ArrayList<>();
		if (points == null)
			return res;
		if (Math.abs(radius) < 1.e-8d) {
			System.err.println(
					"[JPLOTS] SUT.HOG.-ALGORITHM (CIRCLE): radius of cutting circle has to be larger than 1.0e-8");
			res.add(points);
			return res;
		}
//		double sig = radius < 0d ? -1d : 1d;
		double crit = radius * Math.abs(radius);
//		double eps = Math.abs(crit) * 1.0e-10d;
		List<JDPoint> new_points = new ArrayList<>();
		List<double[]> cutids = new ArrayList<>();
		double[] null_dir = { 0d, 0d };
		for (int s = points.length - 1, e = 0; e < points.length; s = e++) {
			double vs = points[s].dist2(center);
			double ve = points[e].dist2(center);
			if (ve <= crit) {
				if (vs > crit) {
					new_points.add(center.circleCrossBetween(points[s], points[e], radius));
					if (new_points.size() > 1)
						checkCutC(new_points, cutids, center, radius, null_dir);
				}
				new_points.add(points[e]);
			} else if (vs <= crit) {
				new_points.add(center.circleCrossBetween(points[s], points[e], radius));
			}
		}
		if (new_points.size() > 1) {
			new_points.add(new_points.get(0));
			checkCutC(new_points, cutids, center, radius, null_dir);
			new_points.remove(new_points.size() - 1);
		}

		int pn = new_points.size();
		for (int c = cutids.size() - 1; c >= 0; c--) {
			double[] ci = cutids.get(c);
			if (ci[0] < 0d)
				continue;
			double vs = ci[1], ve = ci[2];
			boolean overlaps = false;
			for (double[] cj : cutids) {
				if (Math.abs(cj[0] - ci[0]) < 1.0e-10d)
					continue;
				double fm = (cj[1] - vs) / (ve - vs);
				if (0d < fm && fm < 1d) {
					overlaps = true;
					break;
				}
				double fn = (cj[2] - vs) / (ve - vs);
				if (0d < fn && fn < 1d) {
					overlaps = true;
					break;
				}
			}
			if (!overlaps)
				cutids.remove(c);
		}

		if (cutids.size() == 0) {
			insert_circlepoints(new_points, center, radius, angle_resolution_degree, null_dir);
			res.add(new_points.toArray(new JDPoint[0]));
			return res;
		}
		double[] ci0 = cutids.get(0);
		double[] ciS = ci0;
		List<JDPoint> subpoly = new ArrayList<>();
		for (int c = 1; c < cutids.size(); c++) {
			subpoly.clear();
			double[] ci1 = cutids.get(c);
			int s = (int) (Math.abs(ci0[0]) - 0.5d);
			int e = (int) (Math.abs(ci1[0]) - 0.5d);
			int l = (e - s + 2 * pn) % pn;
			for (int n = 0; n < l; n++)
				subpoly.add(new_points.get((n + s) % pn));
			insert_circlepoints(subpoly, center, radius, angle_resolution_degree, null_dir);
			res.add(subpoly.toArray(new JDPoint[0]));
			ci0 = ci1;
		}
		subpoly.clear();
		int s = (int) (Math.abs(ci0[0]) - 0.5d);
		int e = (int) (Math.abs(ciS[0]) - 0.5d);
		int l = (e - s + 2 * pn) % pn;
		for (int n = 0; n < l; n++)
			subpoly.add(new_points.get((n + s) % pn));
		insert_circlepoints(subpoly, center, radius, angle_resolution_degree, null_dir);
		res.add(subpoly.toArray(new JDPoint[0]));

		return res;
	}
	
	
	
	

	private static class Intersection {
	    public int p1s,p1e, p2s,p2e;
	    public JDPoint p;
	    public int n1, n2;
	    
	    public Intersection(int p1i, int p1j, int p2i, int p2j, JDPoint isct, int n1, int n2) {
	        this.p1s = p1i; this.p1e = p1j;
	        this.p2s = p2i; this.p2e = p2j;
	        this.p   = isct.copy();
	        this.n1  = n1;
	        this.n2  = n2;
	    }
	    
	    @Override
	    public String toString() {
	    	return "Isect@"+String.format("%08X",hashCode())+"["+p1s+"/"+p1e+"|"+p2s+"/"+p2e+", "+p+" -> "+n1+"/"+n2+"]";
	    }
	}
}
