package jplots.helper;

import java.util.ArrayList;
import java.util.List;

import jplots.maths.JDPoint;

public class GeometryTools {

	/**
	 * calculates the signed area of a polygon
	 * @param points  sequence of corner points of the polygon
	 * @return        negative area when corners are counterclockwise in right-handed coordinate system
	 */
	public static double area(JDPoint... points) {
		if(points==null) return Double.NaN;
		if(points.length<2) return Double.NaN;
		if(points.length==2) return 0d;
		double a = 0d;
		for(int s=points.length-1,e=0; e<points.length; s=e++)
			a += (points[e].x-points[s].x) * (points[e].y+points[s].y);
		return 0.5d * a;
	}

	public static void checkCut(List<JDPoint> p, List<Double> c, double[] normal, double crit, double eps, JDPoint cen) {
		if(p.size()<2) return;
		int k = p.size()-2;
		int l = p.size()-1;
		double uk = p.get(k).x*normal[0] + p.get(k).y*normal[1];// (xy=='x' ? p.get(k).x : p.get(k).y) * sign;
		double ul = p.get(l).x*normal[0] + p.get(l).y*normal[1];// (xy=='x' ? p.get(l).x : p.get(l).y) * sign;
		if( Math.abs(uk-crit)<eps && Math.abs(ul-crit)<eps ) {
			boolean isNegative = area(cen,p.get(k),p.get(l))<0d;
			if(isNegative) c.add( -1d-l );
			else c.add( 1d+l );
		}
	}
	public static void checkCutC(List<JDPoint> p, List<Double> c, JDPoint center, double radius) {
		if(p.size()<2) return;
		int k = p.size()-2;
		int l = p.size()-1;
		double uk = center.dist(p.get(k));
		double ul = center.dist(p.get(l));
		double crit = Math.abs(radius);
		double eps = crit * 1.0e-10d;
		if( Math.abs(uk-crit)<eps && Math.abs(ul-crit)<eps ) {
			boolean isNegative = area(center,p.get(k),p.get(l))<0d;
			if(radius<0d) isNegative = !isNegative;
			if(isNegative) c.add( -1d-l );
			else c.add( 1d+l );
		}
	}
	
	private static void insert_edgepoints(List<JDPoint> sub_polygon, double[] normal, double crit, double eps, List<JDPoint> inserts) {
		for(int s=sub_polygon.size()-1,e=0; s>=0; e=s--) {
			JDPoint ps = sub_polygon.get(s),
					pe = sub_polygon.get(e);
			double  us = ps.x*normal[0] + ps.y*normal[1],
					ue = pe.x*normal[0] + pe.y*normal[1];
			if(Math.abs(us-crit)<eps && Math.abs(ue-crit)<eps) {
				int idx = 0;
				double vs = ps.x*normal[1] - ps.y*normal[0];
				double ve = pe.x*normal[1] - pe.y*normal[0];
				while(idx>=0) {
					idx = -1;
					double maxf = 0d;
					double vc = Double.NaN;
					for(int i=inserts.size()-1; i>=0; i--) {
						JDPoint pi = inserts.get(i);
						double vi = pi.x*normal[1] - pi.y*normal[0];
						double f = (vi-vs) / (ve-vs);
						if(maxf<f && f<1d) {
							maxf = f;
							idx = i;
							vc = vi;
						}
					}
					if(idx>=0) {
//						System.out.println("[GEOTOOLS] insert point "+inserts.get(idx)+" with s="+s+"/"+sub_polygon.size()+
//								" v="+vs+"|"+ve+" f="+maxf);
						sub_polygon.add(s+1, inserts.get(idx));
						ve = vc;//+(vs<ve?-1.0e-8d:-1.0e-8d);
					}
				}
			}
		}
	}
	private static JDPoint circle_line_intersection(JDPoint start_point, JDPoint end_point, JDPoint circle_center, double circle_radius) {
		double dx = end_point.x - start_point.x;
		double dy = end_point.y - start_point.x;
		double dr = dx*dx + dy*dy;
		double di =   (start_point.x-circle_center.x)*(end_point.y-circle_center.y)
					- (end_point.x-circle_center.x)*(start_point.y-circle_center.y);
		double w = Math.sqrt( circle_radius*circle_radius*dr - di*di );
		double u = (dy<0d?-1d:1d) * dx * w;
		double v = Math.abs(dy) * w;
		double  rx = (  di*dy + u ) / dr + circle_center.x,
				ry = ( -di*dx + v ) / dr + circle_center.y;
		double f = (rx-start_point.x)*dx + (ry-start_point.y)*dy;
		if( f < 0d || f > 1d) {
			rx = (  di*dy - u ) / dr + circle_center.x;
			ry = ( -di*dx - v ) / dr + circle_center.y;
			f = (rx-start_point.x)*dx + (ry-start_point.y)*dy;
		}
		f /= dr;
		double rv = Double.isNaN(start_point.value)?end_point.value:Double.isNaN(end_point.value)?start_point.value:
				start_point.value + f*(end_point.value-start_point.value);
		return new JDPoint(rx, ry, rv);
	}
	private static void insert_circlepoints(List<JDPoint> sub_polygon, JDPoint circle_center, double circle_radius, double angle_resolution) {
		boolean clockwise = circle_radius<0d;
		double ar = Math.max( 0.1d, Math.min(15d, angle_resolution));
		double crit = Math.abs(circle_radius);
		double eps = crit * 1.0e-10d;
		for(int s=sub_polygon.size()-1,e=0; s>=0; e=s--) {
			JDPoint ps = sub_polygon.get(s),
					pe = sub_polygon.get(e);
			double  us = ps.dist(circle_center),
					ue = pe.dist(circle_center);
			if(Math.abs(us-crit)<eps && Math.abs(ue-crit)<eps) {
				//TODO
//				int idx = 0;
//				double vs = ps.x*normal[1] - ps.y*normal[0];
//				double ve = pe.x*normal[1] - pe.y*normal[0];
//				while(idx>=0) {
//					idx = -1;
//					double maxf = 0d;
//					double vc = Double.NaN;
//					for(int i=inserts.size()-1; i>=0; i--) {
//						JDPoint pi = inserts.get(i);
//						double vi = pi.x*normal[1] - pi.y*normal[0];
//						double f = (vi-vs) / (ve-vs);
//						if(maxf<f && f<1d) {
//							maxf = f;
//							idx = i;
//							vc = vi;
//						}
//					}
//					if(idx>=0) {
////						System.out.println("[GEOTOOLS] insert point "+inserts.get(idx)+" with s="+s+"/"+sub_polygon.size()+
////								" v="+vs+"|"+ve+" f="+maxf);
//						sub_polygon.add(s+1, inserts.get(idx));
//						ve = vc;//+(vs<ve?-1.0e-8d:-1.0e-8d);
//					}
//				}
			}
		}
	}
	
	public static List<JDPoint[]> SutherlandHodgmanAlgorithm(JDPoint[] points, double[] normal, double crit, double eps, JDPoint reference) {
		return SutherlandHodgmanAlgorithm(points, normal, crit, eps, reference, null);
	}
	public static List<JDPoint[]> SutherlandHodgmanAlgorithm(JDPoint[] points, double[] normal, double crit, double eps, JDPoint reference, List<JDPoint> inserts) {
		List<JDPoint[]> res = new ArrayList<>();
		if(points==null) return res;
		List<JDPoint> new_points = new ArrayList<>();
		List<JDPoint> sub_polygon = new ArrayList<>();
		List<Double> cutids = new ArrayList<>();
		for(int s=points.length-1,e=0; e<points.length; s=e++) {
			double vs = points[s].x*normal[0] + points[s].y*normal[1];
			double ve = points[e].x*normal[0] + points[e].y*normal[1];
			if(ve<=crit) {
				if(vs>crit) {
					double vf = (crit-vs) / (ve-vs);
					new_points.add(points[s].fractionTowards(vf, points[e]));
					cutids.add(0d-new_points.size());
				}
				new_points.add(points[e]);
			} else
			if(vs<=crit) {
				double vf = (crit-vs) / (ve-vs);
				new_points.add(points[s].fractionTowards(vf, points[e]));
				cutids.add(0d+new_points.size());
			}
		}
//		System.out.println("[GEOTOOLS] "+cutids.size()+" cutids");
		if(cutids.size()>0) {
			boolean containsNegative = false;
			for(Double d: cutids) {
				if(d<0d) {
					containsNegative = true;
					break;
				}
			}
			if(containsNegative) {
//				boolean foundOverlap = true;
//				while(foundOverlap) {
//					foundOverlap = false;
//					for(int i=cutids.size()) 
//				}
			} else {
				if(inserts!=null && inserts.size()>0)
					insert_edgepoints(new_points, normal, crit, eps, inserts);
				res.add(new_points.toArray(new JDPoint[0]));
				return res;
			}
			int steps = cutids.size();
			int pl = new_points.size();
			while(cutids.get(0) > 0d && steps>0) {
				double cid = cutids.get(0);
				cutids.add(cid<0d ? cid-pl : cid+pl);
				cutids.remove(0);
				steps--;
			}
			double cid = cutids.get(0);
			cutids.add(cid<0d ? cid-pl : cid+pl);
			int ci1 = (int) (Math.abs(cid) - 0.5d);
			boolean lastPositive = false;
			for(int n=1; n<cutids.size(); n++) {
				cid = cutids.get( n );
				if(cid>0d && lastPositive) continue;
				int ci2 = (int) (Math.abs(cid)-0.5d);
				lastPositive = (cid>0d);
				if(ci2<ci1) ci2 += pl;
				if(ci2-pl>ci1) ci2 -= pl;
				sub_polygon.clear();
				for(int t=0; t<ci2-ci1; t++)
					sub_polygon.add(new_points.get((ci1+t)%pl));
				if(inserts!=null && inserts.size()>0)
					insert_edgepoints(sub_polygon, normal, crit, eps, inserts);
				res.add(sub_polygon.toArray(new JDPoint[0]));
				ci1 = ci2;
			}
		} else {
			res.add(new_points.toArray(new JDPoint[0]));
		}
		return res;
	}
	

	public static List<JDPoint[]> SutherlandHodgmanAlgorithmC(JDPoint[] points, double circle_radius, JDPoint circle_center) {
		return SutherlandHodgmanAlgorithmC(points, circle_radius, circle_center, 0.5d);
	}
	public static List<JDPoint[]> SutherlandHodgmanAlgorithmC(JDPoint[] points, double circle_radius, JDPoint circle_center, double angle_resolution_degree) {
		List<JDPoint[]> res = new ArrayList<>();
		if(points==null) return res;
		if(Math.abs(circle_radius)<1.e-8d) {
			System.err.println("[JPLOTS] SUT.HOG.-ALGORITHM (CIRCLE): radius of cutting circle has to be larger than 1.0e-8");
			res.add(points);
			return res;
		}
		double sig = circle_radius<0d ? -1d : 1d;
		List<JDPoint> new_points = new ArrayList<>();
		List<JDPoint> sub_polygon = new ArrayList<>();
		List<Double> cutids = new ArrayList<>();
		for(int s=points.length-1,e=0; e<points.length; s=e++) {
			double vs = points[s].dist(circle_center) * sig;
			double ve = points[e].dist(circle_center) * sig;
			if(ve<=circle_radius) {
				if(vs>circle_radius) {
					new_points.add(circle_line_intersection(points[s], points[e], circle_center, Math.abs(circle_radius)));
					checkCutC(new_points, cutids, circle_center, circle_radius);
				}
				new_points.add(points[e]);
			} else
			if(vs<=circle_radius) {
				new_points.add(circle_line_intersection(points[s], points[e], circle_center, Math.abs(circle_radius)));
//				checkCut(new_points, cutids, normal, crit, eps, reference);
			}
		}
		if(new_points.size()>0) {
			new_points.add(new_points.get(0));
//			checkCut(new_points, cutids, normal, crit, eps, reference);
			new_points.remove(new_points.size()-1);
		}
//		System.out.println("[GEOTOOLS] "+cutids.size()+" cutids");
		if(cutids.size()>0) {
			boolean containsNegative = false;
			for(Double d: cutids) {
				if(d<0d) containsNegative = true;
			}
			if(containsNegative) {
			} else {
//				if(inserts!=null && inserts.size()>0)
//					insert_edgepoints(new_points, normal, crit, eps, inserts);
				res.add(new_points.toArray(new JDPoint[0]));
				return res;
			}
			int steps = cutids.size();
			int pl = new_points.size();
			while(cutids.get(0) > 0d && steps>0) {
				double cid = cutids.get(0);
				cutids.add(cid<0d ? cid-pl : cid+pl);
				cutids.remove(0);
				steps--;
			}
			double cid = cutids.get(0);
			cutids.add(cid<0d ? cid-pl : cid+pl);
			int ci1 = (int) (Math.abs(cid) - 0.5d);
			boolean lastPositive = false;
			for(int n=1; n<cutids.size(); n++) {
				cid = cutids.get( n );
				if(cid>0d && lastPositive) continue;
				int ci2 = (int) (Math.abs(cid)-0.5d);
				lastPositive = (cid>0d);
				if(ci2<ci1) ci2 += pl;
				if(ci2-pl>ci1) ci2 -= pl;
				sub_polygon.clear();
				for(int t=0; t<ci2-ci1; t++)
					sub_polygon.add(new_points.get((ci1+t)%pl));
//				if(inserts!=null && inserts.size()>0)
//					insert_edgepoints(sub_polygon, normal, crit, eps, inserts);
				res.add(sub_polygon.toArray(new JDPoint[0]));
				ci1 = ci2;
			}
		} else {
			res.add(new_points.toArray(new JDPoint[0]));
		}
		return res;
	}
}
