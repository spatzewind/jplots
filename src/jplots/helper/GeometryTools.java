package jplots.helper;

import java.util.ArrayList;
import java.util.List;

import jplots.maths.JDPoint;

public class GeometryTools {

	/**
	 * calculates the signed area of a polygon
	 * 
	 * @param points sequence of corner points of the polygon
	 * @return negative area when corners are counterclockwise in right-handed
	 *         coordinate system
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
			List<JDPoint> inserts) {
		for (int s = sub_polygon.size() - 1, e = 0; s >= 0; e = s--) {
			JDPoint ps = sub_polygon.get(s), pe = sub_polygon.get(e);
			double us = ps.x * normal[0] + ps.y * normal[1], ue = pe.x * normal[0] + pe.y * normal[1];
			if (Math.abs(us - crit) < eps && Math.abs(ue - crit) < eps) {
				int idx = 0;
				double vs = ps.x * normal[1] - ps.y * normal[0];
				double ve = pe.x * normal[1] - pe.y * normal[0];
				while (idx >= 0) {
					idx = -1;
					double maxf = 0d;
					double vc = Double.NaN;
					for (int i = inserts.size() - 1; i >= 0; i--) {
						JDPoint pi = inserts.get(i);
						double vi = pi.x * normal[1] - pi.y * normal[0];
						double f = (vi - vs) / (ve - vs);
						if (maxf < f && f < 1d) {
							maxf = f;
							idx = i;
							vc = vi;
						}
					}
					if (idx >= 0) {
//						System.out.println("[GEOTOOLS] insert point "+inserts.get(idx)+" with s="+s+"/"+sub_polygon.size()+
//								" v="+vs+"|"+ve+" f="+maxf);
						sub_polygon.add(s + 1, inserts.get(idx));
						ve = vc;// +(vs<ve?-1.0e-8d:-1.0e-8d);
					}
				}
			}
		}
	}

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

	public static List<JDPoint[]> SutherlandHodgmanAlgorithm(JDPoint[] points, double[] normal, double crit,
			double eps) {
		return SutherlandHodgmanAlgorithm(points, normal, crit, eps, null);
	}

	public static List<JDPoint[]> SutherlandHodgmanAlgorithm(JDPoint[] points, double[] normal, double crit, double eps,
			List<JDPoint> inserts) {
		List<JDPoint[]> res = new ArrayList<>();
		if (points == null)
			return res;
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
			if (ci < 0d)
				continue;
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
			if (inserts != null && inserts.size() > 0)
				insert_edgepoints(new_points, normal, crit, eps, inserts);
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
			if (inserts != null && inserts.size() > 0)
				insert_edgepoints(subpoly, normal, crit, eps, inserts);
			res.add(subpoly.toArray(new JDPoint[0]));
			ci0 = ci1;
		}
		subpoly.clear();
		int s = (int) (Math.abs(ci0) - 0.5d);
		int e = (int) (Math.abs(ciS) - 0.5d);
		int l = (e - s + 2 * pn) % pn;
		for (int n = 0; n < l; n++)
			subpoly.add(new_points.get((n + s) % pn));
		if (inserts != null && inserts.size() > 0)
			insert_edgepoints(subpoly, normal, crit, eps, inserts);
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
		double sig = radius < 0d ? -1d : 1d;
		double crit = radius * Math.abs(radius);
		double eps = Math.abs(crit) * 1.0e-10d;
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
}
