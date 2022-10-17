package jplots.transform;

import java.util.ArrayList;
import java.util.List;

import jplots.JAxis;
import jplots.maths.AffineBuilder;
import jplots.maths.JDLine;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;

public class OrthographicJProjection implements JProjection {

	private static List<JDPoint> left_inserts;
	// private static List<JDPoint> right_inserts;
	static {
		left_inserts = new ArrayList<>();
		double f = Math.PI / 360d;
		for (int a = -359; a <= 359; a++) {
			double x = Math.cos(f * a) + 1d;
			double y = Math.sin(f * a);
			double r = 2d / (x * x + y * y);
			left_inserts.add(new JDPoint(1d, r * y));
		}
	}

	private double cenlat, cenlon, updir;
	private double[][] rotmat_p2x, rotmat_x2p;

	public OrthographicJProjection(double center_longitude, double center_latitude, double upward_direction,
			boolean in_degree) {
		double fac = in_degree ? JPlotMath.DEG_TO_RAD : 1d;
		cenlat = center_latitude * fac;
		cenlon = center_longitude * fac;
		updir = upward_direction * fac;
		calcMatrices();
	}

	@Override
	public void setCentralLatitude(double latitude, boolean in_degree) {
		double fac = in_degree ? JPlotMath.DEG_TO_RAD : 1d;
		cenlat = latitude * fac;
		calcMatrices();
	}

	@Override
	public void setCentralLongitude(double longitude, boolean in_degree) {
		double fac = in_degree ? JPlotMath.DEG_TO_RAD : 1d;
		cenlon = longitude * fac;
		calcMatrices();
	}

	private void calcMatrices() {
		double ca = Math.cos(cenlat), sa = Math.sin(cenlat);
		double co = Math.cos(cenlon), so = Math.sin(cenlon);
		double cu = Math.cos(updir), su = Math.sin(updir);
		// Ru * Ra * Ro
		// / 1 0 0 \ /+ca 0 +sa\ /+co +so 0 \
		// | 0 +cu +su| * | 0 1 0 | * |-so +co 0 |
		// \ 0 -su +cu/ \-sa 0 +ca/ \ 0 0 1 /
		//
		// / 1 0 0 \ /+caco +caso +sa \
		// | 0 +cu +su| * | -so +co 0 |
		// \ 0 -su +cu/ \-saco -saso +ca /
		//
		// / +caco +caso +sa \
		// |-cuso-susaco +cuco-susaso +suca |
		// \+suso-cusaco -suco-cusaso +cuca /
		rotmat_p2x = new double[][] { { ca * co, ca * so, sa },
				{ -cu * so - su * sa * co, cu * co - su * sa * so, su * ca },
				{ +su * so - cu * sa * co, -su * co - cu * sa * so, cu * ca } };
		// (Ru * Ra * Ro)'
		// Ro' * Ra' * Ro'
		// /+co -so 0 \ /+ca 0 -sa\ / 1 0 0 \
		// |+so +co 0 | * | 0 1 0 | * | 0 +cu -su|
		// \ 0 0 1 / \+sa 0 +ca/ \ 0 +su +cu/
		//
		// /+co -so 0 \ / +ca -sasu -sacu \
		// |+so +co 0 | * | 0 +cu -su |
		// \ 0 0 1 / \ +sa +casu +cacu /
		//
		// / +coca -cosasu-socu -cosacu+sosu\
		// | +soca -sosasu+cocu -sosacu-cosu|
		// \ +sa +casu +cacu /
		rotmat_x2p = new double[][] { { ca * co, -so * cu - co * sa * su, so * su - co * sa * cu },
				{ so * ca, co * cu - so * sa * su, -co * su - so * sa * cu }, { sa, ca * su, ca * cu } };
	}

	public double[][] getRotmatP2X() {
		return rotmat_p2x;
	}

	@Override
	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree, boolean cut_outside) {
		if (!(Double.isFinite(x) && Double.isFinite(y)))
			return new double[] { Double.NaN, Double.NaN };
		double fac = output_in_degree ? JPlotMath.RAD_TO_DEG : 1d;
		double j = x / EARTH_RADIUS_MEAN;
		double k = y / EARTH_RADIUS_MEAN;
		double i = j * j + k * k;
		if (i > 1) {
//			j /= i;
//			k /= i;
			i = cut_outside ? Double.NaN : -Math.sqrt(Math.max(0d, 1d - j * j - k * k));
		} else {
			i = Math.sqrt(Math.max(0d, 1d - i));
		}
		double i2 = rotmat_x2p[0][0] * i + rotmat_x2p[0][1] * j + rotmat_x2p[0][2] * k;
		double j2 = rotmat_x2p[1][0] * i + rotmat_x2p[1][1] * j + rotmat_x2p[1][2] * k;
		double k2 = rotmat_x2p[2][0] * i + rotmat_x2p[2][1] * j + rotmat_x2p[2][2] * k;
		return new double[] { fac * Math.acos(i2 / (Math.sqrt(i2 * i2 + j2 * j2) + 1e-24d)) * (j2 < 0d ? -1d : 1d),
				fac * Math.asin(Math.max(-1, Math.min(1, k2))) };
	}

	@Override
	public double[] fromLATLONtoPROJ(double u, double v, boolean input_in_degree, boolean cut_outside) {
		if (!(Double.isFinite(u) && Double.isFinite(v)))
			return new double[] { Double.NaN, Double.NaN };
		double fac = input_in_degree ? JPlotMath.DEG_TO_RAD : 1d;
		double i = Math.cos(fac * v) * Math.cos(fac * u);
		double j = Math.cos(fac * v) * Math.sin(fac * u);
		double k = Math.sin(fac * v);

		double x = rotmat_p2x[1][0] * i + rotmat_p2x[1][1] * j + rotmat_p2x[1][2] * k;
		double y = rotmat_p2x[2][0] * i + rotmat_p2x[2][1] * j + rotmat_p2x[2][2] * k;
		double z = rotmat_p2x[0][0] * i + rotmat_p2x[0][1] * j + rotmat_p2x[0][2] * k;
		double r = EARTH_RADIUS_MEAN;
		if (z < 0d) {
			r = EARTH_RADIUS_MEAN / (x * x + y * y);
			if (r > 1.0e20d || cut_outside)
				r = Double.NaN;
		}
		return new double[] { r * x, r * y };
	}

	@Override
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree) {
		double[] xy = fromLATLONtoPROJ(u, v, input_in_degree, false);
		return tissotFromProj(xy[0], xy[1]);
	}

	@Override
	public double[] tissotFromProj(double x, double y) {
		double i = x / EARTH_RADIUS_MEAN;
		double j = y / EARTH_RADIUS_MEAN;
		double w2 = i * i + j * j;
		if (w2 > 1)
			return new double[] { Double.NaN, Double.NaN, Double.NaN, Double.NaN };
		double smin = Math.sqrt(1d - w2);
		double smax = 1d;
		return new double[] { -smax * j, smax * i, smin * i, smin * j };
	}

	@Override
	public List<JDLine> splitByMapBorder(JDLine line) {
		return line.intersectsCircle(new JDPoint(0,0), EARTH_RADIUS_MEAN);
	}

	@Override
	public List<JDPolygon> splitByMapBorder(JDPolygon poly) {
//		List<JDPolygon> res = new ArrayList<>();
//		JDPoint[] p = new JDPoint[poly.c.length];
//		for (int i = 0; i < poly.c.length; i++) {
//			double x = poly.c[i].x / EARTH_RADIUS_MEAN + 1d;
//			double y = poly.c[i].y / EARTH_RADIUS_MEAN;
//			double r = 2d / (x * x + y * y);
////			if(r>1.0e18d) {
////				x = x<0d ? -5.0e-9d : 5.0e-9d;
////				y = y<0d ? -5.0e-9d : 5.0e-9d;
////				r = 2d / ( x*x + y*y );
////			}
//			p[i] = new JDPoint(r * x, -r * y);
//		}
//		List<JDPoint[]> polys = GeometryTools.SutherlandHodgmanAlgorithm(p, new double[] { -1d, 0d }, -1d, 1.0e-10d,
//				left_inserts);
//		for (JDPoint[] pnts : polys) {
//			if (pnts.length < 3)
//				continue;
//			for (int i = 0; i < pnts.length; i++) {
//				double x = pnts[i].x, y = -pnts[i].y;
//				double r = 2d / (x * x + y * y);
//				x = r * x - 1d;
//				y = r * y;
//				pnts[i].x = x * EARTH_RADIUS_MEAN;
//				pnts[i].y = y * EARTH_RADIUS_MEAN;
//			}
//			res.add(new JDPolygon(pnts));
//		}
////		List<JDPoint[]> polygons = GeometryTools.SutherlandHodgmanAlgorithmC(poly.c, EARTH_RADIUS_MEAN, new JDPoint(0d,0d));
////		for(JDPoint[] pnts: polygons)
////			res.add(new JDPolygon(pnts));
////		System.out.println("[ORTHOGRAPHIC-JPROJ.] return "+res.size()+" polygon(s)");
//		return res;
		return poly.intersectsCircle(new JDPoint(0d,0d), EARTH_RADIUS_MEAN);
	}

	@Override
	public double[] defaultMapExtend() {
		return new double[] { -EARTH_RADIUS_MEAN, EARTH_RADIUS_MEAN, -EARTH_RADIUS_MEAN, EARTH_RADIUS_MEAN };
	}

	@Override
	public void drawBorder(JAxis ax, JGroupShape s) {
		int[] p = ax.getSize();
		double[] r = ax.getRange();
		float x1 = (float) (p[0]
				+ p[2] * Math.max(0d, Math.min(1d, JPlotMath.map(-EARTH_RADIUS_MEAN, r[0], r[1], 0d, 1d))));
		float y1 = (float) (p[1] + p[3] * Math.max(0d, Math.min(1d, JPlotMath.map(0d, r[2], r[3], 0d, 1d))));
		for (int i = -359; i <= 360; i++) {
			double bx = EARTH_RADIUS_MEAN * Math.cos(0.5d * i * JPlotMath.DEG_TO_RAD);
			double by = EARTH_RADIUS_MEAN * Math.sin(0.5d * i * JPlotMath.DEG_TO_RAD);
			float x2 = (float) (p[0] + p[2] * Math.max(0d, Math.min(1d, JPlotMath.map(bx, r[0], r[1], 0d, 1d))));
			float y2 = (float) (p[1] + p[3] * Math.max(0d, Math.min(1d, JPlotMath.map(by, r[2], r[3], 1d, 0d))));
			s.addChild(new JLineShape(3f, 0xff000000, x1, y1, x2, y2));
			x1 = x2;
			y1 = y2;
		}
	}

	@Override
	public void addGrid(JAxis ax, JGroupShape s) {
		// determin range of used longitude and latitude (3rd option: shifted longitude
		// range)
		double minLn3 = 360d, maxLn3 = 0d;
		double minLon = 180d, maxLon = -180d;
		double minLat = 90d, maxLat = -90d;
		double pp1 = 0, pp2 = 0, pp3 = 0, pp4 = 0;
		double[] extent = ax.getRange();
		
		for (int k = 0; k <= 50; k++) {
			double kx = extent[0] + 0.02d * k * (extent[1] - extent[0]);
			double ky = extent[2] + 0.02d * k * (extent[3] - extent[2]);
			double[] xy1 = fromPROJtoLATLON(kx, extent[2], true, false);
			double[] xy2 = fromPROJtoLATLON(kx, extent[3], true, false);
			double[] xy3 = fromPROJtoLATLON(extent[0], ky, true, false);
			double[] xy4 = fromPROJtoLATLON(extent[1], ky, true, false);
			double latI = Math.min(Math.min(xy1[1], xy2[1]), Math.min(xy3[1], xy4[1]));
			double latA = Math.max(Math.max(xy1[1], xy2[1]), Math.max(xy3[1], xy4[1]));
			if (latI < minLat) minLat = latI;
			if (latA > maxLat) maxLat = latA;
			double lonI = Math.min(Math.min(xy1[0], xy2[0]), Math.min(xy3[0], xy4[0]));
			double lonA = Math.max(Math.max(xy1[0], xy2[0]), Math.max(xy3[0], xy4[0]));
			if (lonI < minLon) minLon = lonI;
			if (lonA > maxLon) maxLon = lonA;
			pp1 = xy1[0];
			pp2 = xy2[0];
			pp3 = xy3[0];
			pp4 = xy4[0];
			if (xy1[0] < 0d)
				pp1 += 360d;
			if (xy2[0] < 0d)
				pp2 += 360d;
			if (xy3[0] < 0d)
				pp3 += 360d;
			if (xy4[0] < 0d)
				pp4 += 360d;
			lonI = Math.min(Math.min(pp1, pp2), Math.min(pp3, pp4));
			lonA = Math.max(Math.max(pp1, pp2), Math.max(pp3, pp4));
			if (lonI < minLn3)
				minLn3 = lonI;
			if (lonA > maxLn3)
				maxLn3 = lonA;
			pp1 = xy1[0];
			pp2 = xy2[0];
			pp3 = xy3[0];
			pp4 = xy4[0];
		}
		if (ax.getPlot().isDebug())
			System.out.println("[DEBUG]   first estimate in Orthog.Proj.:  lon={" + minLon + " ... " + maxLon + "}  lat={" + minLat + " ... " + maxLat + "}");
		// longitude range spanes more than 180 degrees, one of the poles is inside the
		// view
		double[] ce = fromPROJtoLATLON(0.5d * (extent[0] + extent[1]), 0.5d * (extent[2] + extent[3]), true, false);
		if (ce[0] < -90d || ce[0] > 90d) {
			minLon = minLn3;
			maxLon = maxLn3;
		}
		if (maxLon - minLon > 180.1d) {
			double[] uv = fromLATLONtoPROJ(0d, 90d, true, false);
			if (extent[0] <= uv[0] && uv[0] <= extent[1] && extent[2] <= uv[1] && uv[1] <= extent[3]) {
				maxLat = 90d;
				minLon = -180d;
				maxLon = 180d;
			}
			uv = fromLATLONtoPROJ(0d, -90d, true, false);
			if (extent[0] <= uv[0] && uv[0] <= extent[1] && extent[2] <= uv[1] && uv[1] <= extent[3]) {
				minLat = -90d;
				minLon = -180d;
				maxLon = 180d;
			}
		}
		// actually plot grid, if it is needed
		if (ax.getPlot().isDebug())
			System.out.println("[DEBUG]   second estimate in Orthog.Proj.: lon={" + minLon + " ... " + maxLon + "}  lat={" + minLat + " ... " + maxLat + "}");
		double deltaO = 0.2d*(maxLon-minLon); // 66.5%tile, to estimate needed spacing
 		double spcO = 180d;
		for (int k = 0; k < SPACING.length && spcO > 179.9d; k++)
			if (deltaO >= SPACING[k])
				spcO = SPACING[k];
		if (spcO > 179.9d)
			spcO = SPACING[SPACING.length - 1];
		int minO = (int) (minLon / spcO - 0.5d);
		int maxO = (int) (maxLon / spcO + 0.5d);
		double deltaA = (maxLat - minLat) / 3d;
		double spcA = 180d;
		for (int k = 0; k < SPACING.length && spcA > 179.9d; k++)
			if (deltaA >= SPACING[k])
				spcA = SPACING[k];
		if (spcA > 179.9d)
			spcA = SPACING[SPACING.length - 1];
		int minA = (int) (minLat / spcA - 0.5d);
		int maxA = (int) (maxLat / spcA + 0.5d);
		if (ax.getPlot().isDebug()) {
			System.out.println(
					"[DEBUG]   found spacing {lon:" + spcO + "/" + deltaO + ", lat:" + spcA + "/" + deltaA + "}");
			System.out.println(
					"[DEBUG]   grid of Orthog.Proj.:  id={" + "lon: " + minO + "..." + maxO + ", lat: "+ minA + "..." + maxA + "}");
		}
		List<JDLine> gridlines = new ArrayList<>();
		if(ax.isXGridVisible()) { //longitudes
			JDPoint[] p = new JDPoint[91];
			for (int o = minO; o <= maxO; o++) {
				double glon = o * spcO;
				double miA = minLat, maA = maxLat;
				if(minLat<-89.9999d && o%3!=0) miA += spcA;
				if(maxLat>89.9999d && o%3!=0)  maA -= spcA;
				for (int a = 0; a <= 90; a++) {
					double glat = miA + a * (maA - miA) / 90d;
					double[] uv1 = fromLATLONtoPROJ(glon, glat, true, false);
					p[a] = new JDPoint(uv1[0], uv1[1]);
				}
				gridlines.addAll(splitByMapBorder(new JDLine(p.clone())));
			}
		}
		if (ax.isYGridVisible()) { //latitudes
			JDPoint[] p = new JDPoint[37];
			for (int a = minA; a <= maxA; a++) {
				double glat = a * spcA;
				for (int o = 0; o <= 36; o++) {
					double glon = minLon + o * (maxLon - minLon) / 36d;
					double[] uv1 = fromLATLONtoPROJ(glon, glat, true, false);
					p[o] = new JDPoint(uv1[0], uv1[1]);
				}
				gridlines.addAll(splitByMapBorder(new JDLine(p.clone())));
			}
		}
		int[] p = ax.getSize();
		double xs = p[2] / (extent[1]-extent[0]), ys = p[3] / (extent[3]-extent[2]);
		AffineBuilder affine = new AffineBuilder().scale(1d,-1d).translate(-extent[0], -extent[2]).scale(xs, ys).translate(p[0], p[1]);
		double[][] affmat = affine.getMatrix();
		for(JDLine gl: gridlines)
			for(JDLine l: gl.affine(affmat).intersectsAABB(p[0], p[0]+p[2], p[1], p[1]+p[3]))
				s.addChild(new JLineShape(2f, 0xff999999, l.getCoordsAsFloats()));
	}
}
