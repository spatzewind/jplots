package jplots.transform;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jplots.JAxis;
import jplots.maths.JDLine;
import jplots.maths.JDPolygon;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;

public class StereographicJProjection implements JProjection {

	private double cenlat, cenlon, updir, earth_scale;
	private double[][] rotmat_p2x, rotmat_x2p;

	public StereographicJProjection(double center_longitude, double center_latitude, double upward_direction,
			boolean in_degree) {
		double fac = in_degree ? JPlotMath.DEG_TO_RAD : 1d;
		cenlat = center_latitude * fac;
		cenlon = center_longitude * fac;
		updir = upward_direction * fac;
		earth_scale = EARTH_RADIUS_MEAN * Math.PI;
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

	@Override
	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree, boolean cut_outside) {
		if (!(Double.isFinite(x) && Double.isFinite(y)))
			return new double[] { Double.NaN, Double.NaN };
		double fac = output_in_degree ? JPlotMath.RAD_TO_DEG : 1d;
		double us = x / earth_scale;
		double vs = y / earth_scale;
		double ws2 = us * us + vs * vs;
		double w = (1 - ws2) / (1 + ws2);
		double u = us * (1 + w);
		double v = vs * (1 + w);
		double i = rotmat_x2p[0][0] * w + rotmat_x2p[0][1] * u + rotmat_x2p[0][2] * v;
		double j = rotmat_x2p[1][0] * w + rotmat_x2p[1][1] * u + rotmat_x2p[1][2] * v;
		double k = rotmat_x2p[2][0] * w + rotmat_x2p[2][1] * u + rotmat_x2p[2][2] * v;
		return new double[] { fac * Math.acos(i / (Math.sqrt(i * i + j * j) + 1e-24d)) * (j < 0d ? -1d : 1d),
				fac * Math.asin(Math.max(-1, Math.min(1, k))) };
	}

	@Override
	public double[] fromLATLONtoPROJ(double longitude, double latitude, boolean input_in_degree, boolean cut_outside) {
		if (!(Double.isFinite(longitude) && Double.isFinite(latitude)))
			return new double[] { Double.NaN, Double.NaN };
		double fac = input_in_degree ? JPlotMath.DEG_TO_RAD : 1d;
		double i = Math.cos(fac * latitude) * Math.cos(fac * longitude);
		double j = Math.cos(fac * latitude) * Math.sin(fac * longitude);
		double k = Math.sin(fac * latitude);
		double w = 1d + rotmat_p2x[0][0] * i + rotmat_p2x[0][1] * j + rotmat_p2x[0][2] * k;
		if (w < 0.01d)
			return new double[] { Double.NaN, Double.NaN };
		w = earth_scale / w;
		return new double[] { w * (rotmat_p2x[1][0] * i + rotmat_p2x[1][1] * j + rotmat_p2x[1][2] * k),
				w * (rotmat_p2x[2][0] * i + rotmat_p2x[2][1] * j + rotmat_p2x[2][2] * k) };
	}

	@Override
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree) {
		double[] xy = fromLATLONtoPROJ(u, v, input_in_degree, false);
		return tissotFromProj(xy[0], xy[1]);
	}

	@Override
	public double[] tissotFromProj(double x, double y) {
		double i = x / earth_scale;
		double j = y / earth_scale;
		// r = 2*R*tan(pi/4-phi/2)
		// dr/dphi = 2*R*(1+tan²(pi/4-phi/2))*(-1/2)
		double f = 0.5d * (1d + i * i + j * j);
		return new double[] { f, 0d, 0d, f };
	}

	@Override
	public List<JDLine> splitByMapBorder(JDLine line) {
		List<JDLine> res = new ArrayList<>();
		res.add(line);
		return res;
	}

	@Override
	public List<JDPolygon> splitByMapBorder(JDPolygon poly) {
		List<JDPolygon> res = new ArrayList<>();
		res.add(poly);
		return res;
	}

	@Override
	public double[] defaultMapExtend() {
		return new double[] { -earth_scale, earth_scale, -earth_scale, earth_scale };
	}

	@Override
	public void drawBorder(JAxis ax, JGroupShape s) {
		int[] p = ax.getSize();
		s.addChild(new JLineShape(3f, 0xff000000, p[0], p[1], p[0] + p[2], p[1]));
		s.addChild(new JLineShape(3f, 0xff000000, p[0], p[1] + p[3], p[0] + p[2], p[1] + p[3]));
		s.addChild(new JLineShape(3f, 0xff000000, p[0], p[1], p[0], p[1] + p[3]));
		s.addChild(new JLineShape(3f, 0xff000000, p[0] + p[2], p[1], p[0] + p[2], p[1] + p[3]));
	}

	@Override
	public void addGrid(JAxis ax, JGroupShape s) {
		// determin range of used longitude and latitude (3rd option: shifted longitude
		// range)
		double minLn3 = 360d, maxLn3 = 0d;
		double minLon = 180d, maxLon = -180d;
		double minLat = 90d, maxLat = -90d;
		List<Double> lonlist = new ArrayList<>();
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
			if (latI < minLat)
				minLat = latI;
			if (latA > maxLat)
				maxLat = latA;
			double lonI = Math.min(Math.min(xy1[0], xy2[0]), Math.min(xy3[0], xy4[0]));
			double lonA = Math.max(Math.max(xy1[0], xy2[0]), Math.max(xy3[0], xy4[0]));
			if (lonI < minLon)
				minLon = lonI;
			if (lonA > maxLon)
				maxLon = lonA;
			if (k > 0) {
				pp1 = Math.abs(pp1 - xy1[0]);
				lonlist.add(Math.min(pp1, 360d - pp1));
				pp2 = Math.abs(pp2 - xy2[0]);
				lonlist.add(Math.min(pp2, 360d - pp2));
				pp3 = Math.abs(pp3 - xy3[0]);
				lonlist.add(Math.min(pp3, 360d - pp3));
				pp4 = Math.abs(pp4 - xy4[0]);
				lonlist.add(Math.min(pp4, 360d - pp4));
			}
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
			System.out.println("[DEBUG]   first estimate in Stereog.Proj.:  lon={" + minLon + " ... " + maxLon
					+ "}  lat={" + minLat + " ... " + maxLat + "}");
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
		Double[] lonarr = lonlist.toArray(new Double[0]);
		Arrays.sort(lonarr);
		if (ax.getPlot().isDebug())
			System.out.println("[DEBUG]   second estimate in Stereog.Proj.: lon={" + minLon + " ... " + maxLon
					+ "}  lat={" + minLat + " ... " + maxLat + "}");
		double deltaO = 50d * lonarr[130]; // 66.5%tile, to estimate needed spacing
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
			System.out.println("[DEBUG]   grid of Stereog.Proj.:  id={" + "lon: " + minO + "..." + maxO + ", lat: "
					+ minA + "..." + maxA + "}");
		}
		JPlotShape.noFill();
		JPlotShape.stroke(0xff999999);
		JPlotShape.strokeWeight(2f);
		if (ax.isXGridVisible()) { //longitudes
			for (int o = minO; o <= maxO; o++) {
				double glon = o * spcO;
				double miA = minLat, maA = maxLat;
				if(minLat<-89.9999d && o%3!=0) miA += spcA;
				if(maxLat>89.9999d && o%3!=0)  maA -= spcA;
				double[] uv0 = fromLATLONtoPROJ(glon, miA, true, false);
				for (int a = 1; a <= 90; a++) {
					double glat = miA + a * (maA - miA) / 90d;
					double[] uv1 = fromLATLONtoPROJ(glon, glat, true, false);
					drawLine(ax, s, uv0, uv1, extent);
					uv0[0] = uv1[0];
					uv0[1] = uv1[1];
				}
			}
		}
		if (ax.isYGridVisible()) { //latitudes
			for (int a = minA; a <= maxA; a++) {
				double glat = a * spcA;
				double[] uv0 = fromLATLONtoPROJ(minLon, glat, true, false);
				for (int o = 1; o <= 36; o++) {
					double glon = minLon + o * (maxLon - minLon) / 36d;
					double[] uv1 = fromLATLONtoPROJ(glon, glat, true, false);
					drawLine(ax, s, uv0, uv1, extent);
					uv0[0] = uv1[0];
					uv0[1] = uv1[1];
				}
			}
		}
	}
	
	private void drawLine(JAxis ax, JGroupShape s, double[] st, double[] en, double[] ext) {
		// for first result, go without cutting
		int[] p = ax.getSize();
		double[] r = ax.getRange();
		double xs = p[2] / (r[1] - r[0]), ys = p[3] / (r[3] - r[2]);
		double x1 = p[0] + xs * (st[0] - r[0]), x2 = p[0] + xs * (en[0] - r[0]);
		double y1 = p[1] + ys * (r[3] - st[1]), y2 = p[1] + ys * (r[3] - en[1]);

		boolean stIn = (p[0] <= x1 && x1 <= p[0] + p[2] && p[1] <= y1 && y1 <= p[1] + p[3]);
		boolean enIn = (p[0] <= x2 && x2 <= p[0] + p[2] && p[1] <= y2 && y2 <= p[1] + p[3]);
		// System.out.println(" draw line from "+x1+"|"+y1+" to "+x2+"|"+y2);
		if (!stIn && !enIn)
			return;
		if (stIn && enIn) {
			s.addChild(new JLineShape(JPlotShape.strokeWeight, JPlotShape.strokeColour, (float) x1, (float) y1, (float) x2, (float) y2));
			return;
		}
		if (!stIn) {
			if (x1 < p[0]) {
				y1 += (y2 - y1) * (p[0] - x1) / (x2 - x1);
				x1 = p[0];
			}
			if (x1 > p[0] + p[2]) {
				y1 += (y2 - y1) * (p[0] + p[2] - x1) / (x2 - x1);
				x1 = p[0] + p[2];
			}
			if (y1 < p[1]) {
				x1 += (x2 - x1) * (p[1] - y1) / (y2 - y1);
				y1 = p[1];
			}
			if (y1 > p[1] + p[3]) {
				x1 += (x2 - x1) * (p[1] + p[3] - y1) / (y2 - y1);
				y1 = p[1] + p[3];
			}
		}
		if (!enIn) {
			if (x2 < p[0]) {
				y2 += (y1 - y2) * (p[0] - x2) / (x1 - x2);
				x2 = p[0];
			}
			if (x2 > p[0] + p[2]) {
				y2 += (y1 - y2) * (p[0] + p[2] - x2) / (x1 - x2);
				x2 = p[0] + p[2];
			}
			if (y2 < p[1]) {
				x2 += (x1 - x2) * (p[1] - y2) / (y1 - y2);
				y2 = p[1];
			}
			if (y2 > p[1] + p[3]) {
				x2 += (x1 - x2) * (p[1] + p[3] - y2) / (y1 - y2);
				y2 = p[1] + p[3];
			}
		}
		if (Math.abs(x2 - x1) < 0.0001d && Math.abs(y2 - y1) < 0.0001d)
			return;
		s.addChild(new JLineShape(JPlotShape.strokeWeight, JPlotShape.strokeColour, (float) x1, (float) y1, (float) x2, (float) y2));
	}
}
