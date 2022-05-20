package jplots.transform;

import java.util.ArrayList;
import java.util.List;

import jplots.JAxis;
import jplots.maths.JDPolygon;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;

public class RectangleJProjection implements JProjection {

	private double xs, xe, ys, ye;

	public RectangleJProjection(double xstart, double xend, double ystart, double yend) {
		xs = xstart;
		xe = xend;
		ys = ystart;
		ye = yend;
	}

	@Override
	public void setCentralLatitude(double latitude, boolean in_degree) {
	}

	@Override
	public void setCentralLongitude(double longitude, boolean in_degree) {
	}

	@Override
	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree) {
		if (Double.isNaN(x) || Double.isNaN(y) || x < xs || x > xe)
			return new double[] { Double.NaN, Double.NaN };
		if (y < ys || y > ye)
			return new double[] { Double.NaN, Double.NaN };
		return new double[] { (x - xs) / (xe - xs), (y - ys) / (ye - ys) };
	}

	@Override
	public double[] fromLATLONtoPROJ(double u, double v, boolean input_in_degree) {
		if (Double.isNaN(u) || Double.isNaN(v) || u < 0d || u > 1d)
			return new double[] { Double.NaN, Double.NaN };
		if (v < 0d || v > 1d)
			return new double[] { Double.NaN, Double.NaN };
		return new double[] { xs + u * (xe - xs), ys + v * (ye - ys) };
	}

	@Override
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree) {
		return new double[] { (xe - xs) / (ye - ys), 0d, 0d, (ye - ys) / (xe - xs) };
	}

	@Override
	public double[] tissotFromProj(double x, double y) {
		return new double[] { (xe - xs) / (ye - ys), 0d, 0d, (ye - ys) / (xe - xs) };
	}

	@Override
	public List<JDPolygon> splitByMapBorder(JDPolygon poly) {
		List<JDPolygon> res = new ArrayList<>();
		res.add(poly);
		return res;
	}

	@Override
	public double[] defaultMapExtend() {
		return new double[] { xs, xe, ys, ye };
	}

	@Override
	public void drawBorder(JAxis ax, JGroupShape s) {
		JPlotShape.stroke(0xff000000);
		JPlotShape.strokeWeight(3f);
		int[] p = ax.getSize();
		s.addChild(new JLineShape(p[0], p[1], p[0] + p[2], p[1]));
		s.addChild(new JLineShape(p[0], p[1] + p[3], p[0] + p[2], p[1] + p[3]));
		s.addChild(new JLineShape(p[0], p[1], p[0], p[1] + p[3]));
		s.addChild(new JLineShape(p[0] + p[2], p[1], p[0] + p[2], p[1] + p[3]));
	}

	@Override
	public void addGrid(JAxis ax, JGroupShape s) {
		// TODO implement GRID
	}
}
