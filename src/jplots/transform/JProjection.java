package jplots.transform;

import jplots.JAxis;
import jplots.shapes.JGroupShape;

public interface JProjection {

	public final static double EARTH_RADIUS_AEQU = 6378137.000d; //m
	public final static double EARTH_RADIUS_MEAN = 6371000.685d; //m
	public final static double EARTH_RADIUS_POL  = 6356752.000d; //m
	public final static double EARTH_FLATTENING  = 1d/298.257222101d;

	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree);
	public double[] fromLATLONtoPROJ(double u, double v, boolean input_in_degree);
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree);
	public double[] tissotFromProj(double x, double y);
	public double[] defaultMapExtend();
	public void drawBorder(JAxis ax, JGroupShape s);
}
