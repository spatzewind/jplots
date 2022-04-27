package jplots.transform;

import java.util.List;

import jplots.JAxis;
import jplots.maths.JDPolygon;
import jplots.shapes.JGroupShape;

public interface JProjection {

	public final static double EARTH_RADIUS_AEQU = 6378137.000d; //m
	public final static double EARTH_RADIUS_MEAN = 6371000.685d; //m
	public final static double EARTH_RADIUS_POL  = 6356752.000d; //m
	public final static double EARTH_FLATTENING  = 1d/298.257222101d;

	public void setCentralLatitude(double latitude, boolean in_degree);
	public void setCentralLongitude(double longitude, boolean in_degree);
	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree);
	public double[] fromLATLONtoPROJ(double u, double v, boolean input_in_degree);
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree);
	public double[] tissotFromProj(double x, double y);
	public List<JDPolygon> splitByMapBorder(JDPolygon poly);
	public double[] defaultMapExtend();
	public void drawBorder(JAxis ax, JGroupShape s);
}
