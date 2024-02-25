package jplots.transform;

import java.util.List;

import jplots.JPlotConstants;
import jplots.axes.JAxis;
import jplots.maths.JDLine;
import jplots.maths.JDPolygon;
import jplots.shapes.JGroupShape;

public interface JProjection extends JPlotConstants {

	final static double[] SPACING = { 30d, 10d, 5d, 2d, 1d, 0.5d, 1d / 3d, 1d / 6d, 1d / 12d, 1d / 30d, 1d / 60d };

	void setCentralLatitude(double latitude, boolean in_degree);

	void setCentralLongitude(double longitude, boolean in_degree);

	double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree, boolean cut_outside);

	double[] fromLATLONtoPROJ(double u, double v, boolean input_in_degree, boolean cut_outside);

	double[] tissotFromLatLon(double u, double v, boolean input_in_degree);

	double[] tissotFromProj(double x, double y);

	List<JDLine> splitByMapBorder(JDLine line);

	List<JDPolygon> splitByMapBorder(JDPolygon poly);

	double[] defaultMapExtend();

	void drawBorder(JAxis ax, JGroupShape s);

	void addGrid(JAxis ax, JGroupShape s);
}
