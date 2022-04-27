package jplots.transform;

import java.util.ArrayList;
import java.util.List;

import jplots.JAxis;
import jplots.maths.JDPolygon;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;

public class PlateCarreeJProjection implements JProjection {
	
	private double full_circle, half_circle, quarter_circle;
	
	public PlateCarreeJProjection(boolean in_degree) {
		full_circle    = in_degree ? 360d : 2d*Math.PI;
		half_circle    = in_degree ? 180d : Math.PI;
		quarter_circle = in_degree ? 90d : 0.5d*Math.PI;
	}
	
	@Override
	public void setCentralLatitude(double latitude, boolean in_degree) { }
	
	@Override
	public void setCentralLongitude(double longitude, boolean in_degree) { }
	
	@Override
	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree) {
		if(!(Double.isFinite(x) && Double.isFinite(y)))
			return new double[] {Double.NaN, Double.NaN};
		double xc = x;
		while(xc>half_circle) xc -= full_circle;
		while(xc<-half_circle) xc += full_circle;
		double yc = y<-quarter_circle ? Double.NaN : y>quarter_circle ? Double.NaN : y;
		double f = (output_in_degree ? 180d : Math.PI) / half_circle;
		return new double[] {xc*f,yc*f};
	}

	@Override
	public double[] fromLATLONtoPROJ(double u, double v, boolean input_in_degree) {
		if(!(Double.isFinite(u) && Double.isFinite(v)))
			return new double[] {Double.NaN, Double.NaN};
		double f = half_circle / (input_in_degree ? 180d : Math.PI);
		double x = u*f, y = v*f;
		while(x>half_circle) x -= full_circle;
		while(x<-half_circle) x += full_circle;
		y = y<-quarter_circle ? Double.NaN : y>quarter_circle ? Double.NaN : y;
		if(Double.isNaN(y)) x = Double.NaN;
		return new double[] {x,y};
	}
	
	@Override
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree) {
		double cv = Math.cos(input_in_degree ? v*JPlotMath.DEG_TO_RAD : v);
		double h = 1d;
		double k = 1d / cv;
		return new double[] {k,0d, 0d,h};
	}
	
	@Override
	public double[] tissotFromProj(double x, double y) {
		double[]xy = fromPROJtoLATLON(x, y, false);
		return tissotFromLatLon(xy[0], xy[1], false);
	}
	
	@Override
	public List<JDPolygon> splitByMapBorder(JDPolygon poly) {
		List<JDPolygon> res = new ArrayList<>();
		res.add(poly);
		return res;
	}
	
	@Override
	public double[] defaultMapExtend() {
		return new double[] {-180d,180d,-90d,90d};
	}
	
	@Override
	public void drawBorder(JAxis ax, JGroupShape s) {
		JPlotShape.stroke(0xff000000); JPlotShape.strokeWeight(3f);
		int[] p = ax.getSize();
		double[] r = ax.getRange();
		float[] e = {
			(float)(p[0]+p[2]*Math.max(0d, Math.min(1d, JPlotMath.map(-180d,r[0],r[1],0d,1d)))),
			(float)(p[0]+p[2]*Math.max(0d, Math.min(1d, JPlotMath.map( 180d,r[0],r[1],0d,1d)))),
			(float)(p[1]+p[3]*Math.max(0d, Math.min(1d, JPlotMath.map( -90d,r[2],r[3],0d,1d)))),
			(float)(p[1]+p[3]*Math.max(0d, Math.min(1d, JPlotMath.map(  90d,r[2],r[3],0d,1d))))
		};
		s.addChild(new JLineShape(e[0],e[3],e[1],e[3]));
		s.addChild(new JLineShape(e[0],e[2],e[1],e[2]));
		s.addChild(new JLineShape(e[0],e[3],e[0],e[2]));
		s.addChild(new JLineShape(e[1],e[3],e[1],e[2]));
	}
}
