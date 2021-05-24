package pplots.transform;

import pplots.PAxis;
import pplots.PPlotMath;
import pplots.shapes.PGroupShape;
import pplots.shapes.PLineShape;
import pplots.shapes.PPlotShape;

public class PPlateCarreeProjection implements PProjection {
	
	private double full_circle, half_circle, quarter_circle;
	
	public PPlateCarreeProjection(boolean in_degree) {
		full_circle    = in_degree ? 360d : 2d*Math.PI;
		half_circle    = in_degree ? 180d : Math.PI;
		quarter_circle = in_degree ? 90d : 0.5d*Math.PI;
	}

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
		double cv = Math.cos(input_in_degree ? v*PPlotMath.DEG_TO_RAD : v);
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
	public double[] defaultMapExtend() {
		return new double[] {-180d,180d,-90d,90d};
	}
	
	@Override
	public void drawBorder(PAxis ax, PGroupShape s) {
		PPlotShape.stroke(0xff000000); PPlotShape.strokeWeight(3f);
		int[] p = ax.getSize();
		s.addChild(new PLineShape(p[0],     p[1],     p[0]+p[2],p[1]));
		s.addChild(new PLineShape(p[0],     p[1]+p[3],p[0]+p[2],p[1]+p[3]));
		s.addChild(new PLineShape(p[0],     p[1],     p[0],     p[1]+p[3]));
		s.addChild(new PLineShape(p[0]+p[2],p[1],     p[0]+p[2],p[1]+p[3]));
	}
}
