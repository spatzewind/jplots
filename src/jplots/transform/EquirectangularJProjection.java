package jplots.transform;

import jplots.JAxis;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;

public class EquirectangularJProjection implements JProjection {
	
	private double cenlat_d,cenlon_d;
	private double cenlat_r,cenlon_r;
	private double radaequ,radpol;
	
	public EquirectangularJProjection(double center_longitude, double center_latitude, boolean in_degree) {
		double fac = in_degree ? JPlotMath.DEG_TO_RAD : 1d;
		cenlat_r = center_latitude * fac;
		cenlon_r = center_longitude * fac;
		fac = in_degree ? 1d : JPlotMath.RAD_TO_DEG;
		cenlat_d = center_latitude * fac;
		cenlon_d = center_longitude * fac;
		radaequ = EARTH_RADIUS_AEQU*Math.PI;
		radpol  = EARTH_RADIUS_POL*Math.PI;
	}

	@Override
	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree) {
		if(!(Double.isFinite(x) && Double.isFinite(y)))
			return new double[] {Double.NaN, Double.NaN};
		double u,v;
		if(output_in_degree) {
			u = cenlon_d + 180d*x/radaequ;
			if(u>180d) u-=360d; if(u<-180d) u+=360d;
			v = cenlat_d + 180d*x/radpol;
			if(v>90d) { v = Double.NaN; u = Double.NaN; }
			if(v<-90d) { v = Double.NaN; u = Double.NaN; }
		} else {
			u = cenlon_r + Math.PI*x/radaequ;
			if(u>Math.PI) u-=2*Math.PI; if(u<-Math.PI) u+=2*Math.PI;
			v = cenlat_r + Math.PI*x/radpol;
			if(v>0.5d*Math.PI) { v = Double.NaN; u = Double.NaN; }
			if(v<-0.5d*Math.PI) { v = Double.NaN; u = Double.NaN; }
		}
		return new double[] {u,v};
	}

	@Override
	public double[] fromLATLONtoPROJ(double u, double v, boolean input_in_degree) {
		if(!(Double.isFinite(u) && Double.isFinite(v)))
			return new double[] {Double.NaN, Double.NaN};
		double x,y;
		if(input_in_degree) {
			x = radaequ * (u-cenlon_d) / 180d;
			y = radpol  * (v-cenlat_d) / 180d;
		} else {
			x = radaequ * (u-cenlon_r) / Math.PI;
			y = radpol  * (v-cenlat_r) / Math.PI;
		}
		if(x<-radaequ) x += 2*radaequ;
		if(x>radaequ)  x -= 2*radaequ;
		return new double[] {x,y};
	}
	
	@Override
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree) {
		double cv = Math.cos(input_in_degree ? v*JPlotMath.DEG_TO_RAD :v);
		double h = 1d;
		double k = 1d / cv;
		return new double[] {k,0d, 0d,h};
	}
	
	@Override
	public double[] tissotFromProj(double x, double y) {
		double[] xy = fromPROJtoLATLON(x, y, false);
		return tissotFromLatLon(xy[0], xy[1], false);
	}
	
	@Override
	public double[] defaultMapExtend() {
		return new double[] {-radaequ,radaequ,-0.5d*radpol,0.5d*radpol};
	}
	
	@Override
	public void drawBorder(JAxis ax, JGroupShape s) {
		JPlotShape.stroke(0xff000000); JPlotShape.strokeWeight(3f);
		int[] p = ax.getSize();
		s.addChild(new JLineShape(p[0],     p[1],     p[0]+p[2],p[1]));
		s.addChild(new JLineShape(p[0],     p[1]+p[3],p[0]+p[2],p[1]+p[3]));
		s.addChild(new JLineShape(p[0],     p[1],     p[0],     p[1]+p[3]));
		s.addChild(new JLineShape(p[0]+p[2],p[1],     p[0]+p[2],p[1]+p[3]));
	}
}
