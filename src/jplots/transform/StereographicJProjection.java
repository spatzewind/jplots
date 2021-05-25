package jplots.transform;

import jplots.JAxis;
import jplots.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;

public class StereographicJProjection implements JProjection {
	
	private double cenlat,cenlon, updir, earth_scale;
	private double[][] rotmat_p2x, rotmat_x2p;
	
	public StereographicJProjection(double center_longitude, double center_latitude, double upward_direction, boolean in_degree) {
		double fac = in_degree ? JPlotMath.DEG_TO_RAD : 1d;
		cenlat = center_latitude * fac;
		cenlon = center_longitude * fac;
		updir = upward_direction * fac;
		earth_scale = EARTH_RADIUS_MEAN * Math.PI;
		calcMatrices();
	}

	private void calcMatrices() {
		double ca = Math.cos(cenlat), sa = Math.sin(cenlat);
		double co = Math.cos(cenlon), so = Math.sin(cenlon);
		double cu = Math.cos(updir),  su = Math.sin(updir);
		// Ru * Ra * Ro
		//   / 1   0   0 \     /+ca  0 +sa\     /+co +so  0 \
		//   | 0  +cu +su|  *  | 0   1  0 |  *  |-so +co  0 |
		//   \ 0  -su +cu/     \-sa  0 +ca/     \ 0   0   1 /
		//
		//   / 1   0   0 \     /+caco +caso  +sa \
		//   | 0  +cu +su|  *  | -so   +co     0 |
		//   \ 0  -su +cu/     \-saco -saso  +ca /
		//
		//   /    +caco        +caso          +sa   \
		//   |-cuso-susaco +cuco-susaso     +suca   |
		//   \+suso-cusaco -suco-cusaso     +cuca   /
		rotmat_p2x = new double[][] {
			{        ca*co,             ca*so,       sa },
			{-cu*so -su*sa*co,  cu*co -su*sa*so,   su*ca},
			{+su*so -cu*sa*co, -su*co -cu*sa*so,   cu*ca}
		};
		// (Ru * Ra * Ro)'
		// Ro' * Ra' * Ro'
		//   /+co -so  0 \     /+ca  0 -sa\     / 1   0   0 \
		//   |+so +co  0 |  *  | 0   1  0 |  *  | 0  +cu -su|
		//   \ 0   0   1 /     \+sa  0 +ca/     \ 0  +su +cu/
		//
		//   /+co -so  0 \     / +ca  -sasu -sacu \
		//   |+so +co  0 |  *  |  0    +cu   -su  |
		//   \ 0   0   1 /     \ +sa  +casu +cacu /
		//
		//   /    +coca    -cosasu-socu -cosacu+sosu\
		//   |    +soca    -sosasu+cocu -sosacu-cosu|
		//   \     +sa         +casu        +cacu   /
		rotmat_x2p = new double[][] {
			{ca*co, -so*cu -co*sa*su,  so*su -co*sa*cu},
			{so*ca,  co*cu -so*sa*su, -co*su -so*sa*cu},
			{  sa,           ca*su,            ca*cu  }
		};
	}

	@Override
	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree) {
		if(!(Double.isFinite(x) && Double.isFinite(y)))
			return new double[] {Double.NaN, Double.NaN};
		double fac = output_in_degree ? JPlotMath.RAD_TO_DEG : 1d;
		double us = x/earth_scale;
		double vs = y/earth_scale;
		double ws2 = us*us+vs*vs;
		double w = (1-ws2)/(1+ws2);
		double u = us*(1+w);
		double v = vs*(1+w);
		double i = rotmat_x2p[0][0]*w + rotmat_x2p[0][1]*u + rotmat_x2p[0][2]*v;
		double j = rotmat_x2p[1][0]*w + rotmat_x2p[1][1]*u + rotmat_x2p[1][2]*v;
		double k = rotmat_x2p[2][0]*w + rotmat_x2p[2][1]*u + rotmat_x2p[2][2]*v;
		return new double[] {
				fac*Math.acos(i/(Math.sqrt(i*i+j*j)+1e-24d))*(j<0d?-1d:1d),
				fac*Math.asin(Math.max(-1,Math.min(1,k)))
		};
	}

	@Override
	public double[] fromLATLONtoPROJ(double u, double v, boolean input_in_degree) {
		if(!(Double.isFinite(u) && Double.isFinite(v)))
			return new double[] {Double.NaN, Double.NaN};
		double fac = input_in_degree ? JPlotMath.DEG_TO_RAD : 1d;
		double i = Math.cos(fac*v) * Math.cos(fac*u);
		double j = Math.cos(fac*v) * Math.sin(fac*u);
		double k = Math.sin(fac*v);
		double w = 1d + rotmat_p2x[0][0]*i+rotmat_p2x[0][1]*j+rotmat_p2x[0][2]*k;
		if(w<0.01d)
			return new double[] {Double.NaN, Double.NaN};
		w = earth_scale / w;
		return new double[] {
				w*(rotmat_p2x[1][0]*i+rotmat_p2x[1][1]*j+rotmat_p2x[1][2]*k),
				w*(rotmat_p2x[2][0]*i+rotmat_p2x[2][1]*j+rotmat_p2x[2][2]*k)
		};
	}
	
	@Override
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree) {
		double[] xy = fromLATLONtoPROJ(u, v, input_in_degree);
		return tissotFromProj(xy[0], xy[1]);
	}
	
	@Override
	public double[] tissotFromProj(double x, double y) {
		double i = x/earth_scale;
		double j = y/earth_scale;
		// r = 2*R*tan(pi/4-phi/2)
		// dr/dphi = 2*R*(1+tan²(pi/4-phi/2))*(-1/2)
		double f = 0.5d*(1d + i*i + j*j);
		return new double[] {f,0d, 0d,f};
	}
	
	@Override
	public double[] defaultMapExtend() {
		return new double[] {-earth_scale,earth_scale,-earth_scale,earth_scale};
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
