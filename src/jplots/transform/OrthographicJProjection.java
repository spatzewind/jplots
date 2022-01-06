package jplots.transform;

import jplots.JAxis;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;

public class OrthographicJProjection implements JProjection {
	
	private double cenlat,cenlon, updir;
	private double[][] rotmat_p2x, rotmat_x2p;
	
	public OrthographicJProjection(double center_longitude, double center_latitude, double upward_direction, boolean in_degree) {
		double fac = in_degree ? JPlotMath.DEG_TO_RAD : 1d;
		cenlat = center_latitude * fac;
		cenlon = center_longitude * fac;
		updir = upward_direction * fac;
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
	
	public double[][] getRotmatP2X() {
		return rotmat_p2x;
	}

	@Override
	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree) {
		if(!(Double.isFinite(x) && Double.isFinite(y)))
			return new double[] {Double.NaN, Double.NaN};
		double fac = output_in_degree ? JPlotMath.RAD_TO_DEG : 1d;
		double j = x/EARTH_RADIUS_MEAN;
		double k = y/EARTH_RADIUS_MEAN;
		double i = j*j+k*k;
		if(i>1)
			return new double[] {Double.NaN, Double.NaN};
		i = Math.sqrt(Math.max(0, 1-i));
		double i2 = rotmat_x2p[0][0]*i + rotmat_x2p[0][1]*j + rotmat_x2p[0][2]*k;
		double j2 = rotmat_x2p[1][0]*i + rotmat_x2p[1][1]*j + rotmat_x2p[1][2]*k;
		double k2 = rotmat_x2p[2][0]*i + rotmat_x2p[2][1]*j + rotmat_x2p[2][2]*k;
		return new double[] {
				fac*Math.acos(i2/(Math.sqrt(i2*i2+j2*j2)+1e-24d))*(j2<0d?-1d:1d),
				fac*Math.asin(Math.max(-1,Math.min(1,k2)))
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
		if(rotmat_p2x[0][0]*i+rotmat_p2x[0][1]*j+rotmat_p2x[0][2]*k<0d)
			return new double[] {Double.NaN, Double.NaN};
		return new double[] {
				EARTH_RADIUS_MEAN*(rotmat_p2x[1][0]*i+rotmat_p2x[1][1]*j+rotmat_p2x[1][2]*k),
				EARTH_RADIUS_MEAN*(rotmat_p2x[2][0]*i+rotmat_p2x[2][1]*j+rotmat_p2x[2][2]*k)
		};
	}
	
	@Override
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree) {
		double[] xy = fromLATLONtoPROJ(u, v, input_in_degree);
		return tissotFromProj(xy[0], xy[1]);
	}
	
	@Override
	public double[] tissotFromProj(double x, double y) {
		double i = x/EARTH_RADIUS_MEAN;
		double j = y/EARTH_RADIUS_MEAN;
		double w2 = i*i + j*j;
		if(w2>1)
			return new double[] {Double.NaN,Double.NaN, Double.NaN,Double.NaN};
		double smin = Math.sqrt(1d-w2);
		double smax = 1d;
		return new double[] {-smax*j,smax*i, smin*i,smin*j};
	}
	
	@Override
	public double[] defaultMapExtend() {
		return new double[] {-EARTH_RADIUS_MEAN,EARTH_RADIUS_MEAN,-EARTH_RADIUS_MEAN,EARTH_RADIUS_MEAN};
	}
	
	@Override
	public void drawBorder(JAxis ax, JGroupShape s) {
		JPlotShape.stroke(0xff000000); JPlotShape.strokeWeight(3f);
		int[] p    = ax.getSize();
		double[] r = ax.getRange();
		float x1 = (float) (p[0]+p[2]*Math.max(0d, Math.min(1d, JPlotMath.map(-EARTH_RADIUS_MEAN,r[0],r[1],0d,1d))));
		float y1 = (float) (p[1]+p[3]*Math.max(0d, Math.min(1d, JPlotMath.map(                0d,r[2],r[3],0d,1d))));
		for(int i=-35; i<=36; i++) {
			double bx = EARTH_RADIUS_MEAN*Math.cos(5*i*JPlotMath.DEG_TO_RAD);
			double by = EARTH_RADIUS_MEAN*Math.sin(5*i*JPlotMath.DEG_TO_RAD);
			float x2 = (float) (p[0]+p[2]*Math.max(0d, Math.min(1d, JPlotMath.map(bx,r[0],r[1],0d,1d))));
			float y2 = (float) (p[1]+p[3]*Math.max(0d, Math.min(1d, JPlotMath.map(by,r[2],r[3],1d,0d))));
			s.addChild(new JLineShape(x1, y1, x2, y2));
			x1 = x2;
			y1 = y2;
		}
	}
}
