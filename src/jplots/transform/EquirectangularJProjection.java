package jplots.transform;

import java.util.ArrayList;
import java.util.List;

import jplots.JAxis;
import jplots.helper.GeometryTools;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
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
	public List<JDPolygon> splitByMapBorder(JDPolygon poly) {
		//double x_crit = 1.9d * radaequ;
		int check_dirs = 0;
		JDPoint[] points = new JDPoint[poly.c.length];
		List<JDPolygon> res = new ArrayList<>();
		double last_x=0d, mw = 2d;
		double minx = Double.POSITIVE_INFINITY, maxx = Double.NEGATIVE_INFINITY;
		double miny = Double.POSITIVE_INFINITY, maxy = Double.NEGATIVE_INFINITY;
		for(int i=0; i<poly.c.length; i++) {
			double next_x = poly.c[i].x / radaequ + 1d;
			double next_y = poly.c[i].y / radpol;
			if(i>0) {
				if(Math.abs(next_x - last_x) > Math.abs(next_x - last_x + mw)) {
					check_dirs += 1;
				} else
				if(Math.abs(next_x - last_x) > Math.abs(next_x - last_x - mw)) {
					check_dirs -= 1;
				}
			}
			double x = next_x + check_dirs*mw;
			double y = next_y;
			points[i] = new JDPoint(x, next_y);
			last_x = next_x;
			if(x < minx) minx = x;
			if(x > maxx) maxx = x;
			if(y < miny) miny = y;
			if(y > maxy) maxy = y;
		}
		double next_x = poly.c[0].x/radaequ + 1d;
		if(Math.abs(next_x - last_x) > Math.abs(next_x - last_x + mw)) {
			check_dirs += 1;
		} else
		if(Math.abs(next_x - last_x) > Math.abs(next_x - last_x - mw)) {
			check_dirs -= 1;
		}
		
		if(check_dirs!=0) {
			System.err.println("Impossible regular polygon, cannot wrap around dateline and map borders!");
			//res.add(poly);
			return res;
		}
		minx *= 0.5d;
		maxx *= 0.5d;
//		System.out.println("[EQUIRECTANGULAR-JPROJ.] x={"+minx+" ... "+maxx+"}  y={"+miny+" ... "+maxy+"}");
		int min_off = (int) minx - (minx<0d ? 1 : 0);
		int max_off = (int) maxx - (maxx<0d ? 1 : 0);
		if(min_off==0 && max_off==0) {
			res.add(poly);
			return res;
		}
		List<JDPoint[]> temp = new ArrayList<>();
		List<JDPoint[]> temp2 = new ArrayList<>();
		double[] left_normal = { -1d, 0d };
		double[] right_normal = { 1d, 0d };
		for(int o=min_off; o<=max_off; o++) {
			double left = 2d*o;
			double right = left + 2d;
			JDPoint ref = new JDPoint(left+1d, 0d);
			temp.clear();
			temp.addAll(GeometryTools.SutherlandHodgmanAlgorithm(points, left_normal, -left, 1.0e-10d, ref));
//			System.out.println("    -> temp size = "+temp.size());
			temp2.clear();
			for(JDPoint[] p: temp)
				temp2.addAll(GeometryTools.SutherlandHodgmanAlgorithm(p, right_normal, right, 1.0e-10d, ref));
//			System.out.println("    -> temp2 size = "+temp2.size());
			for(JDPoint[] p: temp2) {
				JDPoint[] p2 = new JDPoint[p.length];
				for(int i=0; i<p.length; i++)
					p2[i] = new JDPoint( (p[i].x-left-1d)*radaequ, p[i].y*radpol );
				res.add(new JDPolygon(p2));
			}
		}
//		System.out.println("[EQUIRECTANGULAR-JPROJ.] offset={"+min_off+" ... "+max_off+"} found "+res.size()+" polygon(s)");
		return res;
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
