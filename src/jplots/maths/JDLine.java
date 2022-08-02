package jplots.maths;

import java.util.ArrayList;
import java.util.List;

public class JDLine {
	
	private JDPoint[] points;

	public JDLine(JDPoint... p) {
		points = p.clone();
	}
	
	public JDPoint[] getPoints() {
		return points;
	}
	
	public float[] getCoordsAsFloats() {
		return JPlotMath.toFloatArray1D(getCoords());
	}
	public double[] getCoords() {
		double[] coords = new double[points.length*2];
		for(int i=0; i<points.length; i++) {
			coords[ 2*i ] = points[i].x();
			coords[2*i+1] = points[i].y();
		}
		return coords;
	}
	
	public JDLine affine(double[][] transformationMatrix) {
		for(JDPoint p: points) p.affine(transformationMatrix);
		return this;
	}
	
	public boolean join(JDLine other) {
		return join(other, 0.0001d);
	}
	public boolean join(JDLine other, double delta) {
		JDPoint[] p1 = points.clone();
		int len1 = p1.length;
		JDPoint[] p2 = other.getPoints();
		int len2 = p2.length;
//		if(p1[p1.length-1].equals(p2[0], delta)) {
//			points = new JDPoint[len1+len2-1];
//			for(int i=0; i<len2; i++) points[i+len1-1] = p2[i];
//			for(int i=0; i<len1; i++) points[i]        = p1[i];
//			return true;
//		}
		if(p1[0].equals(p2[p2.length-1], delta)) {
			points = new JDPoint[len1+len2-1];
			for(int i=0; i<len2; i++) points[i]        = p2[i];
			for(int i=0; i<len1; i++) points[i+len2-1] = p1[i];
			return true;
		}
//		if(p1[0].equals(p2[0], delta)) {
//			points = new JDPoint[len1+len2-1];
//			for(int i=0; i<len2; i++) points[i]             = p2[i];
//			for(int i=0; i<len1; i++) points[len1+len2-i-2] = p1[i];
//			return true;
//		}
//		if(p1[p1.length-1].equals(p2[p2.length-1], delta)) {
//			points = new JDPoint[len1+len2-1];
//			for(int i=0; i<len2; i++) points[len1+len2-i-2] = p2[i];
//			for(int i=0; i<len1; i++) points[i]             = p1[i];
//			return true;
//		}
		return false;
	}
	
	public List<JDLine> intersectsCircle(JDPoint center, double radius) {
		List<JDLine> res = new ArrayList<>();
		double er2 = radius * radius;
		List<JDPoint> np = new ArrayList<>();
		double d1 = points[0].x*points[0].x + points[0].y*points[0].y - er2;
		if(d1<=0d) np.add(points[0].copy());
		for(int i=1; i<points.length; i++) {
			double d0 = points[i].x*points[i].x + points[i].y*points[i].y - er2;
			if(d0<=0d) {
				if(d1>0d) {
					JDPoint c = center.circleCrossBetween(points[i], points[i-1], radius);
					if(c!=null) np.add(c);
//					np.add(points[i].fractionTowards(d0/(d0-d1), points[i-1]));
				}
				np.add(points[i].copy());
				if(i+1==points.length) {
					res.add(new JDLine(np.toArray(new JDPoint[0])));
					np.clear();
				}
			} else
			if(d1<=0d) {
				JDPoint c = center.circleCrossBetween(points[i-1], points[i], radius);
				if(c!=null) np.add(c);
//				np.add(points[i-1].fractionTowards(d1/(d1-d0), points[i]));
				res.add(new JDLine(np.toArray(new JDPoint[0])));
				np.clear();
			}
			d1 = d0;
		}
		
		// TODO Auto-generated method stub
		return res;
	}
	public List<JDLine> intersectsAABB(double left, double right, double top, double bottom) {
		double xl = Math.min(left, right);
		double xr = Math.max(left, right);
		double yt = Math.min(top, bottom);
		double yb = Math.max(top, bottom);
		double dx = 1d / (xr-xl), dy = 1d / (yb-yt);
		List<JDPoint> np = new ArrayList<>();
		List<JDPoint[]> sres = new ArrayList<>();
		//check left & right
		double d0 = (points[0].x-xl)*dx;
		//d0 = 0d;
		if(0<=d0 && d0<=1d) np.add(points[0].copy());
		for(int i=1; i<points.length; i++) {
			double di = (points[i].x-xl)*dx;
			double f0 = (0d-di)/(d0-di), f1 = (1d-di)/(d0-di);
			if(0<=di && di<=1d) {
				if(d0<0d) np.add(points[i].fractionTowards(f0, points[i-1]));
				if(d0>1d) np.add(points[i].fractionTowards(f1, points[i-1]));
				np.add(points[i].copy());
				if(i+1==points.length) {
					sres.add(np.toArray(new JDPoint[0]));
					np.clear();
				}
			} else
			if(0d<=d0 && d0<=1d) {
				if(di<0d) np.add(points[i].fractionTowards(f0, points[i-1]));
				if(di>1d) np.add(points[i].fractionTowards(f1, points[i-1]));
				sres.add(np.toArray(new JDPoint[0]));
				np.clear();
			} else
			if(d0*di<0d) {
				np.add(points[i].fractionTowards(f0, points[i-1]));
				np.add(points[i].fractionTowards(f1, points[i-1]));
				sres.add(np.toArray(new JDPoint[0]));
				np.clear();
			}
			d0 = di;
		}
		//check top & bottom
		List<JDLine> res = new ArrayList<>();
		for(JDPoint[] pnts: sres) {
			d0 = (pnts[0].y-yt)*dy;
			if(0<=d0 && d0<=1d) np.add(pnts[0].copy());
			for(int i=1; i<pnts.length; i++) {
				double di = (pnts[i].y-yt)*dy;
				double f0 = (0d-di)/(d0-di), f1 = (1d-di)/(d0-di);
				if(0<=di && di<=1d) {
					if(d0<0d) np.add(pnts[i].fractionTowards(f0, pnts[i-1]));
					if(d0>1d) np.add(pnts[i].fractionTowards(f1, pnts[i-1]));
					np.add(pnts[i].copy());
					if(i+1==pnts.length) {
						res.add(new JDLine(np.toArray(new JDPoint[0])));
						np.clear();
					}
				} else
				if(0d<=d0 && d0<=1d) {
					if(di<0d) np.add(pnts[i].fractionTowards(f0, pnts[i-1]));
					if(di>1d) np.add(pnts[i].fractionTowards(f1, pnts[i-1]));
					if(di<0d || 1d<di) {
						res.add(new JDLine(np.toArray(new JDPoint[0])));
						np.clear();
					}
				} else
				if(d0*di<0d) {
					np.add(pnts[i].fractionTowards(f0, pnts[i-1]));
					np.add(pnts[i].fractionTowards(f1, pnts[i-1]));
					res.add(new JDLine(np.toArray(new JDPoint[0])));
					np.clear();
				}
				d0 = di;
			}
		}
		return res;
	}
	
	public List<JDLine> cutByHalfPlane(double nx, double ny, double no) {
		List<JDLine> res = new ArrayList<>();
		List<JDPoint> np = new ArrayList<>();
		double d1 = points[0].x*nx + points[0].y*ny - no;
		if(d1>=0d) np.add(points[0].copy());
		for(int i=1; i<points.length; i++) {
			double d0 = points[i].x*nx + points[i].y*ny - no;
			if(d0>=0d) {
				if(d1<0d) np.add(points[i].fractionTowards(d0/(d0-d1), points[i-1]));
				np.add(points[i].copy());
				if(i+1==points.length) {
					res.add(new JDLine(np.toArray(new JDPoint[0])));
					np.clear();
				}
			} else
			if(d1>=0d) {
				np.add(points[i-1].fractionTowards(d1/(d1-d0), points[i]));
				res.add(new JDLine(np.toArray(new JDPoint[0])));
				np.clear();
			}
			d1 = d0;
		}
		return res;
	}
}
