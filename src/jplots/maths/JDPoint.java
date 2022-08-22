package jplots.maths;

import org.locationtech.jts.geom.Coordinate;

public class JDPoint implements Comparable<JDPoint> {

	public double x, y, value;
	public int idx, lev;
	private Integer hash = null;

	public JDPoint(Coordinate c) {
		x = c.x;
		y = c.y;
		value = c.z;
	}

	public JDPoint(double _x, double _y) {
		x = _x;
		y = _y;
		value = Double.NaN;
	}

	public JDPoint(double _x, double _y, double _v) {
		x = _x;
		y = _y;
		value = _v;
	}

	public double x() {
		return x;
	}

	public double y() {
		return y;
	}

	public double value() {
		return value;
	}

	public JDPoint copy() {
		return new JDPoint(this.x + 0d, this.y + 0d, this.value + 0d);
	}

	@Override
	public int compareTo(JDPoint p) {
		boolean fxsy = Math.abs(p.x-this.x) >= Math.abs(p.y-this.y);
		if(fxsy)
			return Double.compare(this.x, p.x);
		else
			return Double.compare(this.y, p.y);
	}

	@Override
	public String toString() {
		return "p[" + x + ", " + y + "]";
	}

	@Override
	public int hashCode() {
		if (hash != null) {
			return hash;
		}
		return hash = hash(x, y);
	}

	private static int hash(double x, double y) {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(x);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(y);
		return prime * result + (int) (temp ^ (temp >>> 32));
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if ((obj == null) || (getClass() != obj.getClass())) {
			return false;
		}
		
		JDPoint a = this;
		JDPoint b = (JDPoint) obj;
		if (a.x != b.x || a.y != b.y) {
			return false;
		}
		return true;
	}
	public boolean equals(Object obj, double tolerance) {
		if (this == obj) {
			return true;
		}
		if ((obj == null) || (getClass() != obj.getClass())) {
			return false;
		}

		JDPoint a = this;
		JDPoint b = (JDPoint) obj;
		if (Math.abs(a.x - b.x) > tolerance || Math.abs(a.y - b.y) > tolerance) {
			return false;
		}
		return true;
	}

	public JDPoint affine(double[][] transformationMatrix) {
		double a = x * transformationMatrix[0][0] + y * transformationMatrix[0][1] + transformationMatrix[0][2];
		double b = x * transformationMatrix[1][0] + y * transformationMatrix[1][1] + transformationMatrix[1][2];
		x = a;
		y = b;
		return this;
	}

	public double dist(JDPoint point) {
		return Math.sqrt(dist2(point));
	}
	public double dist2(JDPoint point) {
		return (this.x - point.x) * (this.x - point.x) + (this.y - point.y) * (this.y - point.y);
	}

	public JDPoint fractionTowards(double fraction, JDPoint towards) {
		if(towards==null) return null;
		return new JDPoint(
				this.x + fraction * (towards.x - this.x),
				this.y + fraction * (towards.y - this.y),
				Double.isNaN(this.value) ? towards.value :
				Double.isNaN(towards.value) ? this.value :
				this.value + fraction * (towards.value - this.value));
	}

	/**
	 * calculates the (firs) intersection between a line segment from p1 to p2 and a circle around this point with given radius
	 * 
	 * @param p1      start point of the line segment
	 * @param p2      end point of the line segment
	 * @param radius  the radius of the circle
	 * @return        intersection point or null if the intersections is not between p1 and p2 or no intersection exists
	 */
	public JDPoint circleCrossBetween(JDPoint p1, JDPoint p2, double radius) {
		double dx = p2.x - p1.x;
		double dy = p2.y - p1.y;
		double dr2 = dx * dx + dy * dy;
		double D = (this.x-p1.x)*dy - (this.y-p1.y)*dx;
		double w = dr2 * radius * radius - D * D;
		if(w<0d) return null;
		w = Math.sqrt(w);
		D = (this.x-p1.x)*dx + (this.y-p1.y)*dy;
		double f = D-w;
		if(f<0d) f = D+w;
		if(f*f>dr2) return null;
		double nx = f*dx/dr2;
		double ny = f*dy/dr2;
		return new JDPoint(p1.x + nx, p1.y + ny);
	}
	
}
