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

	public static int hash(double x, double y) {
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
		return new JDPoint(this.x + fraction * (towards.x - this.x), this.y + fraction * (towards.y - this.y),
				Double.isNaN(this.value) ? towards.value
						: Double.isNaN(towards.value) ? this.value
								: this.value + fraction * (towards.value - this.value));
	}

	public JDPoint circleCrossBetween(JDPoint p1, JDPoint p2, double radius) {
		double dx = p2.x - p1.x;
		double dy = p2.y - p1.y;
		double dr = dx * dx + dy * dy;
		double D = (p1.x - this.x) * (p2.y - this.y) - (p2.x - this.x) * (p1.y - this.y);
		double w = dr * radius * radius - D * D;
		if (w < 0d)
			return null;
		w = Math.sqrt(w) * (dy < 0d ? -1d : 1d);
		double nx = (D * dy + w * dx) / dr;
		double ny = (-D * dx + w * dy) / dr;
		double f = (nx * dx + ny * dy) / dr;
		if (f < 0d || f > 1d) {
			nx = (D * dy + w * dx) / dr;
			ny = (-D * dx + w * dy) / dr;
		}
		return new JDPoint(this.x + nx, this.y + ny);
	}

//	double dist2(jpoint v) {
//		return (v.x-this.x)*(v.x-this.x) + (v.y-this.y)*(v.y-this.y);
//	}
//	double[] dir_n(jpoint v) {
//		double dx = this.x - v.x;
//		double dy = this.y - v.y;
//		double dr = Math.sqrt(dx*dx+dy*dy);
//		return new double[] {dx/dr, dy/dr};
//	}
//	void calcCntLevel(double[] contourIntervals) {
//		if(Double.isNaN(value))
//			return;
//		lev = 0;
//		for(int k=0; k<contourIntervals.length; k++)
//			if(contourIntervals[k]<=value)
//				lev++;
//		if(contourIntervals[contourIntervals.length-1]<value)
//			lev++;
//	}

}
