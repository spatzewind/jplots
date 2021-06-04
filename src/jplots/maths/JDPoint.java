package jplots.maths;

public class JDPoint implements Comparable<JDPoint> {

	public double x,y,value;
	public int idx, lev;
	private Integer hash = null;
	
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

	@Override
	public int compareTo(JDPoint p) {
		return this.x != p.x ? Double.compare(this.x, p.x) : Double.compare(this.y, p.y);
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
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}

		JDPoint a = this;
		JDPoint b = (JDPoint) obj;
		if (a.x != b.x || a.y != b.y) {
			return false;
		}
		return true;
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
