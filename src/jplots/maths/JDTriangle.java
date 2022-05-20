package jplots.maths;

import org.locationtech.jts.geom.Coordinate;

import jplots.helper.GeometryTools;

public class JDTriangle {

	public int idx;
	public int lev;
	public double[] x;
	public double[] y;
	public double[] value;
	public JDEdge ab;
	public JDEdge bc;
	public JDEdge ca;

	private Integer hash = null;

	public JDTriangle(Coordinate[] c) {
		this(new JDPoint(c[0]), new JDPoint(c[1]), new JDPoint(c[2]), 0);
	}

	public JDTriangle(JDPoint a, JDPoint b, JDPoint c) {
		this(a, b, c, 0);
	}

	public JDTriangle(JDPoint a, JDPoint b, JDPoint c, int index) {
		JDPoint[] tmp = { a, b, c };
		// Arrays.sort(tmp);
		this.idx = index;
		x = new double[] { tmp[0].x, tmp[1].x, tmp[2].x };
		y = new double[] { tmp[0].y, tmp[1].y, tmp[2].y };
		value = new double[] { tmp[0].value, tmp[1].value, tmp[2].value };
	}

	public void edges(JDEdge ab, JDEdge bc, JDEdge ca) {
		JDPoint a = getA();
		JDPoint b = getB();
		JDPoint c = getC();
		this.ab = ab.equals(a, b) ? ab : bc.equals(a, b) ? bc : ca;
		this.bc = ab.equals(b, c) ? ab : bc.equals(b, c) ? bc : ca;
		this.ca = ab.equals(c, a) ? ab : bc.equals(c, a) ? bc : ca;
	}

	public void reverse_orientation() {
		double t = x[2];
		x[2] = x[1];
		x[1] = t;
		t = y[2];
		y[2] = y[1];
		y[1] = t;
		t = value[2];
		value[2] = value[1];
		value[1] = t;
		JDPoint a = getA();
		JDPoint b = getB();
		JDPoint c = getC();
		this.ab = new JDEdge(a, b);
		this.bc = new JDEdge(b, c);
		this.ca = new JDEdge(c, a);
	}

	public JDPoint getA() {
		return new JDPoint(x[0], y[0], value[0]);
	}

	public JDPoint getB() {
		return new JDPoint(x[1], y[1], value[1]);
	}

	public JDPoint getC() {
		return new JDPoint(x[2], y[2], value[2]);
	}

	public JDPoint[] getCorners() {
		return new JDPoint[] { getA(), getB(), getC() };
	}

	public JDTriangle copy() {
		return new JDTriangle(getA(), getB(), getC(), idx);
	}

	public JDPolygon toPolygon() {
		return new JDPolygon(getCorners(), idx);
	}

	public void contourIntervalLevel(double[] intervals) {
		JDPoint[] tmpP = getCorners();
		int[] tmpL = { -1, -1, -1 };
		for (int c = 0; c < 3; c++) {
			if (Double.isNaN(tmpP[c].value))
				continue;
			tmpL[c] = 0;
			for (double interval : intervals)
				if (tmpP[c].value > interval)
					tmpL[c]++;
		}
		lev = tmpL[0] + tmpL[1] + tmpL[2];
		if (lev >= 0) {
			lev /= 3;
		} else {
			lev = -1;
		}
	}

	@Override
	public String toString() {
		return "t[(" + x[0] + "," + y[0] + "," + value[0] + ") - (" + x[1] + "," + y[1] + "," + value[1] + ") - ("
				+ x[2] + "," + y[2] + "," + value[2] + ")]";
	}

	@Override
	public int hashCode() {
		if (hash != null) {
			return hash;
		}
		final int prime = 31;
		int hash = 1;
		hash = prime * hash + getA().hashCode();
		hash = prime * hash + getB().hashCode();
		hash = prime * hash + getC().hashCode();
		return hash;
	}

//	public static int hash(JDPoint a, JDPoint b, JDPoint c) {
//		final int prime = 31;
//		int hash = 1;
//		hash = prime * hash + a.hashCode();
//		hash = prime * hash + b.hashCode();
//		hash = prime * hash + c.hashCode();
//		return hash;
//	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if ((obj == null) || (getClass() != obj.getClass()))
			return false;

		JDPoint[] A = this.getCorners();
		JDPoint[] B = ((JDTriangle) obj).getCorners();

		if (A[0].equals(B[0])) {
			return (A[1].equals(B[1]) && A[2].equals(B[2])) || (A[1].equals(B[2]) && A[2].equals(B[1]));
		} else if (A[0].equals(B[1])) {
			return (A[1].equals(B[0]) && A[2].equals(B[2])) || (A[1].equals(B[2]) && A[2].equals(B[0]));
		} else if (A[0].equals(B[2])) {
			return (A[1].equals(B[0]) && A[2].equals(B[1])) || (A[1].equals(B[1]) && A[2].equals(B[0]));
		}

		return false;
	}

	public JDTriangle affine(double[][] tm) {
		for (int i = 0; i < 3; i++) {
			double xx = x[i] * tm[0][0] + y[i] * tm[0][1] + tm[0][2];
			double yy = x[i] * tm[1][0] + y[i] * tm[1][1] + tm[1][2];
			x[i] = xx;
			y[i] = yy;
		}
		if (tm[0][0] * tm[1][1] - tm[0][1] * tm[1][0] < 0d)
			reverse_orientation();
		return this;
	}

	public JDPolygon intersectsAABB(double le, double to, double ri, double bt) {
		if (le > ri)
			return intersectsAABB(ri, to, le, bt);
		if (to > bt)
			return intersectsAABB(le, bt, ri, to);
		// outer check
		double minX = Math.min(x[0], Math.min(x[1], x[2]));
		double maxX = Math.max(x[0], Math.max(x[1], x[2]));
		double minY = Math.min(y[0], Math.min(y[1], y[2]));
		double maxY = Math.max(y[0], Math.max(y[1], y[2]));
		if (maxX < le || minX > ri || maxY < to || minY > bt)
			return null;
		// inner check
		if (le <= minX && maxX <= ri && to <= minY && maxY <= bt)
			return toPolygon();
		// recalc triangle, if it intersects with vertical borders
		double x1 = x[0], y1 = y[0], v1 = value[0];
		double x2 = x[1], y2 = y[1], v2 = value[1];
		double x3 = x[2], y3 = y[2], v3 = value[2];
		double x4 = Double.NaN, y4 = Double.NaN, v4 = Double.NaN;
		double x5 = Double.NaN, y5 = Double.NaN, v5 = Double.NaN;
		double x6 = Double.NaN, y6 = Double.NaN, v6 = Double.NaN;
		double x7 = Double.NaN, y7 = Double.NaN, v7 = Double.NaN;
		int count = 3;
		if (minX < le) {
			int xcode = (x1 < le ? 1 : 0) | (x2 < le ? 2 : 0) | (x3 < le ? 4 : 0);
			switch (xcode) {
			case 1:
				count = 4;
				double xf12 = JPlotMath.map(le, x1, x2, 0d, 1d), xf13 = JPlotMath.map(le, x1, x3, 0d, 1d);
				x4 = JPlotMath.lerp(x1, x3, xf13);
				y4 = JPlotMath.lerp(y1, y3, xf13);
				v4 = JPlotMath.nlerp(v1, v3, xf13);
				x1 = JPlotMath.lerp(x1, x2, xf12);
				y1 = JPlotMath.lerp(y1, y2, xf12);
				v1 = JPlotMath.nlerp(v1, v2, xf12);
				break;
			case 2:
				count = 4;
				double xf21 = JPlotMath.map(le, x2, x1, 0d, 1d), xf23 = JPlotMath.map(le, x2, x3, 0d, 1d);
				x4 = x1;
				y4 = y1;
				v4 = v1;
				x1 = JPlotMath.lerp(x2, x1, xf21);
				y1 = JPlotMath.lerp(y2, y1, xf21);
				v1 = JPlotMath.nlerp(v2, v1, xf21);
				x2 = JPlotMath.lerp(x2, x3, xf23);
				y2 = JPlotMath.lerp(y2, y3, xf23);
				v2 = JPlotMath.nlerp(v2, v3, xf23);
				break;
			case 4:
				count = 4;
				double xf31 = JPlotMath.map(le, x3, x1, 0d, 1d), xf32 = JPlotMath.map(le, x3, x2, 0d, 1d);
				x4 = JPlotMath.lerp(x3, x1, xf31);
				y4 = JPlotMath.lerp(y3, y1, xf31);
				v4 = JPlotMath.nlerp(v3, v1, xf31);
				x3 = JPlotMath.lerp(x3, x2, xf32);
				y3 = JPlotMath.lerp(y3, y1, xf32);
				v3 = JPlotMath.nlerp(v3, v1, xf32);
				break;
			case 3:
				count = 3;
				double fx32 = JPlotMath.map(le, x3, x2, 0d, 1d), fx31 = JPlotMath.map(le, x3, x1, 0d, 1d);
				x1 = JPlotMath.lerp(x3, x1, fx31);
				y1 = JPlotMath.lerp(y3, y1, fx31);
				v1 = JPlotMath.nlerp(v3, v1, fx31);
				x2 = JPlotMath.lerp(x3, x2, fx32);
				y2 = JPlotMath.lerp(y3, y2, fx32);
				v2 = JPlotMath.nlerp(v3, v2, fx32);
				break;
			case 5:
				count = 3;
				double fx21 = JPlotMath.map(le, x2, x1, 0d, 1d), fx23 = JPlotMath.map(le, x2, x3, 0d, 1d);
				x1 = JPlotMath.lerp(x2, x1, fx21);
				y1 = JPlotMath.lerp(y2, y1, fx21);
				v1 = JPlotMath.nlerp(v2, v1, fx21);
				x3 = JPlotMath.lerp(x2, x3, fx23);
				y3 = JPlotMath.lerp(y2, y3, fx23);
				v3 = JPlotMath.nlerp(v2, v3, fx23);
				break;
			case 6:
				count = 3;
				double fx12 = JPlotMath.map(le, x1, x2, 0d, 1d), fx13 = JPlotMath.map(le, x1, x3, 0d, 1d);
				x2 = JPlotMath.lerp(x1, x2, fx12);
				y2 = JPlotMath.lerp(y1, y2, fx12);
				v2 = JPlotMath.nlerp(v1, v2, fx12);
				x3 = JPlotMath.lerp(x1, x3, fx13);
				y3 = JPlotMath.lerp(y1, y3, fx13);
				v3 = JPlotMath.nlerp(v1, v3, fx13);
				break;
			default:
				count = 3;
				break;
			}
		}
		if (maxX > ri) {
			int xcode = (x1 > ri ? 1 : 0) | (x2 > ri ? 2 : 0) | (x3 > ri ? 4 : 0);
			if (count == 4)
				xcode |= (x4 > ri ? 8 : 0);
			if (count == 3) {
				switch (xcode) {
				case 1:
					count = 4;
					double xf12 = JPlotMath.map(ri, x1, x2, 0d, 1d), xf13 = JPlotMath.map(ri, x1, x3, 0d, 1d);
					x4 = JPlotMath.lerp(x1, x3, xf13);
					y4 = JPlotMath.lerp(y1, y3, xf13);
					v4 = JPlotMath.nlerp(v1, v3, xf13);
					x1 = JPlotMath.lerp(x1, x2, xf12);
					y1 = JPlotMath.lerp(y1, y2, xf12);
					v1 = JPlotMath.nlerp(v1, v2, xf12);
					break;
				case 2:
					count = 4;
					double xf21 = JPlotMath.map(ri, x2, x1, 0d, 1d), xf23 = JPlotMath.map(ri, x2, x3, 0d, 1d);
					x4 = x1;
					y4 = y1;
					v4 = v1;
					x1 = JPlotMath.lerp(x2, x1, xf21);
					y1 = JPlotMath.lerp(y2, y1, xf21);
					v1 = JPlotMath.nlerp(v2, v1, xf21);
					x2 = JPlotMath.lerp(x2, x3, xf23);
					y2 = JPlotMath.lerp(y2, y3, xf23);
					v2 = JPlotMath.nlerp(v2, v3, xf23);
					break;
				case 4:
					count = 4;
					double xf31 = JPlotMath.map(ri, x3, x1, 0d, 1d), xf32 = JPlotMath.map(ri, x3, x2, 0d, 1d);
					x4 = JPlotMath.lerp(x3, x1, xf31);
					y4 = JPlotMath.lerp(y3, y1, xf31);
					v4 = JPlotMath.nlerp(v3, v1, xf31);
					x3 = JPlotMath.lerp(x3, x2, xf32);
					y3 = JPlotMath.lerp(y3, y1, xf32);
					v3 = JPlotMath.nlerp(v3, v1, xf32);
					break;
				case 3:
					count = 3;
					double fx32 = JPlotMath.map(ri, x3, x2, 0d, 1d), fx31 = JPlotMath.map(ri, x3, x1, 0d, 1d);
					x1 = JPlotMath.lerp(x3, x1, fx31);
					y1 = JPlotMath.lerp(y3, y1, fx31);
					v1 = JPlotMath.nlerp(v3, v1, fx31);
					x2 = JPlotMath.lerp(x3, x2, fx32);
					y2 = JPlotMath.lerp(y3, y2, fx32);
					v2 = JPlotMath.nlerp(v3, v2, fx32);
					break;
				case 5:
					count = 3;
					double fx21 = JPlotMath.map(ri, x2, x1, 0d, 1d), fx23 = JPlotMath.map(ri, x2, x3, 0d, 1d);
					x1 = JPlotMath.lerp(x2, x1, fx21);
					y1 = JPlotMath.lerp(y2, y1, fx21);
					v1 = JPlotMath.nlerp(v2, v1, fx21);
					x3 = JPlotMath.lerp(x2, x3, fx23);
					y3 = JPlotMath.lerp(y2, y3, fx23);
					v3 = JPlotMath.nlerp(v2, v3, fx23);
					break;
				case 6:
					count = 3;
					double fx12 = JPlotMath.map(ri, x1, x2, 0d, 1d), fx13 = JPlotMath.map(ri, x1, x3, 0d, 1d);
					x2 = JPlotMath.lerp(x1, x2, fx12);
					y2 = JPlotMath.lerp(y1, y2, fx12);
					v2 = JPlotMath.nlerp(v1, v2, fx12);
					x3 = JPlotMath.lerp(x1, x3, fx13);
					y3 = JPlotMath.lerp(y1, y3, fx13);
					v3 = JPlotMath.nlerp(v1, v3, fx13);
					break;
				default:
					count = 3;
					break;
				}
			} else if (count == 4) {
				switch (xcode) {
				case 1:
					count = 5;
					double xf12 = JPlotMath.map(ri, x1, x2, 0d, 1d), xf14 = JPlotMath.map(ri, x1, x4, 0d, 1d);
					x5 = JPlotMath.lerp(x1, x4, xf14);
					y5 = JPlotMath.lerp(y1, y4, xf14);
					v5 = JPlotMath.nlerp(v1, v4, xf14);
					x1 = JPlotMath.lerp(x1, x2, xf12);
					y1 = JPlotMath.lerp(y1, y2, xf12);
					v1 = JPlotMath.nlerp(v1, v2, xf12);
					break;
				case 2:
					count = 5;
					double xf21 = JPlotMath.map(ri, x2, x1, 0d, 1d), xf23 = JPlotMath.map(ri, x2, x3, 0d, 1d);
					x5 = x1;
					y5 = y1;
					v5 = v1;
					x1 = JPlotMath.lerp(x2, x1, xf21);
					y1 = JPlotMath.lerp(y2, y1, xf21);
					v1 = JPlotMath.nlerp(v2, v1, xf21);
					x2 = JPlotMath.lerp(x2, x3, xf23);
					y2 = JPlotMath.lerp(y2, y3, xf23);
					v2 = JPlotMath.nlerp(v2, v3, xf23);
					break;
				case 4:
					count = 5;
					double xf32 = JPlotMath.map(ri, x3, x2, 0d, 1d), xf34 = JPlotMath.map(ri, x3, x4, 0d, 1d);
					x5 = x4;
					y5 = y4;
					v5 = v4;
					x4 = JPlotMath.lerp(x3, x4, xf34);
					y4 = JPlotMath.lerp(y3, y4, xf34);
					v4 = JPlotMath.nlerp(v3, v4, xf34);
					x3 = JPlotMath.lerp(x3, x2, xf32);
					y3 = JPlotMath.lerp(y3, y2, xf32);
					v3 = JPlotMath.nlerp(v3, v2, xf32);
					break;
				case 8:
					count = 5;
					double xf41 = JPlotMath.map(ri, x4, x1, 0d, 1d), xf43 = JPlotMath.map(ri, x4, x3, 0d, 1d);
					x5 = JPlotMath.lerp(x4, x1, xf41);
					y5 = JPlotMath.lerp(y4, y1, xf41);
					v5 = JPlotMath.nlerp(v4, v1, xf41);
					x4 = JPlotMath.lerp(x4, x3, xf43);
					y4 = JPlotMath.lerp(y4, y3, xf43);
					v4 = JPlotMath.nlerp(v4, v3, xf43);
					break;
				case 3:
					count = 4;
					double fx14 = JPlotMath.map(ri, x1, x4, 0d, 1d), fx23 = JPlotMath.map(ri, x2, x3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x4, fx14);
					y1 = JPlotMath.lerp(y1, y4, fx14);
					v1 = JPlotMath.nlerp(v1, v4, fx14);
					x2 = JPlotMath.lerp(x2, x3, fx23);
					y2 = JPlotMath.lerp(y2, y3, fx23);
					v2 = JPlotMath.nlerp(v2, v3, fx23);
					break;
				case 6:
					count = 4;
					double fx21 = JPlotMath.map(ri, x2, x1, 0d, 1d), fx34 = JPlotMath.map(ri, x3, x4, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, fx21);
					y2 = JPlotMath.lerp(y2, y1, fx21);
					v2 = JPlotMath.nlerp(v2, v1, fx21);
					x3 = JPlotMath.lerp(x3, x4, fx34);
					y3 = JPlotMath.lerp(y3, y4, fx34);
					v3 = JPlotMath.nlerp(v3, v4, fx34);
					break;
				case 9:
					count = 4;
					double fx12 = JPlotMath.map(ri, x1, x2, 0d, 1d), fx43 = JPlotMath.map(ri, x4, x3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, fx12);
					y1 = JPlotMath.lerp(y1, y2, fx12);
					v1 = JPlotMath.nlerp(v1, v2, fx12);
					x4 = JPlotMath.lerp(x4, x3, fx43);
					y4 = JPlotMath.lerp(y4, y3, fx43);
					v4 = JPlotMath.nlerp(v4, v3, fx43);
					break;
				case 12:
					count = 4;
					double fx41 = JPlotMath.map(ri, x4, x1, 0d, 1d), fx32 = JPlotMath.map(ri, x3, x2, 0d, 1d);
					x3 = JPlotMath.lerp(x3, x2, fx32);
					y3 = JPlotMath.lerp(y3, y2, fx32);
					v3 = JPlotMath.nlerp(v3, v2, fx32);
					x4 = JPlotMath.lerp(x4, x1, fx41);
					y4 = JPlotMath.lerp(y4, y1, fx41);
					v4 = JPlotMath.nlerp(v4, v1, fx41);
					break;
				default:
					count = 4;
					break;
				}
			}
		}
		minY = Math.min(y1, Math.min(y2, y3));
		if (count > 3 && minY > y4)
			minY = y4;
		if (count > 4 && minY > y5)
			minY = y5;
		maxY = Math.max(y1, Math.max(y2, y3));
		if (count > 3 && maxY < y4)
			maxY = y4;
		if (count > 4 && maxY < y5)
			maxY = y5;
		if (maxY <= to || minY >= bt)
			return null;

		// recalc triangle/quad/pentagon, if it intersects with horizontal borders
		if (minY < to) {
			int ycode = (y1 < to ? 1 : 0) | (y2 < to ? 2 : 0) | (y3 < to ? 4 : 0);
			if (count > 3)
				ycode |= (y4 < to ? 8 : 0);
			if (count > 4)
				ycode |= (y5 < to ? 16 : 0);
			if (count == 3) {
				switch (ycode) {
				case 1:
					count = 4;
					double yy12 = JPlotMath.map(to, y1, y2, 0d, 1d), yy13 = JPlotMath.map(to, y1, y3, 0d, 1d);
					x4 = JPlotMath.lerp(x1, x3, yy13);
					y4 = JPlotMath.lerp(y1, y3, yy13);
					v4 = JPlotMath.nlerp(v1, v3, yy13);
					x1 = JPlotMath.lerp(x1, x2, yy12);
					y1 = JPlotMath.lerp(y1, y2, yy12);
					v1 = JPlotMath.nlerp(v1, v2, yy12);
					break;
				case 2:
					count = 4;
					double yy21 = JPlotMath.map(to, y2, y1, 0d, 1d), yy23 = JPlotMath.map(to, y2, y3, 0d, 1d);
					x4 = x1;
					y4 = y1;
					v4 = v1;
					x1 = JPlotMath.lerp(x2, x1, yy21);
					y1 = JPlotMath.lerp(y2, y1, yy21);
					v1 = JPlotMath.nlerp(v2, v1, yy21);
					x2 = JPlotMath.lerp(x2, x3, yy23);
					y2 = JPlotMath.lerp(y2, y3, yy23);
					v2 = JPlotMath.nlerp(v2, v3, yy23);
					break;
				case 4:
					count = 4;
					double yy31 = JPlotMath.map(to, y3, y1, 0d, 1d), yy32 = JPlotMath.map(to, y3, y2, 0d, 1d);
					x4 = JPlotMath.lerp(x3, x1, yy31);
					y4 = JPlotMath.lerp(y3, y1, yy31);
					v4 = JPlotMath.nlerp(v3, v1, yy31);
					x3 = JPlotMath.lerp(x3, x2, yy32);
					y3 = JPlotMath.lerp(y3, y2, yy32);
					v3 = JPlotMath.nlerp(v3, v2, yy32);
					break;
				case 3:
					count = 3;
					double ff32 = JPlotMath.map(to, y3, y2, 0d, 1d), ff31 = JPlotMath.map(to, y3, y1, 0d, 1d);
					x1 = JPlotMath.lerp(x3, x1, ff31);
					y1 = JPlotMath.lerp(y3, y1, ff31);
					v1 = JPlotMath.nlerp(v3, v1, ff31);
					x2 = JPlotMath.lerp(x3, x2, ff32);
					y2 = JPlotMath.lerp(y3, y2, ff32);
					v2 = JPlotMath.nlerp(v3, v2, ff32);
					break;
				case 5:
					count = 3;
					double ff21 = JPlotMath.map(to, y2, y1, 0d, 1d), ff23 = JPlotMath.map(to, y2, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x2, x1, ff21);
					y1 = JPlotMath.lerp(y2, y1, ff21);
					v1 = JPlotMath.nlerp(v2, v1, ff21);
					x3 = JPlotMath.lerp(x2, x3, ff23);
					y3 = JPlotMath.lerp(y2, y3, ff23);
					v3 = JPlotMath.nlerp(v2, v3, ff23);
					break;
				case 6:
					count = 3;
					double ff12 = JPlotMath.map(to, y1, y2, 0d, 1d), ff13 = JPlotMath.map(to, y1, y3, 0d, 1d);
					x2 = JPlotMath.lerp(x1, x2, ff12);
					y2 = JPlotMath.lerp(y1, y2, ff12);
					v2 = JPlotMath.nlerp(v1, v2, ff12);
					x3 = JPlotMath.lerp(x1, x3, ff13);
					y3 = JPlotMath.lerp(y1, y3, ff13);
					v3 = JPlotMath.nlerp(v1, v3, ff13);
					break;
				default:
					count = 3;
					break;
				}
			} else if (count == 4) {
				switch (ycode) {
				case 1:
					count = 5;
					double yf12 = JPlotMath.map(to, y1, y2, 0d, 1d), yf14 = JPlotMath.map(to, y1, y4, 0d, 1d);
					x5 = JPlotMath.lerp(x1, x4, yf14);
					y5 = JPlotMath.lerp(y1, y4, yf14);
					v5 = JPlotMath.nlerp(v1, v4, yf14);
					x1 = JPlotMath.lerp(x1, x2, yf12);
					y1 = JPlotMath.lerp(y1, y2, yf12);
					v1 = JPlotMath.nlerp(v1, v2, yf12);
					break;
				case 2:
					count = 5;
					double yf21 = JPlotMath.map(to, y2, y1, 0d, 1d), yf23 = JPlotMath.map(to, y2, y3, 0d, 1d);
					x5 = x1;
					y5 = y1;
					v5 = v1;
					x1 = JPlotMath.lerp(x2, x1, yf21);
					y1 = JPlotMath.lerp(y2, y1, yf21);
					v1 = JPlotMath.nlerp(v2, v1, yf21);
					x2 = JPlotMath.lerp(x2, x3, yf23);
					y2 = JPlotMath.lerp(y2, y3, yf23);
					v2 = JPlotMath.nlerp(v2, v3, yf23);
					break;
				case 4:
					count = 5;
					double yf32 = JPlotMath.map(to, y3, y2, 0d, 1d), yf34 = JPlotMath.map(to, y3, y4, 0d, 1d);
					x5 = x4;
					y5 = y4;
					v5 = v4;
					x4 = JPlotMath.lerp(x3, x4, yf34);
					y4 = JPlotMath.lerp(y3, y4, yf34);
					v4 = JPlotMath.nlerp(v3, v4, yf34);
					x3 = JPlotMath.lerp(x3, x2, yf32);
					y3 = JPlotMath.lerp(y3, y2, yf32);
					v3 = JPlotMath.nlerp(v3, v2, yf32);
					break;
				case 8:
					count = 5;
					double yf41 = JPlotMath.map(to, y4, y1, 0d, 1d), yf43 = JPlotMath.map(to, y4, y3, 0d, 1d);
					x5 = JPlotMath.lerp(x4, x1, yf41);
					y5 = JPlotMath.lerp(y4, y1, yf41);
					v5 = JPlotMath.nlerp(v4, v1, yf41);
					x4 = JPlotMath.lerp(x4, x3, yf43);
					y4 = JPlotMath.lerp(y4, y3, yf43);
					v4 = JPlotMath.nlerp(v4, v3, yf43);
					break;
				case 3:
					count = 4;
					double fy14 = JPlotMath.map(to, y1, y4, 0d, 1d), fy23 = JPlotMath.map(to, y2, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x4, fy14);
					y1 = JPlotMath.lerp(y1, y4, fy14);
					v1 = JPlotMath.nlerp(v1, v4, fy14);
					x2 = JPlotMath.lerp(x2, x3, fy23);
					y2 = JPlotMath.lerp(y2, y3, fy23);
					v2 = JPlotMath.nlerp(v2, v3, fy23);
					break;
				case 6:
					count = 4;
					double fy21 = JPlotMath.map(to, y2, y1, 0d, 1d), fy34 = JPlotMath.map(to, y3, y4, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, fy21);
					y2 = JPlotMath.lerp(y2, y1, fy21);
					v2 = JPlotMath.nlerp(v2, v1, fy21);
					x3 = JPlotMath.lerp(x3, x4, fy34);
					y3 = JPlotMath.lerp(y3, y4, fy34);
					v3 = JPlotMath.nlerp(v3, v4, fy34);
					break;
				case 9:
					count = 4;
					double fy12 = JPlotMath.map(to, y1, y2, 0d, 1d), fy43 = JPlotMath.map(to, y4, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, fy12);
					y1 = JPlotMath.lerp(y1, y2, fy12);
					v1 = JPlotMath.nlerp(v1, v2, fy12);
					x4 = JPlotMath.lerp(x4, x3, fy43);
					y4 = JPlotMath.lerp(y4, y3, fy43);
					v4 = JPlotMath.nlerp(v4, v3, fy43);
					break;
				case 12:
					count = 4;
					double fy41 = JPlotMath.map(to, y4, y1, 0d, 1d), fy32 = JPlotMath.map(to, y3, y2, 0d, 1d);
					x3 = JPlotMath.lerp(x3, x2, fy32);
					y3 = JPlotMath.lerp(y3, y2, fy32);
					v3 = JPlotMath.nlerp(v3, v2, fy32);
					x4 = JPlotMath.lerp(x4, x1, fy41);
					y4 = JPlotMath.lerp(y4, y1, fy41);
					v4 = JPlotMath.nlerp(v4, v1, fy41);
					break;
				case 7:
					count = 3;
					double yy41 = JPlotMath.map(to, y4, y1, 0d, 1d), yy43 = JPlotMath.map(to, y4, y3, 0d, 1d);
					x2 = JPlotMath.lerp(x4, x1, yy41);
					y2 = JPlotMath.lerp(y4, y1, yy41);
					v2 = JPlotMath.nlerp(v4, v1, yy41);
					x3 = JPlotMath.lerp(x4, x3, yy43);
					y3 = JPlotMath.lerp(y4, y3, yy43);
					v3 = JPlotMath.nlerp(v4, v3, yy43);
					x1 = x4;
					y1 = y4;
					v1 = v4;
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					break;
				case 11:
					count = 3;
					double yy32 = JPlotMath.map(to, y3, y2, 0d, 1d), yy34 = JPlotMath.map(to, y3, y4, 0d, 1d);
					x1 = JPlotMath.lerp(x3, x4, yy34);
					y1 = JPlotMath.lerp(y3, y4, yy34);
					v1 = JPlotMath.nlerp(v3, v4, yy34);
					x2 = JPlotMath.lerp(x3, x2, yy32);
					y2 = JPlotMath.lerp(y3, y2, yy32);
					v2 = JPlotMath.nlerp(v3, v2, yy32);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					break;
				case 13:
					count = 3;
					double yy21 = JPlotMath.map(to, y2, y1, 0d, 1d), yy23 = JPlotMath.map(to, y2, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x2, x1, yy21);
					y1 = JPlotMath.lerp(y2, y1, yy21);
					v1 = JPlotMath.nlerp(v2, v1, yy21);
					x3 = JPlotMath.lerp(x2, x3, yy23);
					y3 = JPlotMath.lerp(y2, y3, yy23);
					v3 = JPlotMath.nlerp(v2, v3, yy23);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					break;
				case 14:
					count = 3;
					double yy12 = JPlotMath.map(to, y1, y2, 0d, 1d), yy14 = JPlotMath.map(to, y1, y4, 0d, 1d);
					x2 = JPlotMath.lerp(x1, x2, yy12);
					y2 = JPlotMath.lerp(y1, y2, yy12);
					v2 = JPlotMath.nlerp(v1, v2, yy12);
					x3 = JPlotMath.lerp(x1, x4, yy14);
					y3 = JPlotMath.lerp(y1, y4, yy14);
					v3 = JPlotMath.nlerp(v1, v4, yy14);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					break;
				default:
					count = 4;
					break;
				}
			} else if (count == 5) {
				switch (ycode) {
				case 1:
					count = 6;
					double yf12 = JPlotMath.map(to, y1, y2, 0d, 1d), yf15 = JPlotMath.map(to, y1, y5, 0d, 1d);
					x6 = JPlotMath.lerp(x1, x5, yf15);
					y6 = JPlotMath.lerp(y1, y5, yf15);
					v6 = JPlotMath.nlerp(v1, v5, yf15);
					x1 = JPlotMath.lerp(x1, x2, yf12);
					y1 = JPlotMath.lerp(y1, y2, yf12);
					v1 = JPlotMath.nlerp(v1, v2, yf12);
					break;
				case 2:
					count = 6;
					double yf21 = JPlotMath.map(to, y2, y1, 0d, 1d), yf23 = JPlotMath.map(to, y2, y3, 0d, 1d);
					x6 = x1;
					y6 = y1;
					v6 = v1;
					x1 = JPlotMath.lerp(x2, x1, yf21);
					y1 = JPlotMath.lerp(y2, y1, yf21);
					v1 = JPlotMath.nlerp(v2, v1, yf21);
					x2 = JPlotMath.lerp(x2, x3, yf23);
					y2 = JPlotMath.lerp(y2, y3, yf23);
					v2 = JPlotMath.nlerp(v2, v3, yf23);
					break;
				case 4:
					count = 6;
					double yf32 = JPlotMath.map(to, y3, y2, 0d, 1d), yf34 = JPlotMath.map(to, y3, y4, 0d, 1d);
					x6 = x1;
					y6 = y1;
					v6 = v1;
					x1 = x2;
					y1 = y2;
					v1 = v2;
					x2 = JPlotMath.lerp(x3, x2, yf32);
					y2 = JPlotMath.lerp(y3, y2, yf32);
					v2 = JPlotMath.nlerp(v3, v2, yf32);
					x3 = JPlotMath.lerp(x3, x4, yf34);
					y3 = JPlotMath.lerp(y3, y4, yf34);
					v3 = JPlotMath.nlerp(v3, v4, yf34);
					break;
				case 8:
					count = 6;
					double yf43 = JPlotMath.map(to, y4, y3, 0d, 1d), yf45 = JPlotMath.map(to, y4, y5, 0d, 1d);
					x6 = x5;
					y6 = y5;
					v6 = v5;
					x5 = JPlotMath.lerp(x4, x5, yf45);
					y5 = JPlotMath.lerp(y4, y5, yf45);
					v5 = JPlotMath.nlerp(v4, v5, yf45);
					x4 = JPlotMath.lerp(x4, x3, yf43);
					y4 = JPlotMath.lerp(y4, y3, yf43);
					v4 = JPlotMath.nlerp(v4, v3, yf43);
					break;
				case 16:
					count = 6;
					double yf51 = JPlotMath.map(to, y5, y1, 0d, 1d), yf54 = JPlotMath.map(to, y5, y4, 0d, 1d);
					x6 = JPlotMath.lerp(x5, x1, yf51);
					y6 = JPlotMath.lerp(y5, y1, yf51);
					v6 = JPlotMath.nlerp(v5, v1, yf51);
					x5 = JPlotMath.lerp(x5, x4, yf54);
					y5 = JPlotMath.lerp(y5, y4, yf54);
					v5 = JPlotMath.nlerp(v5, v4, yf54);
					break;
				case 3:
					count = 5;
					double fy15 = JPlotMath.map(to, y1, y5, 0d, 1d), fy23 = JPlotMath.map(to, y2, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x5, fy15);
					y1 = JPlotMath.lerp(y1, y5, fy15);
					v1 = JPlotMath.nlerp(v1, v5, fy15);
					x2 = JPlotMath.lerp(x2, x3, fy23);
					y2 = JPlotMath.lerp(y2, y3, fy23);
					v2 = JPlotMath.nlerp(v2, v3, fy23);
					break;
				case 6:
					count = 5;
					double fy21 = JPlotMath.map(to, y2, y1, 0d, 1d), fy34 = JPlotMath.map(to, y3, y4, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, fy21);
					y2 = JPlotMath.lerp(y2, y1, fy21);
					v2 = JPlotMath.nlerp(v2, v1, fy21);
					x3 = JPlotMath.lerp(x3, x4, fy34);
					y3 = JPlotMath.lerp(y3, y4, fy34);
					v3 = JPlotMath.nlerp(v3, v4, fy34);
					break;
				case 12:
					count = 5;
					double fy32 = JPlotMath.map(to, y3, y2, 0d, 1d), fy45 = JPlotMath.map(to, y4, y5, 0d, 1d);
					x3 = JPlotMath.lerp(x3, x2, fy32);
					y3 = JPlotMath.lerp(y3, y2, fy32);
					v3 = JPlotMath.nlerp(v3, v2, fy32);
					x4 = JPlotMath.lerp(x4, x5, fy45);
					y4 = JPlotMath.lerp(y4, y5, fy45);
					v4 = JPlotMath.nlerp(v4, v5, fy45);
					break;
				case 17:
					count = 5;
					double fy12 = JPlotMath.map(to, y1, y2, 0d, 1d), fy54 = JPlotMath.map(to, y5, y4, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, fy12);
					y1 = JPlotMath.lerp(y1, y2, fy12);
					v1 = JPlotMath.nlerp(v1, v2, fy12);
					x5 = JPlotMath.lerp(x5, x4, fy54);
					y5 = JPlotMath.lerp(y5, y4, fy54);
					v5 = JPlotMath.nlerp(v5, v4, fy54);
					break;
				case 24:
					count = 5;
					double fy43 = JPlotMath.map(to, y4, y3, 0d, 1d), fy51 = JPlotMath.map(to, y5, y1, 0d, 1d);
					x4 = JPlotMath.lerp(x4, x3, fy43);
					y4 = JPlotMath.lerp(y4, y3, fy43);
					v4 = JPlotMath.nlerp(v4, v3, fy43);
					x5 = JPlotMath.lerp(x5, x1, fy51);
					y5 = JPlotMath.lerp(y5, y1, fy51);
					v5 = JPlotMath.nlerp(v5, v1, fy51);
					break;
				case 7:
					count = 4;
					double ff15 = JPlotMath.map(to, y1, y5, 0d, 1d), ff34 = JPlotMath.map(to, y3, y4, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x5, ff15);
					y1 = JPlotMath.lerp(y1, y5, ff15);
					v1 = JPlotMath.nlerp(v1, v5, ff15);
					x2 = JPlotMath.lerp(x3, x4, ff34);
					y2 = JPlotMath.lerp(y3, y4, ff34);
					v2 = JPlotMath.nlerp(v3, v4, ff34);
					x3 = x4;
					y3 = y4;
					v3 = v4;
					x4 = x5;
					y4 = y5;
					v4 = v5;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 14:
					count = 4;
					double ff21 = JPlotMath.map(to, y2, y1, 0d, 1d), ff45 = JPlotMath.map(to, y4, y5, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, ff21);
					y2 = JPlotMath.lerp(y2, y1, ff21);
					v2 = JPlotMath.nlerp(v2, v1, ff21);
					x3 = JPlotMath.lerp(x4, x5, ff45);
					y3 = JPlotMath.lerp(y4, y5, ff45);
					v3 = JPlotMath.nlerp(v4, v5, ff45);
					x4 = x5;
					y4 = y5;
					v4 = v5;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 19:
					count = 4;
					double ff23 = JPlotMath.map(to, y2, y3, 0d, 1d), ff54 = JPlotMath.map(to, y5, y4, 0d, 1d);
					x1 = JPlotMath.lerp(x5, x4, ff54);
					y1 = JPlotMath.lerp(y5, y4, ff54);
					v1 = JPlotMath.nlerp(v5, v4, ff54);
					x2 = JPlotMath.lerp(x2, x3, ff23);
					y2 = JPlotMath.lerp(y2, y3, ff23);
					v2 = JPlotMath.nlerp(v2, v3, ff23);
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 25:
					count = 4;
					double ff12 = JPlotMath.map(to, y1, y2, 0d, 1d), ff43 = JPlotMath.map(to, y4, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, ff12);
					y1 = JPlotMath.lerp(y1, y2, ff12);
					v1 = JPlotMath.nlerp(v1, v2, ff12);
					x4 = JPlotMath.lerp(x4, x3, ff43);
					y4 = JPlotMath.lerp(y4, y3, ff43);
					v4 = JPlotMath.nlerp(v4, v3, ff43);
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 28:
					count = 4;
					double ff32 = JPlotMath.map(to, y3, y2, 0d, 1d), ff51 = JPlotMath.map(to, y5, y1, 0d, 1d);
					x3 = JPlotMath.lerp(x3, x2, ff32);
					y3 = JPlotMath.lerp(y3, y2, ff32);
					v3 = JPlotMath.nlerp(v3, v2, ff32);
					x4 = JPlotMath.lerp(x5, x1, ff51);
					y4 = JPlotMath.lerp(y5, y1, ff51);
					v4 = JPlotMath.nlerp(v5, v1, ff51);
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 15:
					count = 3;
					double yy15 = JPlotMath.map(to, y1, y5, 0d, 1d), yy45 = JPlotMath.map(to, y4, y5, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x5, yy15);
					y1 = JPlotMath.lerp(y1, y5, yy15);
					v1 = JPlotMath.nlerp(v1, v5, yy15);
					x2 = JPlotMath.lerp(x4, x5, yy45);
					y2 = JPlotMath.lerp(y4, y5, yy45);
					v2 = JPlotMath.nlerp(v4, v5, yy45);
					x3 = x5;
					y3 = y5;
					v3 = v5;
					x4 = Double.NaN;
					y4 = Double.NaN;
					v5 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 23:
					count = 3;
					double yy34 = JPlotMath.map(to, y3, y4, 0d, 1d), yy54 = JPlotMath.map(to, y5, y4, 0d, 1d);
					x1 = x4;
					y1 = y4;
					v1 = v4;
					x2 = JPlotMath.lerp(x5, x4, yy54);
					y2 = JPlotMath.lerp(y5, y4, yy54);
					v2 = JPlotMath.nlerp(v5, v4, yy54);
					x3 = JPlotMath.lerp(x3, x4, yy34);
					y3 = JPlotMath.lerp(y3, y4, yy34);
					v3 = JPlotMath.nlerp(v3, v4, yy34);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 27:
					count = 3;
					double yy23 = JPlotMath.map(to, y2, y3, 0d, 1d), yy43 = JPlotMath.map(to, y4, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x4, x3, yy43);
					y1 = JPlotMath.lerp(y4, y3, yy43);
					v1 = JPlotMath.nlerp(v4, v3, yy43);
					x2 = JPlotMath.lerp(x2, x3, yy23);
					y2 = JPlotMath.lerp(y2, y3, yy23);
					v2 = JPlotMath.nlerp(v2, v3, yy23);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 29:
					count = 3;
					double yy12 = JPlotMath.map(to, y1, y2, 0d, 1d), yy32 = JPlotMath.map(to, y3, y2, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, yy12);
					y1 = JPlotMath.lerp(y1, y2, yy12);
					v1 = JPlotMath.nlerp(v1, v2, yy12);
					x3 = JPlotMath.lerp(x3, x2, yy32);
					y3 = JPlotMath.lerp(y3, y2, yy32);
					v3 = JPlotMath.nlerp(v3, v2, yy32);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 30:
					count = 3;
					double yy21 = JPlotMath.map(to, y2, y1, 0d, 1d), yy51 = JPlotMath.map(to, y5, y1, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, yy21);
					y2 = JPlotMath.lerp(y2, y1, yy21);
					v2 = JPlotMath.nlerp(v2, v1, yy21);
					x3 = JPlotMath.lerp(x5, x1, yy51);
					y3 = JPlotMath.lerp(y5, y1, yy51);
					v3 = JPlotMath.nlerp(v5, v1, yy51);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				default:
					count = 5;
					break;
				}
			}
		}
		if (maxY > bt) {
			int ycode = (y1 > bt ? 1 : 0) | (y2 > bt ? 2 : 0) | (y3 > bt ? 4 : 0);
			if (count > 3)
				ycode |= (y4 > bt ? 8 : 0);
			if (count > 4)
				ycode |= (y5 > bt ? 16 : 0);
			if (count > 5)
				ycode |= (y6 > bt ? 32 : 0);
			if (count == 3) {
				switch (ycode) {
				case 1:
					count = 4;
					double yy12 = JPlotMath.map(bt, y1, y2, 0d, 1d), yy13 = JPlotMath.map(bt, y1, y3, 0d, 1d);
					x4 = JPlotMath.lerp(x1, x3, yy13);
					y4 = JPlotMath.lerp(y1, y3, yy13);
					v4 = JPlotMath.nlerp(v1, v3, yy13);
					x1 = JPlotMath.lerp(x1, x2, yy12);
					y1 = JPlotMath.lerp(y1, y2, yy12);
					v1 = JPlotMath.nlerp(v1, v2, yy12);
					break;
				case 2:
					count = 4;
					double yy21 = JPlotMath.map(bt, y2, y1, 0d, 1d), yy23 = JPlotMath.map(bt, y2, y3, 0d, 1d);
					x4 = x1;
					y4 = y1;
					v4 = v1;
					x1 = JPlotMath.lerp(x2, x1, yy21);
					y1 = JPlotMath.lerp(y2, y1, yy21);
					v1 = JPlotMath.nlerp(v2, v1, yy21);
					x2 = JPlotMath.lerp(x2, x3, yy23);
					y2 = JPlotMath.lerp(y2, y3, yy23);
					v2 = JPlotMath.nlerp(v2, v3, yy23);
					break;
				case 4:
					count = 4;
					double yy31 = JPlotMath.map(bt, y3, y1, 0d, 1d), yy32 = JPlotMath.map(bt, y3, y2, 0d, 1d);
					x4 = JPlotMath.lerp(x3, x1, yy31);
					y4 = JPlotMath.lerp(y3, y1, yy31);
					v4 = JPlotMath.nlerp(v3, v1, yy31);
					x3 = JPlotMath.lerp(x3, x2, yy32);
					y3 = JPlotMath.lerp(y3, y2, yy32);
					v3 = JPlotMath.nlerp(v3, v2, yy32);
					break;
				case 3:
					count = 3;
					double ff32 = JPlotMath.map(bt, y3, y2, 0d, 1d), ff31 = JPlotMath.map(bt, y3, y1, 0d, 1d);
					x1 = JPlotMath.lerp(x3, x1, ff31);
					y1 = JPlotMath.lerp(y3, y1, ff31);
					v1 = JPlotMath.nlerp(v3, v1, ff31);
					x2 = JPlotMath.lerp(x3, x2, ff32);
					y2 = JPlotMath.lerp(y3, y2, ff32);
					v2 = JPlotMath.nlerp(v3, v2, ff32);
					break;
				case 5:
					count = 3;
					double ff21 = JPlotMath.map(bt, y2, y1, 0d, 1d), ff23 = JPlotMath.map(bt, y2, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x2, x1, ff21);
					y1 = JPlotMath.lerp(y2, y1, ff21);
					v1 = JPlotMath.nlerp(v2, v1, ff21);
					x3 = JPlotMath.lerp(x2, x3, ff23);
					y3 = JPlotMath.lerp(y2, y3, ff23);
					v3 = JPlotMath.nlerp(v2, v3, ff23);
					break;
				case 6:
					count = 3;
					double ff12 = JPlotMath.map(bt, y1, y2, 0d, 1d), ff13 = JPlotMath.map(bt, y1, y3, 0d, 1d);
					x2 = JPlotMath.lerp(x1, x2, ff12);
					y2 = JPlotMath.lerp(y1, y2, ff12);
					v2 = JPlotMath.nlerp(v1, v2, ff12);
					x3 = JPlotMath.lerp(x1, x3, ff13);
					y3 = JPlotMath.lerp(y1, y3, ff13);
					v3 = JPlotMath.nlerp(v1, v3, ff13);
					break;
				default:
					count = 3;
					break;
				}
			} else if (count == 4) {
				switch (ycode) {
				case 1:
					count = 5;
					double yf12 = JPlotMath.map(bt, y1, y2, 0d, 1d), yf14 = JPlotMath.map(bt, y1, y4, 0d, 1d);
					x5 = JPlotMath.lerp(x1, x4, yf14);
					y5 = JPlotMath.lerp(y1, y4, yf14);
					v5 = JPlotMath.nlerp(v1, v4, yf14);
					x1 = JPlotMath.lerp(x1, x2, yf12);
					y1 = JPlotMath.lerp(y1, y2, yf12);
					v1 = JPlotMath.nlerp(v1, v2, yf12);
					break;
				case 2:
					count = 5;
					double yf21 = JPlotMath.map(bt, y2, y1, 0d, 1d), yf23 = JPlotMath.map(bt, y2, y3, 0d, 1d);
					x5 = x1;
					y5 = y1;
					v5 = v1;
					x1 = JPlotMath.lerp(x2, x1, yf21);
					y1 = JPlotMath.lerp(y2, y1, yf21);
					v1 = JPlotMath.nlerp(v2, v1, yf21);
					x2 = JPlotMath.lerp(x2, x3, yf23);
					y2 = JPlotMath.lerp(y2, y3, yf23);
					v2 = JPlotMath.nlerp(v2, v3, yf23);
					break;
				case 4:
					count = 5;
					double yf32 = JPlotMath.map(bt, y3, y2, 0d, 1d), yf34 = JPlotMath.map(bt, y3, y4, 0d, 1d);
					x5 = x4;
					y5 = y4;
					v5 = v4;
					x4 = JPlotMath.lerp(x3, x4, yf34);
					y4 = JPlotMath.lerp(y3, y4, yf34);
					v4 = JPlotMath.nlerp(v3, v4, yf34);
					x3 = JPlotMath.lerp(x3, x2, yf32);
					y3 = JPlotMath.lerp(y3, y2, yf32);
					v3 = JPlotMath.nlerp(v3, v2, yf32);
					break;
				case 8:
					count = 5;
					double yf41 = JPlotMath.map(bt, y4, y1, 0d, 1d), yf43 = JPlotMath.map(bt, y4, y3, 0d, 1d);
					x5 = JPlotMath.lerp(x4, x1, yf41);
					y5 = JPlotMath.lerp(y4, y1, yf41);
					v5 = JPlotMath.nlerp(v4, v1, yf41);
					x4 = JPlotMath.lerp(x4, x3, yf43);
					y4 = JPlotMath.lerp(y4, y3, yf43);
					v4 = JPlotMath.nlerp(v4, v3, yf43);
					break;
				case 3:
					count = 4;
					double fy14 = JPlotMath.map(bt, y1, y4, 0d, 1d), fy23 = JPlotMath.map(bt, y2, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x4, fy14);
					y1 = JPlotMath.lerp(y1, y4, fy14);
					v1 = JPlotMath.nlerp(v1, v4, fy14);
					x2 = JPlotMath.lerp(x2, x3, fy23);
					y2 = JPlotMath.lerp(y2, y3, fy23);
					v2 = JPlotMath.nlerp(v2, v3, fy23);
					break;
				case 6:
					count = 4;
					double fy21 = JPlotMath.map(bt, y2, y1, 0d, 1d), fy34 = JPlotMath.map(bt, y3, y4, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, fy21);
					y2 = JPlotMath.lerp(y2, y1, fy21);
					v2 = JPlotMath.nlerp(v2, v1, fy21);
					x3 = JPlotMath.lerp(x3, x4, fy34);
					y3 = JPlotMath.lerp(y3, y4, fy34);
					v3 = JPlotMath.nlerp(v3, v4, fy34);
					break;
				case 9:
					count = 4;
					double fy12 = JPlotMath.map(bt, y1, y2, 0d, 1d), fy43 = JPlotMath.map(bt, y4, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, fy12);
					y1 = JPlotMath.lerp(y1, y2, fy12);
					v1 = JPlotMath.nlerp(v1, v2, fy12);
					x4 = JPlotMath.lerp(x4, x3, fy43);
					y4 = JPlotMath.lerp(y4, y3, fy43);
					v4 = JPlotMath.nlerp(v4, v3, fy43);
					break;
				case 12:
					count = 4;
					double fy41 = JPlotMath.map(bt, y4, y1, 0d, 1d), fy32 = JPlotMath.map(bt, y3, y2, 0d, 1d);
					x3 = JPlotMath.lerp(x3, x2, fy32);
					y3 = JPlotMath.lerp(y3, y2, fy32);
					v3 = JPlotMath.nlerp(v3, v2, fy32);
					x4 = JPlotMath.lerp(x4, x1, fy41);
					y4 = JPlotMath.lerp(y4, y1, fy41);
					v4 = JPlotMath.nlerp(v4, v1, fy41);
					break;
				case 7:
					count = 3;
					double yy41 = JPlotMath.map(bt, y4, y1, 0d, 1d), yy43 = JPlotMath.map(bt, y4, y3, 0d, 1d);
					x2 = JPlotMath.lerp(x4, x1, yy41);
					y2 = JPlotMath.lerp(y4, y1, yy41);
					v2 = JPlotMath.nlerp(v4, v1, yy41);
					x3 = JPlotMath.lerp(x4, x3, yy43);
					y3 = JPlotMath.lerp(y4, y3, yy43);
					v3 = JPlotMath.nlerp(v4, v3, yy43);
					x1 = x4;
					y1 = y4;
					v1 = v4;
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					break;
				case 11:
					count = 3;
					double yy32 = JPlotMath.map(bt, y3, y2, 0d, 1d), yy34 = JPlotMath.map(bt, y3, y4, 0d, 1d);
					x1 = JPlotMath.lerp(x3, x4, yy34);
					y1 = JPlotMath.lerp(y3, y4, yy34);
					v1 = JPlotMath.nlerp(v3, v4, yy34);
					x2 = JPlotMath.lerp(x3, x2, yy32);
					y2 = JPlotMath.lerp(y3, y2, yy32);
					v2 = JPlotMath.nlerp(v3, v2, yy32);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					break;
				case 13:
					count = 3;
					double yy21 = JPlotMath.map(bt, y2, y1, 0d, 1d), yy23 = JPlotMath.map(bt, y2, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x2, x1, yy21);
					y1 = JPlotMath.lerp(y2, y1, yy21);
					v1 = JPlotMath.nlerp(v2, v1, yy21);
					x3 = JPlotMath.lerp(x2, x3, yy23);
					y3 = JPlotMath.lerp(y2, y3, yy23);
					v3 = JPlotMath.nlerp(v2, v3, yy23);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					break;
				case 14:
					count = 3;
					double yy12 = JPlotMath.map(bt, y1, y2, 0d, 1d), yy14 = JPlotMath.map(bt, y1, y4, 0d, 1d);
					x2 = JPlotMath.lerp(x1, x2, yy12);
					y2 = JPlotMath.lerp(y1, y2, yy12);
					v2 = JPlotMath.nlerp(v1, v2, yy12);
					x3 = JPlotMath.lerp(x1, x4, yy14);
					y3 = JPlotMath.lerp(y1, y4, yy14);
					v3 = JPlotMath.nlerp(v1, v4, yy14);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					break;
				default:
					count = 4;
					break;
				}
			} else if (count == 5) {
				switch (ycode) {
				case 1:
					count = 6;
					double yf12 = JPlotMath.map(bt, y1, y2, 0d, 1d), yf15 = JPlotMath.map(bt, y1, y5, 0d, 1d);
					x6 = JPlotMath.lerp(x1, x5, yf15);
					y6 = JPlotMath.lerp(y1, y5, yf15);
					v6 = JPlotMath.nlerp(v1, v5, yf15);
					x1 = JPlotMath.lerp(x1, x2, yf12);
					y1 = JPlotMath.lerp(y1, y2, yf12);
					v1 = JPlotMath.nlerp(v1, v2, yf12);
					break;
				case 2:
					count = 6;
					double yf21 = JPlotMath.map(bt, y2, y1, 0d, 1d), yf23 = JPlotMath.map(bt, y2, y3, 0d, 1d);
					x6 = x1;
					y6 = y1;
					v6 = v1;
					x1 = JPlotMath.lerp(x2, x1, yf21);
					y1 = JPlotMath.lerp(y2, y1, yf21);
					v1 = JPlotMath.nlerp(v2, v1, yf21);
					x2 = JPlotMath.lerp(x2, x3, yf23);
					y2 = JPlotMath.lerp(y2, y3, yf23);
					v2 = JPlotMath.nlerp(v2, v3, yf23);
					break;
				case 4:
					count = 6;
					double yf32 = JPlotMath.map(bt, y3, y2, 0d, 1d), yf34 = JPlotMath.map(bt, y3, y4, 0d, 1d);
					x6 = x1;
					y6 = y1;
					v6 = v1;
					x1 = x2;
					y1 = y2;
					v1 = v2;
					x2 = JPlotMath.lerp(x3, x2, yf32);
					y2 = JPlotMath.lerp(y3, y2, yf32);
					v2 = JPlotMath.nlerp(v3, v2, yf32);
					x3 = JPlotMath.lerp(x3, x4, yf34);
					y3 = JPlotMath.lerp(y3, y4, yf34);
					v3 = JPlotMath.nlerp(v3, v4, yf34);
					break;
				case 8:
					count = 6;
					double yf43 = JPlotMath.map(bt, y4, y3, 0d, 1d), yf45 = JPlotMath.map(bt, y4, y5, 0d, 1d);
					x6 = x5;
					y6 = y5;
					v6 = v5;
					x5 = JPlotMath.lerp(x4, x5, yf45);
					y5 = JPlotMath.lerp(y4, y5, yf45);
					v5 = JPlotMath.nlerp(v4, v5, yf45);
					x4 = JPlotMath.lerp(x4, x3, yf43);
					y4 = JPlotMath.lerp(y4, y3, yf43);
					v4 = JPlotMath.nlerp(v4, v3, yf43);
					break;
				case 16:
					count = 6;
					double yf51 = JPlotMath.map(bt, y5, y1, 0d, 1d), yf54 = JPlotMath.map(bt, y5, y4, 0d, 1d);
					x6 = JPlotMath.lerp(x5, x1, yf51);
					y6 = JPlotMath.lerp(y5, y1, yf51);
					v6 = JPlotMath.nlerp(v5, v1, yf51);
					x5 = JPlotMath.lerp(x5, x4, yf54);
					y5 = JPlotMath.lerp(y5, y4, yf54);
					v5 = JPlotMath.nlerp(v5, v4, yf54);
					break;
				case 3:
					count = 5;
					double fy15 = JPlotMath.map(bt, y1, y5, 0d, 1d), fy23 = JPlotMath.map(bt, y2, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x5, fy15);
					y1 = JPlotMath.lerp(y1, y5, fy15);
					v1 = JPlotMath.nlerp(v1, v5, fy15);
					x2 = JPlotMath.lerp(x2, x3, fy23);
					y2 = JPlotMath.lerp(y2, y3, fy23);
					v2 = JPlotMath.nlerp(v2, v3, fy23);
					break;
				case 6:
					count = 5;
					double fy21 = JPlotMath.map(bt, y2, y1, 0d, 1d), fy34 = JPlotMath.map(bt, y3, y4, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, fy21);
					y2 = JPlotMath.lerp(y2, y1, fy21);
					v2 = JPlotMath.nlerp(v2, v1, fy21);
					x3 = JPlotMath.lerp(x3, x4, fy34);
					y3 = JPlotMath.lerp(y3, y4, fy34);
					v3 = JPlotMath.nlerp(v3, v4, fy34);
					break;
				case 12:
					count = 5;
					double fy32 = JPlotMath.map(bt, y3, y2, 0d, 1d), fy45 = JPlotMath.map(bt, y4, y5, 0d, 1d);
					x3 = JPlotMath.lerp(x3, x2, fy32);
					y3 = JPlotMath.lerp(y3, y2, fy32);
					v3 = JPlotMath.nlerp(v3, v2, fy32);
					x4 = JPlotMath.lerp(x4, x5, fy45);
					y4 = JPlotMath.lerp(y4, y5, fy45);
					v4 = JPlotMath.nlerp(v4, v5, fy45);
					break;
				case 17:
					count = 5;
					double fy12 = JPlotMath.map(bt, y1, y2, 0d, 1d), fy54 = JPlotMath.map(bt, y5, y4, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, fy12);
					y1 = JPlotMath.lerp(y1, y2, fy12);
					v1 = JPlotMath.nlerp(v1, v2, fy12);
					x5 = JPlotMath.lerp(x5, x4, fy54);
					y5 = JPlotMath.lerp(y5, y4, fy54);
					v5 = JPlotMath.nlerp(v5, v4, fy54);
					break;
				case 24:
					count = 5;
					double fy43 = JPlotMath.map(bt, y4, y3, 0d, 1d), fy51 = JPlotMath.map(bt, y5, y1, 0d, 1d);
					x4 = JPlotMath.lerp(x4, x3, fy43);
					y4 = JPlotMath.lerp(y4, y3, fy43);
					v4 = JPlotMath.nlerp(v4, v3, fy43);
					x5 = JPlotMath.lerp(x5, x1, fy51);
					y5 = JPlotMath.lerp(y5, y1, fy51);
					v5 = JPlotMath.nlerp(v5, v1, fy51);
					break;
				case 7:
					count = 4;
					double ff15 = JPlotMath.map(bt, y1, y5, 0d, 1d), ff34 = JPlotMath.map(bt, y3, y4, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x5, ff15);
					y1 = JPlotMath.lerp(y1, y5, ff15);
					v1 = JPlotMath.nlerp(v1, v5, ff15);
					x2 = JPlotMath.lerp(x3, x4, ff34);
					y2 = JPlotMath.lerp(y3, y4, ff34);
					v2 = JPlotMath.nlerp(v3, v4, ff34);
					x3 = x4;
					y3 = y4;
					v3 = v4;
					x4 = x5;
					y4 = y5;
					v4 = v5;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 14:
					count = 4;
					double ff21 = JPlotMath.map(bt, y2, y1, 0d, 1d), ff45 = JPlotMath.map(bt, y4, y5, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, ff21);
					y2 = JPlotMath.lerp(y2, y1, ff21);
					v2 = JPlotMath.nlerp(v2, v1, ff21);
					x3 = JPlotMath.lerp(x4, x5, ff45);
					y3 = JPlotMath.lerp(y4, y5, ff45);
					v3 = JPlotMath.nlerp(v4, v5, ff45);
					x4 = x5;
					y4 = y5;
					v4 = v5;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 19:
					count = 4;
					double ff23 = JPlotMath.map(bt, y2, y3, 0d, 1d), ff54 = JPlotMath.map(bt, y5, y4, 0d, 1d);
					x1 = JPlotMath.lerp(x5, x4, ff54);
					y1 = JPlotMath.lerp(y5, y4, ff54);
					v1 = JPlotMath.nlerp(v5, v4, ff54);
					x2 = JPlotMath.lerp(x2, x3, ff23);
					y2 = JPlotMath.lerp(y2, y3, ff23);
					v2 = JPlotMath.nlerp(v2, v3, ff23);
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 25:
					count = 4;
					double ff12 = JPlotMath.map(bt, y1, y2, 0d, 1d), ff43 = JPlotMath.map(bt, y4, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, ff12);
					y1 = JPlotMath.lerp(y1, y2, ff12);
					v1 = JPlotMath.nlerp(v1, v2, ff12);
					x4 = JPlotMath.lerp(x4, x3, ff43);
					y4 = JPlotMath.lerp(y4, y3, ff43);
					v4 = JPlotMath.nlerp(v4, v3, ff43);
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 28:
					count = 4;
					double ff32 = JPlotMath.map(bt, y3, y2, 0d, 1d), ff51 = JPlotMath.map(bt, y5, y1, 0d, 1d);
					x3 = JPlotMath.lerp(x3, x2, ff32);
					y3 = JPlotMath.lerp(y3, y2, ff32);
					v3 = JPlotMath.nlerp(v3, v2, ff32);
					x4 = JPlotMath.lerp(x5, x1, ff51);
					y4 = JPlotMath.lerp(y5, y1, ff51);
					v4 = JPlotMath.nlerp(v5, v1, ff51);
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 15:
					count = 3;
					double yy15 = JPlotMath.map(bt, y1, y5, 0d, 1d), yy45 = JPlotMath.map(bt, y4, y5, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x5, yy15);
					y1 = JPlotMath.lerp(y1, y5, yy15);
					v1 = JPlotMath.nlerp(v1, v5, yy15);
					x2 = JPlotMath.lerp(x4, x5, yy45);
					y2 = JPlotMath.lerp(y4, y5, yy45);
					v2 = JPlotMath.nlerp(v4, v5, yy45);
					x3 = x5;
					y3 = y5;
					v3 = v5;
					x4 = Double.NaN;
					y4 = Double.NaN;
					v5 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 23:
					count = 3;
					double yy34 = JPlotMath.map(bt, y3, y4, 0d, 1d), yy54 = JPlotMath.map(bt, y5, y4, 0d, 1d);
					x1 = x4;
					y1 = y4;
					v1 = v4;
					x2 = JPlotMath.lerp(x5, x4, yy54);
					y2 = JPlotMath.lerp(y5, y4, yy54);
					v2 = JPlotMath.nlerp(v5, v4, yy54);
					x3 = JPlotMath.lerp(x3, x4, yy34);
					y3 = JPlotMath.lerp(y3, y4, yy34);
					v3 = JPlotMath.nlerp(v3, v4, yy34);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 27:
					count = 3;
					double yy23 = JPlotMath.map(bt, y2, y3, 0d, 1d), yy43 = JPlotMath.map(bt, y4, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x4, x3, yy43);
					y1 = JPlotMath.lerp(y4, y3, yy43);
					v1 = JPlotMath.nlerp(v4, v3, yy43);
					x2 = JPlotMath.lerp(x2, x3, yy23);
					y2 = JPlotMath.lerp(y2, y3, yy23);
					v2 = JPlotMath.nlerp(v2, v3, yy23);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 29:
					count = 3;
					double yy12 = JPlotMath.map(bt, y1, y2, 0d, 1d), yy32 = JPlotMath.map(bt, y3, y2, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, yy12);
					y1 = JPlotMath.lerp(y1, y2, yy12);
					v1 = JPlotMath.nlerp(v1, v2, yy12);
					x3 = JPlotMath.lerp(x3, x2, yy32);
					y3 = JPlotMath.lerp(y3, y2, yy32);
					v3 = JPlotMath.nlerp(v3, v2, yy32);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				case 30:
					count = 3;
					double yy21 = JPlotMath.map(bt, y2, y1, 0d, 1d), yy51 = JPlotMath.map(bt, y5, y1, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, yy21);
					y2 = JPlotMath.lerp(y2, y1, yy21);
					v2 = JPlotMath.nlerp(v2, v1, yy21);
					x3 = JPlotMath.lerp(x5, x1, yy51);
					y3 = JPlotMath.lerp(y5, y1, yy51);
					v3 = JPlotMath.nlerp(v5, v1, yy51);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					break;
				default:
					count = 5;
					break;
				}
			} else if (count == 6) {
				switch (ycode) {
				case 1:
					count = 7;
					double yf12 = JPlotMath.map(bt, y1, y2, 0d, 1d), yf16 = JPlotMath.map(bt, y1, y6, 0d, 1d);
					x7 = JPlotMath.lerp(x1, x6, yf16);
					y7 = JPlotMath.lerp(y1, y6, yf16);
					v7 = JPlotMath.nlerp(v1, v6, yf16);
					x1 = JPlotMath.lerp(x1, x2, yf12);
					y1 = JPlotMath.lerp(y1, y2, yf12);
					v1 = JPlotMath.nlerp(v1, v2, yf12);
					break;
				case 2:
					count = 7;
					double yf21 = JPlotMath.map(bt, y2, y1, 0d, 1d), yf23 = JPlotMath.map(bt, y2, y3, 0d, 1d);
					x7 = x1;
					y7 = y1;
					v7 = v1;
					x1 = JPlotMath.lerp(x2, x1, yf21);
					y1 = JPlotMath.lerp(y2, y1, yf21);
					v1 = JPlotMath.nlerp(v2, v1, yf21);
					x2 = JPlotMath.lerp(x2, x3, yf23);
					y2 = JPlotMath.lerp(y2, y3, yf23);
					v2 = JPlotMath.nlerp(v2, v3, yf23);
					break;
				case 4:
					count = 7;
					double yf32 = JPlotMath.map(bt, y3, y2, 0d, 1d), yf34 = JPlotMath.map(bt, y3, y4, 0d, 1d);
					x7 = x1;
					y7 = y1;
					v7 = v1;
					x1 = x2;
					y1 = y2;
					v1 = v2;
					x2 = JPlotMath.lerp(x3, x2, yf32);
					y2 = JPlotMath.lerp(y3, y2, yf32);
					v2 = JPlotMath.nlerp(v3, v2, yf32);
					x3 = JPlotMath.lerp(x3, x4, yf34);
					y3 = JPlotMath.lerp(y3, y4, yf34);
					v3 = JPlotMath.nlerp(v3, v4, yf34);
					break;
				case 8:
					count = 7;
					double yf43 = JPlotMath.map(bt, y4, y3, 0d, 1d), yf45 = JPlotMath.map(bt, y4, y5, 0d, 1d);
					x7 = x6;
					y7 = y6;
					v7 = v6;
					x6 = x5;
					y6 = y5;
					v6 = v5;
					x5 = JPlotMath.lerp(x4, x5, yf45);
					y5 = JPlotMath.lerp(y4, y5, yf45);
					v5 = JPlotMath.nlerp(v4, v5, yf45);
					x4 = JPlotMath.lerp(x4, x3, yf43);
					y4 = JPlotMath.lerp(y4, y3, yf43);
					v4 = JPlotMath.nlerp(v4, v3, yf43);
					break;
				case 16:
					count = 7;
					double yf54 = JPlotMath.map(bt, y5, y4, 0d, 1d), yf56 = JPlotMath.map(bt, y5, y6, 0d, 1d);
					x7 = x6;
					y7 = y6;
					v7 = v6;
					x6 = JPlotMath.lerp(x5, x6, yf56);
					y6 = JPlotMath.lerp(y5, y6, yf56);
					v6 = JPlotMath.nlerp(v5, v6, yf56);
					x5 = JPlotMath.lerp(x5, x4, yf54);
					y5 = JPlotMath.lerp(y5, y4, yf54);
					v5 = JPlotMath.nlerp(v5, v4, yf54);
					break;
				case 32:
					count = 7;
					double yf61 = JPlotMath.map(bt, y6, y1, 0d, 1d), yf65 = JPlotMath.map(bt, y6, y5, 0d, 1d);
					x7 = JPlotMath.lerp(x6, x1, yf61);
					y7 = JPlotMath.lerp(y6, y1, yf61);
					v7 = JPlotMath.nlerp(v6, v1, yf61);
					x6 = JPlotMath.lerp(x6, x5, yf65);
					y6 = JPlotMath.lerp(y6, y5, yf65);
					v6 = JPlotMath.nlerp(v6, v5, yf65);
					break;
				case 3:
					count = 6;
					double fy16 = JPlotMath.map(bt, y1, y6, 0d, 1d), fy23 = JPlotMath.map(bt, y2, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x6, fy16);
					y1 = JPlotMath.lerp(y1, y6, fy16);
					v1 = JPlotMath.nlerp(v1, v6, fy16);
					x2 = JPlotMath.lerp(x2, x3, fy23);
					y2 = JPlotMath.lerp(y2, y3, fy23);
					v2 = JPlotMath.nlerp(v2, v3, fy23);
					break;
				case 6:
					count = 6;
					double fy21 = JPlotMath.map(bt, y2, y1, 0d, 1d), fy34 = JPlotMath.map(bt, y3, y4, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, fy21);
					y2 = JPlotMath.lerp(y2, y1, fy21);
					v2 = JPlotMath.nlerp(v2, v1, fy21);
					x3 = JPlotMath.lerp(x3, x4, fy34);
					y3 = JPlotMath.lerp(y3, y4, fy34);
					v3 = JPlotMath.nlerp(v3, v4, fy34);
					break;
				case 12:
					count = 6;
					double fy32 = JPlotMath.map(bt, y3, y2, 0d, 1d), fy45 = JPlotMath.map(bt, y4, y5, 0d, 1d);
					x3 = JPlotMath.lerp(x3, x2, fy32);
					y3 = JPlotMath.lerp(y3, y2, fy32);
					v3 = JPlotMath.nlerp(v3, v2, fy32);
					x4 = JPlotMath.lerp(x4, x5, fy45);
					y4 = JPlotMath.lerp(y4, y5, fy45);
					v4 = JPlotMath.nlerp(v4, v5, fy45);
					break;
				case 24:
					count = 6;
					double fy43 = JPlotMath.map(bt, y4, y3, 0d, 1d), fy56 = JPlotMath.map(bt, y5, y6, 0d, 1d);
					x4 = JPlotMath.lerp(x4, x3, fy43);
					y4 = JPlotMath.lerp(y4, y3, fy43);
					v4 = JPlotMath.nlerp(v4, v3, fy43);
					x5 = JPlotMath.lerp(x5, x6, fy56);
					y5 = JPlotMath.lerp(y5, y6, fy56);
					v5 = JPlotMath.nlerp(v5, v6, fy56);
					break;
				case 48:
					count = 6;
					double fy54 = JPlotMath.map(bt, y5, y4, 0d, 1d), fy61 = JPlotMath.map(bt, y6, y1, 0d, 1d);
					x5 = JPlotMath.lerp(x5, x4, fy54);
					y5 = JPlotMath.lerp(y5, y4, fy54);
					v5 = JPlotMath.nlerp(v5, v4, fy54);
					x6 = JPlotMath.lerp(x6, x1, fy61);
					y6 = JPlotMath.lerp(y6, y1, fy61);
					v6 = JPlotMath.nlerp(v6, v1, fy61);
					break;
				case 33:
					count = 6;
					double fy12 = JPlotMath.map(bt, y1, y2, 0d, 1d), fy65 = JPlotMath.map(bt, y6, y5, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, fy12);
					y1 = JPlotMath.lerp(y1, y2, fy12);
					v1 = JPlotMath.nlerp(v1, v2, fy12);
					x6 = JPlotMath.lerp(x6, x5, fy65);
					y6 = JPlotMath.lerp(y6, y5, fy65);
					v6 = JPlotMath.nlerp(v6, v5, fy65);
					break;
				case 7:
					count = 5;
					double yy16 = JPlotMath.map(bt, y1, y6, 0d, 1d), yy34 = JPlotMath.map(bt, y3, y4, 0d, 1d);
					x2 = JPlotMath.lerp(x1, x6, yy16);
					y2 = JPlotMath.lerp(y1, y6, yy16);
					v2 = JPlotMath.nlerp(v1, v6, yy16);
					x3 = JPlotMath.lerp(x3, x4, yy34);
					y3 = JPlotMath.lerp(y3, y4, yy34);
					v3 = JPlotMath.nlerp(v3, v4, yy34);
					x1 = x6;
					y1 = y6;
					v1 = v6;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 14:
					count = 5;
					double yy21 = JPlotMath.map(bt, y2, y1, 0d, 1d), yy45 = JPlotMath.map(bt, y4, y5, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, yy21);
					y2 = JPlotMath.lerp(y2, y1, yy21);
					v2 = JPlotMath.nlerp(v2, v1, yy21);
					x3 = JPlotMath.lerp(x4, x5, yy45);
					y3 = JPlotMath.lerp(y4, y5, yy45);
					v3 = JPlotMath.nlerp(v4, v5, yy45);
					x4 = x5;
					y4 = y5;
					v4 = v5;
					x5 = x6;
					y5 = y6;
					v5 = v6;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 28:
					count = 5;
					double yy32 = JPlotMath.map(bt, y3, y2, 0d, 1d), yy56 = JPlotMath.map(bt, y5, y6, 0d, 1d);
					x3 = JPlotMath.lerp(x3, x2, yy32);
					y3 = JPlotMath.lerp(y3, y2, yy32);
					v3 = JPlotMath.nlerp(v3, v2, yy32);
					x4 = JPlotMath.lerp(x5, x6, yy56);
					y4 = JPlotMath.lerp(y5, y6, yy56);
					v4 = JPlotMath.nlerp(v5, v6, yy56);
					x5 = x6;
					y5 = y6;
					v5 = v6;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 36:
					count = 5;
					double yy23 = JPlotMath.map(bt, y2, y3, 0d, 1d), yy65 = JPlotMath.map(bt, y6, y5, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x3, yy23);
					y2 = JPlotMath.lerp(y2, y3, yy23);
					v2 = JPlotMath.nlerp(v2, v3, yy23);
					x1 = JPlotMath.lerp(x6, x5, yy65);
					y1 = JPlotMath.lerp(y6, y5, yy65);
					v1 = JPlotMath.nlerp(v6, v5, yy65);
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 49:
					count = 5;
					double yy12 = JPlotMath.map(bt, y1, y2, 0d, 1d), yy54 = JPlotMath.map(bt, y5, y4, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, yy12);
					y1 = JPlotMath.lerp(y1, y2, yy12);
					v1 = JPlotMath.nlerp(v1, v2, yy12);
					x5 = JPlotMath.lerp(x5, x4, yy54);
					y5 = JPlotMath.lerp(y5, y4, yy54);
					v5 = JPlotMath.nlerp(v5, v4, yy54);
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 56:
					count = 5;
					double yy43 = JPlotMath.map(bt, y4, y3, 0d, 1d), yy61 = JPlotMath.map(bt, y6, y1, 0d, 1d);
					x4 = JPlotMath.lerp(x4, x3, yy43);
					y4 = JPlotMath.lerp(y4, y3, yy43);
					v4 = JPlotMath.nlerp(v4, v3, yy43);
					x5 = JPlotMath.lerp(x6, x1, yy61);
					y5 = JPlotMath.lerp(y6, y1, yy61);
					v5 = JPlotMath.nlerp(v6, v1, yy61);
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 30:
					count = 4;
					double ff21 = JPlotMath.map(bt, y2, y1, 0d, 1d), ff56 = JPlotMath.map(bt, y5, y6, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, ff21);
					y2 = JPlotMath.lerp(y2, y1, ff21);
					v2 = JPlotMath.nlerp(v2, v1, ff21);
					x3 = JPlotMath.lerp(x5, x6, ff56);
					y3 = JPlotMath.lerp(y5, y6, ff56);
					v3 = JPlotMath.nlerp(v5, v6, ff56);
					x4 = x6;
					y4 = y6;
					v4 = v6;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 15:
					count = 4;
					double ff16 = JPlotMath.map(bt, y1, y6, 0d, 1d), ff45 = JPlotMath.map(bt, y4, y5, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x6, ff16);
					y1 = JPlotMath.lerp(y1, y6, ff16);
					v1 = JPlotMath.nlerp(v1, v6, ff16);
					x2 = JPlotMath.lerp(x4, x5, ff45);
					y2 = JPlotMath.lerp(y4, y5, ff45);
					v2 = JPlotMath.nlerp(v4, v5, ff45);
					x3 = x5;
					y3 = y5;
					v3 = v5;
					x4 = x6;
					y4 = y6;
					v4 = v6;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 39:
					count = 4;
					double ff34 = JPlotMath.map(bt, y3, y4, 0d, 1d), ff65 = JPlotMath.map(bt, y6, y5, 0d, 1d);
					x1 = x5;
					y1 = y5;
					v1 = v5;
					x2 = JPlotMath.lerp(x6, x5, ff65);
					y2 = JPlotMath.lerp(y6, y5, ff65);
					v2 = JPlotMath.nlerp(v6, v5, ff65);
					x3 = JPlotMath.lerp(x3, x4, ff34);
					y3 = JPlotMath.lerp(y3, y4, ff34);
					v3 = JPlotMath.nlerp(v3, v4, ff34);
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 51:
					count = 4;
					double ff23 = JPlotMath.map(bt, y2, y3, 0d, 1d), ff54 = JPlotMath.map(bt, y5, y4, 0d, 1d);
					x1 = JPlotMath.lerp(x5, x4, ff54);
					y1 = JPlotMath.lerp(y5, y4, ff54);
					v1 = JPlotMath.nlerp(v5, v4, ff54);
					x2 = JPlotMath.lerp(x2, x3, ff23);
					y2 = JPlotMath.lerp(y2, y3, ff23);
					v2 = JPlotMath.nlerp(v2, v3, ff23);
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 57:
					count = 4;
					double ff12 = JPlotMath.map(bt, y1, y2, 0d, 1d), ff43 = JPlotMath.map(bt, y4, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, ff12);
					y1 = JPlotMath.lerp(y1, y2, ff12);
					v1 = JPlotMath.nlerp(v1, v2, ff12);
					x4 = JPlotMath.lerp(x4, x3, ff43);
					y4 = JPlotMath.lerp(y4, y3, ff43);
					v4 = JPlotMath.nlerp(v4, v3, ff43);
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 60:
					count = 4;
					double ff61 = JPlotMath.map(bt, y6, y1, 0d, 1d), ff32 = JPlotMath.map(bt, y3, y2, 0d, 1d);
					x3 = JPlotMath.lerp(x3, x2, ff32);
					y3 = JPlotMath.lerp(y3, y2, ff32);
					v3 = JPlotMath.nlerp(v3, v2, ff32);
					x4 = JPlotMath.lerp(x6, x1, ff61);
					y4 = JPlotMath.lerp(y6, y1, ff61);
					v4 = JPlotMath.nlerp(v6, v1, ff61);
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 31:
					count = 3;
					double jj16 = JPlotMath.map(bt, y1, y6, 0d, 1d), jj56 = JPlotMath.map(bt, y5, y6, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x6, jj16);
					y1 = JPlotMath.lerp(y1, y6, jj16);
					v1 = JPlotMath.nlerp(v1, v6, jj16);
					x2 = JPlotMath.lerp(x5, x6, jj56);
					y2 = JPlotMath.lerp(y5, y6, jj56);
					v2 = JPlotMath.nlerp(v5, v6, jj56);
					x3 = x6;
					y3 = y6;
					v3 = v6;
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 47:
					count = 3;
					double jj45 = JPlotMath.map(bt, y4, y5, 0d, 1d), jj65 = JPlotMath.map(bt, y6, y5, 0d, 1d);
					x1 = JPlotMath.lerp(x4, x5, jj45);
					y1 = JPlotMath.lerp(y4, y5, jj45);
					v1 = JPlotMath.nlerp(v4, v5, jj45);
					x3 = JPlotMath.lerp(x6, x5, jj65);
					y3 = JPlotMath.lerp(y6, y5, jj65);
					v3 = JPlotMath.nlerp(v6, v5, jj65);
					x2 = x5;
					y2 = y5;
					v2 = v5;
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 55:
					count = 3;
					double jj34 = JPlotMath.map(bt, y3, y4, 0d, 1d), jj54 = JPlotMath.map(bt, y5, y4, 0d, 1d);
					x2 = JPlotMath.lerp(x5, x4, jj54);
					y2 = JPlotMath.lerp(y5, y4, jj54);
					v2 = JPlotMath.nlerp(v5, v4, jj54);
					x3 = JPlotMath.lerp(x3, x4, jj34);
					y3 = JPlotMath.lerp(y3, y4, jj34);
					v3 = JPlotMath.nlerp(v3, v4, jj34);
					x1 = x4;
					y1 = y4;
					v1 = v4;
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 59:
					count = 3;
					double jj23 = JPlotMath.map(bt, y2, y3, 0d, 1d), jj43 = JPlotMath.map(bt, y4, y3, 0d, 1d);
					x1 = JPlotMath.lerp(x4, x3, jj43);
					y1 = JPlotMath.lerp(y4, y3, jj43);
					v1 = JPlotMath.nlerp(v4, v3, jj43);
					x2 = JPlotMath.lerp(x2, x3, jj23);
					y2 = JPlotMath.lerp(y2, y3, jj23);
					v2 = JPlotMath.nlerp(v2, v3, jj23);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 61:
					count = 3;
					double jj12 = JPlotMath.map(bt, y1, y2, 0d, 1d), jj32 = JPlotMath.map(bt, y3, y2, 0d, 1d);
					x1 = JPlotMath.lerp(x1, x2, jj12);
					y1 = JPlotMath.lerp(y1, y2, jj12);
					v1 = JPlotMath.nlerp(v1, v2, jj12);
					x3 = JPlotMath.lerp(x3, x2, jj32);
					y3 = JPlotMath.lerp(y3, y2, jj32);
					v3 = JPlotMath.nlerp(v3, v2, jj32);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				case 62:
					count = 3;
					double jj21 = JPlotMath.map(bt, y2, y1, 0d, 1d), jj61 = JPlotMath.map(bt, y6, y1, 0d, 1d);
					x2 = JPlotMath.lerp(x2, x1, jj21);
					y2 = JPlotMath.lerp(y2, y1, jj21);
					v2 = JPlotMath.nlerp(v2, v1, jj21);
					x3 = JPlotMath.lerp(x6, x1, jj61);
					y3 = JPlotMath.lerp(y6, y1, jj61);
					v3 = JPlotMath.nlerp(v6, v1, jj61);
					x4 = Double.NaN;
					y4 = Double.NaN;
					v4 = Double.NaN;
					x5 = Double.NaN;
					y5 = Double.NaN;
					v5 = Double.NaN;
					x6 = Double.NaN;
					y6 = Double.NaN;
					v6 = Double.NaN;
					break;
				default:
					count = 6;
					break;
				}
			}
		}
		switch (count) {
		case 3:
			return new JDPolygon(new JDPoint(x1, y1, v1), new JDPoint(x2, y2, v2), new JDPoint(x3, y3, v3));
		case 4:
			return new JDPolygon(new JDPoint(x1, y1, v1), new JDPoint(x2, y2, v2), new JDPoint(x3, y3, v3),
					new JDPoint(x4, y4, v4));
		case 5:
			return new JDPolygon(new JDPoint(x1, y1, v1), new JDPoint(x2, y2, v2), new JDPoint(x3, y3, v3),
					new JDPoint(x4, y4, v4), new JDPoint(x5, y5, v5));
		case 6:
			return new JDPolygon(new JDPoint(x1, y1, v1), new JDPoint(x2, y2, v2), new JDPoint(x3, y3, v3),
					new JDPoint(x4, y4, v4), new JDPoint(x5, y5, v5), new JDPoint(x6, y6, v6));
		case 7:
			return new JDPolygon(new JDPoint(x1, y1, v1), new JDPoint(x2, y2, v2), new JDPoint(x3, y3, v3),
					new JDPoint(x4, y4, v4), new JDPoint(x5, y5, v5), new JDPoint(x6, y6, v6), new JDPoint(x7, y7, v7));
		default:
			return null;
		}
	}

	public double area() {
		return GeometryTools.area(getCorners());
	}
}
