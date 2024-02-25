package jplots.maths;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;

import jplots.helper.GeometryTools;

public class JDQuad {

	public int idx;
	public int lev;
	public double[] x;
	public double[] y;
	public double[] value;
	public JDEdge ab;
	public JDEdge bc;
	public JDEdge ca;
	private boolean parallelogram, parallelU, parallelV;
	private double xmin,xmax, ymin,ymax, mx,my, ux,uy, vx,vy;
	public double[][] k;
	
	private Integer hash = null;

	public JDQuad(Coordinate[] c) {
		this(new JDPoint(c[0]), new JDPoint(c[1]), new JDPoint(c[2]), new JDPoint(c[3]), 0);
	}

	public JDQuad(JDPoint a, JDPoint b, JDPoint c, JDPoint d) {
		this(a, b, c, d, 0);
	}

	public JDQuad(JDPoint a, JDPoint b, JDPoint c, JDPoint d, int index) {
		this.idx = index;
		x = new double[] { a.x, b.x, c.x, d.x };
		y = new double[] { a.y, b.y, c.y, d.y };
		value = new double[] { a.value, b.value, c.value, d.value };
		calcBarycenterHelper();
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
		double t = x[3];
		x[3] = x[2];
		x[2] = x[1];
		x[1] = t;
		t = y[3];
		y[3] = y[2];
		y[2] = y[1];
		y[1] = t;
		t = value[3];
		value[3] = value[2];
		value[2] = value[1];
		value[1] = t;
		calcBarycenterHelper();
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
	public JDPoint getD() {
		return new JDPoint(x[3], y[3], value[3]);
	}

	public JDPoint[] getCorners() {
		return new JDPoint[] { getA(), getB(), getC(), getD() };
	}

	public JDQuad copy() {
		return new JDQuad(getA(), getB(), getC(), getD(), idx);
	}

	public JDPolygon toPolygon() {
		return new JDPolygon(getCorners(), idx);
	}
	
//	public JDQuad[] withoutNaNCorners() {
//		JDQuad sub1 = null;
//		JDQuad sub2 = null;
//		int nancode = (Double.isNaN(value[0])?1:0) | (Double.isNaN(value[1])?2:0) | (Double.isNaN(value[2])?4:0);
//		JDPoint ab = getA().fractionTowards(0.5d, getB());
//		JDPoint bc = getB().fractionTowards(0.5d, getC());
//		JDPoint ca = getC().fractionTowards(0.5d, getA());
//		switch(nancode) {
//			case 0: break;
//			case 1: ab.value = value[1]; ca.value = value[2];
//					sub1 = new JDQuad(getB(),getC(),ca); sub2 = new JDQuad(getB(),ca,ab);
//					break;
//			case 2: ab.value = value[0]; bc.value = value[2];
//					sub1 = new JDQuad(getC(),getA(),ab); sub2 = new JDQuad(getC(),ab,bc);
//					break;
//			case 3: bc.value = value[2]; ca.value = value[2];
//					sub1 = new JDQuad(getC(),ca,bc);
//					break;
//			case 4: bc.value = value[1]; ca.value = value[0];
//					sub1 = new JDQuad(getA(),getB(),bc); sub2 = new JDQuad(getA(),bc,ca);
//					break;
//			case 5: ab.value = value[1]; bc.value = value[1];
//					sub1 = new JDQuad(getB(),bc,ab);
//					break;
//			case 6: ab.value = value[0]; ca.value = value[0];
//					sub1 = new JDQuad(getA(),ab,ca);
//					break;
//			default:
//			case 7:	return new JDQuad[0];
//		}
//		if(sub1==null)
//			return new JDQuad[] {this};
//		
//		if(sub2==null)
//			return new JDQuad[] {sub1};
//		
//		return new JDQuad[] {sub1, sub2};
//	}

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
		JDPoint[] B = ((JDQuad) obj).getCorners();

		if (A[0].equals(B[0])) {
			return (A[1].equals(B[1]) && A[2].equals(B[2])) || (A[1].equals(B[2]) && A[2].equals(B[1]));
		} else if (A[0].equals(B[1])) {
			return (A[1].equals(B[0]) && A[2].equals(B[2])) || (A[1].equals(B[2]) && A[2].equals(B[0]));
		} else if (A[0].equals(B[2])) {
			return (A[1].equals(B[0]) && A[2].equals(B[1])) || (A[1].equals(B[1]) && A[2].equals(B[0]));
		}

		return false;
	}
	
	public JDQuad affine(double[][] tm) {
		for (int i = 0; i < 3; i++) {
			double xx = x[i] * tm[0][0] + y[i] * tm[0][1] + tm[0][2];
			double yy = x[i] * tm[1][0] + y[i] * tm[1][1] + tm[1][2];
			x[i] = xx;
			y[i] = yy;
		}
		if(tm[0][0]*tm[1][1]-tm[0][1]*tm[1][0]<0d)
			reverse_orientation();
		else
			calcBarycenterHelper();
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
	public JDPoint[] intersectsLine(JDPoint p1, JDPoint p2) {
		return intersectsLine(p1, p2, 0.00000001d);
	}
	public JDPoint[] intersectsLine(JDPoint p1, JDPoint p2, double tol) {
		JDPoint out1=null, out2=null;
		JDPoint pab = GeometryTools.segmentIntersection(getA(),getB(), p1,p2, tol);
		JDPoint pbc = GeometryTools.segmentIntersection(getB(),getC(), p1,p2, tol);
		JDPoint pcd = GeometryTools.segmentIntersection(getC(),getD(), p1,p2, tol);
		JDPoint pda = GeometryTools.segmentIntersection(getD(),getA(), p1,p2, tol);
		double rel1=Double.POSITIVE_INFINITY, rel2=Double.NEGATIVE_INFINITY;
		if(pab!=null) { double rab = (pab.x-p1.x)*(p2.x-p1.x) + (pab.y-p1.y)*(p2.y-p1.y); 
			if(rab<rel1) { rel1 = rab; out1 = pab; }
			if(rab>rel2) { rel2 = rab; out2 = pab; } }
		if(pbc!=null) { double rbc = (pbc.x-p1.x)*(p2.x-p1.x) + (pbc.y-p1.y)*(p2.y-p1.y); 
			if(rbc<rel1) { rel1 = rbc; out1 = pbc; }
			if(rbc>rel2) { rel2 = rbc; out2 = pbc; } }
		if(pcd!=null) { double rcd = (pcd.x-p1.x)*(p2.x-p1.x) + (pcd.y-p1.y)*(p2.y-p1.y); 
			if(rcd<rel1) { rel1 = rcd; out1 = pcd; }
			if(rcd>rel2) { rel2 = rcd; out2 = pcd; } }
		if(pda!=null) { double rda = (pda.x-p1.x)*(p2.x-p1.x) + (pda.y-p1.y)*(p2.y-p1.y); 
			if(rda<rel1) { rel1 = rda; out1 = pda; }
			if(rda>rel2) { rel2 = rda; out2 = pda; } }
		if(out1==null || out2==null)
			return null;
		if(out1.equals(out2))
			return null;
		return new JDPoint[] {out1,out2};
		
//		double[] w1 = barycentricCoords(p1);
//		double[] w2 = barycentricCoords(p2);
//		double[] zo = { 0d, 1d };
//		// w1->0, w2->1, find 0->t1, 1->
//		double[] t1 = JPlotMath.map(zo, w1[0], w2[0], 0d, 1d);
//		if(t1[0]>t1[1]) t1 = new double[] {t1[1], t1[0]};
//		if(Math.abs(w1[0]-w2[0])<tol) { if(w1[0]<-tol || w1[0]>1d+tol) return null;
//			t1[0] = Double.NEGATIVE_INFINITY; t1[1] = Double.POSITIVE_INFINITY; }
//		double[] t2 = JPlotMath.map(zo, w1[1], w2[1], 0d, 1d);
//		if(t2[0]>t2[1]) t2 = new double[] {t2[1], t2[0]};
//		if(Math.abs(w1[1]-w2[1])<tol) { if(w1[1]<-tol || w1[1]>1d+tol) return null;
//			t2[0] = Double.NEGATIVE_INFINITY; t2[1] = Double.POSITIVE_INFINITY; }
//		double[] t3 = JPlotMath.map(zo, w1[2], w2[2], 0d, 1d);
//		if(t3[0]>t3[1]) t3 = new double[] {t3[1], t3[0]};
//		if(Math.abs(w1[2]-w2[2])<tol) { if(w1[2]<-tol || w1[2]>1d+tol) return null;
//			t3[0] = Double.NEGATIVE_INFINITY; t3[1] = Double.POSITIVE_INFINITY; }
//		double ts = Math.max(Math.max(0d, t1[0]), Math.max(t2[0], t3[0]));
//		double te = Math.min(Math.min(1d, t1[1]), Math.min(t2[1], t3[1]));
//		if(te<ts) return null;
//		return new JDPoint[] {
//				p1.fractionTowards(ts, p2),
//				p1.fractionTowards(te, p2)
//		};
	}
	
	public boolean contains(JDPoint p) {
		return contains(p, 0.0001d);
	}
	public boolean contains(JDPoint p, double delta) {
		if(p==null) return false;
		if(p.x<xmin || p.x>xmax) return false;
		if(p.y<ymin || p.y>ymax) return false;
		double kk = 1d + delta;
		double[] w = barycentricCoords(p);
		return Math.abs(w[0])<=kk && Math.abs(w[1])<=kk;
	}
	public JDPoint closestPointInside(JDPoint p) {
		if(contains(p)) return p;
		JDPoint a = closestOnEdge(x[0],y[0],value[0], x[1],y[1],value[1], p.x,p.y);
		JDPoint b = closestOnEdge(x[1],y[1],value[1], x[2],y[2],value[2], p.x,p.y);
		JDPoint c = closestOnEdge(x[2],y[2],value[2], x[3],y[3],value[3], p.x,p.y);
		JDPoint d = closestOnEdge(x[3],y[3],value[3], x[0],y[0],value[0], p.x,p.y);
		double ad = p.dist2(a), bd = p.dist2(b), cd = p.dist2(c), dd = p.dist2(d);
		if(ad<bd && ad<cd && ad<dd) return a;
		if(bd<cd && bd<dd) return b;
		if(cd<dd) return c;
		return d;
	}
	public double[] barycentricCoords(JDPoint p) {
		return barycentricCoords(p, 0.0001d);
	}
	public double[] barycentricCoords(JDPoint p, double delta) {
		/* Quad:
		         A---B
		         |   |
		         D---C
		 */
		double u = Double.NaN, v = Double.NaN;
		if(p==null) return new double[] {u,v};
		if(parallelU) {
			v = (p.x-mx)*vx + (p.y-my)*vy;
			if(parallelV) {
				u = (p.x-mx)*ux + (p.y-my)*uy;
			} else {
				double xl = 0.5d*(x[0]*(1d-v)+x[3]*(1d+v));
				double yl = 0.5d*(y[0]*(1d-v)+y[3]*(1d+v));
				double xr = 0.5d*(x[1]*(1d-v)+x[2]*(1d+v));
				double yr = 0.5d*(y[1]*(1d-v)+y[2]*(1d+v));
				double am = (p.x-xl)*(xr-xl) + (p.y-yl)*(yr-yl);
				double rl = (xr-xl)*(xr-xl) + (yr-yl)*(yr-yl);
				u = 2*am/rl - 1d;
			}
		} else
		if(parallelV) {
			u = (p.x-mx)*ux + (p.y-my)*uy;
			double xt = 0.5d*(x[0]*(1d-u)+x[1]*(1d+u));
			double yt = 0.5d*(y[0]*(1d-u)+y[1]*(1d+u));
			double xb = 0.5d*(x[3]*(1d-u)+x[2]*(1d+u));
			double yb = 0.5d*(y[3]*(1d-u)+y[2]*(1d+u));
			double am = (p.x-xt)*(xb-xt) + (p.y-yt)*(yb-yt);
			double rl =  (xb-xt)*(xb-xt) +  (yb-yt)*(yb-yt);
			v = 2*am/rl - 1d;
		}
		else {
			double q0 = (2*p.x+k[0][0])*k[1][2]                 - (2*p.y+k[1][0])*k[0][2];
			double q1 = (2*p.x+k[0][0])*k[1][3]+k[0][1]*k[1][2] - (2*p.y+k[1][0])*k[0][3]-k[1][1]*k[0][2];
			double q2 =                         k[0][1]*k[1][3] -                         k[1][1]*k[0][3];
			v = -q0 / q1;
			if(Math.abs(q2)>1e-32d) {
				double qs = Math.sqrt(q1*q1 - 4d*q2*q0);
				v = 0.5d*(-q1+qs)/q2;
				if(Math.abs(v)>1d) v = 0.5d*(-q1-qs)/q2;
			}
			u = 2 * (2*p.x+k[0][0]+k[0][1]*v) / (k[0][2]+k[0][3]*v) - 1d;
		}
		return new double[] {u,v};
	}
	public double valueAt(JDPoint p) {
		double[] w = barycentricCoords(p);
		double u = w[0];
		double v = w[1];
		double vsum = 0d;
		double wsum = 0d;
		int nancode = (Double.isNaN(value[0])?1:0) | (Double.isNaN(value[1])?2:0) | (Double.isNaN(value[2])?4:0) | (Double.isNaN(value[3])?8:0);
		switch(nancode) {
			default:
			case  0: break;
			case  1: if(-u-v>1d) return Double.NaN; break;
			case  2: if( u-v>1d) return Double.NaN; break;
			case  3: if(   v<0d) return Double.NaN; break;
			case  4: if( u+v>1d) return Double.NaN; break;
			case  5: if(-u-v>1d || u+v>1d) return Double.NaN; break;
			case  6: if( u  >0d) return Double.NaN; break;
			case  7: if(-u+v<1d) return Double.NaN; break;
			case  8: if(-u+v>1d) return Double.NaN; break;
			case  9: if( u  <0d) return Double.NaN; break;
			case 10: if(-u+v>1d || u-v>1d) return Double.NaN; break;
			case 11: if( u+v<1d) return Double.NaN; break;
			case 12: if(   v>0d) return Double.NaN; break;
			case 13: if( u-v<1d) return Double.NaN; break;
			case 14: if(-u-v<1d) return Double.NaN; break;
			case 15: return Double.NaN;
		}
		if(!Double.isNaN(value[0])) { vsum += value[0]*(1-u)*(1-v); wsum += (1-u)*(1-v); }
		if(!Double.isNaN(value[1])) { vsum += value[1]*(1+u)*(1-v); wsum += (1+u)*(1-v); }
		if(!Double.isNaN(value[2])) { vsum += value[2]*(1+u)*(1+v); wsum += (1+u)*(1+v); }
		if(!Double.isNaN(value[3])) { vsum += value[3]*(1-u)*(1+v); wsum += (1-u)*(1+v); }
		if(wsum<0.00000001d) return Double.NaN;
		return vsum/wsum;
	}
	public double[] gradientAt(JDPoint p) {
		double[] w = barycentricCoords(p);
		double u = w[0];
		double v = w[1];
		double vsum = 0d, wsum = 0d;
		double gu = 0d, gv = 0d, wu = 0d, wv = 0d;
		if(!Double.isNaN(value[0])) { vsum += value[0]*(1-u)*(1-v); wsum += (1-u)*(1-v); gu -= value[0]*(1-v); wu -= 1-v; gv -= value[0]*(1-u); wv -= 1-u; }
		if(!Double.isNaN(value[1])) { vsum += value[1]*(1+u)*(1-v); wsum += (1+u)*(1-v); gu += value[1]*(1-v); wu += 1-v; gv -= value[1]*(1+u); wv -= 1+u; }
		if(!Double.isNaN(value[2])) { vsum += value[2]*(1+u)*(1+v); wsum += (1+u)*(1+v); gu += value[2]*(1+v); wu += 1+v; gv += value[2]*(1+u); wv += 1+u; }
		if(!Double.isNaN(value[3])) { vsum += value[3]*(1-u)*(1+v); wsum += (1-u)*(1+v); gu -= value[3]*(1+v); wu -= 1+v; gv += value[3]*(1-u); wv += 1-u; }
		if(wsum<0.00000001d) return new double[] {Double.NaN, Double.NaN};
		gu = (gu*wsum-vsum*wu) / (wsum*wsum);
		gv = (gv*wsum-vsum*wv) / (wsum*wsum);
		return new double[] {gu, gv};
	}
	public double[] refineUV(double uf, double vf, double lev) {
		int nancode = (Double.isNaN(value[0])?1:0) | (Double.isNaN(value[1])?2:0) | (Double.isNaN(value[2])?4:0) | (Double.isNaN(value[3])?8:0);
		double u = uf, v = vf, o = 0d;
		if(nancode == 15) return new double[] {u,v};
		boolean stickToNaNBorder = false;
		switch(nancode) {
			default: break;
			case  1: stickToNaNBorder = (Math.abs(-u-v-1d)<0.0001d); break;
			case  2: stickToNaNBorder = (Math.abs( u-v-1d)<0.0001d); break;
			case  4: stickToNaNBorder = (Math.abs( u+v-1d)<0.0001d); break;
			case  8: stickToNaNBorder = (Math.abs(-u+v-1d)<0.0001d); break;
			case  5: stickToNaNBorder = (Math.abs(-u-v-1d)<0.0001d || Math.abs( u+v-1d)<0.0001d); break;
			case 10: stickToNaNBorder = (Math.abs( u-v-1d)<0.0001d || Math.abs(-u+v-1d)<0.0001d); break;
			case  3: stickToNaNBorder = (Math.abs(v)<0.0001d); break;
			case  6: stickToNaNBorder = (Math.abs(u)<0.0001d); break;
			case  9: stickToNaNBorder = (Math.abs(u)<0.0001d); break;
			case 12: stickToNaNBorder = (Math.abs(v)<0.0001d); break;
		}
		boolean stickToUBorder = Math.abs(u)>0.9999d;
		boolean stickToVBorder = Math.abs(v)>0.9999d;
		
		int ii = nancode==0 ? 4 : 10;
		for(int i=0; i<ii; i++) {
			double vsum = 0d, wsum = 0d;
			double gu = 0d, gv = 0d, wu = 0d, wv = 0d;
			if(!Double.isNaN(value[0])) { vsum += value[0]*(1-u)*(1-v); wsum += (1-u)*(1-v); gu -= value[0]*(1-v); wu -= 1-v; gv -= value[0]*(1-u); wv -= 1-u; }
			if(!Double.isNaN(value[1])) { vsum += value[1]*(1+u)*(1-v); wsum += (1+u)*(1-v); gu += value[1]*(1-v); wu += 1-v; gv -= value[1]*(1+u); wv -= 1+u; }
			if(!Double.isNaN(value[2])) { vsum += value[2]*(1+u)*(1+v); wsum += (1+u)*(1+v); gu += value[2]*(1+v); wu += 1+v; gv += value[2]*(1+u); wv += 1+u; }
			if(!Double.isNaN(value[3])) { vsum += value[3]*(1-u)*(1+v); wsum += (1-u)*(1+v); gu -= value[3]*(1+v); wu -= 1+v; gv += value[3]*(1-u); wv += 1-u; }
			if(wsum<0.00000001d) return new double[] {Double.NaN, Double.NaN};
			gu = (gu*wsum-vsum*wu) / (wsum*wsum);
			gv = (gv*wsum-vsum*wv) / (wsum*wsum);
			vsum /= wsum;
			double gr = 0.9d/(gu*gu+gv*gv); //reduce overshooting
			if(gr>1000000d) break;
			u += (lev-vsum)*gu*gr;
			v += (lev-vsum)*gv*gr;
			switch(nancode) {
				default:
				case  0: break;
				case  1: if(-u-v>1d || stickToNaNBorder) { o=0.5d*(-u-v-1d); u+=o; v+=o; } break;
				case  2: if( u-v>1d || stickToNaNBorder) { o=0.5d*( u-v-1d); u-=o; v+=o; } break;
				case  3: if(   v<0d || stickToNaNBorder) { v=0.0d; } break;
				case  4: if( u+v>1d || stickToNaNBorder) { o=0.5d*( u+v-1d); u-=o; v-=o; } break;
				case  5: if(-u-v>1d || (stickToNaNBorder && -u-v>0d)) { o=0.5d*(-u-v-1d); u+=o; v+=o; }
					     if( u+v>1d || (stickToNaNBorder &&  u+v>0d)) { o=0.5d*( u+v-1d); u-=o; v-=o; } break;
				case  6: if( u  >0d || stickToNaNBorder) { u=0.0d; } break;
				case  8: if(-u+v>1d || stickToNaNBorder) { o=0.5d*(-u+v-1d); u+=o; v-=o; } break;
				case  9: if( u  <0d || stickToNaNBorder) { u=0.0d; } break;
				case 10: if( u-v>1d || (stickToNaNBorder &&  u-v>0d)) { o=0.5d*( u-v-1d); u-=o; v+=o; }
					     if(-u+v>1d || (stickToNaNBorder &&  v-u>0d)) { o=0.5d*(-u+v-1d); u+=o; v-=o; } break;
				case 12: if(   v>0d || stickToNaNBorder) { v=0.0d; } break;
			}
			if(u<-1d || (stickToUBorder && u<0d)) u = -1d;
			if(u>1d  || (stickToUBorder && u>0d)) u = 1d;
			if(v<-1d || (stickToVBorder && v<0d)) v = -1d;
			if(v>1d  || (stickToVBorder && v>0d)) v = 1d;
		}
		return new double[] {u,v};
	}
	public JDPoint pointFromUV(double[] uv) {
		return pointFromUV(uv[0], uv[1]);
	}
	public JDPoint pointFromUV(double u, double v) {
		return new JDPoint(
				0.25d*(x[0]*(1d-u)*(1d-v)     + x[1]*(1d+u)*(1d-v)     + x[3]*(1d-u)*(1d+v)     + x[2]*(1d+u)*(1d+v)),
				0.25d*(y[0]*(1d-u)*(1d-v)     + y[1]*(1d+u)*(1d-v)     + y[3]*(1d-u)*(1d+v)     + y[2]*(1d+u)*(1d+v)),
				0.25d*(value[0]*(1d-u)*(1d-v) + value[1]*(1d+u)*(1d-v) + value[3]*(1d-u)*(1d+v) + value[2]*(1d+u)*(1d+v))
		);
	}
	public JDPolygon[] getLevelRangePolygons(double low, double upp) {
		int nancode = (Double.isNaN(value[0])?1:0) | (Double.isNaN(value[1])?2:0) | (Double.isNaN(value[2])?4:0) | (Double.isNaN(value[3])?8:0);
		int nancnt  = (Double.isNaN(value[0])?1:0) + (Double.isNaN(value[1])?1:0) + (Double.isNaN(value[2])?1:0) + (Double.isNaN(value[3])?1:0);
		int[] order = { 0, 1, 2, 3 };
		switch(nancode) {
			default:
			case  0:
			case  8:
			case 10:
			case 12:
			case 14: order[0]=0; order[1]=1; order[2]=2; order[3]=3; break;
			case  1:
			case  5:
			case  9:
			case 13: rotateIndecesBy(1); break;
			case  2:
			case  3:
			case 11: rotateIndecesBy(2); break;
			case  4:
			case  6:
			case  7: rotateIndecesBy(3); break;
			case 15: return null;
		}
		switch(nancnt) {
			default:
			case 4: return null;
			case 3: if(value[0]<low || value[0]>upp) return null;
					return new JDPolygon[] {new JDPolygon(
							new JDPoint(x[0],y[0]),
							new JDPoint(0.5d*(x[3]+x[0]),0.5d*(y[3]+y[0])),
							new JDPoint(0.5d*(x[0]+x[1]),0.5d*(y[0]+y[1])))};
			case 2: if(nancode==5 || nancode==10) {
					int levcode2 = (value[0]<low?0x10:value[0]>upp?0x30:0x20) | (value[2]<low?0x01:value[2]>upp?0x03:0x02);
					if(levcode2==0x11 || levcode2==0x33) return null;
					JDPoint a2 = new JDPoint(x[0],y[0]);
					JDPoint b2 = new JDPoint(x[2],y[2]);
					double low12 = (low-value[0])/(value[2]-value[0]);
					double upp12 = (upp-value[0])/(value[2]-value[0]);
					List<JDPoint> pnts2 = new ArrayList<JDPoint>();
					switch(levcode2) {
						default: break;
						case 0x12: addCurve(pnts2, low12, low12-1d, low12-1d, low12, low, false); pnts2.add(b2); break;
						case 0x13: addCurve(pnts2, low12, low12-1d, low12-1d, low12, low, false); addCurve(pnts2, upp12-1d, upp12, upp12, upp12-1d, upp, false); break;
						case 0x21: addCurve(pnts2, low12-1d, low12, low12, low12-1d, low, false); pnts2.add(a2); break;
						case 0x22: pnts2.add(a2); pnts2.add(new JDPoint(0.5d*(x[3]+x[0]),0.5d*(y[3]+y[0]))); pnts2.add(new JDPoint(0.5d*(x[2]+x[3]),0.5d*(y[2]+y[3])));
							pnts2.add(b2); pnts2.add(new JDPoint(0.5d*(x[1]+x[2]),0.5d*(y[1]+y[2]))); pnts2.add(new JDPoint(0.5d*(x[0]+x[1]),0.5d*(y[0]+y[1]))); break;
						case 0x23: addCurve(pnts2, upp12-1d, upp12, upp12, upp12-1d, upp, false); pnts2.add(a2); break;
						case 0x31: addCurve(pnts2, upp12, upp12-1d, upp12-1d, upp12, upp, false); addCurve(pnts2, low12-1d, low12, low12, low12-1d, low, false); break;
						case 0x32: addCurve(pnts2, upp12, upp12-1d, upp12-1d, upp12, upp, false); pnts2.add(b2); break;
					}
					return new JDPolygon[] { new JDPolygon(pnts2.toArray(new JDPoint[0])) };
				}else {
					int levcode2 = (value[0]<low?0x10:value[0]>upp?0x30:0x20) | (value[1]<low?0x01:value[1]>upp?0x03:0x02);
					JDPoint a2 = new JDPoint(x[0],y[0]);
					JDPoint b2 = new JDPoint(x[1],y[1]);
					JDPoint c2 = new JDPoint(0.5d*(x[1]+x[2]),0.5d*(y[1]+y[2]));
					JDPoint d2 = new JDPoint(0.5d*(x[0]+x[3]),0.5d*(y[0]+y[3]));
					double low12 = (low-value[0])/(value[1]-value[0]);
					JDPoint ab2low = new JDPoint(a2.x*(1-low12)+b2.x*low12, a2.y*(1-low12)+b2.y*low12);
					JDPoint dc2low = new JDPoint(d2.x*(1-low12)+c2.x*low12, d2.y*(1-low12)+c2.y*low12);
					double upp12 = (upp-value[0])/(value[1]-value[0]);
					JDPoint ab2upp = new JDPoint(a2.x*(1-upp12)+b2.x*upp12, a2.y*(1-upp12)+b2.y*upp12);
					JDPoint dc2upp = new JDPoint(d2.x*(1-upp12)+c2.x*upp12, d2.y*(1-upp12)+c2.y*upp12);
					switch(levcode2) {
						default: return null;
						case 0x11: return null;
						case 0x12: return new JDPolygon[] {new JDPolygon(b2,ab2low,dc2low,c2)};
						case 0x13: return new JDPolygon[] {new JDPolygon(ab2upp,ab2low,dc2low,dc2upp)};
						case 0x21: return new JDPolygon[] {new JDPolygon(ab2low,a2,d2,dc2low)};
						case 0x22: return new JDPolygon[] {new JDPolygon(d2,c2,b2,a2)};
						case 0x23: return new JDPolygon[] {new JDPolygon(d2,dc2upp,ab2upp,a2)};
						case 0x31: return new JDPolygon[] {new JDPolygon(ab2low,ab2upp,dc2upp,dc2low)};
						case 0x32: return new JDPolygon[] {new JDPolygon(dc2upp,c2,b2,ab2upp)};
						case 0x33: return null;
					}
				}
			case 1: int levcode1 = (value[0]<low?0x100:value[0]>upp?0x300:0x200) | (value[1]<low?0x010:value[1]>upp?0x030:0x020) | (value[2]<low?0x001:value[2]>upp?0x003:0x002);
				if(levcode1==0x111 || levcode1==0x333) return null;
				JDPoint a1 = new JDPoint(x[0],y[0]), b1 = new JDPoint(x[1],y[1]), c1 = new JDPoint(x[2],y[2]);
				JDPoint ad1 = new JDPoint(0.5d*(x[0]+x[3]),0.5d*(y[0]+y[3])), cd1 = new JDPoint(0.5d*(x[2]+x[3]),0.5d*(y[2]+y[3]));
				double ab1low = 2*(low-value[0])/(value[1]-value[0])-1d, ab1upp = 2*(upp-value[0])/(value[1]-value[0])-1d;
				double bc1low = 2*(low-value[1])/(value[2]-value[1])-1d, bc1upp = 2*(upp-value[1])/(value[2]-value[1])-1d;
				double ac1low = (low-value[0])/(value[2]-value[0]), ac1upp = (upp-value[0])/(value[2]-value[0]);
				List<JDPoint> pnts1 = new ArrayList<JDPoint>();
				switch(levcode1) {
					default: break;
					case 0x112: addCurve(pnts1, 1d, bc1low-1d, ac1low-1d, ac1low, low, false); pnts1.add(cd1); pnts1.add(c1); break;
					case 0x113: addCurve(pnts1, 1d, bc1low, ac1low-1d, ac1low, low, false); addCurve(pnts1, ac1upp-1d, ac1upp, 1d, bc1upp, upp, false); break;
					case 0x121: addCurve(pnts1, ab1low, -1d, 1d, bc1low, low, false); pnts1.add(b1); break;
					case 0x122: addCurve(pnts1, ab1low, -1d, ac1low-1d, ac1low, low, false); pnts1.add(cd1); pnts1.add(c1); pnts1.add(b1); break;
					case 0x123: addCurve(pnts1, ab1low, -1d, ac1low-1d, ac1low, low, false); addCurve(pnts1, ac1upp-1d, ac1upp, 1d, bc1upp, upp, false); pnts1.add(b1); break;
					case 0x131: addCurve(pnts1, ab1low, -1d, 1d, bc1low, low, false); addCurve(pnts1, 1d, bc1upp, ab1upp, -1d, upp, false); break;
					case 0x132: addCurve(pnts1, ab1low, -1d, ac1low-1d, ac1low, low, false); pnts1.add(cd1); pnts1.add(c1); addCurve(pnts1, 1d, bc1upp, ab1upp, -1d, upp, false); break;
					case 0x133: addCurve(pnts1, ab1low, -1d, ac1low-1d, ac1low, low, false); addCurve(pnts1, ac1upp-1d, ac1upp, ab1upp, -1d, upp, false); break;
					case 0x211: addCurve(pnts1, ac1low-1d, ac1low, ab1low, -1d, low, false); pnts1.add(a1); pnts1.add(ad1); break;
					case 0x212: addCurve(pnts1, 1d, bc1low, ab1low, -1d, low, false); pnts1.add(a1); pnts1.add(ad1); pnts1.add(cd1); pnts1.add(c1); break;
					case 0x213: addCurve(pnts1, 1d, bc1low, ab1low, -1d, low, false); pnts1.add(a1); pnts1.add(ad1); addCurve(pnts1, ac1upp-1d, ac1upp, 1d, bc1upp, upp, false); break;
					case 0x221: addCurve(pnts1, ac1low-1d, ac1low, 1d, bc1low, low, false); pnts1.add(b1); pnts1.add(a1); pnts1.add(ad1); break;
					case 0x222: pnts1.add(ad1); pnts1.add(cd1); pnts1.add(c1); pnts1.add(b1); pnts1.add(a1); break;
					case 0x223: addCurve(pnts1, ac1upp-1d, ac1upp, 1d, bc1upp, upp, false); pnts1.add(b1); pnts1.add(a1); pnts1.add(ad1); break;
					case 0x231: addCurve(pnts1, ac1low-1d, ac1low, 1d, bc1low, low, false); addCurve(pnts1, 1d, bc1upp, ab1upp, -1d, upp, false); pnts1.add(a1); pnts1.add(ad1); break;
					case 0x232: addCurve(pnts1, 1d, bc1upp, ab1upp, -1d, upp, false); pnts1.add(a1); pnts1.add(ad1); pnts1.add(cd1); pnts1.add(c1); break;
					case 0x233: addCurve(pnts1, ac1upp-1d, ac1upp, ab1upp, -1d, upp, false); pnts1.add(a1); pnts1.add(ad1); break;
					case 0x311: addCurve(pnts1, ac1low-1d, ac1low, ab1low, -1d, low, false); addCurve(pnts1, ab1upp, -1d, ac1upp-1d, ac1upp, upp, false); break;
					case 0x312: addCurve(pnts1, 1d, bc1low, ab1low, -1d, low, false); addCurve(pnts1, ab1upp, -1d, ac1upp-1d, ac1upp, upp, false); pnts1.add(cd1); pnts1.add(c1); break;
					case 0x313: addCurve(pnts1, 1d, bc1low, ab1low, -1d, low, false); addCurve(pnts1, ab1upp, -1d, 1d, bc1upp, upp, false); break;
					case 0x321: addCurve(pnts1, ac1low-1d, ac1low, 1d, bc1low, low, false); pnts1.add(b1); addCurve(pnts1, ab1upp, -1d, ac1upp-1d, ac1upp, upp, false); break;
					case 0x322: addCurve(pnts1, ab1upp, -1d, ac1upp-1d, ac1upp, upp, false); pnts1.add(cd1); pnts1.add(c1); pnts1.add(b1); break;
					case 0x323: addCurve(pnts1, ab1upp, -1d, 1d, bc1upp, upp, false); pnts1.add(b1); break;
					case 0x331: addCurve(pnts1, ac1low-1d, ac1low, 1d, bc1low, low, false); addCurve(pnts1, 1d, bc1upp, ac1upp-1d, ac1upp, upp, false); break;
					case 0x332: addCurve(pnts1, 1d, bc1upp, ac1upp-1d, ac1upp, upp, false); pnts1.add(cd1); pnts1.add(c1); break;
				}
				return new JDPolygon[] { new JDPolygon(pnts1.toArray(new JDPoint[0])) };
			case 0: int levcode0 = (value[0]<low?0x1000:value[0]>upp?0x3000:0x2000) |
					               (value[1]<low?0x0100:value[1]>upp?0x0300:0x0200) |
					               (value[2]<low?0x0010:value[2]>upp?0x0030:0x0020) |
					               (value[3]<low?0x0001:value[3]>upp?0x0003:0x0002);
				if(levcode0==0x1111 || levcode0==0x3333)
					return null;
				List<JDPoint> pntsA = new ArrayList<JDPoint>();
				List<JDPoint> pntsB = new ArrayList<JDPoint>();
				double	ab0low = 2*(low-value[0])/(value[1]-value[0])-1d, bc0low = 2*(low-value[1])/(value[2]-value[1])-1d,
						dc0low = 2*(low-value[3])/(value[2]-value[3])-1d, ad0low = 2*(low-value[0])/(value[3]-value[0])-1d;
				double	ab0upp = 2*(upp-value[0])/(value[1]-value[0])-1d, bc0upp = 2*(upp-value[1])/(value[2]-value[1])-1d,
						dc0upp = 2*(upp-value[3])/(value[2]-value[3])-1d, ad0upp = 2*(upp-value[0])/(value[3]-value[0])-1d;
				JDPoint a0 = new JDPoint(x[0],y[0]), b0 = new JDPoint(x[1],y[1]), c0 = new JDPoint(x[2],y[2]), d0 = new JDPoint(x[3],y[3]);
				double us = (value[3]-value[0]+value[2]-value[1])/(value[3]-value[0]-value[2]+value[1]);
				double vs = (value[1]-value[0]+value[2]-value[3])/(value[1]-value[0]-value[2]+value[3]);
				double ms = 0.25d*(value[0]*(1-us)*(1-vs)+value[1]*(1+us)*(1-vs)+value[2]*(1+us)*(1+vs)+value[3]*(1-us)*(1+vs));
				switch(levcode0) {
					case 0x1112: addCurve(pntsA, dc0low, 1d, -1d, ad0low, low, false); pntsA.add(d0); break;
					case 0x1113: addCurve(pntsA, dc0low, 1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false); break;
					case 0x1121: addCurve(pntsA, 1d, bc0low, dc0low, 1d, low, false); pntsA.add(c0); break;
					case 0x1122: addCurve(pntsA, 1d, bc0low, -1d, ad0low, low, false); pntsA.add(d0); pntsA.add(c0); break;
					case 0x1123: addCurve(pntsA, 1d, bc0low, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false); pntsA.add(c0); break;
					case 0x1131: addCurve(pntsA, 1d, bc0low, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false);  break;
					case 0x1132: addCurve(pntsA, 1d, bc0low, -1d, ad0low, low, false); pntsA.add(d0); addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false); break;
					case 0x1133: addCurve(pntsA, 1d, bc0low, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, 1d, bc0upp, upp, false); break;
					case 0x1211: addCurve(pntsA, ab0low, -1d, 1d, bc0low, low, false); pntsA.add(b0); break;
					case 0x1212: if(ms<low) {
							addCurve(pntsA, dc0low, 1d, -1d, ad0low, low, false); pntsA.add(d0);
							addCurve(pntsB, ab0low, -1d, 1d, bc0low, low, false); pntsB.add(b0);
						} else {
							addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); pntsA.add(d0);
							addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); pntsA.add(b0);
						} break;
					case 0x1213: if(ms<low) {
							addCurve(pntsA, dc0low, 1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false);
							addCurve(pntsB, ab0low, -1d, 1d, bc0low, low, false); pntsB.add(b0);
						} else {
							addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false);
							addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); pntsA.add(b0);
						} break;
					case 0x1221: addCurve(pntsA, ab0low, -1d, dc0low, 1d, low, false); pntsA.add(c0); pntsA.add(b0); break;
					case 0x1222: addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); pntsA.add(d0); pntsA.add(c0); pntsA.add(b0); break;
					case 0x1223: addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false); pntsA.add(c0); pntsA.add(b0); break;
					case 0x1231: addCurve(pntsA, ab0low, -1d, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false); pntsA.add(b0); break;
					case 0x1232: addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); pntsA.add(d0); addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false); pntsA.add(b0); break;
					case 0x1233: addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, 1d, bc0upp, upp, false); pntsA.add(b0); break;
					case 0x1311: addCurve(pntsA, ab0low, -1d, 1d, bc0low, low, false); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false); break;
					case 0x1312: if(ms<low) {
							addCurve(pntsA, ab0low, -1d, 1d, bc0low, low, false); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false);
							addCurve(pntsB, dc0low, 1d, -1d, ad0low, low, false); pntsB.add(d0);
						} else {
							addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); pntsA.add(d0);
							addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false);
						} break;
					case 0x1313: if(ms<low) {
							addCurve(pntsA, dc0low, 1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false);
							addCurve(pntsB, ab0low, -1d, 1d, bc0low, low, false); addCurve(pntsB, 1d, bc0upp, ab0upp, -1d, upp, false);
						} else if(ms>upp) {
							addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, ab0upp, -1d, upp, false);
							addCurve(pntsB, dc0low, 1d, 1d, bc0low, low, false); addCurve(pntsB, 1d, bc0upp, dc0upp, 1d, upp, false);
						} else {
							addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false);
							addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false);
						} break;
					case 0x1321: addCurve(pntsA, ab0low, -1d, dc0low, 1d, low, false); pntsA.add(c0); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false); break;
					case 0x1322: addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); pntsA.add(d0); pntsA.add(c0); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false); break;
					case 0x1323: if(ms>upp) {
							addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, ab0upp, -1d, upp, false);
							addCurve(pntsB, 1d, bc0upp, dc0upp, 1d, upp, false); pntsB.add(c0);
						} else {
							addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false);
							pntsA.add(c0); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false);
						} break;
					case 0x1331: addCurve(pntsA, ab0low, -1d, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, ab0upp, -1d, upp, false); break;
					case 0x1332: addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); pntsA.add(d0); addCurve(pntsA, dc0upp, 1d, ab0upp, -1d, upp, false); break;
					case 0x1333: addCurve(pntsA, ab0low, -1d, -1d, ad0low, low, false); addCurve(pntsA, -1d, ad0upp, ab0upp, -1d, upp, false); break;
					case 0x2111: addCurve(pntsA, -1d, ad0low, ab0low, -1d, low, false); pntsA.add(a0); break;
					case 0x2112: addCurve(pntsA, dc0low, 1d, ab0low, -1d, low, false); pntsA.add(a0); pntsA.add(d0); break;
					case 0x2113: addCurve(pntsA, dc0low, 1d, ab0low, -1d, low, false); pntsA.add(a0); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false); break;
					case 0x2121: if(ms<low) {
							addCurve(pntsA, -1d, ad0low, ab0low, -1d, low, false); pntsA.add(a0);
							addCurve(pntsB, 1d, bc0low, dc0low, 1d, low, false); pntsB.add(c0);
						} else {
							addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); pntsA.add(c0);
							addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); pntsA.add(a0);
						} break;
					case 0x2122: addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); pntsA.add(a0); pntsA.add(d0); pntsA.add(c0); break;
					case 0x2123: addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); pntsA.add(a0); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false); pntsA.add(c0); break;
					case 0x2131: if(ms<low) {
							addCurve(pntsA, -1d, ad0low, ab0low, -1d, low, false); pntsA.add(a0);
							addCurve(pntsB, 1d, bc0low, dc0low, 1d, low, false); addCurve(pntsB, dc0upp, 1d, 1d, bc0upp, upp, false);
						} else {
							addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false);
							addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); pntsA.add(a0);
						} break;
					case 0x2132: addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); pntsA.add(a0); pntsA.add(d0); addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false); break;
					case 0x2133: addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); pntsA.add(a0); addCurve(pntsA, -1d, ad0upp, 1d, bc0upp, upp, false); break;
					case 0x2211: addCurve(pntsA, -1d, ad0low, 1d, bc0low, low, false); pntsA.add(b0); pntsA.add(a0); break;
					case 0x2212: addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); pntsA.add(b0); pntsA.add(a0); pntsA.add(d0); break;
					case 0x2213: addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); pntsA.add(b0); pntsA.add(a0); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false); break;
					case 0x2221: addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); pntsA.add(c0); pntsA.add(b0); pntsA.add(a0); break;
					case 0x2222: pntsA.add(d0); pntsA.add(c0); pntsA.add(b0); pntsA.add(a0); break;
					case 0x2223: addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false); pntsA.add(c0); pntsA.add(b0); pntsA.add(a0); break;
					case 0x2231: addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false); pntsA.add(b0); pntsA.add(a0); break;
					case 0x2232: addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false); pntsA.add(b0); pntsA.add(a0); pntsA.add(d0); break;
					case 0x2233: addCurve(pntsA, -1d, ad0upp, 1d, bc0upp, upp, false); pntsA.add(b0); pntsA.add(a0); break;
					case 0x2311: addCurve(pntsA, -1d, ad0low, 1d, bc0low, low, false); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false); pntsA.add(a0); break;
					case 0x2312: addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false); pntsA.add(a0); pntsA.add(d0); break;
					case 0x2313: if(ms>upp) {
							addCurve(pntsA, -1d, ad0upp, ab0upp, -1d, upp, false); pntsA.add(a0);
							addCurve(pntsB, dc0low, 1d, 1d, bc0low, low, false); addCurve(pntsB, 1d, bc0upp, dc0upp, 1d, upp, false);
						} else {
							pntsA.add(a0); addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false);
							addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false);
						} break;
					case 0x2321: addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); pntsA.add(c0); addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false); pntsA.add(a0); break;
					case 0x2322: addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false); pntsA.add(a0); pntsA.add(d0); pntsA.add(c0); break;
					case 0x2323: if(ms>upp) {
							addCurve(pntsA, -1d, ad0upp, ab0upp, -1d, upp, false); pntsA.add(a0);
							addCurve(pntsB, 1d, bc0upp, dc0upp, 1d, upp, false); pntsB.add(c0);
						} else {
							addCurve(pntsA, -1d, ad0upp, dc0upp, 1d, upp, false); pntsA.add(c0);
							addCurve(pntsA, 1d, bc0upp, ab0upp, -1d, upp, false); pntsA.add(a0);
						} break;
					case 0x2331: addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, ab0upp, -1d, upp, false); pntsA.add(a0); break;
					case 0x2332: addCurve(pntsA, dc0upp, 1d, ab0upp, -1d, upp, false); pntsA.add(a0); pntsA.add(d0); break;
					case 0x2333: addCurve(pntsA, -1d, ad0upp, ab0upp, -1d, upp, false); pntsA.add(a0); break;
					case 0x3111: addCurve(pntsA, -1d, ad0low, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false); break;
					case 0x3112: addCurve(pntsA, dc0low, 1d, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false); pntsA.add(d0); break;
					case 0x3113: addCurve(pntsA, dc0low, 1d, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, dc0upp, 1d, upp, false); break;
					case 0x3121: if(ms<low) {
							addCurve(pntsA, -1d, ad0low, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false);
							addCurve(pntsB, 1d, bc0low, dc0low, 1d, low, false); pntsB.add(c0);
						} else {
							addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); pntsA.add(c0);
							addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false);
						} break;
					case 0x3122: addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false); pntsA.add(d0); pntsA.add(c0); break;
					case 0x3123: addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, dc0upp, 1d, upp, false); pntsA.add(c0); break;
					case 0x3131: if(ms<low) {
							addCurve(pntsA, -1d, ad0low, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false);
							addCurve(pntsB, 1d, bc0low, dc0low, 1d, low, false); addCurve(pntsB, dc0upp, 1d, 1d, bc0upp, upp, false);
						} else if(ms>upp) {
							addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, -1d, ad0upp, upp, false);
							addCurve(pntsB, 1d, bc0low, ab0low, -1d, low, false); addCurve(pntsB, ab0upp, -1d, 1d, bc0upp, upp, false);
						} else {
							addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false);
							addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false);
						} break;
					case 0x3132: if(ms>upp) {
							addCurve(pntsA, dc0upp, 1d, -1d, ad0upp, upp, false); pntsA.add(d0);
							addCurve(pntsB, 1d, bc0low, ab0low, -1d, low, false); addCurve(pntsB, ab0upp, -1d, 1d, bc0upp, upp, false);
						} else {
							pntsA.add(d0); addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false);
							addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false);
						} break;
					case 0x3133: addCurve(pntsA, 1d, bc0low, ab0low, -1d, low, false); addCurve(pntsA, ab0upp, -1d, 1d, bc0upp, upp, false); break;
					case 0x3211: addCurve(pntsA, -1d, ad0low, 1d, bc0low, low, false); pntsA.add(b0); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false); break;
					case 0x3212: addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); pntsA.add(b0); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false); pntsA.add(d0); break;
					case 0x3213: addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); pntsA.add(b0); addCurve(pntsA, ab0upp, -1d, dc0upp, 1d, upp, false); break;
					case 0x3221: addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); pntsA.add(c0); pntsA.add(b0); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false); break;
					case 0x3222: addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false); pntsA.add(d0); pntsA.add(c0); pntsA.add(b0); break;
					case 0x3223: addCurve(pntsA, ab0upp, -1d, dc0upp, 1d, upp, false); pntsA.add(c0); pntsA.add(b0); break;
					case 0x3231: if(ms>upp) {
							addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, -1d, ad0upp, upp, false);
							addCurve(pntsB, ab0upp, -1d, 1d, bc0upp, upp, false); pntsB.add(b0);
						} else {
							addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false);
							pntsA.add(b0); addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false);
						} break;
					case 0x3232: if(ms>upp) {
							addCurve(pntsA, dc0upp, 1d, -1d, ad0upp, upp, false); pntsA.add(d0);
							addCurve(pntsB, ab0upp, -1d, 1d, bc0upp, upp, false); pntsB.add(b0);
						} else {
							addCurve(pntsA, ab0upp, -1d, -1d, ad0upp, upp, false); pntsA.add(d0);
							addCurve(pntsA, dc0upp, 1d, 1d, bc0upp, upp, false); pntsA.add(b0);
						} break;
					case 0x3233: addCurve(pntsA, ab0upp, -1d, 1d, bc0upp, upp, false); pntsA.add(b0); break;
					case 0x3311: addCurve(pntsA, -1d, ad0low, 1d,  bc0low, low, false); addCurve(pntsA, 1d, bc0upp, -1d, ad0upp, upp, false); break;
					case 0x3312: addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); addCurve(pntsA, 1d, bc0upp, -1d, ad0upp, upp, false); pntsA.add(d0); break;
					case 0x3313: addCurve(pntsA, dc0low, 1d, 1d, bc0low, low, false); addCurve(pntsA, 1d, bc0upp, dc0upp, 1d, upp, false); break;
					case 0x3321: addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); pntsA.add(c0); addCurve(pntsA, 1d, bc0upp, -1d, ad0upp, upp, false); break;
					case 0x3322: addCurve(pntsA, 1d, bc0upp, -1d, ad0upp, upp, false); pntsA.add(d0); pntsA.add(c0); break;
					case 0x3323: addCurve(pntsA, 1d, bc0upp, dc0upp, 1d, upp, false); pntsA.add(c0); break;
					case 0x3331: addCurve(pntsA, -1d, ad0low, dc0low, 1d, low, false); addCurve(pntsA, dc0upp, 1d, -1d, ad0upp, upp, false);
					case 0x3332: addCurve(pntsA, dc0upp, 1d, -1d, ad0upp, upp, false); pntsA.add(d0); break;
					default: return null;
				}
				if(pntsB.size()==0) {
					if(pntsA.size()==0) return null;
					return new JDPolygon[] { new JDPolygon(pntsA.toArray(new JDPoint[0])) };
				}
				if(pntsA.size()==0)
					return new JDPolygon[] { new JDPolygon(pntsB.toArray(new JDPoint[0])) };
				return new JDPolygon[] {
					new JDPolygon(pntsA.toArray(new JDPoint[0])),
					new JDPolygon(pntsB.toArray(new JDPoint[0]))
				};
		}
	}
	public JDLine[] getLevelRangeLines(JDPoint s, JDPoint e, double low, double upp) {
		double dx = e.x-s.x, dy = e.y-s.y;
		double dr2 = 1d / (dx*dx+dy*dy);
		double dr = Math.sqrt(dr2);
		JDPoint[] lineEnds = intersectsLine(s, e);
		if(lineEnds==null) return null;
		if(lineEnds[0]==null || lineEnds[1]==null) return null;
		dx *= dr; dy *= dr;
		//now we have a linesegment, which goes through our quad
//		return new JDLine[] { new JDLine( lineEnds ) };
		JDPoint[] lp = new JDPoint[5]; //9?
		for(int i=0; i<lp.length; i++) {
			lp[i] = lineEnds[0].fractionTowards(0.25d*i, lineEnds[1]);
			lp[i].value = valueAt(lp[i]);
		}
		JDPoint[] ends = new JDPoint[18];
		int j=0, prev=Double.isNaN(lp[0].value)?0:lp[0].value<low?1:lp[0].value>upp?3:2;
		if(prev==2) { j = 1; ends[0] = lp[0]; }
		for(int i=1; i<lp.length; i++) {
			prev = (prev<<4) & 0xff;
			prev |= Double.isNaN(lp[i].value)?0:lp[i].value<low?1:lp[i].value>upp?3:2;
			double fl = (low-lp[i-1].value) / (lp[i].value-lp[i-1].value);
			double fu = (upp-lp[i-1].value) / (lp[i].value-lp[i-1].value);
			switch(prev) {
				default:
				case 0x11:
				case 0x22:
				case 0x33:
					break;
				case 0x02: ends[j] = lp[i]; j++; break;
				case 0x12: ends[j] = lp[i-1].fractionTowards(fl, lp[i]); j++; break;
				case 0x13: ends[j] = lp[i-1].fractionTowards(fl, lp[i]); ends[j+1] = lp[i-1].fractionTowards(fu, lp[i]); j+=2; break;
				case 0x20: if(j>0) { ends[j] = lp[i-1]; j++; } break;
				case 0x21: ends[j] = lp[i-1].fractionTowards(fl, lp[i]); j++; break;
				case 0x23: ends[j] = lp[i-1].fractionTowards(fu, lp[i]); j++; break;
				case 0x31: ends[j] = lp[i-1].fractionTowards(fu, lp[i]); ends[j+1] = lp[i-1].fractionTowards(fl, lp[i]); j+=2; break;
				case 0x32: ends[j] = lp[i-1].fractionTowards(fu, lp[i]); j++; break;
			}
		}
		if((prev&0xf) == 2) { ends[j] = lp[lp.length-1]; j++; }
		switch(j) {
			default:
				return null;
			case 2:
				return new JDLine[] { new JDLine(ends[0], ends[1]) };
			case 4:
				return new JDLine[] { new JDLine(ends[0], ends[1]), new JDLine(ends[2], ends[3]) };
			case 6:
				return new JDLine[] { new JDLine(ends[0], ends[1]), new JDLine(ends[2], ends[3]), new JDLine(ends[4], ends[5]) };
		}
		
////		System.out.println("Found intersection: "+lineEnds[0]+"("+vs+"), "+lineEnds[1]+"("+ve+")"); 
//		//estimate via parabola, if this linesegment is splitted or not:
//		// vs =             c
//		// vm = a/4 + b/2 + c       4(vm-c) =  a + 2b
//		// ve = a   + b   + c       2(ve-c) = 2a + 2b
//		double pa =  2*vs - 4*vm + 2*ve;
//		double pb = -3*vs + 4*vm -   ve;
//		double pc =    vs;
////		System.out.println("  --> estimated parabola:  y = "+pa+"*x² + "+pb+"*x + "+pc);
//		if(Math.abs(pa)<0.00000001d) {
//			if(Math.max(vs, ve)<low) return null;
//			if(Math.min(vs, ve)>upp) return null;
//			double vl = Math.max(0d, Math.min(1d, (low-vs)/(ve-vs)));
//			double vu = Math.max(0d, Math.min(1d, (upp-vs)/(ve-vs)));
////			System.out.println("  --> found segment {"+vl+" ... "+vu+"}");
//			if(vl<vu) return new JDLine[] { new JDLine( lineEnds[0].fractionTowards(vl, lineEnds[1]), lineEnds[0].fractionTowards(vu, lineEnds[1]) ) };
//			if(vl>vu) return new JDLine[] { new JDLine( lineEnds[0].fractionTowards(vu, lineEnds[1]), lineEnds[0].fractionTowards(vl, lineEnds[1]) ) };
//			return null;
//		}
//		//midnight-formula:  y = [-b+-sqrt(b²-4ac)]/(2a)
//		
//		double slow = Math.sqrt(pb*pb - 4d*pa*(pc-low));
//		double supp = Math.sqrt(pb*pb - 4d*pa*(pc-upp));
//		if(Double.isNaN(slow)) {
//			if(Double.isNaN(supp)) return null;
//			//only the lower tip of the parabola is in the desired range
//			double ui = Math.max(0d, Math.min(1d, 0.5d*(-pb-supp)/pa));
//			double ua = Math.max(0d, Math.min(1d, 0.5d*(-pb+supp)/pa));
//			if(pa<0d) { double t = ui; ui = ua; ua = t; }
////			System.out.println("  --> found segment {"+ui+" ... "+ua+"}");
//			if(ui>=ua) return null;
//			return new JDLine[] { new JDLine(lineEnds[0].fractionTowards(ui, lineEnds[1]), lineEnds[0].fractionTowards(ua, lineEnds[1])) };
//		}
//		if(Double.isNaN(supp)) {
//			//only the upper tip of the parabola is in the desired range
//			double ui = Math.max(0d, Math.min(1d, 0.5d*(-pb-slow)/pa));
//			double ua = Math.max(0d, Math.min(1d, 0.5d*(-pb+slow)/pa));
//			if(pa<0d) { double t = ui; ui = ua; ua = t; }
////			System.out.println("  --> found segment {"+ui+" ... "+ua+"}");
//			if(ui>=ua) return null;
//			return new JDLine[] { new JDLine(lineEnds[0].fractionTowards(ui, lineEnds[1]), lineEnds[0].fractionTowards(ua, lineEnds[1])) };
//		}
////		System.out.println("Found intersection: "+lineEnds[0]+"("+vs+"), "+lineEnds[1]+"("+ve+")"); 
//		if(supp<slow) { double temp = slow; slow = supp; supp = temp; };
//		double ui = Math.max(0d, Math.min(1d, 0.5d*(-pb-slow)/pa));
//		double ua = Math.max(0d, Math.min(1d, 0.5d*(-pb+slow)/pa));
//		double vi = Math.max(0d, Math.min(1d, 0.5d*(-pb-supp)/pa));
//		double va = Math.max(0d, Math.min(1d, 0.5d*(-pb+supp)/pa));
//		if(pa<0d) { double t = ui; ui = ua; ua = t; t = vi; vi = va; va = t; }
////		System.out.println("  --> found segment {"+ui+" ... "+ua+"} and {"+vi+" ... "+va+"}");
//		JDLine left=null,right=null;
//		if(vi<ui) left  = new JDLine( lineEnds[0].fractionTowards(vi, lineEnds[1]), lineEnds[0].fractionTowards(ui, lineEnds[1]) );
//		if(ua<va) right = new JDLine( lineEnds[0].fractionTowards(ua, lineEnds[1]), lineEnds[0].fractionTowards(va, lineEnds[1]) );
//		if(left==null) {
//			if(right==null) return null;
//			return new JDLine[] { right };
//		}
//		if(right==null) return new JDLine[] { left };
//		return new JDLine[] { left, right };
	}
	
	public double area() {
		return GeometryTools.area(getCorners());
	}
	public boolean isParallelogram() {
		return parallelogram;
	}
	
	
	private void calcBarycenterHelper() {
		/* Quad:
		         A---B
		         |   |
		         D---C
		 */
		mx = 0.25d*(x[0]+x[1]+x[2]+x[3]);
		my = 0.25d*(y[0]+y[1]+y[2]+y[3]);
		double abx = x[1]-x[0], aby = y[1]-y[0], dcx = x[2]-x[3], dcy = y[2]-y[3];
		parallelU = (Math.abs(abx*dcy - dcx*aby) < 0.0001d * (Math.abs(abx)+Math.abs(aby)+Math.abs(dcx)+Math.abs(dcy)));
		if(parallelU) {
			double dx = abx+dcx, dy = aby+dcy;
			double vn = (0.5d*(x[0]+x[1])-mx)*dy + (my-0.5d*(y[0]+y[1]))*dx;
			double vp = (0.5d*(x[2]+x[3])-mx)*dy + (my-0.5d*(y[2]+y[3]))*dx;
			vx = 2*dy / (vp-vn);
			vy = 2*dx / (vn-vp);
		}
		double adx = x[3]-x[0], ady = y[3]-y[0], bcx = x[2]-x[1], bcy = y[2]-y[1];
		parallelV = (Math.abs(adx*bcy - bcx*ady) < 0.0001d * (Math.abs(adx)+Math.abs(ady)+Math.abs(bcx)+Math.abs(bcy)));
		if(parallelV) {
			double dx = adx+bcx, dy = ady+bcy;
			double un = (0.5d*(x[0]+x[3])-mx)*dy + (my-0.5d*(y[0]+y[3]))*dx;
			double up = (0.5d*(x[1]+x[2])-mx)*dy + (my-0.5d*(y[1]+y[2]))*dx;
			ux = 2*dy / (up-un);
			uy = 2*dx / (un-up);
		}
		parallelogram = parallelU && parallelV;
		if(parallelogram) {
			k = JPlotMath.invert(new double[][] {
				{ 0.5d*(x[1]-x[0]), 0.5d*(x[3]-x[0]), 0.25d*(x[0]+x[1]+x[2]+x[3]) },
				{ 0.5d*(y[1]-y[0]), 0.5d*(y[3]-y[0]), 0.25d*(y[0]+y[1]+y[2]+y[3]) },
				{      0d,               0d,                      1d              }
			});
		} else {
			k = new double[][] {
				{ -x[0]-x[3], x[0]-x[3], x[1]+x[2]-x[0]-x[3], x[2]-x[1]+x[0]-x[3]},
				{ -y[0]-y[3], y[0]-y[3], y[1]+y[2]-y[0]-y[3], y[2]-y[1]+y[0]-y[3]}
			};
		}
		xmin = Math.min(Math.min(x[0], x[1]), Math.min(x[2], x[3]));
		xmax = Math.max(Math.max(x[0], x[1]), Math.max(x[2], x[3]));
		ymin = Math.min(Math.min(y[0], y[1]), Math.min(y[2], y[3]));
		ymax = Math.max(Math.max(y[0], y[1]), Math.max(y[2], y[3]));
	}
	public void addCurve(List<JDPoint> path, double u0, double v0, double u8, double v8, double l, boolean check_emptynes) {
		double[] m0 = { u0, v0 };
		double[] m8 = { u8, v8 };
		m0 = refineUV(u0, v0, l);
		m8 = refineUV(u8, v8, l);
		double[] m4 = refineUV(0.5d*(m0[0]+m8[0]), 0.5d*(m0[1]+m8[1]), l);
		
		double[] m2 = refineUV(0.5d*(m0[0]+m4[0]), 0.5d*(m0[1]+m4[1]), l);
		double[] m6 = refineUV(0.5d*(m4[0]+m8[0]), 0.5d*(m4[1]+m8[1]), l);
		
		double[] m1 = refineUV(0.5d*(m0[0]+m2[0]), 0.5d*(m0[1]+m2[1]), l);
		double[] m3 = refineUV(0.5d*(m2[0]+m4[0]), 0.5d*(m2[1]+m4[1]), l);
		double[] m5 = refineUV(0.5d*(m4[0]+m6[0]), 0.5d*(m4[1]+m6[1]), l);
		double[] m7 = refineUV(0.5d*(m5[0]+m8[0]), 0.5d*(m6[1]+m8[1]), l);
		
		if(path.isEmpty() || !check_emptynes)
			path.add(pointFromUV(m0));
		path.add(pointFromUV(m1));
		path.add(pointFromUV(m2));
		path.add(pointFromUV(m3));
		path.add(pointFromUV(m4));
		path.add(pointFromUV(m5));
		path.add(pointFromUV(m6));
		path.add(pointFromUV(m7));
		path.add(pointFromUV(m8));
	}
	
	private JDPoint closestOnEdge(double x1, double y1, double v1, double x2, double y2, double v2, double px, double py) {
		// (x1 + t*(x2-x1) - px)² + (y1 + t*(y2-y1) - py)² -> minimal
		// d/dt[(a+bt)²] = d/dt[a²+2abt+b²t²] = 2ab + 2b²t
		// ->
		// (x1-px)*(x2-x1)+(x2-x1)²t + (y1-py)*(y2-y1) + (y2-y1)²t = 0
		// [ (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) ] * t = [ (px-x1)*(x2-x1) + (py-y1)*(y2-y1) ]
		if(Double.isNaN(v1) && Double.isNaN(v2)) return null;
		double a = ((px-x1)*(x2-x1)+(py-y1)*(y2-y1)) / ((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		if( a<0d ) a = 0d;
		if(Double.isNaN(v1) && a<0.5001d) a = 0.5001d;
		if( a>1d ) a = 1d;
		if(Double.isNaN(v2) && a>0.4999d) a = 0.4999d;
		return new JDPoint(x1+a*(x2-x1), y1+a*(y2-y1));
	}
	private void rotateIndecesBy(int num) {
		if(num<0) num = (4-num) & 4;
		if(num==0) return;
		if(num==1) {
			double xt = x[0]; x[0] = x[1]; x[1] = x[2]; x[2] = x[3]; x[3] = xt;
			double yt = y[0]; y[0] = y[1]; y[1] = y[2]; y[2] = y[3]; y[3] = yt;
			double vt = value[0]; value[0] = value[1]; value[1] = value[2]; value[2] = value[3]; value[3] = vt;
			calcBarycenterHelper();
		}
		if(num==2) {
			double xt = x[0]; x[0] = x[2]; x[2] = xt; xt = x[1]; x[1] = x[3]; x[3] = xt;
			double yt = y[0]; y[0] = y[2]; y[2] = yt; yt = y[1]; y[1] = y[3]; y[3] = yt;
			double vt = value[0]; value[0] = value[2]; value[2] = vt; vt = value[1]; value[1] = value[3]; value[3] = vt;
			calcBarycenterHelper();
		}
		if(num==3) {
			double xt = x[3]; x[3] = x[2]; x[2] = x[1]; x[1] = x[0]; x[0] = xt;
			double yt = y[3]; y[3] = y[2]; y[2] = y[1]; y[1] = y[0]; y[0] = yt;
			double vt = value[3]; value[3] = value[2]; value[2] = value[1]; value[1] = value[0]; value[0] = vt;
			calcBarycenterHelper();
		}
	}
}
