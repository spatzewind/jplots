package jplots.maths;

public class AffineBuilder {
	double[][] atm;

	public AffineBuilder() {
		atm = new double[][] { { 1d, 0d, 0d }, { 0d, 1d, 0d } };
	}

	public double[][] getMatrix() {
		return atm;
	}
	public double[][] getInverse() {
		double det = 1d / (atm[0][0]*atm[1][1] - atm[0][1]*atm[1][0]);
		return new double[][]{
			{  det*atm[1][1], -det*atm[0][1], det*(atm[0][1]*atm[1][2]-atm[1][1]*atm[0][2]) },
			{ -det*atm[1][0],  det*atm[0][0], det*(atm[1][0]*atm[0][2]-atm[0][0]*atm[1][2]) }
		};
	}
	
	public AffineBuilder translate(double x, double y) {
		atm[0][2] += x;
		atm[1][2] += y;
		return this;
	}
	
	public AffineBuilder rotate(double a) {
		return rotate(a, false);
	}
	public AffineBuilder rotate(double a, boolean input_in_degrees) {
		double v = input_in_degrees ? a*Math.PI/180d : a;
		double c = Math.cos(v), s = Math.sin(v);
		double x,y;
		x = atm[0][0]*c - atm[1][0]*s; y = atm[0][0]*s + atm[1][0]*c;
		atm[0][0] = x; atm[1][0] = y;
		x = atm[0][1]*c - atm[1][1]*s; y = atm[0][1]*s + atm[1][1]*c;
		atm[0][1] = x; atm[1][1] = y;
		x = atm[0][2]*c - atm[1][2]*s; y = atm[0][2]*s + atm[1][2]*c;
		atm[0][2] = x; atm[1][2] = y;
		return this;
	}
	
	public AffineBuilder scale(double s) {
		return scale(s, s);
	}
	public AffineBuilder scale(double x, double y) {
		atm[0][0] *= x;
		atm[0][1] *= x;
		atm[0][2] *= x;
		atm[1][0] *= y;
		atm[1][1] *= y;
		atm[1][2] *= y;
		return this;
	}
	
}
