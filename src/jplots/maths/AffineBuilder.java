package jplots.maths;

public class AffineBuilder {
	double[][] atm;

	public AffineBuilder() {
		atm = new double[][] { { 1d, 0d, 0d }, { 0d, 1d, 0d } };
	}

	public double[][] getMatrix() {
		return atm;
	}

	public AffineBuilder translate(double x, double y) {
		atm[0][2] += x;
		atm[1][2] += y;
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
