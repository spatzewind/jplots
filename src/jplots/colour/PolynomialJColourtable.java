package jplots.colour;

public class PolynomialJColourtable extends JColourtable {

	private int invalidColour;
	private double[] polyRed, polyGreen, polyBlue;

	public PolynomialJColourtable(int nan, double[] red, double[] green, double[] blue) {
		invalidColour = nan;
		polyRed = red;
		polyGreen = green;
		polyBlue = blue;
	}

	@Override
	public int getColour(double percentage) {
		if (Double.isNaN(percentage))
			return invalidColour;
		double p = oclamp(0d, 1d, percentage);
		double	r = polynom(p, polyRed),
				g = polynom(p, polyGreen),
				b = polynom(p, polyBlue);
		return toColour(r, g, b);
	}
	
}
