package jplots.colour;

public class MultiPolynomialJColourtable extends JColourtable {

	private int invalidColour, overflowColour, underflowColour;
	private double[] polySegments;
	private double[][] polyRed, polyGreen, polyBlue;
	
	public MultiPolynomialJColourtable(int nan, int under, int over, double[] segments, double[][] red, double[][] green, double[][] blue) {
		invalidColour = nan;
		underflowColour = under;
		overflowColour = over;
		polySegments = segments;
		polyRed = red;
		polyGreen = green;
		polyBlue = blue;
	}
	
	@Override
	public int getColour(double percentage) {
		if (Double.isNaN(percentage))
			return invalidColour;
		if(percentage<polySegments[0])
			return underflowColour;
		if(percentage>polySegments[polySegments.length-1])
			return overflowColour;
		double r = 0d, g = 0d, b = 0d;
		for(int s=0; s+1<polySegments.length; s++) {
			double[] pr = polyRed[s],
					pg = polyGreen[s],
					pb = polyBlue[s];
			r = polynom(percentage, pr);
			g = polynom(percentage, pg);
			b = polynom(percentage, pb);
		}
		return toColour(r,g,b);
	}
}
