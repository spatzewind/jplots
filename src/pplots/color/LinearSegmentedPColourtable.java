package pplots.color;

public class LinearSegmentedPColourtable extends PColourtable {

	private int overflowColour, underflowColour, invalidColour;
	private double[] linearColourPosition;
	private int[][] linearColourBounds;
	
	LinearSegmentedPColourtable(int under, int over, int nan, double[] pos, int[][] cols) {
		underflowColour = under;
		overflowColour  = over;
		invalidColour   = nan;
		linearColourPosition = new double[pos.length];
		linearColourBounds = new int[pos.length-1][2];
		for(int c=0; c<pos.length; c++) {
			linearColourPosition[c]  = pos[c];
			if(c+1<pos.length) {
				linearColourBounds[c][0] = cols[c][0];
				linearColourBounds[c][1] = cols[c][1];
			}
		}
	}
	public int getColour(double percentage) {
		if(Double.isNaN(percentage))
			return invalidColour;
		if(percentage<0d)
			return underflowColour;
		if(percentage>1d)
			return overflowColour;
		int fircol=0,seccol=0; double rel=0d;
		for(int p=0; p+1<linearColourPosition.length; p++) {
			if(percentage<linearColourPosition[p])
				break;
			fircol = linearColourBounds[p][0];
			seccol = linearColourBounds[p][1];
			rel = (percentage-linearColourPosition[p]) / (linearColourPosition[p+1]-linearColourPosition[p]);
		}
		return colourmix(fircol, seccol, rel);
	}
}
