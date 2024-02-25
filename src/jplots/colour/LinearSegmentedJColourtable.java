package jplots.colour;

public class LinearSegmentedJColourtable extends JColourtable {

	private int overflowColour, underflowColour, invalidColour;
	private double[] linearColourPosition;
	private int[][] linearColourBounds;

	public LinearSegmentedJColourtable(int under, int over, int nan, double[] pos, int[][] cols) {
		underflowColour = under;
		overflowColour = over;
		invalidColour = nan;
		linearColourPosition = new double[pos.length];
		linearColourBounds = new int[pos.length - 1][2];
		double cpmin = Double.MAX_VALUE, cpmax = -Double.MAX_VALUE;
		for (int c = 0; c < pos.length; c++) {
			if (pos[c] < cpmin)
				cpmin = pos[c];
			if (pos[c] > cpmax)
				cpmax = pos[c];
			if (c + 1 < pos.length) {
				linearColourBounds[c][0] = cols[c][0];
				linearColourBounds[c][1] = cols[c][1];
			}
		}
		double cpf = 1d / (cpmax - cpmin);
		for (int c = 0; c < pos.length; c++)
			linearColourPosition[c] = (pos[c] - cpmin) * cpf;
	}

	@Override
	public int getColour(double percentage) {
		if (Double.isNaN(percentage))
			return invalidColour;
		if (percentage < 0d)
			return underflowColour;
		if (percentage > 1d)
			return overflowColour;
		int fircol = 0, seccol = 0;
		double rel = 0d;
		for (int p = 0; p + 1 < linearColourPosition.length; p++) {
			if (percentage < linearColourPosition[p])
				break;
			fircol = linearColourBounds[p][0];
			seccol = linearColourBounds[p][1];
			rel = (percentage - linearColourPosition[p]) / (linearColourPosition[p + 1] - linearColourPosition[p]);
		}
		return colourmix(fircol, seccol, rel);
	}
	
	public void setUnderflowColour(int c) {
		underflowColour = c;
	}
	public void setOverflowColour(int c) {
		overflowColour = c;
	}
	public void setUndefinedColour(int c) {
		invalidColour = c;
	}
}
