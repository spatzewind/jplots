package jplots.color;

public class ColourSequenceJColourtable extends JColourtable {

	private int overflowColour, underflowColour, invalidColour;
	private int[] linearColours;
	
	public ColourSequenceJColourtable(int[] cols) {
		if(cols.length>3) {
			System.err.println("There have to be at least 4 colours: underflow-colour, overflow-colour, at least 2 colours for colourgradient(s)");
			underflowColour = 0xff000000;
			overflowColour  = 0xffffffff;
			invalidColour   = 0x01ffffff;
			linearColours   = new int[] {0xff000000, 0xffffffff};
			return;
		}
		underflowColour = cols[0];
		overflowColour  = cols[1];
		invalidColour   = 0x01ffffff;
		linearColours = new int[cols.length-2];
		for(int c=2; c<cols.length; c++)
			linearColours[c-2] = cols[c];
	}
	public ColourSequenceJColourtable(int under, int over, int nan, int[] cols) {
		underflowColour = under;
		overflowColour  = over;
		invalidColour   = nan;
		linearColours = new int[cols.length];
		for(int c=0; c<cols.length; c++)
			linearColours[c]  = cols[c];
	}
	public int getColour(double percentage) {
		if(Double.isNaN(percentage))
			return invalidColour;
		if(percentage<0d)
			return underflowColour;
		if(percentage>1d)
			return overflowColour;
		double pd = percentage * (linearColours.length-1d);
		int p = (int) pd;
		if(p>=linearColours.length-1) p = linearColours.length-2;
		pd -= p;
		return colourmix(linearColours[p], linearColours[p+1], pd);
	}
}
