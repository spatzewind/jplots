package jplots.colour;

public class ReverseJColourtable extends JColourtable {
	
	private JColourtable baseCT;
	
	public ReverseJColourtable(JColourtable ct) {
		baseCT = ct;
		if(baseCT==null)
			throw new RuntimeException("Cannot create JColortable from null");
	}

	@Override
	public int getColour(double percentage) {
		return baseCT.getColour(1d-percentage);
	}

}
