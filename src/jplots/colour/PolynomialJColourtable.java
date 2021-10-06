package jplots.colour;

public class PolynomialJColourtable extends JColourtable {
	
	private int invalidColour;
	private double[] polyRed, polyGreen, polyBlue;
	
	public PolynomialJColourtable(int nan, double[] red, double[] green, double[] blue) {
		invalidColour = nan;
		polyRed   = red;
		polyGreen = green;
		polyBlue  = blue;
	}
	

	@Override
	public int getColour(double percentage) {
		if(Double.isNaN(percentage)) return invalidColour;
		double r=0d,g=0d,b=0d, p=oclamp(0d, 1d, percentage);
		for(int ri=0; ri<polyRed.length; ri++)
			r += opow(p, ri) * polyRed[ri];
		for(int gi=0; gi<polyGreen.length; gi++)
			g += opow(p, gi) * polyGreen[gi];
		for(int bi=0; bi<polyBlue.length; bi++)
			b += opow(p, bi) * polyBlue[bi];
		return 0xff000000
				| (((int) (0.4999d+255.5d*oclamp(0d,1d,r)))<<16)
				| (((int) (0.4999d+255.5d*oclamp(0d,1d,g)))<< 8)
				| (((int) (0.4999d+255.5d*oclamp(0d,1d,b)))    );
	}
	
	private double opow(double v, int p) {
		if(p==0) return 1d;
		double b = p<0?1d/v:v, res = 1d;
		int pc = p<0 ? -p : p;
		while(pc>0) {
			if(pc%2==1)
				res *= b;
			b *= b;
			pc = (pc>>1);
		}
		return res;
	}
	private double oclamp(double vmin, double vmax, double v) {
		if(v<vmin && v<vmax) return vmin<vmax ? vmin : vmax;
		if(v>vmin && v>vmax) return vmax>vmin ? vmax : vmin;
		return v;
	}
}
