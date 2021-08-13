package jplots.color;

public class PolynomialJColourtable extends JColourtable {
	
	double[] polyRed,
	         polyGreen,
	         polyBlue;
	
	public PolynomialJColourtable(double[] red, double[] green, double[] blue) {
		polyRed   = red;
		polyGreen = green;
		polyBlue  = blue;
	}
	

	@Override
	public int getColour(double percentage) {
		double r=0d,g=0d,b=0d;
		for(int ri=0; ri<polyRed.length; ri++)
			r += opow(percentage, ri) * polyRed[ri];
		for(int gi=0; gi<polyGreen.length; gi++)
			g += opow(percentage, gi) * polyGreen[gi];
		for(int bi=0; bi<polyBlue.length; bi++)
			b += opow(percentage, bi) * polyBlue[bi];
		return 0xff000000
				| (((int) (0.4999d+255.5d*oclamp(0d,1d,r)))<<16)
				| (((int) (0.4999d+255.5d*oclamp(0d,1d,g)))<< 8)
				| (((int) (0.4999d+255.5d*oclamp(0d,1d,b)))    );
	}
	
	private double opow(double v, int p) {
		double b = p<0?1d/v:v, res = 1d;
		int pc = Math.abs(p);
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
