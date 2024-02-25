package jplots.axes;

public abstract class AxisScale {
	public final static int DEFAULT_TICK_COUNT = 10;
	
//	protected final static String ltxMinus = "$\\raisebox{0.15em}{\\scalebox{0.5}[1.0]{\\( - \\)}}$";
	protected final static String ltxMinus = "$\\raisebox{0.15em}{\\text{-}}$";
	
	protected JAxis subplot;
	protected char axis;
	protected int subtickFactor;
	protected double[] ticks;
	protected double[] pos;
	protected String[] tickmarks;
	protected String tickmarkFactor;
	
	public AxisScale(JAxis subplot, char axis) {
		this.subplot = subplot;
		this.axis = axis;
		ticks = new double[0];
		pos = new double[0];
		tickmarks = new String[0];
		tickmarkFactor = "";
		subtickFactor = 1;
	}
	
	public abstract void create(double lower_bound, double upper_bound);
	
	public double       scale(double v) {
		return v; }
	public double[]     scale(double[] v) { double[] res=new double[v.length];
		for(int i=0; i<v.length; i++) res[i]=scale(v[i]);
		return res; }
	public double[][]   scale(double[][] v) { double[][] res=new double[v.length][v[0].length];
		for(int j=0; j<v.length; j++) for(int i=0; i<v[0].length; i++) res[j][i]=scale(v[j][i]);
		return res; }
	public double[][][] scale(double[][][] v) { double[][][] res=new double[v.length][v[0].length][v[0][0].length];
		for(int k=0; k<v.length; k++) for(int j=0; j<v[0].length; j++) for(int i=0; i<v[0][0].length; i++)
			res[k][j][i]=scale(v[k][j][i]);
		return res; }
	public double       invscale(double v) {
		return v; }
	public double[]     invscale(double[] v) { double[] res=new double[v.length];
		for(int i=0; i<v.length; i++) res[i]=invscale(v[i]);
		return res; }
	public double[][]   invscale(double[][] v) { double[][] res=new double[v.length][v[0].length];
		for(int j=0; j<v.length; j++) for(int i=0; i<v[0].length; i++) res[j][i]=invscale(v[j][i]);
		return res; }
	public double[][][] invscale(double[][][] v) { double[][][] res=new double[v.length][v[0].length][v[0][0].length];
		for(int k=0; k<v.length; k++) for(int j=0; j<v[0].length; j++) for(int i=0; i<v[0][0].length; i++)
			res[k][j][i]=invscale(v[k][j][i]);
		return res; }
	
	public double[] getTicks() {
		return ticks;
	}
	public double[] getPos() {
		return pos;
	}
	public int getSubtickFactor() {
		return subtickFactor;
	}
	public void setSubtickFactor(int stf) {
		subtickFactor = stf;
	}
	public String[] getTickmarks() {
		return tickmarks;
	}
	public String getTickmarkFactor() {
		return tickmarkFactor;
	}
	
	protected void debug(String message) {
		if(subplot.getPlot().isDebug()) System.out.println("[DEBUG] "+this.getClass().getSimpleName()+": "+message);
	}
}
