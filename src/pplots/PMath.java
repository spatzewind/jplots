package pplots;

public class PMath {

	public static int imin(int[] arr) {
		int im = Integer.MAX_VALUE;
		for(int i: arr) if(i<im) im = i;
		return im;
	}
	public static int imax(int[] arr) {
		int im = Integer.MIN_VALUE;
		for(int i: arr) if(i>im) im = i;
		return im;
	}
	public static int imin(int[][] arr) {
		int im = Integer.MAX_VALUE;
		for(int[] ia: arr) for(int i: ia) if(i<im) im = i;
		return im;
	}
	public static int imax(int[][] arr) {
		int im = Integer.MIN_VALUE;
		for(int[] ia: arr) for(int i: ia) if(i>im) im = i;
		return im;
	}

	public static float fmin(float[] arr) {
		float fm = Float.POSITIVE_INFINITY;
		for(float f: arr) if(Float.isFinite(f) && f<fm) fm = f;
		return fm;
	}
	public static float fmax(float[] arr) {
		float fm = Float.NEGATIVE_INFINITY;
		for(float f: arr) if(Float.isFinite(f) && f>fm) fm = f;
		return fm;
	}
	public static float fmin(float[][] arr) {
		float fm = Float.POSITIVE_INFINITY;
		for(float[] fa: arr) for(float f: fa) if(Float.isFinite(f) && f<fm) fm = f;
		return fm;
	}
	public static float fmax(float[][] arr) {
		float fm = Float.NEGATIVE_INFINITY;
		for(float[] fa: arr) for(float f: fa) if(Float.isFinite(f) && f>fm) fm = f;
		return fm;
	}
	
	public static double dmin(double[] arr) {
		double dm = Double.POSITIVE_INFINITY;
		for(double d: arr) if(Double.isFinite(d) && d<dm) dm = d;
		return dm;
	}
	public static double dmax(double[] arr) {
		double dm = Double.NEGATIVE_INFINITY;
		for(double d: arr) if(Double.isFinite(d) && d>dm) dm = d;
		return dm;
	}
	public static double dmin(double[][] arr) {
		double dm = Double.POSITIVE_INFINITY;
		for(double[] da: arr) for(double d: da) if(Double.isFinite(d) && d<dm) dm = d;
		return dm;
	}
	public static double dmax(double[][] arr) {
		double dm = Double.NEGATIVE_INFINITY;
		for(double[] da: arr) for(double d: da) if(Double.isFinite(d) && d>dm) dm = d;
		return dm;
	}
}
