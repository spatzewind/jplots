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

	public static float flerp(float val, float in_s, float in_e, float out_s, float out_e) {
		return out_s + (val-in_s)*(out_e-out_s) / (in_e-in_s); }
	public static float[] flerp(float[] arr, float in_s, float in_e, float out_s, float out_e) {
		float convertion = (out_e-out_s) / (in_e-in_s);
		float[] res = new float[arr.length];
		for(int i=0; i<res.length; i++)
			res[i] = out_s + convertion*(arr[i]-in_s);
		return res;
	}
	public static double dlerp(double val, double in_s, double in_e, double out_s, double out_e) {
		return out_s + (val-in_s)*(out_e-out_s) / (in_e-in_s); }
	public static double[] dlerp(double[] arr, double in_s, double in_e, double out_s, double out_e) {
		double convertion = (out_e-out_s) / (in_e-in_s);
		double[] res = new double[arr.length];
		for(int i=0; i<res.length; i++)
			res[i] = out_s + convertion*(arr[i]-in_s);
		return res;
	}

	public static double[] optimalLinearTicks(double vmin, double vmax) {
		double vin = Math.min(vmin,vmax), vax = Math.max(vmin,vmax);
		double p10 = Math.log(vax-vin) / Math.log(10d);
		int p10i = (int) p10 - (p10<0d ? 1 : 0);
		double f10 = 1d;
		p10 = Math.pow(10d, p10i);
		double vIn = vin/p10, vAx = vax/p10;
		while(vAx-vIn<1d) {
			p10i--; p10 /= 10d; vIn *= 10d; vAx *= 10d; }
		if(vAx-vIn<6d) {
			f10 = 0.5d; vIn *= 2d; vAx *= 2d; }
		if(vAx-vIn<6d) {
			f10 = 0.4d; vIn *= 1.25d; vAx *= 1.25d; }
		if(vAx-vIn<6d) {
			f10 = 0.2d; vIn *= 2d; vAx *= 2d; }
		int vS = (int) vIn - (vIn<0d ? 1 : 0), vE = (int) vAx - (vAx<0d ? 1 : 0);
		if(vS>vE) { int vT = vS; vS = vE; vE = vT; }
		double[] ticks = new double[vE+3-vS];
		ticks[0] = p10;
		ticks[1] = f10;
		for(int t=vS; t<=vE; t++)
			ticks[t+2-vS] = t*f10*p10;
		return ticks;
	}
}
