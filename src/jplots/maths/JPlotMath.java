package jplots.maths;

public class JPlotMath {

	public final static double DEG_TO_RAD = Math.PI/180d;
	public final static double RAD_TO_DEG = 180d/Math.PI;

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
		if(fm>0 && Float.isInfinite(fm))
			return Float.NaN;
		return fm;
	}
	public static float fmax(float[] arr) {
		float fm = Float.NEGATIVE_INFINITY;
		for(float f: arr) if(Float.isFinite(f) && f>fm) fm = f;
		if(fm<0 && Float.isInfinite(fm))
			return Float.NaN;
		return fm;
	}
	public static float fmin(float[][] arr) {
		float fm = Float.POSITIVE_INFINITY;
		for(float[] fa: arr) for(float f: fa) if(Float.isFinite(f) && f<fm) fm = f;
		if(fm>0 && Float.isInfinite(fm))
			return Float.NaN;
		return fm;
	}
	public static float fmax(float[][] arr) {
		float fm = Float.NEGATIVE_INFINITY;
		for(float[] fa: arr) for(float f: fa) if(Float.isFinite(f) && f>fm) fm = f;
		if(fm<0 && Float.isInfinite(fm))
			return Float.NaN;
		return fm;
	}
	
	public static double dmin(double[] arr) {
		double dm = Double.POSITIVE_INFINITY;
		for(double d: arr) if(Double.isFinite(d) && d<dm) dm = d;
		if(dm>0 && Double.isInfinite(dm))
			return Double.NaN;
		return dm;
	}
	public static double dmax(double[] arr) {
		double dm = Double.NEGATIVE_INFINITY;
		for(double d: arr) if(Double.isFinite(d) && d>dm) dm = d;
		if(dm<0 && Double.isInfinite(dm))
			return Double.NaN;
		return dm;
	}
	public static double dmin(double[][] arr) {
		double dm = Double.POSITIVE_INFINITY;
		for(double[] da: arr) for(double d: da) if(Double.isFinite(d) && d<dm) dm = d;
		if(dm>0 && Double.isInfinite(dm))
			return Double.NaN;
		return dm;
	}
	public static double dmax(double[][] arr) {
		double dm = Double.NEGATIVE_INFINITY;
		for(double[] da: arr) for(double d: da) if(Double.isFinite(d) && d>dm) dm = d;
		if(dm<0 && Double.isInfinite(dm))
			return Double.NaN;
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
		return optimalLinearTicks(vmin, vmax, 10);
	}
	public static double[] optimalLinearTicks(double vmin, double vmax, int maxTickCount) {
		double vin = Math.min(vmin,vmax), vax = Math.max(vmin,vmax);
		double p10 = Math.log10(Math.max(vax, -vin))/3d;
		int p10i = (int) (p10) - (p10<0d ? 1 : 0);
		if(p10i==-1) p10i = 0;
		p10i *= 3;
		p10 = Math.pow(10d, p10i);
		double f10 = 1d;
		double vIn = vin/p10, vAx = vax/p10;
		while(vAx-vIn>maxTickCount+0.00000001d) {
			f10 *= 10d; vIn /= 10d; vAx /= 10d; }
		while(true) {
			if(vAx-vIn<=0.5d*maxTickCount+0.00000001d) {
				f10 *= 0.5d; vIn *= 2d; vAx *= 2d;
			} else { break; }
			if(vAx-vIn<=0.8d*maxTickCount+0.00000001d) {
				f10 *= 0.8d; vIn *= 1.25d; vAx *= 1.25d;
			} else { break; }
			if(vAx-vIn<=0.5d*maxTickCount+0.00000001d) {
				f10 *= 0.5d; vIn *= 2d; vAx *= 2d;
			} else { break; }
			if(vAx-vIn<=0.5d*maxTickCount+0.00000001d) {
				f10 *= 0.5d; vIn *= 2d; vAx *= 2d;
			} else { break; }
		}
		int vS = (int) vIn - (vIn<0d ? 1 : 0), vE = (int) vAx + (vAx<0d ? 0 : 1);
		if(vS>vE) { int vT = vS; vS = vE; vE = vT; }
		double[] ticks = new double[vE+3-vS];
		ticks[0] = p10;
		ticks[1] = f10;
		for(int t=vS; t<=vE; t++)
			ticks[t+2-vS] = t*f10*p10;
		return ticks;
	}

}
