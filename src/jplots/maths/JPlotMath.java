package jplots.maths;

import processing.core.PApplet;

public class JPlotMath {

	public final static double DEG_TO_RAD = Math.PI / 180d;
	public final static double RAD_TO_DEG = 180d / Math.PI;

	public static int imin(int[] arr) {
		int im = Integer.MAX_VALUE;
		for (int i : arr)
			if (i < im)
				im = i;
		return im;
	}

	public static int imax(int[] arr) {
		int im = Integer.MIN_VALUE;
		for (int i : arr)
			if (i > im)
				im = i;
		return im;
	}

	public static int imin(int[][] arr) {
		int im = Integer.MAX_VALUE;
		for (int[] ia : arr)
			for (int i : ia)
				if (i < im)
					im = i;
		return im;
	}

	public static int imax(int[][] arr) {
		int im = Integer.MIN_VALUE;
		for (int[] ia : arr)
			for (int i : ia)
				if (i > im)
					im = i;
		return im;
	}

	public static float fmin(float[] arr) {
		float fm = Float.POSITIVE_INFINITY;
		for (float f : arr)
			if (Float.isFinite(f) && f < fm)
				fm = f;
		if (fm > 0 && Float.isInfinite(fm))
			return Float.NaN;
		return fm;
	}

	public static float fmax(float[] arr) {
		float fm = Float.NEGATIVE_INFINITY;
		for (float f : arr)
			if (Float.isFinite(f) && f > fm)
				fm = f;
		if (fm < 0 && Float.isInfinite(fm))
			return Float.NaN;
		return fm;
	}

	public static float fmin(float[][] arr) {
		float fm = Float.POSITIVE_INFINITY;
		for (float[] fa : arr)
			for (float f : fa)
				if (Float.isFinite(f) && f < fm)
					fm = f;
		if (fm > 0 && Float.isInfinite(fm))
			return Float.NaN;
		return fm;
	}

	public static float fmax(float[][] arr) {
		float fm = Float.NEGATIVE_INFINITY;
		for (float[] fa : arr)
			for (float f : fa)
				if (Float.isFinite(f) && f > fm)
					fm = f;
		if (fm < 0 && Float.isInfinite(fm))
			return Float.NaN;
		return fm;
	}

	public static double dmin(double[] arr) {
		double dm = Double.POSITIVE_INFINITY;
		for (double d : arr)
			if (Double.isFinite(d) && d < dm)
				dm = d;
		if (dm > 0 && Double.isInfinite(dm))
			return Double.NaN;
		return dm;
	}

	public static double dmax(double[] arr) {
		double dm = Double.NEGATIVE_INFINITY;
		for (double d : arr)
			if (Double.isFinite(d) && d > dm)
				dm = d;
		if (dm < 0 && Double.isInfinite(dm))
			return Double.NaN;
		return dm;
	}

	public static double dmin(double[][] arr) {
		double dm = Double.POSITIVE_INFINITY;
		for (double[] da : arr)
			for (double d : da)
				if (Double.isFinite(d) && d < dm)
					dm = d;
		if (dm > 0 && Double.isInfinite(dm))
			return Double.NaN;
		return dm;
	}

	public static double dmax(double[][] arr) {
		double dm = Double.NEGATIVE_INFINITY;
		for (double[] da : arr)
			for (double d : da)
				if (Double.isFinite(d) && d > dm)
					dm = d;
		if (dm < 0 && Double.isInfinite(dm))
			return Double.NaN;
		return dm;
	}

	public static float map(float val, float in_s, float in_e, float out_s, float out_e) {
		return out_s + (val - in_s) * (out_e - out_s) / (in_e - in_s);
	}

	public static float[] map(float[] arr, float in_s, float in_e, float out_s, float out_e) {
		float convertion = (out_e - out_s) / (in_e - in_s);
		float[] res = new float[arr.length];
		for (int i = 0; i < res.length; i++)
			res[i] = out_s + convertion * (arr[i] - in_s);
		return res;
	}

	public static double map(double val, double in_s, double in_e, double out_s, double out_e) {
		return out_s + (val - in_s) * (out_e - out_s) / (in_e - in_s);
	}

	public static double[] map(double[] arr, double in_s, double in_e, double out_s, double out_e) {
		double convertion = (out_e - out_s) / (in_e - in_s);
		double[] res = new double[arr.length];
		for (int i = 0; i < res.length; i++)
			res[i] = out_s + convertion * (arr[i] - in_s);
		return res;
	}

	public static float lerp(float valA, float valB, float fraction) {
		return valA + fraction * (valB - valA);
	}

	public static float[] lerp(float valA, float valB, float[] fraction) {
		float[] res = new float[fraction.length];
		for (int i = 0; i < fraction.length; i++)
			res[i] = valA + fraction[i] * (valB - valA);
		return res;
	}

	public static double lerp(double valA, double valB, double fraction) {
		return valA + fraction * (valB - valA);
	}

	public static double[] lerp(double valA, double valB, double[] fraction) {
		double[] res = new double[fraction.length];
		for (int i = 0; i < fraction.length; i++)
			res[i] = valA + fraction[i] * (valB - valA);
		return res;
	}

	public static float nlerp(float valA, float valB, float fraction) {
		if (Float.isNaN(valA) && fraction >= 0.4999f)
			return valB;
		if (Float.isNaN(valB) && fraction <= 0.5001f)
			return valA;
		return valA + fraction * (valB - valA);
	}

	public static float[] nlerp(float valA, float valB, float[] fraction) {
		float[] res = new float[fraction.length];
		for (int i = 0; i < fraction.length; i++)
			res[i] = Float.isNaN(valA) && fraction[i] >= 0.4999f ? valB
					: Float.isNaN(valB) && fraction[i] <= 0.5001f ? valA : valA + fraction[i] * (valB - valA);
		return res;
	}

	public static double nlerp(double valA, double valB, double fraction) {
		if (Double.isNaN(valA) && fraction >= 0.4999d)
			return valB;
		if (Double.isNaN(valB) && fraction <= 0.5001d)
			return valA;
		return valA + fraction * (valB - valA);
	}

	public static double[] nlerp(double valA, double valB, double[] fraction) {
		double[] res = new double[fraction.length];
		for (int i = 0; i < fraction.length; i++)
			res[i] = Double.isNaN(valA) && fraction[i] >= 0.4999d ? valB
					: Double.isNaN(valB) && fraction[i] <= 0.5001d ? valA : valA + fraction[i] * (valB - valA);
		return res;
	}

	public static double[] log10(double[] arr) {
		double[] res = new double[arr.length];
		for (int a = 0; a < arr.length; a++)
			res[a] = Math.log10(arr[a]);
		return res;
	}

	public static double[] optimalLinearTicks(double vmin, double vmax) {
		return optimalLinearTicks(vmin, vmax, 10);
	}

	public static double[] optimalLinearTicks(double vmin, double vmax, int maxTickCount) {
		if (vmin == vmax)
			return new double[] { 0d, 0d };
		double vin = Math.min(vmin, vmax), vax = Math.max(vmin, vmax);
		double p10 = Math.log10(Math.max(vax, -vin)) / 3d;
		int p10i = (int) (p10) - (p10 < 0d ? 1 : 0);
		if (p10i == -1)
			p10i = 0;
		p10i *= 3;
		p10 = Math.pow(10d, p10i);
		double f10 = 1d;
		double vIn = vin / p10, vAx = vax / p10;
		while (vAx - vIn > maxTickCount + 0.00000001d) {
			f10 *= 10d;
			vIn /= 10d;
			vAx /= 10d;
		}
		while (true) {
			if (vAx - vIn <= 0.5d * maxTickCount + 0.00000001d) {
				f10 *= 0.5d;
				vIn *= 2d;
				vAx *= 2d;
			} else {
				break;
			}
			if (vAx - vIn <= 0.8d * maxTickCount + 0.00000001d) {
				f10 *= 0.8d;
				vIn *= 1.25d;
				vAx *= 1.25d;
			} else {
				break;
			}
			if (vAx - vIn <= 0.5d * maxTickCount + 0.00000001d) {
				f10 *= 0.5d;
				vIn *= 2d;
				vAx *= 2d;
			} else {
				break;
			}
			if (vAx - vIn <= 0.5d * maxTickCount + 0.00000001d) {
				f10 *= 0.5d;
				vIn *= 2d;
				vAx *= 2d;
			} else {
				break;
			}
		}
		int vS = (int) vIn - (vIn < 0d ? 1 : 0), vE = (int) vAx + (vAx < 0d ? 0 : 1);
		if (vS > vE) {
			int vT = vS;
			vS = vE;
			vE = vT;
		}
		double[] ticks = new double[vE + 3 - vS];
		ticks[0] = p10;
		ticks[1] = f10;
		for (int t = vS; t <= vE; t++)
			ticks[t + 2 - vS] = t * f10 * p10;
		return ticks;
	}

	public static double[] optimalLogarithmicTicks(double vmin, double vmax) {
		return optimalLogarithmicTicks(vmin, vmax, 10);
	}

	public static double[] optimalLogarithmicTicks(double vmin, double vmax, int maxTickCount) {
		if (vmin == vmax)
			return new double[] { 0d, 0d };
		double vin = Math.log10(Math.min(vmin, vmax)), vax = Math.log10(Math.max(vmin, vmax));
		if (Double.isNaN(vin) || Double.isNaN(vax))
			return new double[] { 0d, 0d };
		if (vax - vin < 0.1d)
			return optimalLinearTicks(vmin, vmax, maxTickCount);
		int[] factors = { 99, 49, 19, 9, 5, 2, 1, -2, -5, -10, -20, -50, -100, -200, -500, -1000 };
		int preFac = 9;
		for (int pf : factors) {
			double pfd = pf > 0 ? pf : 1d / (-pf);
			int iin = (int) (pfd * vin) - (vin < 0d ? 1 : 0), iax = (int) (pfd * vax) + (vax > 0d ? 1 : 0);
			double effTickCount = maxTickCount;
			switch (pf) {
			case 2:
				effTickCount *= (1d - Math.log10(5.0d)) / 0.50d;
				break;
			case 5:
				effTickCount *= (1d - Math.log10(8.0d)) / 0.20d;
				break;
			case 9:
				effTickCount *= (1d - Math.log10(9.0d)) / 0.10d;
				break;
			case 19:
				effTickCount *= (1d - Math.log10(9.5d)) / 0.05d;
				break;
			case 49:
				effTickCount *= (1d - Math.log10(9.8d)) / 0.02d;
				break;
			case 99:
				effTickCount *= (1d - Math.log10(9.9d)) / 0.01d;
				break;
			default:
				break;
			}
			preFac = pf;
			// System.out.println("preFac="+preFac+" jin/jax="+iin+"/"+iax+"
			// effTickCount="+effTickCount);
			if (iax - iin < effTickCount)
				break;
		}
		if (preFac > 0) {
			// System.out.println("idiff LESS THAN maxTickCount!");
			int iin = (int) (preFac * vin) - (vin < 0d ? 1 : 0), iax = (int) (preFac * vax) + (vax > 0d ? 1 : 0);
			int idiff = iax - iin;
			double[] ticks = new double[2 + idiff + 1];
			double preFacD = 1d / preFac;
			for (int j = 0; j <= idiff; j++) {
				ticks[j + 2] = Math.pow(10d, Math.floor(preFacD * (iin + j)) - (vin + j < 0 ? 1 : 0));
				if ((iin + j) % preFac != 0)
					switch (preFac) {
					case 2:
						ticks[j + 2] *= 5d;
						break;
					case 5:
						ticks[j + 2] *= 2d * ((iin + j) % 5);
						break;
					case 9:
						ticks[j + 2] *= 1 + (iin + j) % 9;
						break;
					case 19:
						ticks[j + 2] *= 1d + 0.5d * ((iin + j) % 19);
						break;
					case 49:
						ticks[j + 2] *= 1d + 0.2d * ((iin + j) % 49);
						break;
					case 99:
						ticks[j + 2] *= 1d + 0.1d * ((iin + j) % 99);
						break;
					default:
						break;
					}
			}
			ticks[0] = -1d;
			ticks[1] = preFacD;
			return ticks;
		} else {
			preFac = -preFac;
			int iin = (int) (vin / preFac) - (vin < 0d ? 1 : 0), iax = (int) (vax / preFac) + (vax > 0d ? 1 : 0);
			int idiff = iax - iin;
			double[] ticks = new double[2 + idiff + 1];
			iin *= preFac;
			iax *= preFac;
			for (int j = 0; j <= idiff; j++) {
				ticks[j + 2] = Math.pow(10d, iin + j * preFac);
			}
			ticks[0] = -1d;
			ticks[1] = 1d;
			return ticks;
		}
	}

	public static double[] optimalTimeTicks(double vmin, double vmax, String unit, String calendar) {
		return optimalTimeTicks(vmin, vmax, unit, calendar, 10);
	}

	public static double[] optimalTimeTicks(double vmin, double vmax, String unit, String calendar, int maxTickCount) {
		if (vmin == vmax)
			return new double[] { 0d, 0d };
		double vin = Math.min(vmin, vmax), vax = Math.max(vmin, vmax);
		DateTime din = DateTime.fromDouble(vin, unit, calendar);
		DateTime dax = DateTime.fromDouble(vax, unit, calendar);
		int[] tin = din.toDate(calendar);
		int[] tax = dax.toDate(calendar);
		if (tin[2] > 0)
			tin[2]--;
		if (tax[2] > 0)
			tax[2]--;
		long diffA = (long) tax[2] - (long) tin[2];
		if (diffA > 5) {
			if (diffA < maxTickCount) {
				double[] ticks = new double[(int) diffA + 3];
				for (int j = 0; j <= diffA; j++) {
					int y = j + tin[2];
					if (y >= 0)
						y++;
					ticks[j + 2] = new DateTime(y + "-01-01 00:00:00", calendar).toDouble(unit, calendar);
				}
				ticks[0] = 0d;
				ticks[1] = 0d;
				return ticks;
			} else {
				int mtc = (int) Math.min(diffA, 0x000000003fffffffL);
				mtc = Math.min(maxTickCount, 4*mtc+1);
				double[] temp = optimalLinearTicks(tin[2], tax[2], mtc);
				double[] ticks = new double[temp.length];
				for (int j = 2; j < ticks.length; j++) {
					int y = (int) (temp[j]+0.5d) - (temp[j]<-0.5d?1:0);
					if (y >= 0)
						y++;
					ticks[j] = new DateTime(y + "-01-01 00:00:00", calendar).toDouble(unit, calendar);
				}
				ticks[0] = 0d;
				ticks[1] = 0d;
				return ticks;
			}
		}
		int diffM = 12 * ((int) diffA) + (tax[1] - tin[1]);
		if (diffM > 5) {
			if (diffM < maxTickCount) {
				double[] ticks = new double[diffM + 3];
				for (int j = 0; j <= diffM; j++) {
					int m = tin[1] + j - 1;
					int y = tin[2] + m / 12;
					if(y>=0) y++;
					m = (m % 12) + 1;
					ticks[j + 2] = new DateTime(y + "-" + m + "-01 00:00:00", calendar).toDouble(unit, calendar);
				}
				ticks[0] = 0d;
				ticks[1] = 0d;
				return ticks;
			} else {
				double[] temp = optimalLinearTicks(tin[1] - 1, tin[1] + diffM - 1, maxTickCount);
				double[] ticks = new double[temp.length];
				for (int j = 2; j < ticks.length; j++) {
					int tm = (int) temp[j];
					int m = 1 + (tm % 12);
					int y = tin[2] + (tm / 12);
					if(y>=0) y++;
					int d = 1 + (int) ((m == 2 ? 27d : 29d) * (temp[j] - tm));
					ticks[j] = new DateTime(y + "-" + m + "-" + d + " 00:00:00", calendar).toDouble(unit, calendar);
				}
				ticks[0] = 0d;
				ticks[1] = 0d;
				return ticks;
			}
		}
		double diffD = dax.diff(din, "days");
		if (diffD > 5d) {
			double[] temp = optimalLinearTicks(din.getDays() + din.getDayFraction(),
					dax.getDays() + dax.getDayFraction(), maxTickCount);
			double[] ticks = new double[temp.length];
			for (int j = 2; j < ticks.length; j++) {
				long l = (long) temp[j];
				double d = temp[j] - l;
				ticks[j] = new DateTime(l, d).toDouble(unit, calendar);
			}
			ticks[0] = 0d;
			ticks[1] = 0d;
			return ticks;
		}

		// double diffH = dax.diff(din, "hours");
		double diffm = dax.diff(din, "minutes");
		double diffS = dax.diff(din, "seconds");
		String u = diffm > 59.99d ? "hours" : diffS > 59.99d ? "minutes" : "seconds";
		int[] ref = din.toDate(calendar);
		int h = u.equals("hours") ? 0 : ref[3];
		int m = u.equals("seconds") ? ref[4] : 0;
		DateTime refDate = new DateTime(
				ref[2] + "-" + ref[1] + "-" + ref[0] + " " + PApplet.nf(h, 2) + ":" + PApplet.nf(m, 2) + ":00",
				calendar);
		double vin2 = refDate.diff(din, u);
		double vax2 = refDate.diff(dax, u);
		double[] temp = optimalLinearTicks(vin2, vax2, maxTickCount);
		double[] ticks = map(temp, vin2, vax2, vin, vax);
		ticks[0] = 0d;
		ticks[1] = 0d;
		return ticks;
	}

	public static double[] mult(double[][] mat, double[] vec) {
		if(mat[0].length != vec.length)
			throw new IllegalArgumentException("Length of vector has to be the same as inner dimension of matrix!");
		double[] res = new double[mat.length];
		for(int j=0; j<res.length; j++) {
			res[j] = 0d;
			for(int i=0; i<vec.length; i++)
				res[j] += mat[j][i] * vec[i];
		}
		return res;
	}
	
	public static double[][] invert(double[][] mat) {
		if(mat.length!=3)
			throw new IllegalAccessError("Cannot invert non-3x3 matrices!");
		double det = 1d / (mat[0][0]*mat[1][1]*mat[2][2]+mat[0][1]*mat[1][2]*mat[2][0]+mat[0][2]*mat[1][0]*mat[2][1]
							-mat[0][0]*mat[1][2]*mat[2][1]-mat[0][1]*mat[1][0]*mat[2][2]-mat[0][2]*mat[1][1]*mat[2][0]);
		double[][] res = new double[3][3];
		res[0][0] =  (mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1]) * det;
		res[0][1] = -(mat[0][1]*mat[2][2]-mat[0][2]*mat[2][1]) * det;
		res[0][2] =  (mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1]) * det;
		res[1][0] = -(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0]) * det;
		res[1][1] =  (mat[0][0]*mat[2][2]-mat[0][2]*mat[2][0]) * det;
		res[1][2] = -(mat[0][0]*mat[1][2]-mat[0][2]*mat[1][0]) * det;
		res[2][0] =  (mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]) * det;
		res[2][1] = -(mat[0][0]*mat[2][1]-mat[0][1]*mat[2][0]) * det;
		res[2][2] =  (mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0]) * det;
		return res;
	}
	
	public static class DateTime {
		private long days;
		private double dayfraction;

		public DateTime() {
			long ctm = System.currentTimeMillis();
			long d = ctm / 86400000L;
			long s = ctm % 86400000L;
			days = 2440587L + d;
			dayfraction = s / 86400000d;
		}

		public DateTime(long i, double f) {
			days = i;
			dayfraction = f;
		}

		public DateTime(String value, String _cal) {
			String cal = _cal == null ? "proleptic" : _cal;
			DateTime refDate = null;
			if (refDate == null && cal.equalsIgnoreCase("gregorian"))
				refDate = parseGregorian(value.trim());
			if (refDate == null && cal.equalsIgnoreCase("julian"))
				refDate = parseJulian(value.trim());
			if (refDate == null && cal.equalsIgnoreCase("proleptic"))
				refDate = parseProleptic(value.trim());
			if (refDate == null)
				throw new IllegalArgumentException("Unknown calendar '" + _cal + "'");
			this.days = refDate.getDays();
			this.dayfraction = refDate.getDayFraction();
		}

		public long getDays() {
			return days;
		}

		public double getDayFraction() {
			return dayfraction;
		}

		public static DateTime fromDouble(double _in, String _unit, String _cal) {
			String[] uuu = _unit.trim().split("since");
			if (uuu.length < 1 || uuu.length > 2)
				throw new IllegalArgumentException(
						"Illegal unit-string. Expected format: <sec,min,days,...> since <reference dd.mm.yyyy hh:mm:ss>");
			for (int u = 0; u < uuu.length; u++)
				uuu[u] = uuu[u].trim();
			if (uuu.length == 1)
				uuu = new String[] { _unit.trim(), "2000-01-01 00:00:00" };
			String cal = _cal == null ? "proleptic" : _cal;
			DateTime refDate = null;
			if (refDate == null && cal.equalsIgnoreCase("gregorian"))
				refDate = parseGregorian(uuu[1].trim());
			if (refDate == null && cal.equalsIgnoreCase("julian"))
				refDate = parseJulian(uuu[1].trim());
			if (refDate == null && cal.equalsIgnoreCase("proleptic"))
				refDate = parseProleptic(uuu[1].trim());
			if (refDate == null)
				throw new IllegalArgumentException("Unknown calendar '" + _cal + "'");
			return refDate.add(_in, uuu[0]);
		}

		public int[] toDate(String _cal) {
			int calID = -1;
			if (_cal.equalsIgnoreCase("gregorian"))
				calID = 1;
			if (_cal.equalsIgnoreCase("julian"))
				calID = 2;
			if (_cal.equalsIgnoreCase("proleptic"))
				calID = 3;
			if (calID < 0)
				throw new IllegalArgumentException("Unknown calendar: " + _cal);
			int stage = 1000000000;
			int year = Integer.MIN_VALUE;
			long dayDiff = Long.MIN_VALUE;
			DateTime testDate = null;
			while (dayDiff < 0 || stage > 1) {
				if (dayDiff >= 0) {
					year -= stage;
					stage /= 10;
				}
				switch (calID) {
				case 1:
					testDate = parseGregorian(year + "-01-01 00:00:00");
					break;
				case 2:
					testDate = parseJulian(year + "-01-01 00:00:00");
					break;
				case 3:
					testDate = parseProleptic(year + "-01-01 00:00:00");
					break;
				default:
					break;
				}
				dayDiff = testDate.getDays() - this.days;
				if (dayDiff < 0)
					year += stage;
			}
			if (dayDiff > 0)
				year -= stage;
			switch (calID) {
			case 1:
				testDate = parseGregorian(year + "-01-01 00:00:00");
				break;
			case 2:
				testDate = parseJulian(year + "-01-01 00:00:00");
				break;
			case 3:
				testDate = parseProleptic(year + "-01-01 00:00:00");
				break;
			default:
				break;
			}
			dayDiff = this.days - testDate.getDays();
			boolean leap = (year % 4 == 0);
			if (calID == 3 || (calID == 1 && year > 1582))
				leap = (year % 400 == 0 ? true : year % 100 == 0 ? false : year % 4 == 0);
			int[] dm = daycountToDate((int) dayDiff, leap);
			int sec = (int) (86400d * dayfraction);
			int min = sec / 60;
			int hour = min / 60;
			return new int[] { dm[0], dm[1], year, hour, min % 60, sec % 60 };
		}

		public String toDateString(String _cal) {
			return toDateString(_cal, false);
		}

		public String toDateString(String _cal, boolean withAbbr) {
			int[] dd = toDate(_cal);
			int cy = PApplet.year();
			cy -= cy % 10;
			int dy = cy - (cy % 100);
			if (dd[2] >= cy - 80 && dd[2] <= cy + 20 && withAbbr)
				return dd[0] + "." + dd[1] + ".'" + (dd[2] - dy) + "  " + PApplet.nf(dd[3], 2) + ":"
						+ PApplet.nf(dd[4], 2) + ":" + PApplet.nf(dd[5], 2);
			return dd[0] + "." + dd[1] + "." + dd[2] + " " + PApplet.nf(dd[3], 2) + ":" + PApplet.nf(dd[4], 2) + ":"
					+ PApplet.nf(dd[5], 2);
		}

		public String format(String _format, String _cal) {
			int[] dd = toDate(_cal);
			String res = "" + _format + "";
			// System.out.println("res=<"+res+"> and
			// dd={"+dd[0]+","+dd[1]+","+dd[2]+","+dd[3]+","+dd[4]+","+dd[5]+"}");
			if (res.contains("y")) {
				int a = res.indexOf("y");
				int b = res.lastIndexOf("y");
				String ys = "";
				switch (b + 1 - a) {
				case 1:
					int y1 = dd[2] / 10 - (dd[2] < 0 ? 1 : 0);
					ys = "" + (dd[2] - 10 * y1);
					break;
				case 2:
					int y2 = dd[2] / 100 - (dd[2] < 0 ? 1 : 0);
					ys = dd[2] > -99 && dd[2] < 0 ? "" + dd[2] : PApplet.nf(dd[2] - 100 * y2, 2);
					break;
				case 3:
					int y3 = dd[2] / 1000 - (dd[2] < 0 ? 1 : 0);
					ys = dd[2] > -999 && dd[2] < 0 ? "" + dd[2] : PApplet.nf(dd[2] - 1000 * y3, 3);
					break;
				case 4:
					int y4 = dd[2] / 10000 - (dd[2] < 0 ? 1 : 0);
					ys = dd[2] > -9999 && dd[2] < 0 ? "" + dd[2] : PApplet.nf(dd[2] - 10000 * y4, 4);
					break;
				case 5:
					int y5 = dd[2] / 100000 - (dd[2] < 0 ? 1 : 0);
					ys = dd[2] > -99999 && dd[2] < 0 ? "" + dd[2] : PApplet.nf(dd[2] - 100000 * y5, 5);
					break;
				case 6:
					int y6 = dd[2] / 1000000 - (dd[2] < 0 ? 1 : 0);
					ys = dd[2] > -999999 && dd[2] < 0 ? "" + dd[2] : PApplet.nf(dd[2] - 1000000 * y6, 6);
					break;
				default:
					ys = PApplet.nf(dd[2], b + 1 - a);
				}
				res = res.substring(0, a) + ys + res.substring(b + 1);
			}
			if (res.contains("m")) {
				int a = res.indexOf("m");
				int b = res.lastIndexOf("m");
				res = res.substring(0, a) + PApplet.nf(dd[1], b + 1 - a) + res.substring(b + 1);
			}
			if (res.contains("d")) {
				int a = res.indexOf("d");
				int b = res.lastIndexOf("d");
				res = res.substring(0, a) + PApplet.nf(dd[0], b + 1 - a) + res.substring(b + 1);
			}
			if (res.contains("H")) {
				int a = res.indexOf("H");
				int b = res.lastIndexOf("H");
				res = res.substring(0, a) + PApplet.nf(dd[3], b + 1 - a) + res.substring(b + 1);
			}
			if (res.contains("M")) {
				int a = res.indexOf("M");
				int b = res.lastIndexOf("M");
				res = res.substring(0, a) + PApplet.nf(dd[4], b + 1 - a) + res.substring(b + 1);
			}
			if (res.contains("s")) {
				int a = res.indexOf("s");
				int b = res.lastIndexOf("s");
				res = res.substring(0, a) + PApplet.nf(dd[5], b + 1 - a) + res.substring(b + 1);
			}
			return res;
		}

		public double toDouble(String _unit, String _cal) {
			String[] uuu = _unit.trim().split("since");
			if (uuu.length < 1 || uuu.length > 2)
				throw new IllegalArgumentException(
						"Illegal unit-string. Expected format: <sec,min,days,...> since <reference dd.mm.yyyy hh:mm:ss>");
			for (int u = 0; u < uuu.length; u++)
				uuu[u] = uuu[u].trim();
			if (uuu.length == 1)
				uuu = new String[] { _unit.trim(), "1.1.2000 00:00:00" };
			String cal = _cal == null ? "proleptic" : _cal;
			DateTime refDate = null;
			if (refDate == null && cal.equalsIgnoreCase("gregorian"))
				refDate = parseGregorian(uuu[1]);
			if (refDate == null && cal.equalsIgnoreCase("julian"))
				refDate = parseJulian(uuu[1]);
			if (refDate == null && cal.equalsIgnoreCase("proleptic"))
				refDate = parseProleptic(uuu[1]);
			if (refDate == null)
				throw new IllegalArgumentException("Unknown calendar '" + _cal + "'");
			return diff(refDate, uuu[0]);
		}

		public DateTime add(double diff) {
			return this.add(diff, "days");
		}

		public DateTime add(double diff, String _unit) {
			double ddd = diff * timeFactor(_unit);
			long dti = (long) ddd;
			double dtf = ddd - dti;
			days += dti;
			dayfraction += dtf;
			if (dayfraction < 0d) {
				days--;
				dayfraction += 1d;
			}
			if (dayfraction >= 1d) {
				days++;
				dayfraction -= 1d;
			}
			return this;
		}

		public double diff(DateTime dt, String _unit) {
			double dif = (days - dt.getDays()) + (dayfraction - dt.getDayFraction());
			return dif / timeFactor(_unit);
		}

		private static double timeFactor(String _unit) {
			if (_unit.equalsIgnoreCase("nanoseconds") || _unit.equalsIgnoreCase("ns"))
				return 1d / 86400000000d;
			if (_unit.equalsIgnoreCase("milliseconds") || _unit.equalsIgnoreCase("ms"))
				return 1d / 86400000d;
			if (_unit.equalsIgnoreCase("seconds") || _unit.equalsIgnoreCase("sec"))
				return 1d / 86400d;
			if (_unit.equalsIgnoreCase("minutes") || _unit.equalsIgnoreCase("min"))
				return 1d / 1440d;
			if (_unit.equalsIgnoreCase("hours") || _unit.equalsIgnoreCase("h"))
				return 1d / 24d;
			if (_unit.equalsIgnoreCase("days"))
				return 1d;
			if (_unit.equalsIgnoreCase("weaks"))
				return 7d;
			if (_unit.equalsIgnoreCase("years") || _unit.equalsIgnoreCase("a"))
				return 365.2422d;
			if (_unit.equalsIgnoreCase("kiloyears") || _unit.equalsIgnoreCase("ka"))
				return 365242.2d;
			if (_unit.equalsIgnoreCase("megayears") || _unit.equalsIgnoreCase("ma"))
				return 365242200d;
			if (_unit.equalsIgnoreCase("gigayears") || _unit.equalsIgnoreCase("ga"))
				return 365242200000d;
			throw new IllegalArgumentException("Unknown unit '" + _unit + "'");
		}

		private static DateTime parseGregorian(String date) {
			String[] ss = date.split(" ");
			String[] dd = ss[0].substring(ss[0].charAt(0) == '-' ? 1 : 0).split("-");
			int[] ddi = { _int(dd[2]), _int(dd[1]), _int((ss[0].charAt(0) == '-' ? "-" : "") + dd[0]) };
			boolean isJul = ddi[2] < 1582;
			if (ddi[2] == 1582) {
				isJul = ddi[1] < 10;
				if (ddi[1] == 10)
					isJul = ddi[0] < 15;
			}
			return isJul ? parseJulian(date) : parseProleptic(date);
		}

		private static DateTime parseJulian(String date) { // also proleptic julian
			String[] ss = date.split(" ");
			String[] dd = ss[0].substring(ss[0].charAt(0) == '-' ? 1 : 0).split("-"), tt = null;
			if (ss.length > 1)
				tt = ss[1].split(":");
			else
				tt = new String[] { "00", "00", "00" };
			int[] ddi = { _int(dd[2]), _int(dd[1]), _int((ss[0].charAt(0) == '-' ? "-" : "") + dd[0]) };
			int[] tti = { 0, 0, 0 };
			for (int i = 0; i < tt.length; i++)
				tti[i] = _int(tt[i]);
			double df = tti[2] / 86400d + tti[1] / 1440d + tti[0] / 24d;
			int dc = dateToDaycount(ddi, ddi[2] % 4 == 0);
			int y = ddi[2] < 0 ? 4716 + ddi[2] : 4715 + ddi[2];
			int y4 = y / 4 - (y < 0 ? 1 : 0);
			y -= 4 * y4;
			long dl = 1461L * y4 + 365L * y + dc;
			return new DateTime(dl, df);
		}

		private static DateTime parseProleptic(String date) { // proleptic gregorian
			String[] ss = date.split(" ");
			String[] dd = ss[0].substring(ss[0].charAt(0) == '-' ? 1 : 0).split("-"), tt = null;
			if (ss.length > 1)
				tt = ss[1].split(":");
			else
				tt = new String[] { "00", "00", "00" };
			int[] ddi = { _int(dd[2]), _int(dd[1]), _int((ss[0].charAt(0) == '-' ? "-" : "") + dd[0]) };
			int[] tti = { 0, 0, 0 };
			for (int i = 0; i < tt.length; i++)
				tti[i] = _int(tt[i]);
			double df = tti[2] / 86400d + tti[1] / 1440d + tti[0] / 24d;
			boolean leap = ddi[2] % 4 == 0;
			if (leap && ddi[2] % 100 == 0)
				leap = ddi[2] % 400 == 0;
			int dc = dateToDaycount(ddi, leap);
			int y = ddi[2] < 0 ? ddi[2] : ddi[2] - 1;
			int yQ = y / 400 - (y < 0 ? 1 : 0);
			y -= 400 * yQ;
			int yC = y / 100 - (y < 0 ? 1 : 0);
			y -= 100 * yC;
			int y4 = y / 4 - (y < 0 ? 1 : 0);
			y -= 4 * y4;
			long dl = 1721425L + 146097L * yQ + 36524L * yC + 1461L * y4 + 365L * y + dc;
			return new DateTime(dl, df);
		}

		private static int _int(String s) {
			String j = (s.charAt(0) == '-' || s.charAt(0) == '+') ? s.substring(1) : s;
			return Integer.parseInt((s.charAt(0) == '-' ? "-" : "") + j.replaceFirst("^0+(?!$)", ""));
		}

		private static int dateToDaycount(int[] dd, boolean isLeap) {
			switch (dd[1]) {
			case 1:
				return dd[0] - 1;
			case 2:
				return 30 + dd[0];
			case 3:
				return (isLeap ? 59 : 58) + dd[0];
			case 4:
				return (isLeap ? 90 : 89) + dd[0];
			case 5:
				return (isLeap ? 120 : 119) + dd[0];
			case 6:
				return (isLeap ? 151 : 150) + dd[0];
			case 7:
				return (isLeap ? 181 : 180) + dd[0];
			case 8:
				return (isLeap ? 212 : 211) + dd[0];
			case 9:
				return (isLeap ? 243 : 242) + dd[0];
			case 10:
				return (isLeap ? 273 : 272) + dd[0];
			case 11:
				return (isLeap ? 304 : 303) + dd[0];
			case 12:
				return (isLeap ? 334 : 333) + dd[0];
			default:
				return 0;
			}
		}

		private static int[] daycountToDate(int dc, boolean isLeap) {
			int d = dc + 1;
			if (d < 32)
				return new int[] { d, 1 }; // january
			d -= 31;
			// february
			if ((d < 29) || (isLeap && d < 30))
				return new int[] { d, 2 };
			d -= (isLeap ? 29 : 28);
			if (d < 32)
				return new int[] { d, 3 }; // march
			d -= 31;
			if (d < 31)
				return new int[] { d, 4 }; // april
			d -= 30;
			if (d < 32)
				return new int[] { d, 5 }; // mai
			d -= 31;
			if (d < 31)
				return new int[] { d, 6 }; // june
			d -= 30;
			if (d < 32)
				return new int[] { d, 7 }; // july
			d -= 31;
			if (d < 32)
				return new int[] { d, 8 }; // august
			d -= 31;
			if (d < 31)
				return new int[] { d, 9 }; // september
			d -= 30;
			if (d < 32)
				return new int[] { d, 10 }; // october
			d -= 31;
			if (d < 31)
				return new int[] { d, 11 }; // november
			d -= 30;
			if (d < 32)
				return new int[] { d, 12 }; // december
			throw new RuntimeException("Something went wrong during date-conversion <daycountToDate(dc=" + dc
					+ ", isLeap=" + isLeap + ")>");
		}
	}
}
