package jplots.axes;

import java.util.Arrays;

import jplots.maths.JPlotMath;
import jplots.maths.JPlotMath.DateTime;
import processing.core.PApplet;
import processing.core.PConstants;

public class DateTimeScale extends AxisScale {
	
	private String tUnit, tCalendar, tFormat;
	
	public DateTimeScale(JAxis subplot, char axis, String unit, String calendar, String format) {
		super(subplot, axis);
		tUnit = unit;
		tCalendar = calendar;
		tFormat = format;
	}
	
	@Override
	public double scale(double v) { return v; }
	@Override
	public double[] scale(double[] v) { return v; }
	@Override
	public double[][] scale(double[][] v) { return v; }
	@Override
	public double[][][] scale(double[][][] v) { return v; }
	@Override
	public double invscale(double v) { return v; }
	@Override
	public double[] invscale(double[] v) { return v; }
	@Override
	public double[][] invscale(double[][] v) { return v; }
	@Override
	public double[][][] invscale(double[][][] v) { return v; }
	
	@Override
	public void create(double lower_bound, double upper_bound) {
		int[] p = subplot.getSize();
		int ps = axis=='x' ? p[0] : p[1];
		int pe = axis=='x' ? p[0]+p[2] : p[1]+p[3];
		int tickcount = DEFAULT_TICK_COUNT;
		if(axis=='x') {
			createTicks(lower_bound, upper_bound, tickcount);
			createTickmarks();
			double tmlen = 0d;
			subplot.getPlot().getGraphic().textSize(200);
			subplot.getPlot().getGraphic().textAlign(PConstants.LEFT, PConstants.TOP);
			for (int t = 0; t < ticks.length; t++)
				tmlen += subplot.getPlot().getGraphic().textWidth(tickmarks[t]) / 200f;
			tmlen *= subplot.txtsize / (ticks.length);
			debug("mean tickmark length = "+tmlen+" px for available space of "+(pe-ps)+" px");
			tickcount = Math.max(2, (int) ((pe-ps) / (1.2d * tmlen) + 0.99999999d));
			debug("new desired tickcount = "+tickcount);
		} else {
			tickcount = Math.max(2, (int) ((pe-ps)/(1.2*subplot.txtsize)));
			if(tickcount>DEFAULT_TICK_COUNT)
				tickcount = Math.max(DEFAULT_TICK_COUNT, (int) ((pe-ps)/(3.5*subplot.txtsize)));
		}
		createTicks(lower_bound, upper_bound, tickcount);
		createTickmarks();
		pos = JPlotMath.map(ticks, lower_bound, upper_bound, ps, pe);
		
		if(subplot.getPlot().isDebug()) {
			System.out.println("[DEBUG] JAxis-object: "+(""+axis).toUpperCase()+"timeAxis - unit=\""+
					tUnit+"\" calendar=\""+tCalendar+"\" format=\""+tFormat+"\"");
		}
	}
	
	private void createTicks(double min_value, double max_value, int desired_count) {
		debug("min/max: "+min_value+", "+max_value);
		if (min_value == max_value) {
			ticks = new double[0];
			debug("min/max values are equal --> return 0 ticks");
			return;
		}
		double vin = Math.min(min_value, max_value), vax = Math.max(min_value, max_value);
		DateTime din = DateTime.fromDouble(vin, tUnit, tCalendar);
		DateTime dax = DateTime.fromDouble(vax, tUnit, tCalendar);
		int[] tin = din.toDate(tCalendar);
		int[] tax = dax.toDate(tCalendar);
		debug("dates: min="+Arrays.toString(tin)+" max="+Arrays.toString(tax));
		long diffA = (long) tax[2] - (long) tin[2];
		debug("diffA = "+diffA);
		if (diffA > 5) {
			if (diffA < desired_count) {
				ticks = new double[(int) diffA + 3];
				for (int j = 0; j <= diffA; j++) {
					int y = j + tin[2];
					ticks[j] = new DateTime(y + "-01-01 00:00:00", tCalendar).toDouble(tUnit, tCalendar);
				}
			} else {
				int mtc = (int) Math.min(diffA, 0x000000003fffffffL);
				mtc = Math.min(desired_count, 4*mtc+1);
				double[] temp = optimalLinearTicks(tin[2], tax[2], mtc);
				ticks = new double[temp.length];
				for (int j = 0; j < ticks.length; j++) {
					int y = (int) (temp[j]+0.5d) - (temp[j]<-0.5d?1:0);
					ticks[j] = new DateTime(y + "-01-01 00:00:00", tCalendar).toDouble(tUnit, tCalendar);
				}
			}
			debug("diffA>5 --> return "+ticks.length+" ticks");
			return;
		}
		if(tin[2]>0) tin[2]--;
		if(tax[2]>0) tax[2]--;
		int diffM = 12 * ((int) diffA) + (tax[1] - tin[1]);
		debug("diffM = "+diffM);
		if (diffM > 5) {
			if (diffM < desired_count) {
				ticks = new double[diffM + 1];
				for (int j = 0; j <= diffM; j++) {
					int m = tin[1] + j - 1;
					int y = tin[2] + m / 12;
					if(y>=0) y++;
					m = (m % 12) + 1;
					ticks[j] = new DateTime(y + "-" + m + "-01 00:00:00", tCalendar).toDouble(tUnit, tCalendar);
				}
			} else {
				double[] temp = optimalLinearTicks(tin[1] - 1, tin[1] + diffM - 1, desired_count);
				ticks = new double[temp.length];
				for (int j = 0; j < ticks.length; j++) {
					int tm = (int) temp[j];
					int m = 1 + (tm % 12);
					int y = tin[2] + (tm / 12);
					if(y>=0) y++;
					int d = 1 + (int) ((m == 2 ? 27d : 29d) * (temp[j] - tm));
					ticks[j] = new DateTime(y + "-" + m + "-" + d + " 00:00:00", tCalendar).toDouble(tUnit, tCalendar);
				}
			}
			debug("diffM>5 --> return "+ticks.length+" ticks");
			return;
		}
		double diffD = dax.diff(din, "days");
		debug("diffD = "+diffD);
		if (diffD > 5d) {
			double[] temp = optimalLinearTicks(din.getDays() + din.getDayFraction(),
					dax.getDays() + dax.getDayFraction(), desired_count);
			ticks = new double[temp.length];
			for (int j = 0; j < ticks.length; j++) {
				long l = (long) temp[j];
				double d = temp[j] - l;
				ticks[j] = new DateTime(l, d).toDouble(tUnit, tCalendar);
			}
			debug("diffD>5.0 --> return "+ticks.length+" ticks");
			return;
		}

		// double diffH = dax.diff(din, "hours");
		double diffm = dax.diff(din, "minutes");
		double diffS = dax.diff(din, "seconds");
		String u = diffm > 59.99d ? "hours" : diffS > 59.99d ? "minutes" : "seconds";
		int[] ref = din.toDate(tCalendar);
		int h = u.equals("hours") ? 0 : ref[3];
		int m = u.equals("seconds") ? ref[4] : 0;
		DateTime refDate = new DateTime(
				ref[2] + "-" + ref[1] + "-" + ref[0] + " " + PApplet.nf(h, 2) + ":" + PApplet.nf(m, 2) + ":00",
				tCalendar);
		double vin2 = refDate.diff(din, u);
		double vax2 = refDate.diff(dax, u);
		double[] temp = optimalLinearTicks(vin2, vax2, desired_count);
		ticks = JPlotMath.map(temp, vin2, vax2, vin, vax);
	}
	private double[] optimalLinearTicks(double vmin, double vmax, int maxTickCount) {
		double vin = Math.min(vmin, vmax), vax = Math.max(vmin, vmax);
		double p10 = Math.log10(Math.max(vax, -vin)) / 3d - 0.176091259d;
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
		double[] ticks = new double[vE + 1 - vS];
		for (int t = vS; t <= vE; t++)
			ticks[t - vS] = t * f10 * p10;
		return ticks;
	}
	
	private void createTickmarks() {
		tickmarks = new String[ticks.length];
		for(int i=0; i<ticks.length; i++)
			tickmarks[i] = DateTime.fromDouble(ticks[i], tUnit, tCalendar).format(tFormat, tCalendar);
		debug("created "+tickmarks.length+" tickmarks");
	}
	
}
