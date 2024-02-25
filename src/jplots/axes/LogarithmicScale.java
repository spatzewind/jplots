package jplots.axes;

import jplots.JPlot;
import jplots.maths.JPlotMath;
import processing.core.PApplet;
import processing.core.PConstants;

public class LogarithmicScale extends AxisScale {
	
	private boolean linear;
	private int subdiv;
	private double p10, f10;
	
	public LogarithmicScale(JAxis subplot, char axis) {
		super(subplot, axis);
		p10 = 1d;
		f10 = 1d;
		linear = false;
	}
	
	@Override
	public double scale(double v) {
		return Math.log10(v);
	}
	@Override
	public double invscale(double v) {
		return Math.pow(10, v);
	}
	
	@Override
	public void create(double lower_bound, double upper_bound) {
		int[] p = subplot.getSize();
		int ps = axis=='x' ? p[0] : p[1];
		int pe = axis=='x' ? p[0]+p[2] : p[1]+p[3];
		int tickcount = DEFAULT_TICK_COUNT;
		if(axis=='x') {
			subdiv = 1;
			createTicks(lower_bound, upper_bound, tickcount);
			createTickmarks();
			double tmlen = 0d;
			subplot.getPlot().getGraphic().textSize(200);
			subplot.getPlot().getGraphic().textAlign(PConstants.LEFT, PConstants.TOP);
			for (int t = 0; t < ticks.length; t++)
				tmlen += subplot.getPlot().getGraphic().textWidth(tickmarks[t]) / 200f;
			tmlen *= subplot.txtsize / (ticks.length);
			tickcount = Math.max(2, (int) ((pe-ps) / (1.2d * tmlen) + 0.99999999d));
		} else {
			tickcount = Math.max(2, (int) ((pe-ps)/(1.2*subplot.txtsize)));
			if(tickcount>DEFAULT_TICK_COUNT)
				tickcount = Math.max(DEFAULT_TICK_COUNT, (int) ((pe-ps)/(3.5*subplot.txtsize)));
		}
		subdiv = subtickFactor;
		createTicks(lower_bound, upper_bound, tickcount);
		createTickmarks();
		double lmin = Math.log10(lower_bound), lmax = Math.log10(upper_bound);
		pos = new double[ticks.length];
		for(int i=0; i<ticks.length; i++)
			pos[i] = JPlotMath.map(Math.log10(ticks[i]), lmin, lmax, ps, pe);
		
		if(subplot.getPlot().isDebug()) {
			System.out.println("[DEBUG] JAxis-object: "+(""+axis).toUpperCase()+"logAxis - linear section = "+linear);
			System.out.println("[DEBUG] JAxis-object: "+(""+axis).toUpperCase()+"tickfactors={p10: " + p10 + ", f: " + f10+ "}");
		}
	}
	
	private void createTicks(double min_value, double max_value, int desired_count) {
		if (min_value == max_value) {
			ticks = new double[0];
			linear = false;
			return;
		}
		double vin = Math.log10(Math.min(min_value, max_value)), vax = Math.log10(Math.max(min_value, max_value));
		if (Double.isNaN(vin) || Double.isNaN(vax)) {
			ticks = new double[0];
			linear = false;
			return;
		}
		if (vax - vin < 0.1d) {
			createLinear(min_value, max_value, desired_count);
			linear = true;
			return;
		}
		int[] factors = { 99, 49, 19, 9, 5, 2, 1, -2, -5, -10, -20, -50, -100, -200, -500, -1000 };
		int preFac = 9;
		for (int pf : factors) {
			double pfd = pf > 0 ? pf : 1d / (-pf);
			int iin = (int) (pfd * vin) - (vin < 0d ? 1 : 0), iax = (int) (pfd * vax) + (vax > 0d ? 1 : 0);
			double effTickCount = desired_count;
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
			if (iax - iin < effTickCount)
				break;
		}
		if (preFac > 0) {
			int iin = (int) (preFac * vin) - (vin < 0d ? 1 : 0), iax = (int) (preFac * vax) + (vax > 0d ? 1 : 0);
			int idiff = iax - iin;
			ticks = new double[idiff + 1];
			double preFacD = 1d / preFac;
			for (int j = 0; j <= subdiv*idiff; j++) {
				ticks[j] = Math.pow(10d, Math.floor(preFacD * (iin + j)) - (vin + j < 0 ? 1 : 0));
				if ((iin + j) % preFac != 0)
					switch (preFac) {
						case 2:
							ticks[j] *= 5d;
							break;
						case 5:
							ticks[j] *= 2d * ((iin + j) % 5);
							break;
						case 9:
							ticks[j] *= 1 + (iin + j) % 9;
							break;
						case 19:
							ticks[j] *= 1d + 0.5d * ((iin + j) % 19);
							break;
						case 49:
							ticks[j] *= 1d + 0.2d * ((iin + j) % 49);
							break;
						case 99:
							ticks[j] *= 1d + 0.1d * ((iin + j) % 99);
							break;
						default:
							break;
					}
			}
			linear = false;
		} else {
			preFac = -preFac;
			int iin = (int) (vin / preFac) - (vin < 0d ? 1 : 0), iax = (int) (vax / preFac) + (vax > 0d ? 1 : 0);
			int idiff = iax - iin;
			double[] ticks = new double[idiff + 1];
			iin *= preFac;
			iax *= preFac;
			for (int j = 0; j <= idiff; j++) {
				ticks[j] = Math.pow(10d, iin + j * preFac);
			}
			linear = false;
		}
	}
	private void createLinear(double min_value, double max_value, int desired_count) {
		double vin = Math.min(min_value, max_value), vax = Math.max(min_value, max_value);
		p10 = Math.log10(Math.max(vax, -vin)) / 3d - 0.176091259d;
		int p10i = (int) (p10) - (p10 < 0d ? 1 : 0);
		if (p10i == -1)
			p10i = 0;
		p10i *= 3;
		p10 = Math.pow(10d, p10i);
		f10 = 1d;
		double vIn = vin / p10, vAx = vax / p10;
		while (vAx - vIn > desired_count + 0.00000001d) {
			f10 *= 10d;
			vIn /= 10d;
			vAx /= 10d;
		}
		while (true) {
			if (vAx - vIn <= 0.5d * desired_count + 0.00000001d) {
				f10 *= 0.5d;
				vIn *= 2d;
				vAx *= 2d;
			} else {
				break;
			}
			if (vAx - vIn <= 0.8d * desired_count + 0.00000001d) {
				f10 *= 0.8d;
				vIn *= 1.25d;
				vAx *= 1.25d;
			} else {
				break;
			}
			if (vAx - vIn <= 0.5d * desired_count + 0.00000001d) {
				f10 *= 0.5d;
				vIn *= 2d;
				vAx *= 2d;
			} else {
				break;
			}
			if (vAx - vIn <= 0.5d * desired_count + 0.00000001d) {
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
		vS *= subdiv;
		vE *= subdiv;
		ticks = new double[vE + 1 - vS];
		for (int t = vS; t <= vE; t++)
			ticks[t - vS] = t * f10 * p10 / subdiv;
		linear = true;
	}
	
	private void createTickmarks() {
		tickmarks = new String[ticks.length];
		if(linear) {
			double vf = 1d / p10;
			int decimal = (int) (1000d * ticks[1] + 0.5d);
			decimal = decimal % 1000 == 0 ? 0 : decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
			for (int t = 0; t < tickmarks.length; t++) {
				if(decimal==0) tickmarks[t] = ""+(int)(ticks[t]*vf+0.0005d-(ticks[t]*vf<0d?1:0));
				else tickmarks[t] = PApplet.nf((float) (ticks[t] * vf), 0, decimal).replace(",", ".");
			}
			double lvf = Math.log10(p10);
			if(Math.abs(lvf)>2.9d) {
				int ivf = (int) (lvf+0.5d) - (lvf<0d ? -1 : 0);
				tickmarkFactor = JPlot.supportLatex ? "10$^{"+ivf+"}$" : "10^"+ivf;
			}
		} else {
			for(int i=0; i<ticks.length; i++)
				tickmarks[i] = ""+ticks[i]+"";
			tickmarkFactor = "";
		}
	}

}
