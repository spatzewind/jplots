package jplots.axes;

import jplots.JPlot;
import jplots.maths.JPlotMath;
import processing.core.PApplet;
import processing.core.PConstants;

public class LinearScale extends AxisScale {
	
	private int subdiv;
	private double p10, f10;
	
	public LinearScale(JAxis subplot, char axis) {
		super(subplot, axis);
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
			subdiv = 1;
			createTicks(lower_bound, upper_bound, tickcount);
			createTickmarks(false);
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
		createTickmarks(true);
		pos = JPlotMath.map(ticks, lower_bound, upper_bound, ps, pe);
		
		if(subplot.getPlot().isDebug()) {
			System.out.println("[DEBUG] JAxis-object: "+(""+axis).toUpperCase()+"tickfactors={p10: " + p10 + ", f: " + f10+ "}");
		}
	}
	
	private void createTicks(double min_value, double max_value, int desired_count) {
		if (min_value == max_value) {
			p10 = 0d;
			f10 = 0d;
			ticks = new double[0];
			tickmarkFactor = "";
			return;
		}
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
	}
	
	private void createTickmarks(boolean is_final) {
		tickmarks = new String[ticks.length];
		double vf = 1d / p10;
		int decimal = (int) (1000d * f10 + 0.5d);
		decimal = decimal % 1000 == 0 ? 0 : decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
		for (int t = 0; t < tickmarks.length; t++) {
			if(decimal==0) tickmarks[t] = ""+(int)(ticks[t]*vf+0.0005d-(ticks[t]*vf<0d?1:0));
			else tickmarks[t] = PApplet.nf((float) (ticks[t] * vf), 0, decimal).replace(",", ".");
			if(is_final && ticks[t]<0d && JPlot.supportLatex) tickmarks[t] = ltxMinus + tickmarks[t].substring(1);
			if(t%subdiv!=0) tickmarks[t] = "";
		}
		double lvf = Math.log10(p10);
		if(Math.abs(lvf)>2.9d) {
			int ivf = (int) (lvf+0.5d) - (lvf<0d ? -1 : 0);
			tickmarkFactor = JPlot.supportLatex ? "10$^{"+ivf+"}$" : "10^"+ivf;
		}
	}

}
