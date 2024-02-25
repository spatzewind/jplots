package jplots.axes;

import jplots.maths.JPlotMath;

public class StaticTicksScale extends AxisScale {
	
	private double[] tick_values;
	private String[] tick_labels;
	
	public StaticTicksScale(JAxis subplot, char axis, double[] values, String[] labels) {
		super(subplot, axis);
		tick_values = values;
		tick_labels = labels;
	}
	

	@Override
	public double scale(double v) {
		return v;
	}
	@Override
	public double invscale(double v) {
		return v;
	}
	
	@Override
	public void create(double lower_bound, double upper_bound) {
		int[] p = subplot.getSize();
		int ps = axis=='x' ? p[0] : p[1];
		int pe = axis=='x' ? p[0]+p[2] : p[1]+p[3];
		
		ticks = tick_values.clone();
		tickmarks = tick_labels.clone();
		pos = JPlotMath.map(ticks, lower_bound,upper_bound, ps,pe);
	}
}
