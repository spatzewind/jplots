package pplots.transform;

import pplots.PAxis;
import pplots.shapes.PGroupShape;
import pplots.shapes.PLineShape;
import pplots.shapes.PPlotShape;

public class PRectangleProjection implements PProjection {
	
	private double xs,xe,ys,ye;
	
	public PRectangleProjection(double xstart, double xend, double ystart, double end) {
	}

	@Override
	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree) {
		if(Double.isNaN(x) || Double.isNaN(y)) return new double[] {Double.NaN, Double.NaN};
		if(x<xs || x>xe) return new double[] {Double.NaN, Double.NaN};
		if(y<ys || y>ye) return new double[] {Double.NaN, Double.NaN};
		return new double[] {(x-xs)/(xe-xs),(y-ys)/(ye-ys)};
	}

	@Override
	public double[] fromLATLONtoPROJ(double u, double v, boolean input_in_degree) {
		if(Double.isNaN(u) || Double.isNaN(v)) return new double[] {Double.NaN, Double.NaN};
		if(u<0d || u>1d) return new double[] {Double.NaN, Double.NaN};
		if(v<0d || v>1d) return new double[] {Double.NaN, Double.NaN};
		return new double[] {xs+u*(xe-xs),ys+v*(ye-ys)};
	}
	
	@Override
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree) {
		return new double[] {(xe-xs)/(ye-ys),0d, 0d,(ye-ys)/(xe-xs)};
	}
	
	@Override
	public double[] tissotFromProj(double x, double y) {
		return new double[] {(xe-xs)/(ye-ys),0d, 0d,(ye-ys)/(xe-xs)};
	}
	
	@Override
	public double[] defaultMapExtend() {
		return new double[] {xs,xe,ys,ye};
	}
	
	@Override
	public void drawBorder(PAxis ax, PGroupShape s) {
		PPlotShape.stroke(0xff000000); PPlotShape.strokeWeight(3f);
		int[] p = ax.getSize();
		s.addChild(new PLineShape(p[0],     p[1],     p[0]+p[2],p[1]));
		s.addChild(new PLineShape(p[0],     p[1]+p[3],p[0]+p[2],p[1]+p[3]));
		s.addChild(new PLineShape(p[0],     p[1],     p[0],     p[1]+p[3]));
		s.addChild(new PLineShape(p[0]+p[2],p[1],     p[0]+p[2],p[1]+p[3]));
	}
}
