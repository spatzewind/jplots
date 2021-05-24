package pplots.transform;

import pplots.PAxis;
import pplots.shapes.PGroupShape;
import pplots.shapes.PLineShape;
import pplots.shapes.PPlotShape;

public class PIdentityProjection implements PProjection {
	
	public PIdentityProjection() {
	}

	@Override
	public double[] fromPROJtoLATLON(double x, double y, boolean output_in_degree) {
		return new double[] {x,y};
	}

	@Override
	public double[] fromLATLONtoPROJ(double u, double v, boolean input_in_degree) {
		return new double[] {u,v};
	}
	
	@Override
	public double[] tissotFromLatLon(double u, double v, boolean input_in_degree) {
		return new double[] {1d,0d, 0d,1d};
	}
	
	@Override
	public double[] tissotFromProj(double x, double y) {
		return new double[] {1d,0d, 0d,1d};
	}
	
	@Override
	public double[] defaultMapExtend() {
		return new double[] {-1d,1d,-1d,1d};
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
