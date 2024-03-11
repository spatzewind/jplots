package jplots.layer;

import jplots.JPlot;
import jplots.axes.JAxis;
import jplots.axes.LogarithmicScale;
import jplots.maths.AffineBuilder;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLatexShape;
import jplots.shapes.JTextShape;
import processing.core.PGraphics;

public class JTextLayer extends JPlotsLayer {

	private double x, y;
	private int x_align, y_align;
	private double rot;
	private String t_style;
	private boolean useAxis;
	
	public JTextLayer(boolean relative_to_this_axis, String message, double pos_x, double pos_y, double size, int textcolour, int x_align, int y_align,
			double rot, String style) {
		this.label = message;
		this.x = pos_x;
		this.y = pos_y;
		this.useAxis = relative_to_this_axis;
		this.lw = size;
		this.pc = textcolour;
		this.x_align = x_align;
		this.y_align = y_align;
		this.rot = rot;
		this.t_style = "\\text";
		if(style!=null) {
			if(style.equalsIgnoreCase("math")) this.t_style = "";
			if(style.equalsIgnoreCase("bb")) this.t_style = "\\mathbb";
			if(style.equalsIgnoreCase("fraktur")) this.t_style = "\\mathfrak";
		}
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		double qx = this.x;
		double qy = this.y;
		if(useAxis) {
			int[] p = ax.getSize();
			//TODO: deal with JMultiAxis
			boolean isXlog = (ax.getScaleX() instanceof LogarithmicScale);
			boolean isYlog = (ax.getScaleY() instanceof LogarithmicScale);
			double Xin = isXlog ? Math.log10(minX) : minX;
			double Xax = isXlog ? Math.log10(maxX) : maxX;
			double Yin = isYlog ? Math.log10(minY) : minY;
			double Yax = isYlog ? Math.log10(maxY) : maxY;
			double xs = p[2] / (Xax - Xin), ys = p[3] / (Yax - Yin);
			double[] xy = {
					isXlog ? Math.log10(this.x) : this.x,
					isYlog ? Math.log10(this.y) : this.y
			};
			if(ax.isGeoAxis()) {
				xy = inputProj.fromPROJtoLATLON(xy[0], xy[1], false, false);
				xy = ax.getGeoProjection().fromLATLONtoPROJ(xy[0], xy[1], false, false);
			}
			AffineBuilder affine = new AffineBuilder().scale(invertAxisX ? -1d : 1d, invertAxisY ? 1d : -1d)
					.translate(invertAxisX ? Xax : -Xin, invertAxisY ? -Yin : Yax).scale(xs, ys).translate(p[0], p[1]);
			double[][] affmat = affine.getMatrix();
			qx = affmat[0][0]*xy[0] + affmat[0][1]*xy[1] + affmat[0][2];
			qy = affmat[1][0]*xy[0] + affmat[1][1]*xy[1] + affmat[1][2];
		} else {
			int[] p = ax.getPlot().getSize();
			qx *= p[0];
			qy *= p[1];
		}
		if(JPlot.supportLatex) {
			s.addChild(new JLatexShape(label, (float) qx, (float) qy, (float) (ax.getTextSize() * lw), x_align, y_align, pc, (float) rot, t_style));
		} else {
			s.addChild(new JTextShape(label, (float) qx, (float) qy, (float) (ax.getTextSize() * lw), x_align, y_align, pc, (float) rot, t_style));
		}
	}
}
