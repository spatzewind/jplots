package jplots.layer;

import jplots.JAxis;
import jplots.JPlot;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import processing.core.PGraphics;

public class JLineLayer extends JPlotsLayer {
	
	private double value;
	private boolean isHorizontal;
	private int col;
	
	public JLineLayer(float v, char orientation, int colour, float linewidth, String linestyle) {
		value = v;
		isHorizontal = (orientation=='h');
		col = colour;
		lw = linewidth;
		ls = linestyle;
	}
	public JLineLayer(double v, char orientation, int colour, double linewidth, String linestyle) {
		value = v;
		isHorizontal = (orientation=='h');
		col = colour;
		lw = linewidth;
		ls = linestyle;
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		double Xin = ax.isXlogAxis()?Math.log10(minX):minX, Xax = ax.isXlogAxis()?Math.log10(maxX):maxX;
		double Yin = ax.isYlogAxis()?Math.log10(minY):minY, Yax = ax.isYlogAxis()?Math.log10(maxY):maxY;
		JGroupShape xyShape = new JGroupShape();
		JPlotShape.stroke(col); JPlotShape.strokeWeight((float)lw);
		double lln=1d, llf=0d, lpn=0d, lpf=0d, loff = 0d;
		if("-".equals(ls)) { lln=1000*lw; llf=0; lpn=0; lpf=0; }
		if(".".equals(ls)) { lln=0; llf=0; lpn=1*lw; lpf=3*lw; }
		if(",".equals(ls)) { lln=8*lw; llf=7*lw; lpn=0; lpf=0; }
		if(";".equals(ls)) { lln=8*lw; llf=3*lw; lpn=1*lw; lpf=3*lw; }
		double xs = p[2]/(Xax-Xin), ys = p[3]/(Yax-Yin);
		int li = 0;
		double vv,x1,x2,y1,y2;
		if(isHorizontal) {
			vv = ax.isYlogAxis() ? Math.log10(value) : value;
			x1 = p[0]; x2 = p[0]+p[2];
			y1 = p[1]+ys*(invertAxisY ? vv-Yin : Yax-vv);
			y2 = y1;
		} else {
			vv = ax.isXlogAxis() ? Math.log10(value) : value;
			x1 = p[0]+xs*(invertAxisX ? Xax-vv : vv-Xin);
			x2 = x1;
			y1 = p[1]; y2 = p[1]+p[3];
		}
		double dx = x2-x1, dy = y2-y1;
		double l = Math.sqrt(dx*dx+dy*dy);
		dx /= l; dy /= l;
		double lpos = 0d, ldif = 0d;
		while(lpos<l) {
			ldif = 0d;
			switch(li) {
				case 0: if(lln==0d) break; ldif = Math.min(l-lpos, lln-loff); break;
				case 1: if(llf==0d) break; ldif = Math.min(l-lpos, llf-loff); break;
				case 2: if(lpn==0d) break; ldif = Math.min(l-lpos, lpn-loff); break;
				case 3: if(lpf==0d) break; ldif = Math.min(l-lpos, lpf-loff); break;
			}
			float xf1 = (float)(x1+lpos*dx), yf1 = (float)(y1+lpos*dy),
					  xf2 = (float)(x1+(lpos+ldif)*dx), yf2 = (float)(y1+(lpos+ldif)*dy);
			if(xf1<p[0]      && xf2>=p[0])      { yf1 = JPlotMath.map(p[0],      xf1, xf2, yf1, yf2); }
			if(xf1>p[0]+p[2] && xf2<=p[0]+p[2]) { yf1 = JPlotMath.map(p[0]+p[2], xf1, xf2, yf1, yf2); }
			if(xf2<p[0]      && xf1>=p[0])      { yf2 = JPlotMath.map(p[0],      xf1, xf2, yf1, yf2); }
			if(xf2>p[0]+p[2] && xf1<=p[0]+p[2]) { yf2 = JPlotMath.map(p[0]+p[2], xf1, xf2, yf1, yf2); }
			if(yf1<p[1]      && yf2>=p[1])      { xf1 = JPlotMath.map(p[1],      yf1, yf2, xf1, xf2); }
			if(yf1>p[1]+p[3] && yf2<=p[1]+p[3]) { xf1 = JPlotMath.map(p[1]+p[3], yf1, yf2, xf1, xf2); }
			if(yf2<p[1]      && yf1>=p[1])      { xf2 = JPlotMath.map(p[1],      yf1, yf2, xf1, xf2); }
			if(yf2>p[1]+p[3] && yf1<=p[1]+p[3]) { xf2 = JPlotMath.map(p[1]+p[3], yf1, yf2, xf1, xf2); }
			if(xf1>=p[0] && xf1<=p[0]+p[2] && xf2>=p[0] && xf2<=p[0]+p[2] &&
					yf1>=p[1] && yf1<=p[1]+p[3] && yf2>=p[1] && yf2<=p[1]+p[3]) {
				if(li%2==0 && ldif>0d)
					xyShape.addChild(new JLineShape(xf1,yf1,xf2,yf2));
				loff += ldif;
				switch(li) {
					case 0: if(loff>=lln) {loff -= lln; li=1; } break;
					case 1: if(loff>=llf) {loff -= llf; li=2; } break;
					case 2: if(loff>=lpn) {loff -= lpn; li=3; } break;
					case 3: if(loff>=lpf) {loff -= lpf; li=0; } break;
				}
			}
			lpos += ldif;
		}
		s.addChild(xyShape);
	}

}
