package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import jplots.JAxis;
import jplots.JPlot;
import jplots.colour.JColourtable;
import jplots.maths.JPlotMath;
import jplots.shapes.JEllipseShape;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JRectShape;
import jplots.shapes.JTextShape;
import processing.core.PApplet;
import processing.core.PGraphics;

public class JLegend extends JPlotsLayer {

	
	private JAxis srcAxis;
	private boolean isHorizontal, borders;
	private double rts;
	private List<LegendEntry> entries;

	public JLegend(JAxis parent) {
		this(parent,
				(int)(parent.getSize()[0]+1.04d*parent.getSize()[2]+0.5d),parent.getSize()[1],
				(int)(0.08d*parent.getSize()[2]+0.5d),parent.getSize()[3],
				false, 1d);
	}
	public JLegend(JAxis parent, double relTextSize) {
		this(parent,
				(int)(parent.getSize()[0]+1.04d*parent.getSize()[2]+0.5d),parent.getSize()[1],
				(int)(0.08d*parent.getSize()[2]+0.5d),parent.getSize()[3],
				false, relTextSize);
	}
	public JLegend(JAxis parent, int pos_x, int pos_y, int width, int height, boolean horizontal) {
		this(parent, pos_x, pos_y, width, height, horizontal, 1d);
	}
	public JLegend(JAxis parent, int pos_x, int pos_y, int width, int height, boolean horizontal, double relTextSize) {
		srcAxis = parent;
		isHorizontal = horizontal;
		borders = true;
		rts = relTextSize;
		
		entries = new ArrayList<JLegend.LegendEntry>();
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		entries.clear();
		for(JPlotsLayer layer: ax.getLayers()) {
			if(layer instanceof JContourLayer) {
				//maybe do something with contourlayers -> hatching? -> CONTOURHATCHING_LABEL
			}
			//no JImageLayer
			//no other JLegend
			if(layer instanceof JLineLayer)
				entries.add(new LegendEntry(layer.getLabel(), LegendEntry.LINEPLOT_LABEL, layer.getLineColour(), layer.getStyle()));
			//ignore JPlotsLayer, because it is abstract
			if(layer instanceof JScatterLayer)
				entries.add(new LegendEntry(layer.getLabel(), LegendEntry.SCATTERPLOT_LABEL, layer.getLineColour(), layer.getStyle()));
			if(layer instanceof JShapesLayer)
				entries.add(new LegendEntry(layer.getLabel(), LegendEntry.LINEPLOT_LABEL, layer.getLineColour(), layer.getStyle()));
			if(layer instanceof JXYLayer)
				entries.add(new LegendEntry(layer.getLabel(), LegendEntry.LINEPLOT_LABEL, layer.getLineColour(), layer.getStyle()));
		}
		for(int le=entries.size()-1; le>=0; le--)
			if(entries.get(le).getName().length()==0)
				entries.remove(le);
		
		int[] p = ax.getSize();
		double labelWidth = 0d;
		JPlot plot = ax.getPlot();
		plot.getGraphic().textSize(200f);
		for(LegendEntry le: entries) {
			labelWidth = Math.max(labelWidth, plot.getGraphic().textWidth(le.getName())/200d);
		}
		double ts = ax.getTextSize() * rts;
		labelWidth += 3d;
		labelWidth *= ts;
		double toplefX = (p[0]+p[2]-labelWidth-0.5d*ts), toplefY=p[1]+0.5d*ts;
		
		JGroupShape lggs = new JGroupShape();
		lggs.addChild(new JRectShape(
				(float)(toplefX-0.5d*ts), (float)(toplefY-0.5d*ts),
				(float)(toplefX+labelWidth+0.5d*ts), (float)(toplefY+(0.5d+1.5d*entries.size())*ts),
				(float)Math.min(labelWidth, 1.5d*entries.size())*0.1f, 0x3fffffff, 0x3f999999, 2f, true));
		for(int le=0; le<entries.size(); le++)
			entries.get(le).toShape(lggs, (float)toplefX, (float)(toplefY+le*ts*1.5d), (float)ts);
		s.addChild(lggs);
	}

	public class LegendEntry {
		public final static int LINEPLOT_LABEL = 1;
		public final static int SCATTERPLOT_LABEL = 2;
		public final static int CONTOURHATCHING_LABEL = 3;
		
		private String name;
		private int col;
		private int type; //line, marker or hatching
		private String linestyle;
		
		public LegendEntry(String label, int plotType, int colour, String style) {
			name      = label;
			type      = plotType;
			col       = colour;
			linestyle = style;
		}
		
		public String getName() {
			return name;
		}
		public int getColour() {
			return col;
		}
		public int getType() {
			return type;
		}
		public String getStyle() {
			return linestyle;
		}
		
		public void toShape(JGroupShape s, float x, float y, float ts) {
			switch(type) {
				case LINEPLOT_LABEL: drawLine(s,x,y+0.667f*ts, 0.125f*ts); break;
				case SCATTERPLOT_LABEL: drawSymbol(s,x+1.250f*ts,y+0.667f*ts, 0.125f*ts); break;
				case CONTOURHATCHING_LABEL: drawHatching(s,x,y, 0.125f*ts); break;
			}
			JPlotShape.fill(0xff000000);
			s.addChild(new JTextShape(name, x+3*ts, y, ts, PApplet.LEFT, PApplet.TOP, 0xff000000, 0f));
		}
		
		private void drawLine(JGroupShape s, float x, float y, float lw) {
			System.out.println("[JLEGEND] draw Line:   found linestyle \""+linestyle+"\"");
			double lln=1d, llf=0d, lpn=0d, lpf=0d, loff = 0d;
			if("-".equals(linestyle)) { lln=1000*lw; llf=0; lpn=0; lpf=0; }
			if(".".equals(linestyle)) { lln=0; llf=0; lpn=1*lw; lpf=3*lw; }
			if(",".equals(linestyle)) { lln=8*lw; llf=7*lw; lpn=0; lpf=0; }
			if(";".equals(linestyle)) { lln=8*lw; llf=3*lw; lpn=1*lw; lpf=3*lw; }
			int li = 0;
			float x1=x, x2=x+20*lw, y1=y, y2=y;
			double dx = x2-x1, dy = y2-y1;
			double l = Math.sqrt(dx*dx+dy*dy);
			dx /= l; dy /= l;
			double lpos = 0d, ldif = 0d;
			JPlotShape.stroke(col);
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
				//xf1=x1; xf2=x2; yf1=y1; yf2=y2;
				if(li%2==0 && ldif>0d)
					s.addChild(new JLineShape(xf1,yf1,xf2,yf2));
				loff += ldif;
				switch(li) {
					case 0: if(loff>=lln) {loff -= lln; li=1; } break;
					case 1: if(loff>=llf) {loff -= llf; li=2; } break;
					case 2: if(loff>=lpn) {loff -= lpn; li=3; } break;
					case 3: if(loff>=lpf) {loff -= lpf; li=0; } break;
				}
				lpos += ldif;
			}
		}
		private void drawSymbol(JGroupShape s, float x1, float y1, float lw) {
			System.out.println("[JLEGEND] draw Symbol: found symbol \""+linestyle+"\"");
			JPlotShape.stroke(col); JPlotShape.strokeWeight(lw);
			if("o".equals(linestyle) || "@".equals(linestyle)) {
				JPlotShape.fill(col); } else { JPlotShape.noFill(); }
			if("()".equals(linestyle) || "o".equals(linestyle)) {
				s.addChild(new JEllipseShape(x1, y1, 6f*lw, 6f*lw)); }
			if("[]".equals(linestyle) || "@".equals(linestyle)) {
				s.addChild(new JRectShape(x1-2f*lw, y1-2f*lw, x1+2f*lw, y1+2f*lw)); }
			if("x".equals(linestyle)) {
				s.addChild(new JLineShape(x1-2f*lw, y1-2f*lw, x1+2f*lw, y1+2f*lw));
				s.addChild(new JLineShape(x1-2f*lw, y1+2f*lw, x1+2f*lw, y1-2f*lw)); }
			if("+".equals(linestyle)) {
				s.addChild(new JLineShape(x1-3f*lw, y1, x1+3f*lw, y1));
				s.addChild(new JLineShape(x1, y1-3f*lw, x1, y1+3f*lw)); }
			if("<>".equals(linestyle)) {
				s.addChild(new JLineShape(x1, y1-3f*lw, x1+3f*lw, y1));
				s.addChild(new JLineShape(x1+3f*lw, y1, x1, y1+3f*lw));
				s.addChild(new JLineShape(x1, y1+3f*lw, x1-3f*lw, y1));
				s.addChild(new JLineShape(x1-3f*lw, y1, x1, y1-3f*lw)); }
		}
		private void drawHatching(JGroupShape s, float x, float y, float t) {
			
		}
	}
}
