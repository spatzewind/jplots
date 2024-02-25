package jplots.layer;

import java.util.ArrayList;
import java.util.List;

import jplots.JPlot;
import jplots.axes.JSingleAxis;
import jplots.axes.JAxis;
import jplots.axes.JMultiAxis;
import jplots.shapes.JEllipseShape;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLatexShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JRectShape;
import jplots.shapes.JTextShape;
import processing.core.PConstants;
import processing.core.PGraphics;
import processing.core.PShape;

public class JLegend extends JPlotsLayer {

//	private JAxis srcAxis;
//	private boolean isHorizontal, borders;
	private double rts, set_pos_x, set_pos_y;
	private int horpos,vertpos;
	private List<LegendEntry> entries;

	public JLegend(JAxis parent) {
		this(parent, PConstants.RIGHT, PConstants.TOP, false, 1d);
	}
	public JLegend(JAxis parent, double relTextSize) {
		this(parent, PConstants.RIGHT, PConstants.TOP, false, relTextSize);
	}
	public JLegend(JAxis parent, int left_right_positioning, int top_bottom_positioning) {
		this(parent, left_right_positioning, top_bottom_positioning, false, 1d);
	}
	public JLegend(JAxis parent, int left_right_positioning, int top_bottom_positioning, boolean horizontal) {
		this(parent, left_right_positioning, top_bottom_positioning, horizontal, 1d);
	}
	public JLegend(JAxis parent, int left_right_positioning, int top_bottom_positioning, boolean horizontal, double relTextSize) {
		rts = relTextSize;
		horpos = left_right_positioning;
		vertpos = top_bottom_positioning;
		entries = new ArrayList<>();
	}

	public JLegend(JAxis parent, int pos_x, int pos_y, int width, int height, boolean horizontal) {
		this(parent, pos_x, pos_y, width, height, horizontal, 1d);
	}
	public JLegend(JAxis parent, int pos_x, int pos_y, int width, int height, boolean horizontal, double relTextSize) {
		rts = relTextSize;
		
		entries = new ArrayList<>();
	}
	
	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		entries.clear();
		ArrayList<JPlotsLayer> allLayers = new ArrayList<>();
		if(ax instanceof JSingleAxis) allLayers.addAll(((JSingleAxis)ax).getLayers());
		if(ax instanceof JMultiAxis) allLayers.addAll(((JMultiAxis)ax).getLayers());
		for (JPlotsLayer layer : allLayers) {
			if (layer instanceof JContourLayer) {
				// maybe do something with contourlayers -> hatching? -> CONTOURHATCHING_LABEL
			}
			// no JImageLayer
			// no other JLegend
			if (layer instanceof JAxisAlignedLineLayer)
				addEntry(layer.getLabel(), LegendEntry.LINEPLOT_LABEL, layer.getStyle(), layer.getLineColour());
			// ignore JPlotsLayer, because it is abstract
			if (layer instanceof JScatterLayer)
				addEntry(layer.getLabel(), LegendEntry.SCATTERPLOT_LABEL, layer.getStyle(), layer.getLineColour());
			if (layer instanceof JShapesLayer)
				addEntry(layer.getLabel(), LegendEntry.LINEPLOT_LABEL, layer.getStyle(), layer.getLineColour());
			if (layer instanceof JXYLayer)
				addEntry(layer.getLabel(), LegendEntry.LINEPLOT_LABEL, layer.getStyle(), layer.getLineColour());
		}
		for (int le = entries.size() - 1; le >= 0; le--)
			if (entries.get(le).getName().length() == 0)
				entries.remove(le);
		
		int[] p = ax.getSize();
		double labelWidth = 0d;
		JPlot plot = ax.getPlot();
		plot.getGraphic().textSize(200f);
		for (LegendEntry le : entries) {
			if(JPlot.supportLatex) {
				PShape ps = JLatexShape.toPShape(le.getName(), 200, 0xff000000, "\\text");
				labelWidth = Math.max(ps.width/200d, labelWidth);
			} else {
				labelWidth = Math.max(labelWidth, plot.getGraphic().textWidth(le.getName()) / 200d);
			}
		}
		double ts = ax.getTextSize() * rts;
		labelWidth += 3d;
		labelWidth *= ts;
		double toplefX = p[0], toplefY = p[1];
		switch(horpos) {
			case 0:                 toplefX = set_pos_x; break;
			case PConstants.LEFT:   break;
			case PConstants.CENTER: toplefX += 0.5d*(p[2]-labelWidth)-0.5d*ts; break;
			default:
			case PConstants.RIGHT:  toplefX += p[2] - labelWidth - ts; break;
		}
		switch(vertpos) {
			case 0:                 toplefY = set_pos_y; break;
			case PConstants.TOP:    break;
			case PConstants.CENTER: toplefY += 0.5d*(p[3]-1.5d*entries.size()*ts) - 0.5d*ts; break;
			default:
			case PConstants.BOTTOM: toplefY += p[3] - 1.5d*entries.size()*ts - ts; break;
		}
		
		JGroupShape lggs = new JGroupShape();
		lggs.addChild(new JRectShape((float) toplefX, (float) toplefY,
				(float) (toplefX + labelWidth + ts), (float) (toplefY + (1.5d * entries.size()+1d) * ts),
				(float) Math.min(labelWidth, 1.5d * entries.size()) * 0.1f, 0x7fffffff, 0x7f999999, 2f, true));
		toplefX += 0.5d*ts;
		toplefY += 0.5d*ts;
		for (int le = 0; le < entries.size(); le++)
			entries.get(le).toShape(lggs, (float) toplefX, (float) (toplefY + le * ts * 1.5d), (float) ts);
		s.addChild(lggs);
	}
	
	private void addEntry(String label, int type, String style, int colour) {
		boolean exists = false;
		for(int i=entries.size()-1; i>=0; i--) {
			LegendEntry le = entries.get(i);
			if(le.getName().equals(label) && le.getType()==type) {
				if(le instanceof MultiLegendEntry) {
					((MultiLegendEntry) le).addEntry(colour, style);
				} else {
					entries.remove(le);
					MultiLegendEntry mle = new MultiLegendEntry(label, type, le.getColour(), le.getStyle());
					mle.addEntry(colour, style);
					entries.add(mle);
				}
				exists = true;
				break;
			}
		}
		if(!exists)
			entries.add(new LegendEntry(label, type, colour, style));
	}

	private class LegendEntry {
		public final static int LINEPLOT_LABEL = 1;
		public final static int SCATTERPLOT_LABEL = 2;
		public final static int CONTOURHATCHING_LABEL = 3;
		
		protected String name;
		protected int type; // line, marker or hatching
		
		private int col;
		private String linestyle;
		
		public LegendEntry(String label, int plotType, int colour, String style) {
			name = label;
			type = plotType;
			col = colour;
			linestyle = style;
		}
		
		public String getName() {
			return name;
		}
		public int getType() {
			return type;
		}
		public String getStyle() {
			return linestyle;
		}
		public int getColour() {
			return col;
		}
		
		public void toShape(JGroupShape s, float x, float y, float ts) {
			switch (type) {
			case LINEPLOT_LABEL:
				drawLine(s, x, y + 0.667f * ts, 0.125f * ts, linestyle, col);
				break;
			case SCATTERPLOT_LABEL:
				drawSymbol(s, x + 1.250f * ts, y + 0.667f * ts, 0.125f * ts, linestyle, col);
				break;
			case CONTOURHATCHING_LABEL:
				drawHatching(s, x, y, 0.125f * ts);
				break;
			}
			JPlotShape.fill(0xff000000);
			if(JPlot.supportLatex)
				s.addChild(new JLatexShape(name, x+3*ts, y, ts, PConstants.LEFT, PConstants.TOP, 0xff000000, 0f, null));
			else
				s.addChild(new JTextShape(name, x+3*ts, y, ts, PConstants.LEFT, PConstants.TOP, 0xff000000, 0f, null));
		}
		
		protected void drawLine(JGroupShape s, float x, float y, float lw, String style, int color) {
			// System.out.println("[JLEGEND] draw Line: found linestyle \""+linestyle+"\"");
			double lln = 1d, llf = 0d, lpn = 0d, lpf = 0d, loff = 0d;
			if ("-".equals(style)) {
				lln = 1000 * lw;
				llf = 0;
				lpn = 0;
				lpf = 0;
			}
			if (".".equals(style)) {
				lln = 0;
				llf = 0;
				lpn = 1 * lw;
				lpf = 3 * lw;
			}
			if (",".equals(style)) {
				lln = 8 * lw;
				llf = 7 * lw;
				lpn = 0;
				lpf = 0;
			}
			if (";".equals(style)) {
				lln = 8 * lw;
				llf = 3 * lw;
				lpn = 1 * lw;
				lpf = 3 * lw;
			}
			int li = 0;
			float x1 = x, x2 = x + 20 * lw, y1 = y, y2 = y;
			double dx = x2 - x1, dy = y2 - y1;
			double l = Math.sqrt(dx * dx + dy * dy);
			dx /= l;
			dy /= l;
			double lpos = 0d, ldif = 0d;
			while (lpos < l) {
				ldif = 0d;
				switch (li) {
				case 0:
					if (lln == 0d)
						break;
					ldif = Math.min(l - lpos, lln - loff);
					break;
				case 1:
					if (llf == 0d)
						break;
					ldif = Math.min(l - lpos, llf - loff);
					break;
				case 2:
					if (lpn == 0d)
						break;
					ldif = Math.min(l - lpos, lpn - loff);
					break;
				case 3:
					if (lpf == 0d)
						break;
					ldif = Math.min(l - lpos, lpf - loff);
					break;
				}
				float xf1 = (float) (x1 + lpos * dx), yf1 = (float) (y1 + lpos * dy),
						xf2 = (float) (x1 + (lpos + ldif) * dx), yf2 = (float) (y1 + (lpos + ldif) * dy);
				// xf1=x1; xf2=x2; yf1=y1; yf2=y2;
				if (li % 2 == 0 && ldif > 0d)
					s.addChild(new JLineShape(lw, color, xf1, yf1, xf2, yf2));
				loff += ldif;
				switch (li) {
				case 0:
					if (loff >= lln) {
						loff -= lln;
						li = 1;
					}
					break;
				case 1:
					if (loff >= llf) {
						loff -= llf;
						li = 2;
					}
					break;
				case 2:
					if (loff >= lpn) {
						loff -= lpn;
						li = 3;
					}
					break;
				case 3:
					if (loff >= lpf) {
						loff -= lpf;
						li = 0;
					}
					break;
				}
				lpos += ldif;
			}
		}
		protected void drawSymbol(JGroupShape s, float x1, float y1, float lw, String style, int color) {
			System.out.println("[JLEGEND] draw Symbol: found symbol \"" + style + "\"");
			JPlotShape.stroke(color);
			JPlotShape.strokeWeight(lw);
			if ("o".equals(style) || "@".equals(style)) {
				JPlotShape.fill(color);
			} else {
				JPlotShape.noFill();
			}
			if ("()".equals(style) || "o".equals(style)) {
				s.addChild(new JEllipseShape(x1, y1, 6f * lw, 6f * lw));
			}
			if ("[]".equals(style) || "@".equals(style)) {
				s.addChild(new JRectShape(x1 - 2f * lw, y1 - 2f * lw, x1 + 2f * lw, y1 + 2f * lw));
			}
			if ("x".equals(style)) {
				s.addChild(new JLineShape(lw, lc, x1 - 2f * lw, y1 - 2f * lw, x1 + 2f * lw, y1 + 2f * lw));
				s.addChild(new JLineShape(lw, lc, x1 - 2f * lw, y1 + 2f * lw, x1 + 2f * lw, y1 - 2f * lw));
			}
			if ("+".equals(style)) {
				s.addChild(new JLineShape(lw, lc, x1 - 3f * lw, y1, x1 + 3f * lw, y1));
				s.addChild(new JLineShape(lw, lc, x1, y1 - 3f * lw, x1, y1 + 3f * lw));
			}
			if ("<>".equals(style)) {
				s.addChild(new JLineShape(lw, lc, x1, y1 - 3f * lw, x1 + 3f * lw, y1));
				s.addChild(new JLineShape(lw, lc, x1 + 3f * lw, y1, x1, y1 + 3f * lw));
				s.addChild(new JLineShape(lw, lc, x1, y1 + 3f * lw, x1 - 3f * lw, y1));
				s.addChild(new JLineShape(lw, lc, x1 - 3f * lw, y1, x1, y1 - 3f * lw));
			}
		}
		private void drawHatching(JGroupShape s, float x, float y, float t) {
			
		}
	}
	private class MultiLegendEntry extends LegendEntry {
		private int[] cols;
		private String[] linestyles;
		
		public MultiLegendEntry(String label, int plotType, int colour, String style) {
			super(label, plotType, colour, style);
			cols = new int[] {colour};
			linestyles = new String[] {style};
		}
		public void addEntry(int colour, String style) {
			int[] newcols = new int[cols.length+1];
			String[] newstyles = new String[linestyles.length+1];
			for(int i=0; i<cols.length; i++) {
				newcols[i] = cols[i];
				newstyles[i] = linestyles[i];
			}
			newcols[cols.length] = colour;
			newstyles[linestyles.length] = style;
			cols = newcols.clone();
			linestyles = newstyles.clone();
		}
		
		@Override
		public void toShape(JGroupShape s, float x, float y, float ts) {
			for(int i=0; i<cols.length; i++) {
				float frac = (i + 0.5f)/cols.length-0.5f;
				switch(type) {
					case LINEPLOT_LABEL:
						drawLine(s, x, y+(0.667f+frac)*ts, 0.125f*ts, linestyles[i], cols[i]);
						break;
					case SCATTERPLOT_LABEL:
						drawSymbol(s, x+(1.250f+2*frac)*ts, y+0.667f*ts, 0.125f*ts, linestyles[i], cols[i]);
						break;
					default:
						break;
				}
			}
			
			JPlotShape.fill(0xff000000);
			if(JPlot.supportLatex)
				s.addChild(new JLatexShape(name, x+3*ts, y, ts, PConstants.LEFT, PConstants.TOP, 0xff000000, 0f, null));
			else
				s.addChild(new JTextShape(name, x+3*ts, y, ts, PConstants.LEFT, PConstants.TOP, 0xff000000, 0f, null));
		}
	}
}
