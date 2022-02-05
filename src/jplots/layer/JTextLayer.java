package jplots.layer;

import jplots.JAxis;
import jplots.JPlot;
import jplots.shapes.JGroupShape;
import jplots.shapes.JTextShape;
import processing.core.PGraphics;

public class JTextLayer extends JPlotsLayer {
	
	private double x, y;
	private int x_align, y_align;
	private double rot;

	public JTextLayer(String message, double pos_x, double pos_y, double size, int textcolour, int x_align, int y_align, double rot) {
		this.label = message;
		this.x       = pos_x;
		this.y       = pos_y;
		this.lw      = size;
		this.lc      = textcolour;
		this.x_align = x_align;
		this.y_align = y_align;
		this.rot     = rot;
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		
//		entries.clear();
//		for(JPlotsLayer layer: ax.getLayers()) {
//			if(layer instanceof JContourLayer) {
//				//maybe do something with contourlayers -> hatching? -> CONTOURHATCHING_LABEL
//			}
//			//no JImageLayer
//			//no other JLegend
//			if(layer instanceof JLineLayer)
//				entries.add(new LegendEntry(layer.getLabel(), LegendEntry.LINEPLOT_LABEL, layer.getColour(), layer.getStyle()));
//			//ignore JPlotsLayer, because it is abstract
//			if(layer instanceof JScatterLayer)
//				entries.add(new LegendEntry(layer.getLabel(), LegendEntry.SCATTERPLOT_LABEL, layer.getColour(), layer.getStyle()));
//			if(layer instanceof JShapesLayer)
//				entries.add(new LegendEntry(layer.getLabel(), LegendEntry.LINEPLOT_LABEL, layer.getColour(), layer.getStyle()));
//			if(layer instanceof JXYLayer)
//				entries.add(new LegendEntry(layer.getLabel(), LegendEntry.LINEPLOT_LABEL, layer.getColour(), layer.getStyle()));
//		}
//		for(int le=entries.size()-1; le>=0; le--)
//			if(entries.get(le).getName().length()==0)
//				entries.remove(le);
//		
//		int[] p = ax.getSize();
//		double labelWidth = 0d;
//		JPlot plot = ax.getPlot();
//		plot.getGraphic().textSize(200f);
//		for(LegendEntry le: entries) {
//			labelWidth = Math.max(labelWidth, plot.getGraphic().textWidth(le.getName())/200d);
//		}
//		labelWidth += 3d;
//		labelWidth *= ax.getTextSize();
//		double toplefX = (p[0]+p[2]-labelWidth-0.5d*ax.getTextSize()), toplefY=p[1]+0.5d*ax.getTextSize();
		
		JGroupShape lggs = new JGroupShape();
//		lggs.addChild(new JRectShape(
//				(float)(toplefX-0.5d*ax.getTextSize()), (float)(toplefY-0.5d*ax.getTextSize()),
//				(float)(toplefX+labelWidth+0.5d*ax.getTextSize()), (float)(toplefY+(0.5d+1.5d*entries.size())*ax.getTextSize()),
//				(float)Math.min(labelWidth, 1.5d*entries.size())*0.1f, 0x3fffffff, 0x3f999999, 2f, true));
		lggs.addChild(new JTextShape(label, (float)x, (float)y, (float)(ax.getTextSize()*lw), x_align, y_align, lc, (float)rot));
		s.addChild(lggs);
	}
}
