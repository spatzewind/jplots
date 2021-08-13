package jplots.color;

import jplots.JAxis;
import jplots.layer.JContourLayer;
import jplots.layer.JPlotsLayer;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JImageShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JTextShape;
import processing.core.PApplet;
import processing.core.PImage;

public class JColourbar extends JAxis {
	
	private JAxis srcAxis;
	private JColourtable srcColT;
	private boolean isHorizontal, borders;
	private double minC;
	private double maxC;
	private double minV;
	private double maxV;
	private double[] contourLevels;
	private String titleC;
	private PImage img;
	
	public JColourbar(JAxis parent) {
		this(parent,""); }
	public JColourbar(JAxis parent, String colorbar_title) {
		this(parent,
				(int)(parent.getSize()[0]+1.04d*parent.getSize()[2]+0.5d),parent.getSize()[1],
				(int)(0.08d*parent.getSize()[2]+0.5d),parent.getSize()[3], colorbar_title, false);
	}
	public JColourbar(JAxis parent, int pos_x, int pos_y, int width, int height, String colorbar_title, boolean horizontal) {
		super(parent.getPlot(), pos_x, pos_y, width, height);
		srcAxis = parent;
		srcColT = null;
		titleC = colorbar_title;
		isHorizontal = horizontal;
		minC = Double.NaN;
		maxC = Double.NaN;
		contourLevels = new double[0];
		borders = true;
	}

	@Override
	public JGroupShape createPlot(PApplet applet, int w, int h) {
		if(pplot.isDebug()) {
			System.out.println("[DEBUG] JColourbar: begin colourbar\n"+
					"                    x/y/width/height: "+px+"/"+py+"/"+pw+"/"+ph);
		}
		JGroupShape graph = new JGroupShape();
		minV = Double.NaN;
		maxV = Double.NaN;
		if(srcAxis!=null) {
			if(pplot.isDebug())
				System.out.println("[DEBUG] JColourbar: source axis exist.");
			for(JPlotsLayer layer: srcAxis.getLayers())
				if(layer instanceof JContourLayer) {
					if(pplot.isDebug())
						System.out.println("[DEBUG] JColourbar: found JContourLayer: "+layer);
					double[] zr = ((JContourLayer) layer).getZRange();
					minV = zr[0];
					maxV = zr[1];
					contourLevels = ((JContourLayer) layer).getLevels();
					if(Double.isNaN(minC)) minC = minV;
					if(Double.isNaN(maxC)) maxC = maxV;
					srcColT = layer.getColourtable();
					if(pplot.isDebug())
						System.out.println("[DEBUG] JColourbar: set vmin/vmax:  "+minV+"/"+maxV);
					break;
				}
		}
		if(Double.isNaN(minV)) minV = minC;
		if(Double.isNaN(maxV)) maxV = maxC;
		if(srcColT==null || Double.isNaN(minV) || Double.isNaN(maxV)) {
			System.err.println("[ERR] JColourbar: no source for colourbar found!");
			return graph;
		}
		if(pplot.isDebug())
			System.out.println("[DEBUG] JColourbar: use\n"+
					"                    cmin/cmax = "+minC+"/"+maxC+"\n"+
					"                    vmin/vmax = "+minV+"/"+maxV+"\n"+
					"                    colourtable "+srcColT);
		if(isHorizontal) graph.addChild(createXAxis());
		if(!isHorizontal) graph.addChild(createYAxis());
		graph.addChild(createBar(applet));
		if(borders) {
			graph.addChild(new JLineShape(px, py,    px+pw, py   ));
			graph.addChild(new JLineShape(px, py+ph, px+pw, py+ph));
			graph.addChild(new JLineShape(px,    py, px,    py+ph));
			graph.addChild(new JLineShape(px+pw, py, px+pw, py+ph));
		}
		return graph;
	}
	
	private JGroupShape createXAxis() {
		JGroupShape axisgrid = new JGroupShape();
		double[] ticks = JPlotMath.optimalLinearTicks(minC, maxC);
		double[] tcpos = JPlotMath.dlerp(ticks,minC,maxC,px,px+pw);
		String[] tickmark = new String[ticks.length];
		double vf = 1d/(ticks[0]);
		int decimal = (int) (1000d*ticks[1]+0.5d);
		decimal = decimal%100==0 ? 1 : decimal%10==0 ? 2 : 3;
		double tmlen = 0d, tmlc = 0d;
		pplot.getGraphic().textSize(200);
		for(int t=0; t<ticks.length; t++) {
			tickmark[t] = PApplet.nf((float)(ticks[t]*vf),0,decimal);
			if(t>1) {
				tmlen += pplot.getGraphic().textWidth(tickmark[t]) / 200f;
				tmlc += 1d;
			}
		}
		tmlen *= this.txtsize / tmlc;
		int tickcount = Math.max(2, (int) (pw/(1.2d*tmlen)+0.9999d));
		if(pplot.isDebug())
			System.out.println("[DEBUG] JAxis-object: tmlen="+tmlen+" -> tickcount approx. "+tickcount);
		ticks = JPlotMath.optimalLinearTicks(minC, maxC, tickcount);
		tcpos = JPlotMath.dlerp(ticks,minC,maxC,px,px+pw);
		tickmark = new String[ticks.length];
		vf = 1d/(ticks[0]); decimal = (int) (1000d*ticks[1]+0.5d);
		decimal = decimal%100==0 ? 1 : decimal%10==0 ? 2 : 3;
		for(int t=0; t<ticks.length; t++)
			tickmark[t] = PApplet.nf((float)(ticks[t]*vf),0,decimal);
		if(pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for(int t=2; t<ticks.length; t++) {
				tickStr += ", "+PApplet.nf((float)ticks[t],0,2);
				posStr  += ", "+PApplet.nf((float)tcpos[t],0,2);
			}
			System.out.println("[DEBUG] JAxis-object: Xtickfactors={p10: "+ticks[0]+", f: "+ticks[1]+"}");
			System.out.println("[DEBUG] JAxis-object: Xtickval={"+tickStr.substring(2)+"}");
			System.out.println("[DEBUG] JAxis-object: Xtickpos={"+posStr.substring(2)+"}");
		}
		JPlotShape.stroke(0xff000000); JPlotShape.strokeWeight(2f);
		for(int t=2; t<ticks.length; t++)
			if(ticks[t]>=Math.min(minC, maxC) && ticks[t]<=Math.max(minC, maxC)) {
				axisgrid.addChild(new JLineShape((float)tcpos[t],py+ph,(float)tcpos[t],py+1.125f*ph));
				//axisgrid.addChild(ap.createShape(PShape.TEXT, "H", (float)tcpos[t],py-0.1f*ph,(float)tcpos[t],py));
				axisgrid.addChild(new JTextShape(tickmark[t], (float)tcpos[t], py+1.188f*ph, (float)txtsize, PApplet.CENTER, PApplet.TOP, 0xff000000, 0));
			}
		if(titleC.length()>0)
			axisgrid.addChild(new JTextShape(titleC, px+0.5f*pw, py+1.250f*ph+(float)txtsize, (float)(1.1d*txtsize), PApplet.CENTER, PApplet.TOP, 0xff000000, 0));
		return axisgrid;
	}
	private JGroupShape createYAxis() {
		JGroupShape axisgrid = new JGroupShape();
		double[] ticks = JPlotMath.optimalLinearTicks(minC, maxC);
		double[] tcpos = JPlotMath.dlerp(ticks,minC,maxC,py+ph,py);
//		for(int t=0; t<tcpos.length; t++)
//			tcpos[t] = 2*py+ph-tcpos[t];
		if(pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for(int t=2; t<ticks.length; t++) {
				tickStr += ", "+PApplet.nf((float)ticks[t],0,2);
				posStr  += ", "+PApplet.nf((float)tcpos[t],0,2);
			}
			System.out.println("[DEBUG] JAxis-object: Ytickfactors={p10: "+ticks[0]+", f: "+ticks[1]+"}");
			System.out.println("[DEBUG] JAxis-object: Ytickval={"+tickStr.substring(2)+"}");
			System.out.println("[DEBUG] JAxis-object: Ytickpos={"+posStr.substring(2)+"}");
		}
		JPlotShape.stroke(0xff000000); JPlotShape.strokeWeight(2f);
		float tw = 0f;
		double vf = 1d/(ticks[0]);
		int decimal = (int) (1000d*ticks[1]+0.5d);
		decimal = decimal%100==0 ? 1 : decimal%10==0 ? 2 : 3;
		for(int t=2; t<ticks.length; t++)
			if(ticks[t]>=Math.min(minC, maxC) && ticks[t]<=Math.max(minC, maxC)) {
				axisgrid.addChild(new JLineShape(px+pw,(float)tcpos[t],px+1.250f*pw,(float)tcpos[t]));
				String tmstr = PApplet.nf((float)(ticks[t]*vf),0,decimal);
				tw = Math.max(tw, (float)txtsize * pplot.getGraphic().textWidth(tmstr)/pplot.getGraphic().textSize);
				axisgrid.addChild(new JTextShape(tmstr, px+1.375f*pw, (float)tcpos[t], (float)txtsize,
						PApplet.LEFT, PApplet.CENTER, 0xff000000, 0));
			}
		if(titleC.length()>0)
			axisgrid.addChild(new JTextShape(titleC, px+1.500f*pw+tw, py+0.5f*ph, (float)(1.1d*txtsize), PApplet.CENTER, PApplet.TOP, 0xff000000, JPlotShape.ROTATE_COUNTERCLOCKWISE));
		return axisgrid;
	}
	private JImageShape createBar(PApplet applet) {
		if(img==null) {
			img = applet.createImage(pw, ph, PApplet.ARGB);
		} else if(img.width!=pw || img.height!=ph) {
			img = applet.createImage(pw, ph, PApplet.ARGB);
		}
		img.loadPixels();
		double minCI = Double.NaN, maxCI = Double.NaN;
		if(contourLevels.length>0) {
			minCI = contourLevels[0]<1d ? contourLevels[0]*2d : contourLevels[0]-10d;
			maxCI = contourLevels[contourLevels.length-1]>1d ? contourLevels[contourLevels.length-1]*2d : contourLevels[contourLevels.length-1]+10d;
		}
		if(isHorizontal) {
			for(int i=0; i<pw; i++) {
				double v = JPlotMath.dlerp(i, 0, pw, minC, maxC);
				if(contourLevels.length>0) {
					int l = JContourLayer.getLevel(v, contourLevels, -1);
					v = l<1 ? minCI : l>contourLevels.length-1 ? maxCI : 0.5d*(contourLevels[l-1]+contourLevels[l]);
				}
				int c = srcColT.getColour(v, minV, maxV);
				for(int j=0; j<ph; j++)
					img.pixels[j*pw+i] = c;
			}
		} else {
			for(int j=0; j<ph; j++) {
				double v = JPlotMath.dlerp(j, 0, ph, maxC, minC);
				if(contourLevels.length>0) {
					int l = JContourLayer.getLevel(v, contourLevels, -1);
					v = l<1 ? minCI : l>contourLevels.length-1 ? maxCI : 0.5d*(contourLevels[l-1]+contourLevels[l]);
				}
				int c = srcColT.getColour(v, minV, maxV);
				for(int i=0; i<pw; i++)
					img.pixels[j*pw+i] = c;
			}
		}
		img.updatePixels();
		return new JImageShape(img, px, py, pw, ph);
	}
}
