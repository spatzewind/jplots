package jplots.colour;

import java.util.ArrayList;
import java.util.List;

import jplots.JAxis;
import jplots.JPlot;
import jplots.layer.JContourLayer;
import jplots.layer.JContourLayer2D;
import jplots.layer.JHatchLayer;
import jplots.layer.JPlotsLayer;
import jplots.layer.JXYLayer;
import jplots.maths.JDPoint;
import jplots.maths.JDTriangle;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JImageShape;
import jplots.shapes.JLatexShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JTextShape;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PImage;

public class JColourbar extends JAxis {

	private JAxis srcAxis;
	private JColourtable srcColT;
	private boolean isHorizontal, borders, foundPrimary, extendLower, extendUpper;
	private double minC, maxC, minH, maxH, minV, maxV;
	private double[] contourLevels;
	private String titleC;
	private PImage img;
	private List<JHatchLayer> hatchPatterns;
	
	public JColourbar(JAxis parent) {
		this(parent, "");
	}

	public JColourbar(JAxis parent, String colorbar_title) {
		this(parent, (int) (parent.getSize()[0] + 1.04d * parent.getSize()[2] + 0.5d), parent.getSize()[1],
				(int) (0.08d * parent.getSize()[2] + 0.5d), parent.getSize()[3], colorbar_title, false);
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
		extendLower = false;
		extendUpper = false;
		hatchPatterns = new ArrayList<>();
	}

	@Override
	public JGroupShape createPlot(PApplet applet, int w, int h) {
		if (pplot.isDebug()) {
			System.out.println("[DEBUG] JColourbar: begin colourbar\n"+
					"                    x/y/width/height: "+px+"/"+py+"/"+pw+"/"+ph);
		}
		int trihgt = Math.min(pw,ph);
		int xs = px, ys = py, xe = px+pw, ye = py+ph;
		if(isHorizontal) {
			if(extendLower) xs += trihgt;
			if(extendUpper) xe -= trihgt;
		} else {
			if(extendLower) ye -= trihgt;
			if(extendUpper) ys += trihgt;
		}
		JGroupShape graph = new JGroupShape();
		if(srcAxis==null)
			return graph;
		minV = Double.NaN;
		maxV = Double.NaN;
		foundPrimary = false;
		minH = Double.POSITIVE_INFINITY;
		maxH = Double.NEGATIVE_INFINITY;
		if (pplot.isDebug())
			System.out.println("[DEBUG] JColourbar: source axis exist.");
		collectHatchPatterns();
		collectPrimarySource();
		if(!foundPrimary && hatchPatterns.isEmpty()) {
			System.err.println("[ERR] JColourbar: Could not find source to create JColourbar.");
			return graph;
		}
		if(!foundPrimary) {
			minV = minH;
			maxV = maxH;
		}
		if (Double.isNaN(minV) || Double.isNaN(maxV)) {
			System.err.println("[ERR] JColourbar: no source for colourbar found!");
			return graph;
		}

//		System.out.println("Reloaded datarange: [ z={"+minV+" ... "+maxV+"}");
		
		if (pplot.isDebug())
			System.out.println("[DEBUG] JColourbar: use\n"+
					"                    cmin/cmax = "+minC+"/"+maxC+"\n"+
					"                    vmin/vmax = "+minV+"/"+maxV+"\n"+
					"                    colourtable "+srcColT);
		if (isHorizontal)
			graph.addChild(createXAxis());
		if (!isHorizontal)
			graph.addChild(createYAxis());
		JDPoint lt = new JDPoint(xs, ys);
		JDPoint rt = new JDPoint(xe, ys);
		JDPoint lb = new JDPoint(xs, ye);
		JDPoint rb = new JDPoint(xe, ye);
		List<JDTriangle> triangles = new ArrayList<>();
		if(isHorizontal) {
			lt.value = minV; rt.value = maxV;
			lb.value = minV; rb.value = maxV;
		} else {
			lt.value = maxV; rt.value = maxV;
			lb.value = minV; rb.value = minV;
		}
		triangles.add(new JDTriangle(lt,rt,lb));
		triangles.add(new JDTriangle(lb,rt,rb));
		if(extendLower) {
			double low_val = 1.5d*minV - 0.5d*maxV;
			if(isHorizontal) {
				triangles.add(new JDTriangle(new JDPoint(px,py+0.5d*ph,low_val), new JDPoint(xs,ys,low_val), new JDPoint(xs,ye,low_val)));
			} else {
				triangles.add(new JDTriangle(new JDPoint(px+0.5d*pw,py+ph,low_val), new JDPoint(xs,ye,low_val), new JDPoint(xe,ye,low_val)));
			}
		}
		if(extendUpper) {
			double upp_val = 1.5d*maxV - 0.5d*minV;
			if(isHorizontal) {
				triangles.add(new JDTriangle(new JDPoint(px+pw,py+0.5d*ph,upp_val), new JDPoint(xe,ye,upp_val), new JDPoint(xe,ys,upp_val)));
			} else {
				triangles.add(new JDTriangle(new JDPoint(px+0.5d*pw,py,upp_val), new JDPoint(xe,ys,upp_val), new JDPoint(xs,ys,upp_val)));
			}
		}
		//drawPrimarySource(graph);
		drawPrimarySource(applet, graph, triangles);
		drawHatchLayers(graph, triangles);
		
		if (borders) {
			if(isHorizontal) {
				graph.addChild(new JLineShape(3f, 0xff000000, xs, ys, xe, ys));
				graph.addChild(new JLineShape(3f, 0xff000000, xs, ye, xe, ye));
			} else {
				graph.addChild(new JLineShape(3f, 0xff000000, xs, ys, xs, ye));
				graph.addChild(new JLineShape(3f, 0xff000000, xe, ys, xe, ye));
			}
		}
		if(extendLower) {
			if(isHorizontal) {
				graph.addChild(new JLineShape(3f, 0xff000000, px,py+0.5f*ph, xs,ys));
				graph.addChild(new JLineShape(3f, 0xff000000, px,py+0.5f*ph, xs,ye));
			} else {
				graph.addChild(new JLineShape(3f, 0xff000000, px+0.5f*pw,py+ph, xs,ye));
				graph.addChild(new JLineShape(3f, 0xff000000, px+0.5f*pw,py+ph, xe,ye));
			}
		} else {
			if(isHorizontal)
				graph.addChild(new JLineShape(3f, 0xff000000, xs, ys, xs, ye));
			else
				graph.addChild(new JLineShape(3f, 0xff000000, xs, ye, xe, ye));
		}
		if(extendUpper) {
			if(isHorizontal) {
				graph.addChild(new JLineShape(3f, 0xff000000, px+pw,py+0.5f*ph, xe,ys));
				graph.addChild(new JLineShape(3f, 0xff000000, px+pw,py+0.5f*ph, xe,ye));
			} else {
				graph.addChild(new JLineShape(3f, 0xff000000, px+0.5f*pw,py, xs,ys));
				graph.addChild(new JLineShape(3f, 0xff000000, px+0.5f*pw,py, xe,ys));			}
		} else {
			if(isHorizontal)
				graph.addChild(new JLineShape(3f, 0xff000000, xe, ys, xe, ye));
			else
				graph.addChild(new JLineShape(3f, 0xff000000, xs, ys, xe, ys));
		}
		return graph;
	}
	@Override
	public JPlotShape createPlotOnlyAxes(PApplet applet, int w, int h) {
//		System.out.println("[JCOLOURBAR] crearePlotOnlyAxes(...) called.");
//		if (pplot.isDebug()) {
//			System.out.println("[DEBUG] JColourbar: begin colourbar\n"+
//					"                    x/y/width/height: "+px+"/"+py+"/"+pw+"/"+ph);
//		}
//		int trihgt = Math.min(pw,ph);
//		int xs = px, ys = py, xe = px+pw, ye = py+ph;
//		if(isHorizontal) {
//			if(extendLower) xs += trihgt;
//			if(extendUpper) xe -= trihgt;
//		} else {
//			if(extendLower) ye -= trihgt;
//			if(extendUpper) ys += trihgt;
//		}
		JGroupShape graph = new JGroupShape();
//		if(srcAxis==null)
//			return graph;
//		minV = Double.NaN;
//		maxV = Double.NaN;
//		foundPrimary = false;
//		minH = Double.POSITIVE_INFINITY;
//		maxH = Double.NEGATIVE_INFINITY;
//		if (pplot.isDebug())
//			System.out.println("[DEBUG] JColourbar: source axis exist.");
		collectHatchPatterns();
		collectPrimarySource();
//		if(!foundPrimary && hatchPatterns.isEmpty()) {
//			System.err.println("[ERR] JColourbar: Could not find source to create JColourbar.");
//			return graph;
//		}
//		if(!foundPrimary) {
//			minV = minH;
//			maxV = maxH;
//		}
//		if (Double.isNaN(minV) || Double.isNaN(maxV)) {
//			System.err.println("[ERR] JColourbar: no source for colourbar found!");
//			return graph;
//		}
//		System.out.println("Preloaded datarange: [ z={"+minV+" ... "+maxV+"}");
		if (isHorizontal)
			graph.addChild(createXAxis());
		if (!isHorizontal)
			graph.addChild(createYAxis());
		return graph;
	}
	
	public void setExtent(String what) {
		if(what.equalsIgnoreCase("neither")) { extendLower = false; extendUpper = false; }
		if(what.equalsIgnoreCase("lower")) { extendLower = true; extendUpper = false; }
		if(what.equalsIgnoreCase("upper")) { extendLower = false; extendUpper = true; }
		if(what.equalsIgnoreCase("both")) { extendLower = true; extendUpper = true; }
	}
	
	private JGroupShape createXAxis() {
		int trihgt = Math.min(pw,ph);
		int xs = extendLower ? px+trihgt : px;
		int xe = extendUpper ? px+pw-trihgt : px+pw;
		int xw = xe-xs;
		JGroupShape axisgrid = new JGroupShape();
		double[] ticks = JPlotMath.optimalLinearTicks(minC, maxC);
		double[] tcpos = JPlotMath.map(ticks, minC, maxC, xs, xe);
		String[] tickmark = new String[ticks.length];
		double vf = 1d / (ticks[0]);
		int decimal = (int) (1000d * ticks[1] + 0.5d);
		decimal = decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
		double tmlen = 0d, tmlc = 0d;
		pplot.getGraphic().textSize(200);
		for (int t = 0; t < ticks.length; t++) {
			tickmark[t] = PApplet.nf((float) (ticks[t] * vf), 0, decimal);
			if (t > 1) {
				tmlen += pplot.getGraphic().textWidth(tickmark[t]) / 200f;
				tmlc += 1d;
			}
		}
		tmlen *= this.txtsize / tmlc;
		int tickcount = Math.max(2, (int) (xw / (1.2d * tmlen) + 0.9999d));
		if (pplot.isDebug())
			System.out.println("[DEBUG] JAxis-object: tmlen=" + tmlen + " -> tickcount approx. " + tickcount);
		ticks = JPlotMath.optimalLinearTicks(minC, maxC, tickcount);
		tcpos = JPlotMath.map(ticks, minC, maxC, xs, xe);
		tickmark = new String[ticks.length];
		vf = 1d / (ticks[0]);
		decimal = (int) (1000d * ticks[1] + 0.5d);
		decimal = decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
		for (int t = 0; t < ticks.length; t++)
			tickmark[t] = PApplet.nf((float) (ticks[t] * vf), 0, decimal);
		if (pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for (int t = 2; t < ticks.length; t++) {
				tickStr += ", " + PApplet.nf((float) ticks[t], 0, 2);
				posStr += ", " + PApplet.nf((float) tcpos[t], 0, 2);
			}
			System.out.println("[DEBUG] JAxis-object: Xtickfactors={p10: " + ticks[0] + ", f: " + ticks[1] + "}");
			System.out.println("[DEBUG] JAxis-object: Xtickval={" + tickStr.substring(2) + "}");
			System.out.println("[DEBUG] JAxis-object: Xtickpos={" + posStr.substring(2) + "}");
		}
		for (int t = 2; t < ticks.length; t++)
			if (ticks[t] >= Math.min(minC, maxC) && ticks[t] <= Math.max(minC, maxC)) {
				axisgrid.addChild(new JLineShape(2f, 0xff000000, (float) tcpos[t], py + ph, (float) tcpos[t], py + 1.125f * ph));
				// axisgrid.addChild(ap.createShape(PShape.TEXT, "H",
				// (float)tcpos[t],py-0.1f*ph,(float)tcpos[t],py));
				axisgrid.addChild(new JTextShape(tickmark[t], (float) tcpos[t], py + 1.188f * ph, (float) txtsize,
						PConstants.CENTER, PConstants.TOP, 0xff000000, 0, null));
			}
		if (titleC.length() > 0) {
			if(JPlot.supportLatex)
				axisgrid.addChild(new JLatexShape(titleC, xs+0.5f*xw, py+1.250f*ph+(float)txtsize, (float)(1.1d*txtsize), PConstants.CENTER, PConstants.TOP, 0xff000000, 0f, null));
			else
				axisgrid.addChild(new JTextShape(titleC, xs+0.5f*xw, py+1.250f*ph+(float)txtsize, (float)(1.1d*txtsize), PConstants.CENTER, PConstants.TOP, 0xff000000, 0, null));
		}
		return axisgrid;
	}
	private JGroupShape createYAxis() {
		int trihgt = Math.min(pw,ph);
		int ys = extendUpper ? py+trihgt : py;
		int ye = extendLower ? py+ph-trihgt : py+ph;
		int yh = ye-ys;
		JGroupShape axisgrid = new JGroupShape();
		int nticks = (int)((ye-ys) / txtsize + 0.5d);
		double[] ticks = JPlotMath.optimalLinearTicks(minC, maxC, nticks);
		double[] tcpos = JPlotMath.map(ticks, minC, maxC, ye, ys);
//		for(int t=0; t<tcpos.length; t++)
//			tcpos[t] = 2*py+ph-tcpos[t];
		if (pplot.isDebug()) {
			String tickStr = "", posStr = "";
			for (int t = 2; t < ticks.length; t++) {
				tickStr += ", " + PApplet.nf((float) ticks[t], 0, 2);
				posStr += ", " + PApplet.nf((float) tcpos[t], 0, 2);
			}
			System.out.println("[DEBUG] JAxis-object: Ytickfactors={p10: " + ticks[0] + ", f: " + ticks[1] + "}");
			System.out.println("[DEBUG] JAxis-object: Ytickval={" + tickStr.substring(2) + "}");
			System.out.println("[DEBUG] JAxis-object: Ytickpos={" + posStr.substring(2) + "}");
		}
		float tw = 0f;
		double vf = 1d / (ticks[0]);
		int decimal = (int) (1000d * ticks[1] + 0.5d);
		decimal = decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
		for (int t = 2; t < ticks.length; t++)
			if (ticks[t] >= Math.min(minC, maxC) && ticks[t] <= Math.max(minC, maxC)) {
				axisgrid.addChild(new JLineShape(2f, 0xff000000, px + pw, (float) tcpos[t], px + 1.250f * pw, (float) tcpos[t]));
				String tmstr = PApplet.nf((float) (ticks[t] * vf), 0, decimal);
				tw = Math.max(tw, (float) txtsize * pplot.getGraphic().textWidth(tmstr) / pplot.getGraphic().textSize);
				axisgrid.addChild(new JTextShape(tmstr, px + 1.375f * pw, (float) tcpos[t], (float) txtsize,
						PConstants.LEFT, PConstants.CENTER, 0xff000000, 0, null));
			}
		if (titleC.length() > 0) {
			if(JPlot.supportLatex)
				axisgrid.addChild(new JLatexShape(titleC, px+1.500f*pw+tw, ys+0.5f*yh, (float)(1.1d*txtsize), PConstants.CENTER, PConstants.TOP, 0xff000000, JPlotShape.ROTATE_COUNTERCLOCKWISE, null));
			else
				axisgrid.addChild(new JTextShape(titleC, px+1.500f*pw+tw, ys+0.5f*yh, (float)(1.1d*txtsize), PConstants.CENTER, PConstants.TOP, 0xff000000, JPlotShape.ROTATE_COUNTERCLOCKWISE, null));
		}
		return axisgrid;
	}
	
	private void collectPrimarySource() {
		JContourLayer contours = null;
		JContourLayer2D contours2 = null;
		JXYLayer xys = null;
		for (JPlotsLayer layer : srcAxis.getLayers()) {
			if (layer instanceof JContourLayer) {
				contours = (JContourLayer) layer;
				srcColT = contours.getColourtable();
				if(srcColT==null) {
					contours = null;
				} else {
					break;
				}
			}
			if (layer instanceof JContourLayer2D) {
				contours2 = (JContourLayer2D) layer;
				srcColT = contours2.getColourtable();
				if(srcColT==null) {
					contours2 = null;
				} else {
					break;
				}
			}
			if (layer instanceof JXYLayer) {
				xys = (JXYLayer) layer;
				srcColT = xys.getColourtable();
				if(srcColT==null) {
					xys = null;
				} else {
					break;
				}
			}
		}
		if(contours!=null) {
			if (pplot.isDebug())
				System.out.println("[DEBUG] JColourbar: found JContourLayer: " + contours);
			double[] zr = contours.getZRange();
			minV = zr[0];
			maxV = zr[1];
			contourLevels = contours.getLevels();
			minC= JPlotMath.dmin(contourLevels);
			maxC= JPlotMath.dmax(contourLevels);
			if(minC<minV) minV = minC;
			if(maxC>maxV) maxV = maxC;
			if (Double.isNaN(minC))
				minC = minV;
			if (Double.isNaN(maxC))
				maxC = maxV;
			if (pplot.isDebug())
				System.out.println("[DEBUG] JColourbar: set vmin/vmax:  " + minV + "/" + maxV);
			foundPrimary = true;
			return;
		}
		if(contours2!=null) {
			if (pplot.isDebug())
				System.out.println("[DEBUG] JColourbar: found JContourLayer: " + contours2);
			double[] zr = contours2.getZRange();
			minV = zr[0];
			maxV = zr[1];
			contourLevels = contours2.getLevels();
			minC= JPlotMath.dmin(contourLevels);
			maxC= JPlotMath.dmax(contourLevels);
			if(minC<minV) minV = minC;
			if(maxC>maxV) maxV = maxC;
			if (Double.isNaN(minC))
				minC = minV;
			if (Double.isNaN(maxC))
				maxC = maxV;
			if (pplot.isDebug())
				System.out.println("[DEBUG] JColourbar: set vmin/vmax:  " + minV + "/" + maxV);
			foundPrimary = true;
			return;
		}
		if(xys!=null) {
			if (pplot.isDebug())
				System.out.println("[DEBUG] JColourbar: found JXYLayer: " + xys);
			double[] zr = xys.getZRange();
			minV = zr[0];
			maxV = zr[1];
			contourLevels = new double[1001];
			for(int i=0; i<=1000; i++)
				contourLevels[i] = zr[0] + (zr[1]-zr[0])*i/1000d;
			if (Double.isNaN(minC))
				minC = minV;
			if (Double.isNaN(maxC))
				maxC = maxV;
			if (pplot.isDebug())
				System.out.println("[DEBUG] JColourbar: set vmin/vmax:  " + minV + "/" + maxV);
			foundPrimary = true;
			return;
		}
		foundPrimary = false;
	}
	private void collectHatchPatterns() {
		for(JPlotsLayer layer : srcAxis.getLayers())
			if (layer instanceof JHatchLayer) {
				if (pplot.isDebug())
					System.out.println("[DEBUG] JColourbar: found JHatchLayer: " + layer);
				JHatchLayer hl = (JHatchLayer) layer;
				double[] zr = hl.getZRange();
				minH = Math.min(minH, Math.max(zr[0], -1000000d));
				maxH = Math.max(maxH, Math.min(zr[1],  1000000d));
//				hatchPatterns.add(
//						hl.getPattern()+";;;;"+
//						zr[0]+";;;;"+
//						zr[1]+";;;;"+
//						hl.getLineWidth()+";;;;"+
//						hl.getLineColour()
//				);
				hatchPatterns.add(hl);
			}
	}
	
	
	private int getColor(double value) {
		if(Double.isNaN(value)) return srcColT.getColour(Double.NaN);
		if(value<contourLevels[0]) return srcColT.getColour(-1d);
		int cl = contourLevels.length-1;
		if(value>contourLevels[cl]) return srcColT.getColour(2d);
		for(int i=1; i<=cl; i++)
			if(value<=contourLevels[i])
				return srcColT.getColour((i-0.5d)/cl);
		return srcColT.getColour(Double.NaN);
	}
	private void drawPrimarySource(PApplet applet, JGroupShape s, List<JDTriangle> triangles) {
		if (img == null) {
			img = applet.createImage(pw, ph, PConstants.ARGB);
		} else if (img.width != pw || img.height != ph) {
			img = applet.createImage(pw, ph, PConstants.ARGB);
		}
		img.loadPixels();
		for(int p=0; p<img.pixels.length; p++)
			img.pixels[p] = 0x00999999;
		img.updatePixels();
		
		if(srcColT==null)
			return; // new JImageShape(img, px, py, pw, ph);
		
		for(JDTriangle tri: triangles) {
			double	txi = Math.min(Math.min(tri.x[0],tri.x[1]), tri.x[2]),
					txa = Math.max(Math.max(tri.x[0],tri.x[1]), tri.x[2]);
			double	tyi = Math.min(Math.min(tri.y[0],tri.y[1]), tri.y[2]),
					tya = Math.max(Math.max(tri.y[0],tri.y[1]), tri.y[2]);
			int ixs =  Math.max((int) txi    - (txi < 0 ? 1 : 0), px),
				ixe = -Math.max((int) (-txa) - (txa > 0 ? 1 : 0), 1 - px - pw);
			int iys =  Math.max((int) tyi    - (tyi < 0 ? 1 : 0), py),
				iye = -Math.max((int) (-tya) - (tya > 0 ? 1 : 0), 1 - py - ph);
			if ((ixe < ixs) || (iye < iys)) continue;
			for (int v = iys; v <= iye; v++) {
				for (int u = ixs; u <= ixe; u++) {
					JDPoint ij = new JDPoint(u+0.5d,v+0.5d);
					if(!tri.contains(ij)) continue;
					img.pixels[(v-py) * pw + u-px] = getColor(tri.valueAt(ij));
				}
			}
		}
		img.updatePixels();
		s.addChild(new JImageShape(img, px, py, pw, ph));
	}
	private void drawHatchLayers(JGroupShape s, List<JDTriangle> triangles) {
		if(foundPrimary) return;
		for(JHatchLayer layer: hatchPatterns) {
			int[] codes = JHatchLayer.patternCodes(layer.getPattern());
			layer.drawPatterns(triangles, s, codes[0], codes[1], this);
		}
	}
}
