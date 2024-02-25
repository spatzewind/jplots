package jplots.axes;

import java.util.ArrayList;
import java.util.List;

import jplots.JPlot;
import jplots.colour.JColourtable;
import jplots.layer.JContourLayer;
import jplots.layer.JContourLayer2D;
import jplots.layer.JHatchLayer;
import jplots.layer.JHatchLayer2D;
import jplots.layer.JLegend;
import jplots.layer.JAxisAlignedLineLayer;
import jplots.layer.JPColourLayer;
import jplots.layer.JPlotsLayer;
import jplots.layer.JScatterLayer;
import jplots.layer.JTextLayer;
import jplots.layer.JXYLayer;
import jplots.maths.JPlotMath;
import jplots.maths.JPlotMath.DateTime;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLatexShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JTextShape;
import jplots.transform.JProjection;
import processing.core.PApplet;
import processing.core.PConstants;

public class JMultiAxis extends JAxis {

	private int activeSubAxis, subdivX, subdivY;
	private LayerList[] layers;
	private boolean[] xRangeFix, yRangeFix, xAxInv, yAxInv, xLog, yLog, xTim, yTim;
	private int[] subPosX, subPosY, subSizX, subSizY;
	private double[] minX, maxX, minY, maxY;
	private String[] titleX, titleY, unitX, unitY;
	private String[] xTimUnit, yTimUnit, xTimCal, yTimCal, xTimFormat, yTimFormat;
	
	public JMultiAxis(JPlot plot, int pos_x, int pos_y, int width, int height, int subDivX, int subDivY) {
		super(plot, pos_x,pos_y, width,height);
		subdivX = subDivX;
		subdivY = subDivY;
		defineSubAxes(true);
	}
	public JMultiAxis(JAxis src_axis, int subDivX, int subDivY) {
		super(src_axis);
		subdivX = subDivX;
		subdivY = subDivY;
		defineSubAxes(true);
	}
	public JMultiAxis(JMultiAxis src_axis) {
		super(src_axis);
		subdivX = src_axis.subdivX;
		subdivY = src_axis.subdivY;
		defineSubAxes(true);
	}
	
	@Override
	public JAxis setPositionAndSize(int pos_x, int pos_y, int width, int height) {
		super.setPositionAndSize(pos_x, pos_y, width, height);
		defineSubAxes(false);
		return this;
	}
	
	
	
	public JMultiAxis addLayer(JPlotsLayer layer, Object... params) {
		return addLayer(activeSubAxis, layer, params);
	}
	public JMultiAxis addLayer(int subaxis, JPlotsLayer layer, Object... params) {
		layers[subaxis].add(layer);
		readParams(subaxis, layer, params);
		activeSubAxis = subaxis;
		return this;
	}
	
	
	
	
	// *********************************
	// **** IMPLEMENTATION ABSTRACT ****
	// *********************************
	
	@Override
	public JMultiAxis copy() {
		return new JMultiAxis(this);
	}
	@Override
	public double[] getRange() {
		if(activeSubAxis<0 || activeSubAxis>=layers.length) return null;
		int sdx=activeSubAxis%subdivX, sdy = activeSubAxis/subdivX;
		return new double[] {minX[sdx],maxX[sdx], minY[sdy], maxY[sdy]};
	}
	@Override
	public int[] getSize() {
		if(activeSubAxis<0 || activeSubAxis>=layers.length) return null;
		int sdx=activeSubAxis%subdivX, sdy = activeSubAxis/subdivX;
		return new int[] {subPosX[sdx],subPosY[sdy], subSizX[sdx], subSizY[sdy]};
	}
//	@Override
//	public boolean isXlogAxis() {
//		for(int sdx=0; sdx<subdivX; sdx++)
//			if(xLog[sdx])
//				return true;
//		return false;
//	}
//	@Override
//	public boolean isYlogAxis() {
//		for(int sdy=0; sdy<subdivY; sdy++)
//			if(yLog[sdy])
//				return true;
//		return false;
//	}
	@Override
	public AxisScale getScaleX() {
		if(subdivX==1)
			return scaleX[0];
		System.err.println("Mulitple x-scales are present, use getScaleX(int column)");
		return null;
	}
	@Override
	public AxisScale getScaleY() {
		if(subdivX==1)
			return scaleY[0];
		System.err.println("Mulitple x-scales are present, use getScaleX(int column)");
		return null;
	}
	
	@Override
	public boolean isGeoAxis() {
		return false;
	}
	@Override
	public void setGeoProjection(JProjection proj) {
		System.err.println("This JMultiAxis-object cannot have a GeoProjection.");
	}
	@Override
	public List<JPlotsLayer> getLayers() {
		List<JPlotsLayer> all = new ArrayList<>();
		for(LayerList ll: layers) all.addAll(ll);
		return all;
	}
	
	
	
	
	
	
	// **********************
	// **** PLOT ROUTINE ****
	// **********************

	@Override
	public JGroupShape createPlot(PApplet applet, int w, int h) {
		JGroupShape graph = new JGroupShape();
		boolean foundError = false;
		for(int subax=0; subax<subdivX; subax++)
			if (xLog[subax] && (minX[subax] < 0d || maxX[subax] < 0d)) {
//					xLog[subax] = false;
				foundError = true;
				System.err.println("found negative values in X-range [x=" + minX[subax] + " ... " + maxX[subax] + "]");
			}
		for(int subax=0; subax<subdivY; subax++)
			if (yLog[subax] && (minY[subax] < 0d || maxY[subax] < 0d)) {
//				yLog[subax] = false;
				foundError = true;
				System.err.println("found negative values in Y-range [y=" + minY + " ... " + maxY + "]");
			}
		if(foundError) return graph;
		if (pplot.isDebug()) {
			int numlayers = 0;
			for(LayerList ll: layers) numlayers += ll.size();
			System.out.println("[DEBUG] JMultiAxis-object: min/max={x:" + minX + "/" + maxX + ", y:" + minY + "/" + maxY
					+ "} with " + numlayers + " layer" + (numlayers > 1 ? "s" : ""));
		}
		for(int subfig=0; subfig<layers.length; subfig++) {
			activeSubAxis = subfig;
			int sdx = subfig%subdivX, sdy = subfig/subdivX;
			for (int l = 0; l < layers[subfig].size(); l++) {
				JPlotsLayer layer = layers[subfig].get(l);
				if (layer instanceof JLegend)
					continue;
				layer.setRange(minX[sdx], maxX[sdx], minY[sdy], maxY[sdy]);
				if (xAxInv[sdx] || yAxInv[sdy])
					layer.invert(xAxInv[sdx] ? (yAxInv[sdy] ? "both" : "x") : "y", true);
				layer.createVectorImg(this, l, graph);
			}
		}
		for(int subfig=0; subfig<layers.length; subfig++) {
			activeSubAxis = subfig;
			int sdx = subfig%subdivX, sdy = subfig/subdivX;
			for (int l = 0; l < layers[subfig].size(); l++) {
				JPlotsLayer layer = layers[subfig].get(l);
				if (layer instanceof JLegend) {
					layer.setRange(minX[sdx], maxX[sdx], minY[sdy], maxY[sdy]);
					if (xAxInv[sdx] || yAxInv[sdy])
						layer.invert(xAxInv[sdx] ? (yAxInv[sdy] ? "both" : "x") : "y", true);
					layer.createVectorImg(this, l, graph);
				}
			}
		}
		if (xAxOn || xGrdOn)
			graph.addChild(createXAxis());
		if (yAxOn || yGrdOn)
			graph.addChild(createYAxis());
		if (xAxOn) {
			for(int sdx=0; sdx<subdivX; sdx++) {
				float sx = subPosX[sdx]; if(sdx>0) sx -= 0.35f*(subPosX[sdx]-subPosX[sdx-1]-subSizX[sdx-1]);
				float ex = subPosX[sdx]+subSizX[sdx]; if(sdx+1<subdivX) ex += 0.35f*(subPosX[sdx+1]-subPosX[sdx]-subSizX[sdx]);
				graph.addChild(new JLineShape(3f, 0xff000000, sx, py,      ex, py));
				graph.addChild(new JLineShape(3f, 0xff000000, sx, py + ph, ex, py + ph));
			}
		}
		if (yAxOn) {
			for(int sdy=0; sdy<subdivY; sdy++) {
				float sy = subPosY[sdy]; if(sdy>0) sy -= 0.35f*(subPosY[sdy]-subPosY[sdy-1]-subSizY[sdy-1]);
				float ey = subPosY[sdy]+subSizY[sdy]; if(sdy+1<subdivY) ey += 0.35f*(subPosY[sdy+1]-subPosY[sdy]-subSizY[sdy]);
				graph.addChild(new JLineShape(3f, 0xff000000, px,      sy, px,      ey));
				graph.addChild(new JLineShape(3f, 0xff000000, px + pw, sy, px + pw, ey));
			}
		}
		if (titleP.length() > 0) {
			if (pplot.isDebug())
				System.out.println("[DEBUG] JAxis: add title \"" + titleP + "\" to graphic.");
			if(JPlot.supportLatex)
				graph.addChild(new JLatexShape(titleP, px+0.5f*pw, py-0.04f*ph, (float)(1.3d*txtsize), PConstants.CENTER, PConstants.BOTTOM, 0xff000000, 0f, null));
			else
				graph.addChild(new JTextShape(titleP, px+0.5f*pw, py-0.04f*ph, (float)(1.3d*txtsize), PConstants.CENTER, PConstants.BOTTOM, 0xff000000, 0f, null));
		}
		graph.addChild(new JLineShape(0f, 0x00999999, 0f,0f, 1f,1f));
		return graph;
	}
	@Override
	public JGroupShape createPlotOnlyAxes(PApplet applet, int w, int h) {
		JGroupShape graph = new JGroupShape();
		boolean foundError = false;
		for(int subax=0; subax<subdivX; subax++)
			if (xLog[subax] && (minX[subax] < 0d || maxX[subax] < 0d)) {
//					xLog[subax] = false;
				foundError = true;
				System.err.println("found negative values in X-range [x=" + minX[subax] + " ... " + maxX[subax] + "]");
			}
		for(int subax=0; subax<subdivY; subax++)
			if (yLog[subax] && (minY[subax] < 0d || maxY[subax] < 0d)) {
//				yLog[subax] = false;
				foundError = true;
				System.err.println("found negative values in Y-range [y=" + minY + " ... " + maxY + "]");
			}
		if(foundError) return graph;
		for(int subax=0; subax<layers.length; subax++) {
			for(JPlotsLayer layer: layers[subax]) {
				if(layer instanceof JLegend) layer.createVectorImg(this, 0, graph);
				if(layer instanceof JTextLayer) layer.createVectorImg(this, 0, graph);
			}
		}
		if (xAxOn || xGrdOn)
			graph.addChild(createXAxis());
		if (yAxOn || yGrdOn)
			graph.addChild(createYAxis());
		if (xAxOn) {
			graph.addChild(new JLineShape(3f, 0xff000000, px, (float)py, px + pw, py));
			graph.addChild(new JLineShape(3f, 0xff000000, px, (float)py + ph, px + pw, py + ph));
		}
		if (yAxOn) {
			graph.addChild(new JLineShape(3f, 0xff000000, px, (float)py, px, py + ph));
			graph.addChild(new JLineShape(3f, 0xff000000, px + pw, (float)py, px + pw, py + ph));
		}
		if (titleP.length() > 0) {
			if (pplot.isDebug())
				System.out.println("[DEBUG] JAxis: add title \"" + titleP + "\" to graphic.");
			if(JPlot.supportLatex) {
				graph.addChild(new JLatexShape(titleP, px + 0.5f * pw, py - 0.04f * ph, (float) (1.3d * txtsize), PConstants.CENTER,
						PConstants.BOTTOM, 0xff000000, 0f, null));
			} else {
				graph.addChild(new JTextShape(titleP, px + 0.5f * pw, py - 0.04f * ph, (float) (1.3d * txtsize), PConstants.CENTER,
						PConstants.BOTTOM, 0xff000000, 0f, null));
			}
		}
		graph.addChild(new JLineShape(0f, 0x00999999, 0f,0f, 1f,1f));
		return graph;
	}
	
	
	
	// ****************
	// **** PUBLIC ****
	// ****************
	
	@Override
	public void clear() {
		super.clear();
		defineSubAxes(false);
		for(LayerList ll: layers) ll.clear();
	}
	
	public void contour(int subaxis, float[] x, float[] y, float[][] z) {
		this.contour(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contour(int subaxis, float[] x, float[] y, float[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, null);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contour(int subaxis, float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), 0, levels);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contour(int subaxis, float[][] x, float[][] y, float[][] z) {
		this.contour(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contour(int subaxis, float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, null);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contour(int subaxis, float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), 0, levels);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contour(int subaxis, double[] x, double[] y, double[][] z) {
		this.contour(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contour(int subaxis, double[] x, double[] y, double[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, null);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contour(int subaxis, double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), 0, levels);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contour(int subaxis, double[][] x, double[][] y, double[][] z) {
		this.contour(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contour(int subaxis, double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, null);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contour(int subaxis, double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), 0, levels);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	
	public void contourf(int subaxis, float[] x, float[] y, float[][] z) {
		this.contourf(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contourf(int subaxis, float[] x, float[] y, float[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels,
				JColourtable.pctables.get("default"), 2.0f, false, true, false);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourf(int subaxis, float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0f, false, true,
				false);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourf(int subaxis, float[][] x, float[][] y, float[][] z) {
		this.contourf(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contourf(int subaxis, float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels,
				JColourtable.pctables.get("default"), 2.0f, false, true, false);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourf(int subaxis, float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0f, false, true,
				false);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourf(int subaxis, double[] x, double[] y, double[][] z) {
		this.contourf(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contourf(int subaxis, double[] x, double[] y, double[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels,
				JColourtable.pctables.get("default"), 2.0d, false, true, false);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourf(int subaxis, double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0d, false, true,
				false);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourf(int subaxis, double[][] x, double[][] y, double[][] z) {
		this.contourf(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contourf(int subaxis, double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels,
				JColourtable.pctables.get("default"), 2.0d, false, true, false);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourf(int subaxis, double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		JPlotsLayer cnl = new JContourLayer(x, y, z, levels, JColourtable.pctables.get("default"), 2.0d, false, true,
				false);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	
	public void contourp(int subaxis, float[] x, float[] y, float[][] z) {
		this.contourf(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contourp(int subaxis, float[] x, float[] y, float[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, null);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourp(int subaxis, float[] x, float[] y, float[][] z, float[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z),JPlotMath.fmax(z),0, levels);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourp(int subaxis, float[][] x, float[][] y, float[][] z) {
		this.contourf(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contourp(int subaxis, float[][] x, float[][] y, float[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z), JPlotMath.fmax(z), levels, null);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourp(int subaxis, float[][] x, float[][] y, float[][] z, float[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.fmin(z),JPlotMath.fmax(z),0, levels);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourp(int subaxis, double[] x, double[] y, double[][] z) {
		this.contourp(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contourp(int subaxis, double[] x, double[] y, double[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, null);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourp(int subaxis, double[] x, double[] y, double[][] z, double[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z),JPlotMath.dmax(z), 0, levels);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourp(int subaxis, double[][] x, double[][] y, double[][] z) {
		this.contourp(subaxis, x, y, z, 10, (Object[]) null);
	}
	public void contourp(int subaxis, double[][] x, double[][] y, double[][] z, int levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z), JPlotMath.dmax(z), levels, null);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	public void contourp(int subaxis, double[][] x, double[][] y, double[][] z, double[] levels, Object... params) {
		JContourLayer2D cnl = new JContourLayer2D(x, y, z, JPlotMath.dmin(z),JPlotMath.dmax(z),0, levels);
		cnl.setColourtable(JColourtable.pctables.get("default"));
		cnl.setLineWidth(2.0f); cnl.lines(false); cnl.setFilled(true); cnl.setPixelFilling(true);
		layers[subaxis].add(cnl);
		readParams(subaxis, cnl, params);
		updateRange(subaxis, cnl);
	}
	
	public void pcolour(int subaxis, float[] x, float[] y, float[][] z, float cmin, float cmax, Object... params) {
		JPColourLayer pcl = new JPColourLayer(x, y, z, cmin, cmax);
		pcl.setColourtable(JColourtable.pctables.get("default"));
		layers[subaxis].add(pcl);
		readParams(subaxis, pcl, params);
		updateRange(subaxis, pcl);
	}
	public void pcolour(int subaxis, float[][] x, float[][] y, float[][] z, float cmin, float cmax, Object... params) {
		JPColourLayer pcl = new JPColourLayer(x, y, z, cmin, cmax);
		pcl.setColourtable(JColourtable.pctables.get("default"));
		layers[subaxis].add(pcl);
		readParams(subaxis, pcl, params);
		updateRange(subaxis, pcl);
	}
	public void pcolour(int subaxis, double[] x, double[] y, double[][] z, double cmin, double cmax, Object... params) {
		JPColourLayer pcl = new JPColourLayer(x, y, z, cmin, cmax);
		pcl.setColourtable(JColourtable.pctables.get("default"));
		layers[subaxis].add(pcl);
		readParams(subaxis, pcl, params);
		updateRange(subaxis, pcl);
	}
	public void pcolour(int subaxis, double[][] x, double[][] y, double[][] z, double cmin, double cmax, Object... params) {
		JPColourLayer pcl = new JPColourLayer(x, y, z, cmin, cmax);
		pcl.setColourtable(JColourtable.pctables.get("default"));
		layers[subaxis].add(pcl);
		readParams(subaxis, pcl, params);
		updateRange(subaxis, pcl);
	}
	
	public void hatch(int subaxis, float[] x, float[] y, float[][] z, float lower, float upper, String pattern) {
		this.hatch(subaxis, x,y,z, lower,upper, pattern, (Object)null);
	}
	public void hatch(int subaxis, float[] x, float[] y, float[][] z, float lower, float upper, String pattern, Object... params) {
		JPlotsLayer hl = new JHatchLayer(x, y, null, null, z, 1f, lower, upper, pattern);
		layers[subaxis].add(hl);
		readParams(subaxis, hl, params);
		updateRange(subaxis, hl);
	}
	public void hatch(int subaxis, float[][] x, float[][] y, float[][] z, float lower, float upper, String pattern) {
		this.hatch(subaxis, x,y,z, lower,upper, pattern, (Object)null);
	}
	public void hatch(int subaxis, float[][] x, float[][] y, float[][] z, float lower, float upper, String pattern, Object... params) {
		JPlotsLayer hl = new JHatchLayer(null, null, x, y, z, 1f, lower, upper, pattern);
		layers[subaxis].add(hl);
		readParams(subaxis, hl, params);
		updateRange(subaxis, hl);
	}
	public void hatch(int subaxis, double[] x, double[] y, double[][] z, double lower, double upper, String pattern) {
		this.hatch(subaxis, x, y, z, lower, upper, pattern, (Object)null);
	}
	public void hatch(int subaxis, double[] x, double[] y, double[][] z, double lower, double upper, String pattern, Object... params) {
		JPlotsLayer hl = new JHatchLayer(x, y, null, null, z, 1d, lower, upper, pattern);
		layers[subaxis].add(hl);
		readParams(subaxis, hl, params);
		updateRange(subaxis, hl);
	}
	public void hatch(int subaxis, double[][] x, double[][] y, double[][] z, double lower, double upper, String pattern) {
		this.hatch(subaxis, x, y, z, lower, upper, pattern, (Object)null);
	}
	public void hatch(int subaxis, double[][] x, double[][] y, double[][] z, double lower, double upper, String pattern, Object... params) {
		JPlotsLayer hl = new JHatchLayer(null, null, x, y, z, 1d, lower, upper, pattern);
		layers[subaxis].add(hl);
		readParams(subaxis, hl, params);
		updateRange(subaxis, hl);
	}
	
	public void plot(int subaxis, float[] x, float[] y) {
		this.plot(subaxis, x, y, 0xff000000, 3f, "-", (Object) null);
	}
	public void plot(int subaxis, float[] x, float[] y, int colour, float linewidth, String linestyle, Object... params) {
		JPlotsLayer xyl = new JXYLayer(x, y, colour, linewidth, linestyle);
		layers[subaxis].add(xyl);
		readParams(subaxis, xyl, params);
		updateRange(subaxis, xyl);
	}
	public void plot(int subaxis, double[] x, double[] y) {
		this.plot(subaxis, x, y, 0xff000000, 3d, "-", (Object) null);
	}
	public void plot(int subaxis, double[] x, double[] y, int colour, double linewidth, String linestyle, Object... params) {
		JPlotsLayer xyl = new JXYLayer(x, y, colour, linewidth, linestyle);
		layers[subaxis].add(xyl);
		readParams(subaxis, xyl, params);
		updateRange(subaxis, xyl);
	}
	
	public void scatter(int subaxis, float[] x, float[] y) {
		this.scatter(subaxis, x, y, 0xff000000, 1f, "c", (Object) null);
	}
	public void scatter(int subaxis, float[] x, float[] y, int colour, float iconsize, String symbol, Object... params) {
		JPlotsLayer scl = new JScatterLayer(x, y, colour, iconsize, symbol);
		layers[subaxis].add(scl);
		readParams(subaxis, scl, params);
		updateRange(subaxis, scl);
	}
	public void scatter(int subaxis, double[] x, double[] y) {
		this.scatter(subaxis, x, y, 0xff000000, 1d, "c", (Object) null);
	}
	public void scatter(int subaxis, double[] x, double[] y, int colour, double iconsize, String symbol, Object... params) {
		JPlotsLayer scl = new JScatterLayer(x, y, colour, iconsize, symbol);
		layers[subaxis].add(scl);
		readParams(subaxis, scl, params);
		updateRange(subaxis, scl);
	}
	
	public void axhline(int y_section, double y) {
		axhline(y_section, y, 0xff000000, 3f, "-");
	}
	public void axhline(int y_section, double y, int colour, double linewidth, String linestyle, Object... params) {
		JPlotsLayer xyl = new JAxisAlignedLineLayer(y, 'h', colour, linewidth, linestyle);
		for(int sdx=0; sdx<subdivX; sdx++) {
			layers[y_section*subdivX+sdx].add(xyl);
			readParams(sdx, xyl, params);
		}
		updateRange(y_section, xyl, "y");
	}
	public void axvline(int x_section, double x) {
		axvline(x_section, x, 0xff000000, 3f, "-");
	}
	public void axvline(int x_section, double x, int colour, double linewidth, String linestyle, Object... params) {
		JPlotsLayer xyl = new JAxisAlignedLineLayer(x, 'v', colour, linewidth, linestyle);
		for(int sdy=0; sdy<subdivY; sdy++) {
			layers[sdy*subdivX+x_section].add(xyl);
			readParams(sdy, xyl, params);
		}
		updateRange(x_section, xyl, "x");
	}
	
	public void legend(int subaxis) {
		JPlotsLayer lgl = new JLegend(this, PConstants.RIGHT, PConstants.TOP, false, 1d);
		layers[subaxis].add(lgl);
	}
	public void legend(int subaxis, double rts) {
		JPlotsLayer lgl = new JLegend(this, PConstants.RIGHT, PConstants.TOP, false, rts);
		layers[subaxis].add(lgl);
	}
	public void legend(int subaxis, int left_right, int top_bottom) {
		JPlotsLayer lgl = new JLegend(this, left_right, top_bottom, false, 1d);
		layers[subaxis].add(lgl);
	}
	public void legend(int subaxis, int left_right, int top_bottom, double rts) {
		JPlotsLayer lgl = new JLegend(this, left_right, top_bottom, false, rts);
		layers[subaxis].add(lgl);
	}
	
	public void setXTitle(int x_section, String xtitle) {
		titleX[x_section] = "";
		if(xtitle!=null) titleX[x_section] = xtitle;
		unitX[x_section]  = "";
	}
	public void setXTitle(int x_section, String xtitle, String xunit) {
		titleX[x_section] = "";
		if(xtitle!=null) titleX[x_section] = xtitle;
		unitX[x_section]  = "";
		if(xunit!=null) unitX[x_section] = xunit;
	}
	public void setYTitle(int y_section, String ytitle) {
		titleY[y_section] = "";
		if(ytitle!=null) titleY[y_section] = ytitle;
		unitY[y_section]  = "";
	}
	public void setYTitle(int y_section, String ytitle, String yunit) {
		titleY[y_section] = "";
		if(ytitle!=null) titleY[y_section] = ytitle;
		unitY[y_section]  = "";
		if(yunit!=null) unitY[y_section] = yunit;
	}
	
	public JMultiAxis setXRange(int x_section, double xmin, double xmax) {
		minX[x_section] = Math.min(xmin,xmax);
		maxX[x_section] = Math.max(xmin,xmax);
		xRangeFix[x_section] = true;
		return this;
	}
	public JMultiAxis setYRange(int y_section, double ymin, double ymax) {
		minY[y_section] = Math.min(ymin,ymax);
		maxY[y_section] = Math.max(ymin,ymax);
		yRangeFix[y_section] = true;
		return this;
	}
	public JMultiAxis setRange(int subaxis, double xmin, double xmax, double ymin, double ymax) {
		setXRange(subaxis%subdivX, xmin, xmax);
		setYRange(subaxis/subdivX, ymin, ymax);
		return this;
	}
	
	public void setLogarithmicAxis(int subaxis, char which) {
		switch (which) {
		case 'b':
			setLogarithmicAxis(subaxis, 'x');
			setLogarithmicAxis(subaxis, 'y');
			break;
		case 'x': int x = subaxis%subdivX;
			xLog[x] = true;
			xTim[x] = false;
			break;
		case 'y': int y = subaxis/subdivX;
			yLog[y] = true;
			yTim[y] = false;
			break;
		default:
			System.err.println("Unknown parameter '" + which + "' for which in <setLogarithmicAxis(subaxis,which)>.");
			break;
		}
	}
	public void setAsTimeAxis(int subaxis, char which, String unit) {
		setAsTimeAxis(subaxis, which, unit, "gregorian", "dd.mm.yyyy");
	}
	public void setAsTimeAxis(int subaxis, char which, String unit, String calendar) {
		setAsTimeAxis(subaxis, which, unit, calendar, "dd.mm.yyyy");
	}
	public void setAsTimeAxis(int subaxis, char which, String unit, String calendar, String format) {
		switch (which) {
		case 'b':
			setAsTimeAxis(subaxis, 'x', unit, calendar, format);
			setAsTimeAxis(subaxis, 'y', unit, calendar, format);
			break;
		case 'x': int x = subaxis%subdivX;
			xLog[x] = false;
			xTim[x] = true;
			xTimUnit[x] = unit;
			xTimCal[x] = calendar;
			if (format != null)
				xTimFormat[x] = format;
			break;
		case 'y': int y = subaxis/subdivX;
			yLog[y] = false;
			yTim[y] = true;
			yTimUnit[y] = unit;
			yTimCal[y] = calendar;
			if (format != null)
				yTimFormat[y] = format;
			break;
		default:
			System.err.println("Unknown parameter '" + which + "' for which in <setLogarithmicAxis(subaxis,which)>.");
			break;
		}
	}
	
	
	
	
	// *****************
	// **** PRIVATE ****
	// *****************
	
	private void updateRange(int sa, JPlotsLayer layer) {
		//TODO
		double[] r = layer.getRange();
		double xmin = r[0], xmax = r[1], ymin = r[2], ymax = r[3];
		if (r[1] - r[0] < 1.e-20) {
			double xm = 0.5d * (r[0] + r[1]);
			double xr = Math.max(1.0e-10d, Math.abs(xm) * 1.0e-10d);
			xmin = xm - xr;
			xmax = xm + xr;
		}
		if (r[3] - r[2] < 1.e-20) {
			double ym = 0.5d * (r[2] + r[3]);
			double yr = Math.max(1.0e-10d, Math.abs(ym) * 1.e-10d);
			ymin = ym - yr;
			ymax = ym + yr;
		}
		int sx = sa%subdivX;
		int sy = sa/subdivX;
		if (layers[sa].size() == 1) {
			if (!xRangeFix[sx]) {
				minX[sx] = xmin;
				maxX[sx] = xmax;
			}
			if (!yRangeFix[sy]) {
				minY[sy] = ymin;
				maxY[sy] = ymax;
			}
		} else {
			if (!xRangeFix[sx]) {
				if (xmin < minX[sx])
					minX[sx] = xmin;
				if (xmax > maxX[sx])
					maxX[sx] = xmax;
			}
			if (!yRangeFix[sy]) {
				if (ymin < minY[sy])
					minY[sy] = ymin;
				if (ymax > maxY[sy])
					maxY[sy] = ymax;
			}
		}
	}
	private void updateRange(int sc, JPlotsLayer layer, String axis) {
		//TODO
		double[] r = layer.getRange();
		double xmin = r[0], xmax = r[1], ymin = r[2], ymax = r[3];
		if (r[1] - r[0] < 1.e-20) {
			double xm = 0.5d * (r[0] + r[1]);
			double xr = Math.max(1.0e-10d, Math.abs(xm) * 1.0e-10d);
			xmin = xm - xr;
			xmax = xm + xr;
		}
		if (r[3] - r[2] < 1.e-20) {
			double ym = 0.5d * (r[2] + r[3]);
			double yr = Math.max(1.0e-10d, Math.abs(ym) * 1.e-10d);
			ymin = ym - yr;
			ymax = ym + yr;
		}
		if(axis.equals("x")) {
			if (!xRangeFix[sc]) {
				if (xmin < minX[sc])
					minX[sc] = xmin;
				if (xmax > maxX[sc])
					maxX[sc] = xmax;
			}
		}
		if(axis.equals("y")) {
			if (!yRangeFix[sc]) {
				if (ymin < minY[sc])
					minY[sc] = ymin;
				if (ymax > maxY[sc])
					maxY[sc] = ymax;
			}
		}
	}
	private void readParams(int subaxis, JPlotsLayer layer, Object... params) {
		if (params == null)
			return;
		int o = 0;
		while (o < params.length) {
			if (params[o] instanceof String) {
				String p = ((String) params[o]).toLowerCase();
				boolean isunread = true;
				if (isunread && ("l".equals(p) || "lines".equals(p)) && o + 1 < params.length) {
					layer.lines((boolean) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("ls".equals(p) || "linestyle".equals(p)) && o + 1 < params.length) {
					if(params[o+1] instanceof String[]) {
						if(layer instanceof JContourLayer)   ((JContourLayer)   layer).setStyles((String[]) params[o+1]);
						if(layer instanceof JContourLayer2D) ((JContourLayer2D) layer).setStyles((String[]) params[o+1]);
					} else
						layer.setStyle((String) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("lw".equals(p) || "linewidth".equals(p)) && o + 1 < params.length) {
					layer.setLineWidth(params[o + 1] instanceof Float ? (float) params[o + 1] : (double) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("lc".equals(p) || "linecolor".equals(p) || "linecolour".equals(p))
						&& o + 1 < params.length) {
					if (params[o + 1] instanceof Integer) {
						layer.setLineColour((int) params[o + 1]);
						o++;
						isunread = false;
					} else if (params[o + 1] instanceof int[]) {
						layer.setLineColours((int[]) params[o + 1]);
						o++;
						isunread = false;
					}
				}
				if (isunread && ("fc".equals(p) || "fillcolor".equals(p) || "fillcolour".equals(p))
						&& o + 1 < params.length) {
					if (params[o + 1] instanceof Integer) {
						layer.setFillColour((int) params[o + 1]);
						o++;
						isunread = false;
					} else if (params[o + 1] instanceof int[]) {
						layer.setFillColours((int[]) params[o + 1]);
						o++;
						isunread = false;
					}
				}
				if (isunread && ("fs".equals(p) || "ts".equals(p) || "fontsize".equals(p) || "textsize".equals(p))
						&& o + 1 < params.length) {
					layer.setLineWidth((double) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("lb".equals(p) || "label".equals(p))) {
					if (params[o + 1] instanceof String) {
						layer.setLabel((String) params[o + 1]);
						o++;
						isunread = false;
					}
				}
				if (isunread && ("ct".equals(p) || "colortable".equals(p) || "colourtable".equals(p))
						&& o + 1 < params.length) {
					layer.setColourtable((JColourtable) params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("z".equals(p)) && o + 1 < params.length) {
					layer.addParallelArray(params[o + 1]);
					o++;
					isunread = false;
				}
				if (isunread && ("zrange".equals(p)) && o + 2 < params.length) {
					layer.setZRange((double)params[o+1], (double)params[o+2]);
					o+=2;
					isunread = false;
				}
				if (isunread && ("density".equals(p) || "pd".equals(p)) && o + 1 < params.length) {
					if(layer instanceof JHatchLayer2D)
						((JHatchLayer2D) layer).setDensity((double)params[o+1]);
					o+=1;
					isunread = false;
				}
				if (isunread && "invertxaxis".equals(p)) {
					layer.invert("x", true);
					xAxInv[subaxis%subdivX] = true;
					isunread = false;
				}
				if (isunread && "invertyaxis".equals(p)) {
					layer.invert("y", true);
					yAxInv[subaxis/subdivX] = true;
					isunread = false;
				}
			} else {
				System.err.println("[ERROR] Cannot interprete param " + o + ": " + params[o]);
			}
			o++;
		}
	}
	

	private JGroupShape createXAxis() {
		JGroupShape axisgrid = new JGroupShape();
		for(int sdx=0; sdx<subdivX-1; sdx++) {
			double xl = subPosX[sdx]+subSizX[sdx];
			double dx = subPosX[sdx+1] - xl;
			if(xDrawSide==TOP || xDrawSide==BOTH) {
				axisgrid.addChild(new JLineShape(3f, 0xff000000, (float)(xl+0.2*dx), py+0.02f*ph, (float)(xl+0.5*dx), py-0.02f*ph));
				axisgrid.addChild(new JLineShape(3f, 0xff000000, (float)(xl+0.5*dx), py+0.02f*ph, (float)(xl+0.8*dx), py-0.02f*ph));
			}
			if(xDrawSide==BOTTOM || xDrawSide==BOTH) {
				axisgrid.addChild(new JLineShape(3f, 0xff000000, (float)(xl+0.2*dx), py+1.02f*ph, (float)(xl+0.5*dx), py+0.98f*ph));
				axisgrid.addChild(new JLineShape(3f, 0xff000000, (float)(xl+0.5*dx), py+1.02f*ph, (float)(xl+0.8*dx), py+0.98f*ph));
			}
		}
		for(int sdx=0; sdx<subdivX; sdx++) {
			int px = subPosX[sdx], pw = subSizX[sdx];
			double Xin = xLog[sdx] ? Math.log10(minX[sdx]) : minX[sdx];
			double Xax = xLog[sdx] ? Math.log10(maxX[sdx]) : maxX[sdx];
			// first estimate of ticks
			double[] oticks = null;
			String[] otickmark = null;
			if (xLog[sdx]) {
				oticks = JPlotMath.optimalLogarithmicTicks(minX[sdx], maxX[sdx]);
				otickmark = new String[oticks.length];
				for (int t = 2; t < otickmark.length; t++) {
					otickmark[t] = oticks[t] + "";
					oticks[t] = Math.log10(oticks[t]);
				}
			} else if (xTim[sdx]) {
				oticks = JPlotMath.optimalTimeTicks(minX[sdx], maxX[sdx], xTimUnit[sdx], xTimCal[sdx]);
				otickmark = new String[oticks.length];
				for (int t = 2; t < otickmark.length; t++)
					otickmark[t] = DateTime.fromDouble(oticks[t], xTimUnit[sdx], xTimCal[sdx]).format(xTimFormat[sdx], xTimCal[sdx]);
			} else {
				oticks = JPlotMath.optimalLinearTicks(minX[sdx], maxX[sdx]);
				double vf = 1d / (oticks[0]);
				int decimal = (int) (1000d * oticks[1] + 0.5d);
				decimal = decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
				otickmark = new String[oticks.length];
				for (int t = 2; t < otickmark.length; t++)
					otickmark[t] = PApplet.nf((float) (oticks[t] * vf), 0, decimal).replace(",", ".");
			}
			otickmark[0] = "";
			otickmark[1] = "";
			double tmlen = 0d;
			pplot.getGraphic().textSize(200);
			pplot.getGraphic().textAlign(PConstants.LEFT, PConstants.TOP);
			// create tickmark strings and calc mean tickmark text width
			for (int t = 2; t < oticks.length; t++) {
				tmlen += pplot.getGraphic().textWidth(otickmark[t]) / 200f;
			}
			tmlen *= this.txtsize / (oticks.length - 2);
			// with upper bound of number of ticks
			int tickcount = Math.max(2, (int) (pw / (1.2d * tmlen) + 0.99999999d));
			if (pplot.isDebug())
				System.out.println("[DEBUG] JAxis-object: tmlen=" + tmlen + " -> tickcount approx. " + tickcount);
			// create new ticks
			double[] ticks = null;
			String[] tickmark = null;
			String tickmarkFactor = "";
			if (xLog[sdx]) {
				ticks = JPlotMath.optimalLogarithmicTicks(minX[sdx], maxX[sdx], tickcount);
				tickmark = new String[ticks.length];
				for (int t = 2; t < tickmark.length; t++) {
					tickmark[t] = ticks[t] + "";
					ticks[t] = Math.log10(ticks[t]);
				}
			} else if (xTim[sdx]) {
				ticks = JPlotMath.optimalTimeTicks(minX[sdx], maxX[sdx], xTimUnit[sdx], xTimCal[sdx], tickcount);
				tickmark = new String[ticks.length];
				for (int t = 2; t < tickmark.length; t++)
					tickmark[t] = DateTime.fromDouble(ticks[t], xTimUnit[sdx], xTimCal[sdx]).format(xTimFormat[sdx], xTimCal[sdx]);
			} else {
				ticks = JPlotMath.optimalLinearTicks(minX[sdx], maxX[sdx], tickcount);
				double vf = 1d / (ticks[0]);
				int decimal = (int) (1000d * ticks[1] + 0.5d);
				decimal = decimal % 1000 == 0 ? 0 : decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
				tickmark = new String[ticks.length];
				for (int t = 2; t < tickmark.length; t++) {
					if(decimal==0) tickmark[t] = ""+(int)(ticks[t]*vf+0.0005d-(ticks[t]*vf<0d?1:0));
					else tickmark[t] = PApplet.nf((float) (ticks[t] * vf), 0, decimal).replace(",", ".");
				}
				double lvf = Math.log10(ticks[0]);
				if(Math.abs(lvf)>2.9d) {
					int ivf = (int) (lvf+0.5d) - (lvf<0d ? -1 : 0);
					tickmarkFactor = JPlot.supportLatex ? "10$^{"+ivf+"}$" : "10^"+ivf;
				}
			}
			tickmark[0] = "";
			tickmark[1] = "";
			double[] tcpos = JPlotMath.map(ticks, Xin, Xax, px, px + pw);
			if (xAxInv[sdx])
				for (int t = 0; t < tcpos.length; t++)
					tcpos[t] = 2 * px + pw - tcpos[t];
			if (pplot.isDebug()) {
				String tickStr = "", posStr = "";
				for (int t = 2; t < ticks.length; t++) {
					tickStr += ", " + tickmark[t];
					posStr += ", " + PApplet.nf((float) tcpos[t], 0, 2).replace(",", ".");
				}
				System.out.println("[DEBUG] JAxis-object: Xtickfactors={p10: " + ticks[0] + ", f: " + ticks[1] + "}");
				System.out.println("[DEBUG] JAxis-object: Xtickval={" + tickStr.substring(2) + "}");
				System.out.println("[DEBUG] JAxis-object: Xtickpos={" + posStr.substring(2) + "}");
			}
			if (xGrdOn) {
				for(int sdy=0; sdy<subdivY; sdy++) {
					float sy = subPosY[sdy]; if(sdy>0) sy -= 0.35f*(subPosY[sdy]-subPosY[sdy-1]-subSizY[sdy-1]);
					float ey = subPosY[sdy]+subSizY[sdy]; if(sdy+1<subdivY) ey += 0.35f*(subPosY[sdy+1]-subPosY[sdy]-subSizY[sdy]);
					for (int t = 2; t < ticks.length; t++)
						if (ticks[t] >= Math.min(Xin, Xax) && ticks[t] <= Math.max(Xin, Xax))
							axisgrid.addChild(new JLineShape(2f, 0xff999999, (float)tcpos[t], sy, (float)tcpos[t], ey));
				}
			}
			if (xAxOn) {
				if (pplot.isDebug())
					System.out.println("[DEBUG] JAxi-object: add x-axis title \""+titleX[sdx]+"\" and unit \""+unitX[sdx]+"\" with text size "+txtsize);
				String txtemp = ""+titleX[sdx];
				if(unitX[sdx].length()>0 || tickmarkFactor.length()>0) txtemp += (titleX[sdx].length()>0?" ":"")+"[";
				if(tickmarkFactor.length()>0) txtemp += tickmarkFactor;
				if(unitX[sdx].length()>0) txtemp += (tickmarkFactor.length()>0?" ":"")+unitX[sdx];
				if(unitX[sdx].length()>0 || tickmarkFactor.length()>0) txtemp += "]";
				//TODO set txtsize back to 1.1d*txtsize
				if(JPlot.supportLatex) {
					if(xDrawSide==TOP || xDrawSide==BOTH) {
						axisgrid.addChild(new JLatexShape(txtemp, px + 0.5f * pw, py - 0.04f * ph - (float) txtsize,
								(float) (1.1d * txtsize), CENTER, BOTTOM, 0xff000000, 0, null));
					}
					if(xDrawSide==BOTTOM || xDrawSide==BOTH) {
						axisgrid.addChild(new JLatexShape(txtemp, px + 0.5f * pw, py + 1.04f * ph + (float) txtsize,
								(float) (1.1d * txtsize), CENTER, TOP, 0xff000000, 0, null));
					}
				} else {
					if(xDrawSide==TOP || xDrawSide==BOTH) {
						axisgrid.addChild(new JTextShape(txtemp, px + 0.5f * pw, py - 0.04f * ph - (float) txtsize,
								(float) (1.1d * txtsize), CENTER, BOTTOM, 0xff000000, 0, null));
					}
					if(xDrawSide==BOTTOM || xDrawSide==BOTH) {
						axisgrid.addChild(new JTextShape(txtemp, px + 0.5f * pw, py + 1.04f * ph + (float) txtsize,
								(float) (1.1d * txtsize), CENTER, TOP, 0xff000000, 0, null));
					}
				}
				if (xTkOn) {
					for (int t = 2; t < ticks.length; t++)
						if (ticks[t] >= Math.min(Xin, Xax) && ticks[t] <= Math.max(Xin, Xax)) {
							if(xDrawSide==TOP || xDrawSide==BOTH) {
								axisgrid.addChild(new JLineShape(2f, 0xff000000, (float) tcpos[t], py, (float) tcpos[t], py-0.02f*ph));
								axisgrid.addChild(new JTextShape(tickmark[t], (float) tcpos[t], py-0.03f*ph,
										(float) txtsize, CENTER, BOTTOM, 0xff000000, 0, null));
							}
							if(xDrawSide==BOTTOM || xDrawSide==BOTH) {
								axisgrid.addChild(new JLineShape(2f, 0xff000000, (float) tcpos[t], py+ph, (float) tcpos[t], py+1.02f*ph));
								axisgrid.addChild(new JTextShape(tickmark[t], (float) tcpos[t], py+1.03f*ph,
										(float) txtsize, CENTER, TOP, 0xff000000, 0, null));
							}
						}
				}
			}
		}
		return axisgrid;
	}
	private JGroupShape createYAxis() {
		JGroupShape axisgrid = new JGroupShape();
		for(int sdy=0; sdy<subdivY-1; sdy++) {
			double yl = subPosY[sdy]+subSizY[sdy];
			double dy = subPosY[sdy+1] - yl;
			if(yDrawSide==LEFT || yDrawSide==BOTH) {
				axisgrid.addChild(new JLineShape(3f, 0xff000000, px-0.02f*pw, (float)(yl+0.2*dy), px+0.02f*pw, (float)(yl+0.5*dy)));
				axisgrid.addChild(new JLineShape(3f, 0xff000000, px-0.02f*pw, (float)(yl+0.5*dy), px+0.02f*pw, (float)(yl+0.8*dy)));
			}
			if(yDrawSide==RIGHT || yDrawSide==BOTH) {
				axisgrid.addChild(new JLineShape(3f, 0xff000000, px+0.98f*pw, (float)(yl+0.2*dy), px+1.02f*pw, (float)(yl+0.5*dy)));
				axisgrid.addChild(new JLineShape(3f, 0xff000000, px+0.98f*pw, (float)(yl+0.5*dy), px+1.02f*pw, (float)(yl+0.8*dy)));
			}
		}
		for(int sdy=0; sdy<subdivY; sdy++) {
			int spy = subPosY[sdy], sph = subSizY[sdy];
			double Yin = yLog[sdy] ? Math.log10(minY[sdy]) : minY[sdy];
			double Yax = yLog[sdy] ? Math.log10(maxY[sdy]) : maxY[sdy];
			double[] ticks = null;
			String[] tickmark = null;
			String   tickmarkFactor = "";
			if (yLog[sdy]) {
				ticks = JPlotMath.optimalLogarithmicTicks(minY[sdy], maxY[sdy]);
				tickmark = new String[ticks.length];
				for (int t = 2; t < ticks.length; t++) {
					tickmark[t] = "" + ticks[t] + "";
					ticks[t] = Math.log10(ticks[t]);
				}
			} else if (yTim[sdy]) {
				ticks = JPlotMath.optimalTimeTicks(minY[sdy], maxY[sdy], yTimUnit[sdy], yTimCal[sdy]);
				tickmark = new String[ticks.length];
				for (int t = 2; t < ticks.length; t++)
					tickmark[t] = DateTime.fromDouble(ticks[t], yTimUnit[sdy], yTimCal[sdy]).format(yTimFormat[sdy], yTimCal[sdy]);
			} else {
				ticks = JPlotMath.optimalLinearTicks(minY[sdy], maxY[sdy]);
				double vf = 1d / (ticks[0]);
				int decimal = (int) (1000d * ticks[1] + 0.5d);
				decimal = decimal % 1000 == 0 ? 0 : decimal % 100 == 0 ? 1 : decimal % 10 == 0 ? 2 : 3;
				tickmark = new String[ticks.length];
				for (int t = 2; t < ticks.length; t++) {
					if(decimal==0) tickmark[t] = ""+(int)(ticks[t]*vf+0.0005d - (ticks[t]*vf<0d?1:0));
					tickmark[t] = PApplet.nf((float) (ticks[t] * vf), 0, decimal).replace(",", ".");
				}
				double lvf = Math.log10(ticks[0]);
				if(Math.abs(lvf)>2.9d) {
					int ivf = (int) lvf - (lvf<0d ? 1 : 0);
					tickmarkFactor = JPlot.supportLatex ? "10$^{"+ivf+"}$" : "10^"+ivf;
				}
			}
			double[] tcpos = JPlotMath.map(ticks, Yin, Yax, spy + sph, spy);
			if (yAxInv[sdy])
				for (int t = 0; t < tcpos.length; t++)
					tcpos[t] = 2 * spy + sph - tcpos[t];
			if (pplot.isDebug()) {
				String tickStr = "", posStr = "";
				for (int t = 2; t < ticks.length; t++) {
					tickStr += ", " + tickmark[t];
					posStr += ", " + PApplet.nf((float) tcpos[t], 0, 2).replace(",", ".");
				}
				System.out.println("[DEBUG] JAxis-object: Ytickfactors={p10: " + ticks[0] + ", f: " + ticks[1] + "}");
				System.out.println("[DEBUG] JAxis-object: Ytickval={" + tickStr.substring(2) + "}");
				System.out.println("[DEBUG] JAxis-object: Ytickpos={" + posStr.substring(2) + "}");
			}
			if (yGrdOn) {
				for(int sdx=0; sdx<subdivX; sdx++) {
					float sx = subPosX[sdx]; if(sdx>0) sx -= 0.35f*(subPosX[sdx]-subPosX[sdx-1]-subSizX[sdx-1]);
					float ex = subPosX[sdx]+subSizX[sdx]; if(sdx+1<subdivX) ex += 0.35f*(subPosX[sdx+1]-subPosX[sdx]-subSizX[sdx]);
					for (int t = 2; t < ticks.length; t++)
						if (ticks[t] >= Math.min(Yin, Yax) && ticks[t] <= Math.max(Yin, Yax))
								axisgrid.addChild(new JLineShape(2f, 0xff999999, sx, (float) tcpos[t], ex, (float) tcpos[t]));
				}
			}
			if (yAxOn) {
				float tw = 0f;
				if (yTkOn) {
					for (int t = 2; t < ticks.length; t++)
						if (ticks[t] >= Math.min(Yin, Yax) && ticks[t] <= Math.max(Yin, Yax)) {
							tw = Math.max(tw, (float) txtsize * pplot.getGraphic().textWidth(tickmark[t]) / pplot.getGraphic().textSize);
							if(yDrawSide==LEFT || yDrawSide==BOTH) {
								axisgrid.addChild(new JLineShape(2f, 0xff000000, px - 0.02f * pw, (float) tcpos[t], px, (float) tcpos[t]));
								axisgrid.addChild(new JTextShape(tickmark[t], px - 0.03f * pw, (float) tcpos[t],
										(float) txtsize, RIGHT, CENTER, 0xff000000, 0, null));
							}
							if(yDrawSide==RIGHT || yDrawSide==BOTH) {
								axisgrid.addChild(new JLineShape(2f, 0xff000000, px+pw, (float) tcpos[t], px+ 1.02f * pw, (float) tcpos[t]));
								axisgrid.addChild(new JTextShape(tickmark[t], px + 1.03f * pw, (float) tcpos[t],
										(float) txtsize, LEFT, CENTER, 0xff000000, 0, null));
							}
						}
				}
				if (pplot.isDebug())
					System.out.println(
							"[DEBUG] JAxi-object: add y-axis title \""+titleY[sdy]+"\" and unit \""+unitY[sdy]+"\" with text size "+txtsize);
				String tytemp = ""+titleY[sdy];
				if(unitY[sdy].length()>0 || tickmarkFactor.length()>0) tytemp += (titleY[sdy].length()>0?" ":"")+"[";
				if(tickmarkFactor.length() > 0) tytemp += tickmarkFactor;
				if(unitY[sdy].length()>0) tytemp += (tickmarkFactor.length()>0?" ":"")+unitY[sdy];
				if(unitY[sdy].length()>0 || tickmarkFactor.length()>0) tytemp += "]";
				if(JPlot.supportLatex) {
					if(yDrawSide==LEFT || yDrawSide==BOTH) {
						axisgrid.addChild(new JLatexShape(tytemp, px - 0.03f * pw - tw, spy + 0.5f * sph, (float) (1.1d * txtsize),
								CENTER, BOTTOM, 0xff000000, ROTATE_COUNTERCLOCKWISE, null));
					}
					if(yDrawSide==RIGHT || yDrawSide==BOTH) {
						axisgrid.addChild(new JLatexShape(tytemp, px + 1.03f * pw + tw, spy + 0.5f * sph, (float) (1.1d * txtsize),
								CENTER, BOTTOM, 0xff000000, ROTATE_CLOCKWISE, null));
					}
				} else {
					if(yDrawSide==LEFT || yDrawSide==BOTH) {
						axisgrid.addChild(new JTextShape(tytemp, px - 0.03f * pw - tw, spy + 0.5f * sph, (float) (1.1d * txtsize),
								CENTER, BOTTOM, 0xff000000, ROTATE_COUNTERCLOCKWISE, null));
					}
					if(yDrawSide==RIGHT || yDrawSide==BOTH) {
						axisgrid.addChild(new JTextShape(tytemp, px + 1.03f * pw + tw, spy + 0.5f * sph, (float) (1.1d * txtsize),
								CENTER, BOTTOM, 0xff000000, ROTATE_CLOCKWISE, null));
					}
				}
			}
		}
		return axisgrid;
	}
	
	private void defineSubAxes(boolean create) {
		if(create) {
			int numSubAxs = subdivX*subdivY;
			layers = new LayerList[numSubAxs];
			for(int i=0; i<layers.length; i++)
				layers[i] = new LayerList();
			xRangeFix = new boolean[subdivX];
			yRangeFix = new boolean[subdivY];
			xAxInv    = new boolean[subdivX];
			yAxInv    = new boolean[subdivY];
			xLog      = new boolean[subdivX];
			yLog      = new boolean[subdivY];
			xTim      = new boolean[subdivX];
			yTim      = new boolean[subdivY];
			minX      = new double[subdivX];
			maxX      = new double[subdivX];
			minY      = new double[subdivY];
			maxY      = new double[subdivY];
			subPosX   = new int[subdivX];
			subPosY   = new int[subdivY];
			subSizX   = new int[subdivX];
			subSizY   = new int[subdivY];
			xTimUnit  = new String[subdivX];
			yTimUnit  = new String[subdivY];
			xTimCal   = new String[subdivX];
			yTimCal   = new String[subdivY];
			xTimFormat = new String[subdivX];
			yTimFormat = new String[subdivY];
			titleX    = new String[subdivX];
			titleY    = new String[subdivY];
			unitX     = new String[subdivX];
			unitY     = new String[subdivY];
		}
		double xw  = 1.0d*pw / (1.1d*subdivX-0.1d);
		double yh  = 1.0d*ph / (1.1d*subdivY-0.1d);
		double xwo = 1.1d*pw / (1.1d*subdivX-0.1d);
		double yho = 1.1d*ph / (1.1d*subdivY-0.1d);
		for(int sdx=0; sdx<subdivX; sdx++) {
			subPosX[sdx] = (int) (xwo*sdx + 0.5d);
			subSizX[sdx] = (int) (xwo*sdx + xw + 0.5d) - subPosX[sdx];
			subPosX[sdx] += px;
			xRangeFix[sdx] = false;
			xAxInv[sdx]    = false;
			xLog[sdx]      = false;
			xTim[sdx]      = false;
			minX[sdx] = -1d;
			maxX[sdx] =  1d;
			xTimUnit[sdx]   = "days since 1850-01-01 00:00:00";
			xTimCal[sdx]    = "gregorian";
			xTimFormat[sdx] = "dd.mm.yyyy";
			titleX[sdx] = "";
			unitX[sdx]  = "";
		}
		for(int sdy=0; sdy<subdivY; sdy++) {
			subPosY[sdy] = (int) (yho*sdy + 0.5d);
			subSizY[sdy] = (int) (yho*sdy + yh + 0.5d) - subPosY[sdy];
			subPosY[sdy] += py;
			yRangeFix[sdy] = false;
			yAxInv[sdy]    = false;
			yLog[sdy]      = false;
			yTim[sdy]      = false;
			minY[sdy] = -1d;
			maxY[sdy] =  1d;
			yTimUnit[sdy]   = "days since 1850-01-01 00:00:00";
			yTimCal[sdy]    = "gregorian";
			yTimFormat[sdy] = "dd.mm.yyyy";
			titleY[sdy] = "";
			unitY[sdy]  = "";
		}
	}
	
	private class LayerList extends ArrayList<JPlotsLayer> {
		private static final long serialVersionUID = -6120624294306257474L;
	}


}
