package pplots.gfx;

import java.util.ArrayList;
import java.util.List;

import processing.core.PGraphics;

public class PGroupShape extends PPlotShape {
	
	private List<PPlotShape> childs;
	
	public PGroupShape() {
		childs = new ArrayList<PPlotShape>();
	}
	
	@Override
	public void draw(PGraphics g) {
		for(PPlotShape ppe: childs)
			ppe.draw(g);
	}
	
	public void addChild(PPlotShape child) {
		childs.add(child);
	}
	public void removeChild(int idx) {
		childs.remove(idx);
	}
	public void removeChild(PPlotShape child) {
		childs.remove(child);
	}
	public void removeAllChildren() {
		childs.clear();
	}

}
