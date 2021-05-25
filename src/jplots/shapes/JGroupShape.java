package jplots.shapes;

import java.util.ArrayList;
import java.util.List;

import processing.core.PGraphics;

public class JGroupShape extends JPlotShape {
	
	private List<JPlotShape> childs;
	
	public JGroupShape() {
		childs = new ArrayList<JPlotShape>();
	}
	
	@Override
	public void draw(PGraphics g) {
		for(JPlotShape ppe: childs)
			ppe.draw(g);
	}
	
	public void addChild(JPlotShape child) {
		childs.add(child);
	}
	public void removeChild(int idx) {
		childs.remove(idx);
	}
	public void removeChild(JPlotShape child) {
		childs.remove(child);
	}
	public void removeAllChildren() {
		childs.clear();
	}

}
