package jplots.shapes;

import java.util.ArrayList;
import java.util.List;

import jplots.JPlot;
import processing.core.PGraphics;

public class JGroupShape extends JPlotShape {

	private List<JPlotShape> childs;

	public JGroupShape() {
		childs = new ArrayList<>();
	}

	@Override
	public void draw(JPlot plot, PGraphics g) {
		for (JPlotShape ppe : childs)
			ppe.draw(plot, g);
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
