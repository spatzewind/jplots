package jplots.shapes;

import java.io.PrintWriter;
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
	
	public int childCount() {
		return childs.size();
	}
	
	@Override
	public void printStack(PrintWriter pw, String off) {
		pw.println(off+this.getClass().getSimpleName()+" {");
		for(JPlotShape child: childs)
			child.printStack(pw, "   "+off);
		pw.println(off+"}");
	}
}
