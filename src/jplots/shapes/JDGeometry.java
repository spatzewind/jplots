package jplots.shapes;

import java.util.ArrayList;
import java.util.List;

import jplots.maths.JDPolygon;

public class JDGeometry {

	private List<JDPolygon> polys;
	
	public JDGeometry() {
		polys = new ArrayList<JDPolygon>();
	}
	public JDGeometry(JDPolygon... polygons) {
		polys = new ArrayList<JDPolygon>();
		for(JDPolygon p: polygons)
			polys.add(p);
	}
	public JDGeometry(List<JDPolygon> polygons) {
		polys = new ArrayList<JDPolygon>(polygons);
	}
	
	public void add(JDPolygon... polygons) {
		for(JDPolygon p: polygons)
			polys.add(p);
	}
	public void add(List<JDPolygon> polygons) {
		polys.addAll(polygons);
	}
	
	public List<JDPolygon> getPolygons() {
		return polys;
	}
}
