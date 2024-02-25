package jplots.shapes;

import java.util.ArrayList;
import java.util.List;

import jplots.helper.GeometryTools;
import jplots.maths.JDPolygon;

public class JDGeometry {

	private List<JDPolygon> polys;
	
	public JDGeometry() {
		polys = new ArrayList<>();
	}

	public JDGeometry(JDPolygon... polygons) {
		polys = new ArrayList<>();
		for (JDPolygon p : polygons)
			addPolygon(p);
	}

	public JDGeometry(List<JDPolygon> polygons) {
		polys = new ArrayList<>(polygons);
	}

	public void add(JDPolygon... polygons) {
		for (JDPolygon p : polygons)
			addPolygon(p);
	}

	public void add(List<JDPolygon> polygons) {
		polys.addAll(polygons);
	}

	public List<JDPolygon> getPolygons() {
		return polys;
	}

	private void addPolygon(JDPolygon p) {
		boolean addToExistingPoly = false;
		for (int pi = polys.size() - 1; pi >= 0 && !addToExistingPoly; pi--) {
			JDPolygon temp = GeometryTools.union(polys.get(pi), p, Double.NaN);
			if(temp==null) continue;
			polys.set(pi, temp);
			addToExistingPoly = true;
		}
		if (!addToExistingPoly)
			polys.add(p);
	}
}
