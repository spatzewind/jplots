package jplots.layer;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.geotools.data.FeatureSource;
import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

import jplots.JAxis;
import jplots.JPlot;
import jplots.maths.JPlotMath;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import processing.core.PGraphics;

public class JShapesLayer extends JPlotsLayer {
	
	private Map<String, String> connect;
	private String shapeType;
	
	public JShapesLayer(String filepath, String type) {
		File f = new File(filepath);
		connect = new HashMap<String, String>();
		connect.put("url", f.getAbsolutePath());
		shapeType = type;
	}
	public JShapesLayer(File filepath, String type) {
		connect = new HashMap<String, String>();
		connect.put("url", filepath.getAbsolutePath());
		shapeType = type;
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		double Xin = ax.isXlogAxis()?Math.log10(minX):minX, Xax = ax.isXlogAxis()?Math.log10(maxX):maxX;
		double Yin = ax.isYlogAxis()?Math.log10(minY):minY, Yax = ax.isYlogAxis()?Math.log10(maxY):maxY; 
		double xs = p[2]/(Xax-Xin), ys = p[3]/(Yax-Yin);
		boolean debug = ax.getPlot().isDebug();
		if(debug)
			System.out.println("[DEBUG] JShapeLayer: begin reading shape file \""+connect.get("url")+"\"");
		try {
			//DataStore dataStore = DataStoreFinder.getDataStore(connect);
			FileDataStore dataStore = FileDataStoreFinder.getDataStore(new File(connect.get("url")));
			if(debug)
				System.out.println("[DEBUG] JShapeLayer: dataStore is <"+dataStore+">");
			if(dataStore==null)
				return;
			String[] typeNames = dataStore.getTypeNames();
			if(typeNames==null)
				return;
			if(debug) {
				String names = "";
				for(String n: typeNames) names += (names.length()==0?"":", ")+n;
				System.out.println("[DEBUG] JShapeLayer: found typeName = <"+names+">");
			}
			if(typeNames.length==0)
				return;
			String typeName = typeNames[0];
			if(debug)
				System.out.println("[DEBUG] JShapeLayer: reading content <"+typeName+">");

			FeatureSource<SimpleFeatureType, SimpleFeature> featureSource = dataStore.getFeatureSource(typeName);
			FeatureCollection<SimpleFeatureType, SimpleFeature> collection = featureSource.getFeatures();
			FeatureIterator<SimpleFeature> iterator = collection.features();

			try {
				while (iterator.hasNext()) {
					SimpleFeature feature = iterator.next();
					//GeometryAttribute sourceGeometry = feature.getDefaultGeometryProperty();
					//System.out.println(sourceGeometry);
					Geometry geom = (Geometry) feature.getDefaultGeometry();
					Coordinate[] points = geom.getCoordinates();
					double[][] coords = new double[points.length][2];
					for(int c=0; c<points.length; c++) {
						coords[c][0] = points[c].x;
						coords[c][1] = points[c].y;
						if(ax.isGeoAxis())
							coords[c] = ax.getGeoProjection().fromLATLONtoPROJ(coords[c][0], coords[c][1], true);
						if(ax.isXlogAxis())
							coords[c][0] = Math.log10(coords[c][0]);
						if(ax.isYlogAxis())
							coords[c][1] = Math.log10(coords[c][1]);
					}
					if(shapeType.equalsIgnoreCase("line")) {
						JPlotShape.noFill(); JPlotShape.stroke(lc); JPlotShape.strokeWeight((float)lw);
						for(int c=1; c<points.length; c++) {
							double x1 = p[0]+xs*(invertAxisX ? Xax-coords[c-1][0] : coords[c-1][0]-Xin);
							double x2 = p[0]+xs*(invertAxisX ? Xax-coords[ c ][0] : coords[ c ][0]-Xin);
							double y1 = p[1]+ys*(invertAxisY ? coords[c-1][1]-Yin : Yax-coords[c-1][1]);
							double y2 = p[1]+ys*(invertAxisY ? coords[ c ][1]-Yin : Yax-coords[ c ][1]);
							if(x1<p[0]      && x2>=p[0])      { y1 = JPlotMath.map(p[0],      x1, x2, y1, y2); }
							if(x1>p[0]+p[2] && x2<=p[0]+p[2]) { y1 = JPlotMath.map(p[0]+p[2], x1, x2, y1, y2); }
							if(x2<p[0]      && x1>=p[0])      { y2 = JPlotMath.map(p[0],      x1, x2, y1, y2); }
							if(x2>p[0]+p[2] && x1<=p[0]+p[2]) { y2 = JPlotMath.map(p[0]+p[2], x1, x2, y1, y2); }
							
							if(y1<p[1]      && y2>=p[1])      { x1 = JPlotMath.map(p[1],      y1, y2, x1, x2); }
							if(y1>p[1]+p[3] && y2<=p[1]+p[3]) { x1 = JPlotMath.map(p[1]+p[3], y1, y2, x1, x2); }
							if(y2<p[1]      && y1>=p[1])      { x2 = JPlotMath.map(p[1],      y1, y2, x1, x2); }
							if(y2>p[1]+p[3] && y1<=p[1]+p[3]) { x2 = JPlotMath.map(p[1]+p[3], y1, y2, x1, x2); }
							if(x1>=p[0] && x1<=p[0]+p[2] && x2>=p[0] && x2<=p[0]+p[2] && y1>=p[1] && y1<=p[1]+p[3] && y2>=p[1] && y2<=p[1]+p[3])
								s.addChild(new JLineShape((float)x1, (float)y1, (float)x2, (float)y2));
						}
					}
				}
			} finally {
				iterator.close();
			}

		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

}
