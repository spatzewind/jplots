package jplots.layer;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.FeatureSource;
import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.geometry.jts.Geometries;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.MultiLineString;
import org.opengis.feature.GeometryAttribute;
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
	
	private File shapeFile;
	private Map<String, String> connect;
	private String shapeType;
	
	public JShapesLayer(String filepath, String type) {
		shapeFile = new File(filepath);
		connect = new HashMap<String, String>();
		connect.put("url", filepath);
		shapeType = type;
	}
	public JShapesLayer(File filepath, String type) {
		shapeFile = filepath;
		connect = new HashMap<String, String>();
		connect.put("url", shapeFile.getAbsolutePath());
		shapeType = type;
	}

	@Override
	public void createRasterImg(JPlot plot, PGraphics g) {
		// TODO Auto-generated method stub

	}

	@Override
	public void createVectorImg(JAxis ax, int layernum, JGroupShape s) {
		int[] p = ax.getSize();
		double xs = p[2]/(maxX-minX), ys = p[3]/(maxY-minY);
		boolean debug = ax.getPlot().isDebug();
		if(debug)
			System.out.println("[DEBUG] JShapeLayer: begin reading shape file \""+connect.get("url")+"\"");
		try {
			//DataStore dataStore = DataStoreFinder.getDataStore(connect);
			FileDataStore dataStore = FileDataStoreFinder.getDataStore(shapeFile);
			System.out.println(dataStore);
			String[] typeNames = dataStore.getTypeNames();
			System.out.println(typeNames);
			String typeName = typeNames[0];
			System.out.println("Reading content " + typeName);

			FeatureSource<SimpleFeatureType, SimpleFeature> featureSource = dataStore.getFeatureSource(typeName);
			FeatureCollection<SimpleFeatureType, SimpleFeature> collection = featureSource.getFeatures();
			FeatureIterator<SimpleFeature> iterator = collection.features();

			try {
				while (iterator.hasNext()) {
					SimpleFeature feature = iterator.next();
					GeometryAttribute sourceGeometry = feature.getDefaultGeometryProperty();
					//System.out.println(sourceGeometry);
					Geometry geom = (Geometry) feature.getDefaultGeometry();
					Coordinate[] points = geom.getCoordinates();
					double[][] coords = new double[points.length][2];
					for(int c=0; c<points.length; c++) {
						coords[c][0] = points[c].x;
						coords[c][1] = points[c].y;
						if(ax.isGeoAxis())
							coords[c] = ax.getGeoProjection().fromLATLONtoPROJ(coords[c][0], coords[c][1], true);
					}
					JPlotShape.noFill(); JPlotShape.stroke(lc); JPlotShape.strokeWeight((float)lw);
					for(int c=1; c<points.length; c++) {
						double x1 = p[0]+xs*(invertAxisX ? maxX-coords[c-1][0] : coords[c-1][0]-minX);
						double x2 = p[0]+xs*(invertAxisX ? maxX-coords[ c ][0] : coords[ c ][0]-minX);
						double y1 = p[1]+ys*(invertAxisY ? coords[c-1][1]-minY : maxY-coords[c-1][1]);
						double y2 = p[1]+ys*(invertAxisY ? coords[ c ][1]-minY : maxY-coords[ c ][1]);
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
			} finally {
				iterator.close();
			}

		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

}
