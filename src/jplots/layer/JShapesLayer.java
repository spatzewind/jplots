package jplots.layer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.geotools.data.FeatureSource;
import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import jplots.JAxis;
import jplots.JPlot;
import jplots.maths.AffineBuilder;
import jplots.maths.JDPoint;
import jplots.maths.JDPolygon;
import jplots.maths.JPlotMath;
import jplots.shapes.JEllipseShape;
import jplots.shapes.JGroupShape;
import jplots.shapes.JLineShape;
import jplots.shapes.JPlotShape;
import jplots.shapes.JPolygonShape;
import processing.core.PGraphics;

public class JShapesLayer extends JPlotsLayer {
	
	private static CoordinateReferenceSystem dest = null;
	static {
		try {
			dest = CRS.decode("EPSG:4326");
		} catch (FactoryException fe) {
			fe.printStackTrace();
		}
	}
	
	private Map<String, String> connect;
	private String shapeType;
	private CoordinateReferenceSystem user_crs;
	
	public JShapesLayer(String filepath, String type) {
		this(new File(filepath), type, null);
	}
	public JShapesLayer(File filepath, String type) {
		this(filepath, type, null);
	}
	public JShapesLayer(File filepath, String type, CoordinateReferenceSystem crs) {
		super();
		pc = 0xff99ccff;
		lc = 0xff336699;
		connect = new HashMap<String, String>();
		connect.put("url", filepath.getAbsolutePath());
		shapeType = type;
		user_crs = crs;
	}
	
	public JShapesLayer(String filepath, String type, int epsg) {
		this(new File(filepath), type, epsg);
	}
	public JShapesLayer(File filepath, String type, int epsg) {
		super();
		pc = 0xff99ccff;
		lc = 0xff336699;
		connect = new HashMap<String, String>();
		connect.put("url", filepath.getAbsolutePath());
		shapeType = type;
		try {
			if(epsg>=0)
				user_crs = CRS.decode("EPSG:"+epsg);
			else
				user_crs = null;
		} catch (NoSuchAuthorityCodeException nsace) {
			nsace.printStackTrace();
		} catch (FactoryException fe) {
			fe.printStackTrace();
		}
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
				int fidx = 0;
				while (iterator.hasNext()) {
					int ppcc = pc, llcc = lc;
					if(!singleFillColour) ppcc = pcs[fidx % pcs.length];
					if(!singleLineColour) llcc = lcs[fidx % lcs.length];
					SimpleFeature feature = iterator.next();
					SimpleFeatureType type = feature.getFeatureType();
					//GeometryAttribute sourceGeometry = feature.getDefaultGeometryProperty();
					//System.out.println(sourceGeometry);
					Geometry geom = (Geometry) feature.getDefaultGeometry();
					if(ax.isGeoAxis()) {
						CoordinateReferenceSystem crs = type.getCoordinateReferenceSystem();
						if(crs==null && user_crs!=null)
							crs = user_crs;
						if(crs==null) {
							System.err.println("No CoordinateReferenceSystem defined!");
							continue;
						} else {
//							if(debug)
//								System.out.println("[DEBUG] JShapeLayer: found transform "+crs);
						}
						MathTransform transform = CRS.findMathTransform(crs, dest, true);
						geom = JTS.transform(geom, transform);
					}
					Coordinate[] points = geom.getCoordinates();
					double[][] coords = new double[points.length][2];
					for(int c=0; c<points.length; c++) {
						coords[c][0] = points[c].x;
						coords[c][1] = points[c].y;
						if(ax.isGeoAxis())
							coords[c] = ax.getGeoProjection().fromLATLONtoPROJ(points[c].y, points[c].x, true);
						if(ax.isXlogAxis())
							coords[c][0] = Math.log10(coords[c][0]);
						if(ax.isYlogAxis())
							coords[c][1] = Math.log10(coords[c][1]);
					}
					if(shapeType.equalsIgnoreCase("point")) {
						JPlotShape.noFill(); JPlotShape.stroke(llcc); JPlotShape.strokeWeight((float)lw);
						for(int c=0; c<points.length; c++) {
							double x1 = p[0]+xs*(invertAxisX ? Xax-coords[c][0] : coords[c][0]-Xin);
							double y1 = p[1]+ys*(invertAxisY ? coords[c][1]-Yin : Yax-coords[c][1]);
							if(x1<p[0] || x1>=p[0]+p[2]) continue;
							if(y1<p[1] || y1>=p[1]+p[3]) continue;
							s.addChild(new JEllipseShape((float)x1, (float)y1, 1f, 1f));
						}
					}
					if(shapeType.equalsIgnoreCase("line")) {
						JPlotShape.noFill(); JPlotShape.stroke(llcc); JPlotShape.strokeWeight((float)lw);
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
							if(x1>=p[0] && x1<=p[0]+p[2] && x2>=p[0] && x2<=p[0]+p[2] &&
									y1>=p[1] && y1<=p[1]+p[3] && y2>=p[1] && y2<=p[1]+p[3])
								s.addChild(new JLineShape((float)x1, (float)y1, (float)x2, (float)y2));
						}
					}
					if(shapeType.equalsIgnoreCase("polygon")) {
						JPlotShape.fill(ppcc); JPlotShape.stroke(llcc); JPlotShape.strokeWeight((float)lw);
//				        double eps  = Math.min(p[2], p[3]) * 1.0e-9d;
//				        double eps2 = eps*0.0001d;
				        AffineBuilder affine = new AffineBuilder()
				        		.scale(invertAxisX?-1d:1d, invertAxisY?1d:-1d)
				        		.translate(invertAxisX?Xax:-Xin, invertAxisY?-Yin:Yax)
				        		.scale(xs, ys)
				        		.translate(p[0], p[1]);
						JDPolygon poly = new JDPolygon(coords);
				        if(debug) System.out.println("[JShapesLayer] found geometry: p["+
				        		poly.c.length+"points,"+poly.area()+"area]");
				        //System.out.print("    ("+(fidx<9?"  ":fidx<99?" ":"")+(fidx+1)+") ");
				        List<JDPolygon> polys = poly.splitByMapBorder(ax);
//				        List<JDPolygon> polys = new ArrayList<>();
//				        polys.add(poly);
				        for(int p_i=polys.size()-1; p_i>=0; p_i--) {
					        polys.get(p_i).affine(affine.getMatrix());
				        	polys.addAll(polys.get(p_i).intersectsAABB(p[0],p[1],  p[0]+p[2],p[1]+p[3]));
				        	polys.remove(p_i);
				        }
				        int pidx = 0;
				        for(JDPolygon pg: polys) {
				        	if(pg.c.length==0) continue;
				        	if(debug)
				        		System.out.println("[JShapesLayer] add polygon p["+pg.c.length+"points,"+pg.area()+"area]");
				        	s.addChild(new JPolygonShape(pg.getCoords(), true, true));
//				        	if(fidx==0 && pidx==0) {
//				        		double minx = Double.POSITIVE_INFINITY, maxx = Double.NEGATIVE_INFINITY;
//				        		for(JDPoint pt: pg.c) {
//				        			if(pt.x < minx) minx = pt.x;
//				        			if(pt.x > maxx) maxx = pt.x;
//				        		}
//				        		System.out.println("     Polygon extent: x={"+minx+" ... "+maxx+"}");
//				        	}
				        	pidx++;
				        }
					}
					fidx++;
				}
			} catch (NoSuchAuthorityCodeException nsace) {
				nsace.printStackTrace();
			} catch (FactoryException fe) {
				fe.printStackTrace();
			} catch (MismatchedDimensionException mde) {
				mde.printStackTrace();
			} catch (TransformException te) {
				te.printStackTrace();
			} finally {
				iterator.close();
			}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

}
