package jplots.helper;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;
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

import jplots.JPlot;

public class FileLoader {

	static Map<String, File> resources = new HashMap<>();
	public static Map<String, Coordinate[][]> shapefiles = new HashMap<>();
	private static CoordinateReferenceSystem dest = null;
	static {
		try {
			dest = CRS.decode("EPSG:4326");
		} catch (FactoryException fe) {
			fe.printStackTrace();
		}
	}

	public static File loadResourceFile(String resource) {
		if (resources.containsKey(resource))
			return resources.get(resource);
		File file = null;
		URL res = JPlot.class.getResource(resource);
		System.out.println("[FILELOADER] URL res = " + res);
		if (res.getProtocol().equals("jar")) {
			System.out.println("[FILELOADER] with JAR-protocol:");
			try {
				String fs = res.toExternalForm();
				InputStream input = JPlot.class.getResourceAsStream(resource);
				file = File.createTempFile("tempfile", fs.substring(fs.lastIndexOf(".")));
				OutputStream out = new FileOutputStream(file);
				int read;
				byte[] bytes = new byte[1024];

				while ((read = input.read(bytes)) != -1) {
					out.write(bytes, 0, read);
				}
				out.close();
				file.deleteOnExit();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		} else {
			System.out.println("[FILELOADER] with normal protocol");
			// this will probably work in your IDE, but not from a JAR
			file = new File(res.getFile());
		}

		if (file != null && !file.exists()) {
			throw new RuntimeException("Error: File " + file + " not found!");
		}

		resources.put(resource, file);
		return file;
	}

	public static String loadResourceShapeFile(String resource, int user_epsg) {
		if (resource == null) {
			System.err.println("cannot read file <null>");
			return null;
		}
		if (shapefiles.containsKey(resource))
			return resource;
		CoordinateReferenceSystem crs = null;
		try {
			crs = CRS.decode("EPSG:" + user_epsg);
		} catch (NoSuchAuthorityCodeException nsace) {
			nsace.printStackTrace();
		} catch (FactoryException fe) {
			fe.printStackTrace();
		}
		return loadResourceShapeFile(resource, crs);
	}

	public static String loadResourceShapeFile(String resource, CoordinateReferenceSystem user_crs) {
		if (resource == null) {
			System.err.println("cannot read file <null>");
			return null;
		}
		if (shapefiles.containsKey(resource))
			return resource;

		int id = Math.max(Math.max(resource.lastIndexOf(".shp"), resource.lastIndexOf(".shx")),
				Math.max(resource.lastIndexOf(".dbf"), resource.lastIndexOf(".prj")));
		if (id < 1)
			id = resource.length();
		String resourceN = resource.substring(0, id);
		File file = null;
		URL res = JPlot.class.getResource(resourceN + ".shp");
		// System.out.println("[FILELOADER] URL res = "+res);
		if (res != null && res.getProtocol().equals("jar")) {
			// System.out.println("[FILELOADER] with JAR-protocol:");
			try {
				InputStream input = JPlot.class.getResourceAsStream(resourceN + ".shp");
				file = File.createTempFile("tempshape", ".shp");
				String fs = file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."));
				OutputStream out = new FileOutputStream(file);
				int read;
				byte[] bytes = new byte[1024];
				// copy *.shp to temp/
				if (input != null) {
					while ((read = input.read(bytes)) != -1) {
						out.write(bytes, 0, read);
					}
					out.close();
					input.close();
					file.deleteOnExit();
				}
				// copy *.shx to temp/
				input = JPlot.class.getResourceAsStream(resourceN + ".shx");
				if (input != null) {
					File filx = new File(fs + ".shx");
					out = new FileOutputStream(filx);
					while ((read = input.read(bytes)) != -1) {
						out.write(bytes, 0, read);
					}
					out.close();
					input.close();
					filx.deleteOnExit();
				}
				// copy *.prj to temp/
				input = JPlot.class.getResourceAsStream(resourceN + ".prj");
				if (input != null) {
					File filj = new File(fs + ".prj");
					out = new FileOutputStream(filj);
					while ((read = input.read(bytes)) != -1) {
						out.write(bytes, 0, read);
					}
					out.close();
					input.close();
					filj.deleteOnExit();
				}
				// copy *.dbf to temp/
				input = JPlot.class.getResourceAsStream(resourceN + ".dbf");
				if (input != null) {
					File filb = new File(fs + ".dbf");
					out = new FileOutputStream(filb);
					while ((read = input.read(bytes)) != -1) {
						out.write(bytes, 0, read);
					}
					out.close();
					input.close();
					filb.deleteOnExit();
				}
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		} else {
			// System.out.println("[FILELOADER] with normal protocol");
			// this will probably work in your IDE, but not from a JAR
			file = new File(resourceN + ".shp");
		}

		if (file != null && !file.exists()) {
			System.err.println("Error: File " + file + " not found!");
			// throw new RuntimeException("Error: File " + file + " not found!");
			return null;
		}

		// System.out.println("[DEBUG] JShapeLayer: begin reading shape file
		// \""+connect.get("url")+"\"");
		FeatureCollection<SimpleFeatureType, SimpleFeature> collection = null;
		FeatureIterator<SimpleFeature> iterator = null;
		try {
			// DataStore dataStore = DataStoreFinder.getDataStore(connect);
			FileDataStore dataStore = FileDataStoreFinder.getDataStore(file);
			if (dataStore == null) {
				System.err.println("Error: Cannot find dataStore for file " + file + "");
				return null;
			}
			// System.out.println("[DEBUG] JShapeLayer: dataStore is <"+dataStore+">");
			String[] typeNames = dataStore.getTypeNames();
			if (typeNames.length == 0) {
				System.err.println("Error: File " + file + " does not contain readable data!");
				return null;
			}
//			if(debug) {
//				String names = "";
//				for(String n: typeNames) names += (names.length()==0?"":", ")+n;
//				System.out.println("[DEBUG] JShapeLayer: found typeName = <"+names+">");
//			}
			// if(typeNames.length==0)
			// return;
			if (typeNames.length > 0) {
				String typeName = typeNames[0];
				// System.out.println("[DEBUG] JShapeLayer: reading content <"+typeName+">");

				FeatureSource<SimpleFeatureType, SimpleFeature> featureSource = dataStore.getFeatureSource(typeName);
				collection = featureSource.getFeatures();
				iterator = collection.features();

				List<Coordinate[]> geometries = new ArrayList<>();
				while (iterator.hasNext()) {
					SimpleFeature feature = iterator.next();
					SimpleFeatureType type = feature.getFeatureType();
					// GeometryAttribute sourceGeometry = feature.getDefaultGeometryProperty();
					// System.out.println(sourceGeometry);
					Geometry geom = (Geometry) feature.getDefaultGeometry();
					CoordinateReferenceSystem crs = type.getCoordinateReferenceSystem();
					if (crs == null && user_crs != null)
						crs = user_crs;
					if (crs == null) {
						System.err.println("No CoordinateReferenceSystem defined!");
						continue;
					} else {
//								if(debug)
//									System.out.println("[DEBUG] JShapeLayer: found transform "+crs);
					}
					MathTransform transform = CRS.findMathTransform(crs, dest, true);
					geom = JTS.transform(geom, transform);
					geometries.add(geom.getCoordinates());
				}
				shapefiles.put(resource, geometries.toArray(new Coordinate[0][]));
			}
		} catch (IOException ioe) {
			ioe.printStackTrace();
			return null;
		} catch (FactoryException fe) {
			fe.printStackTrace();
		} catch (MismatchedDimensionException mde) {
			mde.printStackTrace();
		} catch (TransformException te) {
			te.printStackTrace();
		} finally {
			iterator.close();
		}
		return resource;
	}
}
