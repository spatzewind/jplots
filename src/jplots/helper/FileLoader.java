package jplots.helper;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.math.MathContext;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

import jplots.JPlot;

public class FileLoader {

	static Map<String, File> resources = new HashMap<String, File>();
	
	public static File loadResourceFile(String resource) {
		if(resources.containsKey(resource))
			return resources.get(resource);
		File file = null;
		URL res = JPlot.class.getResource(resource);
		System.out.println("[FILELOADER] URL res = "+res);
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
		    //this will probably work in your IDE, but not from a JAR
		    file = new File(res.getFile());
		}

		if (file != null && !file.exists()) {
		    throw new RuntimeException("Error: File " + file + " not found!");
		}
		
		resources.put(resource, file);
		return file;
	}
	

	public static File loadResourceShapeFile(String resource) {
		if(resource==null) {
			System.err.println("cannot read file <null>");
			return null;
		}
		if(resources.containsKey(resource))
			return resources.get(resource);
		int id = Math.max(
				Math.max(resource.lastIndexOf(".shp"),
						resource.lastIndexOf(".shx")),
				Math.max(resource.lastIndexOf(".dbf"),
						resource.lastIndexOf(".prj"))
				);
		if(id<1) id = resource.length();
		String resourceN = resource.substring(0, id);
		File file = null;
		URL res = JPlot.class.getResource(resourceN+".shp");
		//System.out.println("[FILELOADER] URL res = "+res);
		if (res!=null && res.getProtocol().equals("jar")) {
			//System.out.println("[FILELOADER] with JAR-protocol:");
		    try {
		        InputStream input = JPlot.class.getResourceAsStream(resourceN+".shp");
		        file = File.createTempFile("tempshape", ".shp");
		        String fs = file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."));
		        OutputStream out = new FileOutputStream(file);
		        int read;
		        byte[] bytes = new byte[1024];
		        //copy *.shp to temp/
		        if(input!=null) {
			        while ((read = input.read(bytes)) != -1) {
			            out.write(bytes, 0, read);
			        }
			        out.close();
			        input.close();
			        file.deleteOnExit();
		        }
		        //copy *.shx to temp/
		        input = JPlot.class.getResourceAsStream(resourceN+".shx");
		        if(input!=null) {
			        File filx = new File(fs+".shx");
			        out = new FileOutputStream(filx);
			        while ((read = input.read(bytes)) != -1) {
			            out.write(bytes, 0, read);
			        }
			        out.close();
			        input.close();
			        filx.deleteOnExit();
		        }
		        //copy *.prj to temp/
		        input = JPlot.class.getResourceAsStream(resourceN+".prj");
		        if(input!=null) {
			        File filj = new File(fs+".prj");
			        out = new FileOutputStream(filj);
			        while ((read = input.read(bytes)) != -1) {
			            out.write(bytes, 0, read);
			        }
			        out.close();
			        input.close();
			        filj.deleteOnExit();
		        }
		        //copy *.dbf to temp/
		        input = JPlot.class.getResourceAsStream(resourceN+".dbf");
		        if(input!=null) {
			        File filb = new File(fs+".dbf");
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
			//System.out.println("[FILELOADER] with normal protocol");
		    //this will probably work in your IDE, but not from a JAR
		    file = new File(resourceN+".shp");
		}

		if (file != null && !file.exists()) {
		    throw new RuntimeException("Error: File " + file + " not found!");
		}
		
		resources.put(resource, file);
		return file;
	}
}
