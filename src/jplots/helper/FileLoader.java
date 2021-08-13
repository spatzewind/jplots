package jplots.helper;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
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
		if(resources.containsKey(resource))
			return resources.get(resource);
		File file = null;
		URL res = JPlot.class.getResource(resource+".shp");
		//System.out.println("[FILELOADER] URL res = "+resP);
		if (res.getProtocol().equals("jar")) {
			//System.out.println("[FILELOADER] with JAR-protocol:");
		    try {
		        InputStream input = JPlot.class.getResourceAsStream(resource+".shp");
		        file = File.createTempFile("tempfile", ".shp");
		        OutputStream out = new FileOutputStream(file);
		        int read;
		        byte[] bytes = new byte[1024];

		        while ((read = input.read(bytes)) != -1) {
		            out.write(bytes, 0, read);
		        }
		        out.close();
		        input.close();
		        file.deleteOnExit();
		        String fs = file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf("."));
		        File filx = new File(fs+".shx");
		        input = JPlot.class.getResourceAsStream(resource+".shx");
		        out = new FileOutputStream(filx);
		        while ((read = input.read(bytes)) != -1) {
		            out.write(bytes, 0, read);
		        }
		        out.close();
		        input.close();
		        filx.deleteOnExit();
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
}
