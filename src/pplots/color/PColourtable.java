package pplots.color;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import processing.core.PApplet;

public abstract class PColourtable {
	
	public static List<String> ctpaths;
	public static Map<String,PColourtable> pctables;
	static {
		ctpaths = new ArrayList<String>();
		pctables = new HashMap<String,PColourtable>();
		pctables.put("default", new LinearSegmentedPColourtable(0xff00007f, 0xff7f0000, 0x007f7f7f,
				new double[] {0d, 0.4d,0.6d,1d},
				new int[][] {{0xff0000ff,0xff00ffff},{0xff00ffff,0xffffff00},{0xffaa0000}}));
	}
	
	public static PColourtable load(String path) {
		String end = path.substring(path.lastIndexOf("."));
//		if(".rgb".equals(end))
//			return readRGB(path);
		if(".cpt".equals(end))
			return readCPT(path);
		System.err.println("Cannot determine type of colourtable. Check the existence of\n"+
			"filename extension to be one of following: .cpt, .rgb");
		return null;
	}
//	private static PColourtable readRGB(String path) {
//		String[] lines = PApplet.loadStrings(new File(path));
//		if(lines == null)
//			return null;
//		return null;
//	}
	private static PColourtable readCPT(String path) {
		String[] lines = PApplet.loadStrings(new File(path));
		if(lines==null)
			return null;
		int col_count = 0;
		for(int l=0; l<lines.length; l++) {
			if(lines[l].length()<3)
				continue;
			if(lines[l].startsWith("#"))
				continue;
			if(lines[l].startsWith("F") || lines[l].startsWith("B") || lines[l].startsWith("N"))
				continue;
			col_count++;
		}
		double[] positions = new double[col_count+1];
		int[][] colours = new int[col_count][2];
		int[] specColours = {0xff000000, 0xffffffff, 0x007f7f7f};
		double minpos = Double.POSITIVE_INFINITY, maxpos = Double.NEGATIVE_INFINITY;
		int idx=0;
		for(int l=0; l<lines.length; l++) {
			if(lines[l].length()<3)
				continue;
			if(lines[l].startsWith("#"))
				continue;
			String[] pcol = lines[l].split("\\s*");
			if(lines[l].startsWith("F")) {
				specColours[1] = 0xff000000 | (Integer.parseInt(pcol[1])<<16) | (Integer.parseInt(pcol[2])<<8) | Integer.parseInt(pcol[3]);
				continue;
			}
			if(lines[l].startsWith("B")) {
				specColours[0] = 0xff000000 | (Integer.parseInt(pcol[1])<<16) | (Integer.parseInt(pcol[2])<<8) | Integer.parseInt(pcol[3]);
				continue;
			}
			if(lines[l].startsWith("N")) {
				specColours[2] = 0xff000000 | (Integer.parseInt(pcol[1])<<16) | (Integer.parseInt(pcol[2])<<8) | Integer.parseInt(pcol[3]);
				continue;
			}
			positions[idx]   = Double.parseDouble(pcol[0]);
			positions[idx+1] = Double.parseDouble(pcol[4]);
			colours[idx][0]  = 0xff000000 | (Integer.parseInt(pcol[1])<<16) | (Integer.parseInt(pcol[2])<<8) | Integer.parseInt(pcol[3]);
			colours[idx][2]  = 0xff000000 | (Integer.parseInt(pcol[5])<<16) | (Integer.parseInt(pcol[6])<<8) | Integer.parseInt(pcol[7]);
			if(positions[idx]<minpos)   minpos = positions[idx];
			if(positions[idx]>maxpos)   maxpos = positions[idx];
			if(positions[idx+1]<minpos) minpos = positions[idx+1];
			if(positions[idx+1]>maxpos) maxpos = positions[idx+1];
			idx++;
		}
		double posfac = 1d/(maxpos-minpos);
		for(int p=0; p<positions.length; p++)
			positions[p] = (positions[p]-minpos) * posfac;
		return new LinearSegmentedPColourtable(specColours[0], specColours[1], specColours[2], positions, colours);
	}

	public static int colourmix(int col1, int col2, double frac) {
		int a = (int) (((col1>>24)&255)*(1d-frac) + ((col2>>24)&255)*frac + 0.5d);
		int r = (int) (((col1>>16)&255)*(1d-frac) + ((col2>>16)&255)*frac + 0.5d);
		int g = (int) (((col1>> 8)&255)*(1d-frac) + ((col2>> 8)&255)*frac + 0.5d);
		int b = (int) (( col1     &255)*(1d-frac) + ( col2     &255)*frac + 0.5d);
		return (a<<24) | (r<<16) | (g<<8) | b;
	}

	public abstract int getColour(double percentage);
	public int getColour(double value, double minval, double maxval) {
		return getColour((value-minval)/(maxval-minval));
	}
}
