package jplots.color;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import processing.core.PApplet;

public abstract class JColourtable {

	public static List<String> ctpaths;
	public static Map<String,JColourtable> pctables;
	static {
		ctpaths = new ArrayList<String>();
		pctables = new HashMap<String,JColourtable>();
		pctables.put("default", new LinearSegmentedJColourtable(0xff00007f, 0xff7f0000, 0x007f7f7f,
				new double[] {0d, 0.4d,0.6d,1d},
				new int[][] {{0xff0000ff,0xff00ffff},{0xff00ffff,0xffffff00},{0xffffff00,0xffaa0000}}));
		pctables.put("viridis", new PolynomialJColourtable(
				new double[] {1.00140166562551d, -3.18751400748827d, 3.22471285473439d,
						-35.0504042872711d, 209.112595595727d, -503.543781409822d,
						604.616471526084d, -362.940118305759d, 86.9682071894914d},
				new double[] {0.901263970278695d, -0.528138050670421d, 4.24357591744261d,
						-39.7684145926684d, 141.424665245748d, -279.228342175845d,
						321.682111320121d, -201.373277956237d, 52.6386411613611d},
				new double[] {0.124237625850364d, -1.93834142286999d, 28.8167568387114d,
						-129.179049793897d, 338.195602476375d, -572.509341852231d,
						606.370604063703d, -359.257649372089d, 89.6360378222204d}));
		pctables.put("magma", new PolynomialJColourtable(
				new double[] {0.972643073613039d, 1.05603463643587d, -20.3739565374434d,
						168.267698278739d, -684.7241212886d, 1436.61797121108d,
						-1630.11620086329d, 955.798163369022d, -227.502465906169d},
				new double[] {1.01092264386188d, -3.44600529318444d, 63.2483466404809d,
						-1486.75734318594d, 18375.8378753408d, -134425.885768738d,
						626778.267707128d, -1946246.54540799d, 4116505.34600292d,
						-5953554.33999223d, 5791782.45656149d, -3621885.23038504d,
						1314373.3805903d, -210277.363781746d},
				new double[] {0.699329889616736d, -1.41059876981806d, -30.8509315000795d,
						584.523858664168d, -5699.6946144749d, 34656.0458114193d,
						-140663.485080509d, 397612.489785934d, -794737.594802251d,
						1115363.42713131d, -1070085.93532432d, 665468.225150372d,
						-240922.544885487d, 38456.1195946251d}));
	}
	
	public static JColourtable load(String path) {
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
	private static JColourtable readCPT(String path) {
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
		return new LinearSegmentedJColourtable(specColours[0], specColours[1], specColours[2], positions, colours);
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
