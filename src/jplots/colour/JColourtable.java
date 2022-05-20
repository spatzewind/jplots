package jplots.colour;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import processing.core.PApplet;

public abstract class JColourtable {

	public static List<String> ctpaths;
	public static Map<String, JColourtable> pctables;
	static {
		ctpaths = new ArrayList<>();
		pctables = new HashMap<>();
		pctables.put("default", new LinearSegmentedJColourtable(0xff00007f, 0xff7f0000, 0x007f7f7f,
				new double[] { 0d, 0.4d, 0.6d, 1d },
				new int[][] { { 0xff0000ff, 0xff00ffff }, { 0xff00ffff, 0xffffff00 }, { 0xffffff00, 0xffaa0000 } }));
		pctables.put("bw", new ColourSequenceJColourtable(0xff000000, 0xffffffff, 0x01999999,
				new int[] { 0xff000000, 0xffffffff }));
		pctables.put("viridis",
				new PolynomialJColourtable(0x007f7f7f,
						new double[] { 0.201570821189837d, 0.294167300047315d, -0.0743508175540542d, -25.607056228173d,
								145.511124038256d, -372.632167385426d, 499.145447735656d, -332.805541786635d,
								86.9682079885005d },
						new double[] { -0.0743537330341358d, 2.76377597149641d, -23.2555302184104d, 158.859539937159d,
								-521.355773730161d, 137.988091791251d, 5924.44032042842d, -26971.3599494796d,
								65857.2411114973d, -102985.323379299d, 106023.724659725d, -69657.9652005744d,
								26491.8632343267d, -4436.64922455579d },
						new double[] { 0.256353527047776d, 1.30541834921538d, 3.07385522780448d, -74.2809484151974d,
								711.990552038198d, -5089.62769276751d, 24942.9548613986d, -81885.8282080333d,
								182018.355331301d, -275322.350623966d, 279413.005987244d, -182065.297593211d,
								68803.0314713548d, -11456.4694989453d }));
		pctables.put("magma",
				new PolynomialJColourtable(0x007f7f7f,
						new double[] { -0.00585332096864856d, 0.888089304296543d, -8.58149233240007d, 98.701618420612d,
								-440.878039344941d, 1040.81783197537d, -1339.87847076277d, 881.475578541223d,
								-231.566779140571d },
						new double[] { -0.0206775739990087d, 3.93287870107854d, -179.224656381348d, 3921.70142834599d,
								-44482.9549601978d, 293527.124975206d, -1224970.44155448d, 3402615.37657402d,
								-6458450.81009414d, 8425548.82106851d, -7439600.79067079d, 4251016.63174419d,
								-1419224.47177957d, 210276.136646589d },
						new double[] { 0.0144250061175345d, 1.34703662295192d, 0.398656521459664d, -171.369218273499d,
								4158.2584982316d, -38051.8228068549d, 187013.08257451d, -566957.693531997d,
								1128121.15855736d, -1508193.65828506d, 1347417.47641253d, -773851.231933887d,
								258964.349740317d, -38449.6107947707d }));
	}

	public static JColourtable load(String path) {
		String end = path.substring(path.lastIndexOf("."));
//		if(".rgb".equals(end))
//			return readRGB(path);
		if (".cpt".equals(end))
			return readCPT(path);
		System.err.println("Cannot determine type of colourtable. Check the existence of\n"
				+ "filename extension to be one of following: .cpt, .rgb");
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
		if (lines == null)
			return null;
		int col_count = 0;
		for (String line : lines) {
			if ((line.length() < 3) || line.startsWith("#"))
				continue;
			if (line.startsWith("F") || line.startsWith("B") || line.startsWith("N"))
				continue;
			col_count++;
		}
		double[] positions = new double[col_count + 1];
		int[][] colours = new int[col_count][2];
		int[] specColours = { 0xff000000, 0xffffffff, 0x007f7f7f };
		double minpos = Double.POSITIVE_INFINITY, maxpos = Double.NEGATIVE_INFINITY;
		int idx = 0;
		for (String line : lines) {
			if ((line.length() < 3) || line.startsWith("#"))
				continue;
			String[] pcol = line.split("\\s*");
			if (line.startsWith("F")) {
				specColours[1] = 0xff000000 | (Integer.parseInt(pcol[1]) << 16) | (Integer.parseInt(pcol[2]) << 8)
						| Integer.parseInt(pcol[3]);
				continue;
			}
			if (line.startsWith("B")) {
				specColours[0] = 0xff000000 | (Integer.parseInt(pcol[1]) << 16) | (Integer.parseInt(pcol[2]) << 8)
						| Integer.parseInt(pcol[3]);
				continue;
			}
			if (line.startsWith("N")) {
				specColours[2] = 0xff000000 | (Integer.parseInt(pcol[1]) << 16) | (Integer.parseInt(pcol[2]) << 8)
						| Integer.parseInt(pcol[3]);
				continue;
			}
			positions[idx] = Double.parseDouble(pcol[0]);
			positions[idx + 1] = Double.parseDouble(pcol[4]);
			colours[idx][0] = 0xff000000 | (Integer.parseInt(pcol[1]) << 16) | (Integer.parseInt(pcol[2]) << 8)
					| Integer.parseInt(pcol[3]);
			colours[idx][2] = 0xff000000 | (Integer.parseInt(pcol[5]) << 16) | (Integer.parseInt(pcol[6]) << 8)
					| Integer.parseInt(pcol[7]);
			if (positions[idx] < minpos)
				minpos = positions[idx];
			if (positions[idx] > maxpos)
				maxpos = positions[idx];
			if (positions[idx + 1] < minpos)
				minpos = positions[idx + 1];
			if (positions[idx + 1] > maxpos)
				maxpos = positions[idx + 1];
			idx++;
		}
		double posfac = 1d / (maxpos - minpos);
		for (int p = 0; p < positions.length; p++)
			positions[p] = (positions[p] - minpos) * posfac;
		return new LinearSegmentedJColourtable(specColours[0], specColours[1], specColours[2], positions, colours);
	}

	public static int colourmix(int col1, int col2, double frac) {
		int a = (int) (((col1 >> 24) & 255) * (1d - frac) + ((col2 >> 24) & 255) * frac + 0.5d);
		int r = (int) (((col1 >> 16) & 255) * (1d - frac) + ((col2 >> 16) & 255) * frac + 0.5d);
		int g = (int) (((col1 >> 8) & 255) * (1d - frac) + ((col2 >> 8) & 255) * frac + 0.5d);
		int b = (int) ((col1 & 255) * (1d - frac) + (col2 & 255) * frac + 0.5d);
		return (a << 24) | (r << 16) | (g << 8) | b;
	}

	public abstract int getColour(double percentage);

	public int getColour(double value, double minval, double maxval) {
		return getColour((value - minval) / (maxval - minval));
	}
}
