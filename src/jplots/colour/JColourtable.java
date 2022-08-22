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
		pctables.put("default",
				new LinearSegmentedJColourtable(0xff00007f, 0xff7f0000, 0x007f7f7f, new double[] { 0d, 0.4d, 0.6d, 1d },
						new int[][] { { 0xff0000ff, 0xff00ffff }, { 0xff00ffff, 0xffffff00 }, { 0xffffff00, 0xffaa0000 } }));
		pctables.put("bw",
				new ColourSequenceJColourtable(0xff000000, 0xffffffff, 0x01999999, new int[] { 0xff000000, 0xffffffff }));
		pctables.put("bwr",
				new LinearSegmentedJColourtable(0xff0000aa, 0xffaa0000, 0xff999999, new double[] {0d, 0.5d, 1d},
						new int[][] {{0xff0000ff,0xffffffff},{0xffffffff,0xffff0000}}));
		pctables.put("enmediff",
				new LinearSegmentedJColourtable(0xff00cc00, 0xffffff66, 0x007f7f7f, new double[] {0d,0.25d,0.5d,0.75d,1d},
						new int[][] {{0xff00b200,0xff0000ff},{0xff0000ff,0xffffffff},{0xffffffff,0xffff0000},{0xffff0000,0xffffff00}}));
		pctables.put("hot",
				new LinearSegmentedJColourtable( 0xff030000, 0xffffffff, 0x00999999, new double[] {0d,0.36d,0.75d,1d},
						new int[][] {{0xff030000,0xffff0000},{0xffff0000,0xffffff00},{0xffffff00,0xffffffff}}));
		pctables.put("magma",
				new PolynomialJColourtable(0x007f7f7f,
						new double[] { -0.00585332096864856d, 0.888089304296543d, -8.58149233240007d, 98.701618420612d, -440.878039344941d,
								1040.81783197537d, -1339.87847076277d, 881.475578541223d, -231.566779140571d },
						new double[] { -0.0206775739990087d, 3.93287870107854d, -179.224656381348d, 3921.70142834599d, -44482.9549601978d,
								293527.124975206d, -1224970.44155448d, 3402615.37657402d, -6458450.81009414d, 8425548.82106851d,
								-7439600.79067079d, 4251016.63174419d, -1419224.47177957d, 210276.136646589d },
						new double[] { 0.0144250061175345d, 1.34703662295192d, 0.398656521459664d, -171.369218273499d, 4158.2584982316d,
								-38051.8228068549d, 187013.08257451d, -566957.693531997d, 1128121.15855736d, -1508193.65828506d,
								1347417.47641253d, -773851.231933887d, 258964.349740317d, -38449.6107947707d }));
		pctables.put("ocean",
				new LinearSegmentedJColourtable( 0xff007f00, 0xffffffff, 0x007f7f7f, new double[] {0d,1d/3d,2d/3d,1d},
						new int[][] {{0xff007f00,0xff000055},{0xff000055,0xff007faa},{0xff007faa,0xffffffff} }));
		pctables.put("terrain",
				new LinearSegmentedJColourtable(0xff333399, 0xffffffff, 0x007f7f7f, new double[] {0d,0.15d,0.25d,0.5d,0.75d,1d},
						new int[][] {{0xff333399,0xff0099ff},{0xff0099ff,0xff00cc66},{0xff00cc66,0xffffff99},{0xffffff99,0xff7f5d55},{0xff7f5d55,0xffffffff}}));
		pctables.put("turner",
				new LinearSegmentedJColourtable(0xff00cc00, 0xffffff66, 0xff999999, new double[] {0d,0.25d,0.3125d,0.375d,0.5d,0.625d,0.6875d,0.75d,1d},
						new int[][] {{0xff7f7f7f,0xff002a33},{0xff002a33,0xff00bfff},{0xff00bfff,0xff000000},{0xff000000,0xff00ff00},
									 {0xff00ff00,0xffffffff},{0xffffffff,0xffff0000},{0xffff0000,0xff330000},{0xff330000,0xff7f7f7f}}));
		pctables.put("twilight",
				new MultiPolynomialJColourtable(0x007f7f7f, 0xffe3d9e3, 0xffe3d9e3, new double[] {0d, 0.5d, 1d},
						new double[][] {{0.8857331275784759d, -0.6750071450587711d, -27.80790241691284d, 125.93029705248773d, -83.8675944507122d, -355.8921470269561d, 461.56191185489297d},
										{11.440466403961182d, -34.27635192871094d, -101.29255676269531d, 577.7247314453125d, -951.5903778076172d, 682.9105911254883d, -184.0284824371338d}},
						new double[][] {{0.8485095066162245d, 0.605482972565369d, -28.997725043504033d, 184.94879733538255d, -630.542937216349d, 1018.3341409992427d, -598.6334923543036d},
										{-40.593764543533325d, 354.4781436920166d, -1258.2324981689453d, 2332.84578704834d, -2394.460029602051d, 1299.4271392822266d, -292.61353397369385d}},
						new double[][] {{0.8926008338073075d, -0.934815051696205d, -1.5657150344923139d, 11.383007136173546d, 76.91811736300588d, -436.39043478667736d, 485.70714078843594d},
										{101.16589748859406d, -803.6041326522827d, 2577.301239013672d, -4244.141311645508d, 3761.6750411987305d, -1684.4910850524902d, 292.9793186187744d}}
						));
		pctables.put("twilight_shifted",
				new MultiPolynomialJColourtable(0x007f7f7f, 0xff2f1337, 0xff2f1336, new double[] {0d, 0.5d, 1d},
						new double[][] {{0.18766755521437517d, 0.5608648217021255d, 24.54853762150742d, -232.90330296382308d, 794.9527024403214d, -1092.4720372110605d, 502.31343837827444d},
										{-47.497127175331116d, 397.86163997650146d, -1324.698444366455d, 2281.51904296875d, -2136.766456604004d, 1025.3919792175293d, -195.6248025894165d}},
						new double[][] {{0.07931472603343082d, -0.9689002125669504d, 13.493183587444946d, 31.79733134806156d, -344.68831496685743d, 800.5009777843952d, -612.2958381175995d},
										{-54.75673669576645d, 445.5318088531494d, -1472.903712272644d, 2590.342580795288d, -2568.006576538086d, 1357.7884616851807d, -297.9217383861542d}},
						new double[][] {{0.21766537517271445d, 0.732190937967971d, 42.00476239598356d, -297.3306187223643d, 836.3242402374744d, -1065.1492710709572d, 511.2624908387661d},
										{-42.92372006177902d, 250.30328798294067d, -424.81978034973145d, -48.1270751953125d, 859.4668121337891d, -867.9975528717041d, 274.3111753463745d}}
						));
		pctables.put("viridis",
				new PolynomialJColourtable(0x007f7f7f,
						new double[] { 0.201570821189837d, 0.294167300047315d, -0.0743508175540542d, -25.607056228173d, 145.511124038256d,
								-372.632167385426d, 499.145447735656d, -332.805541786635d, 86.9682079885005d },
						new double[] { -0.0743537330341358d, 2.76377597149641d, -23.2555302184104d, 158.859539937159d, -521.355773730161d,
								137.988091791251d, 5924.44032042842d, -26971.3599494796d, 65857.2411114973d, -102985.323379299d,
								106023.724659725d, -69657.9652005744d, 26491.8632343267d, -4436.64922455579d },
						new double[] { 0.256353527047776d, 1.30541834921538d, 3.07385522780448d, -74.2809484151974d, 711.990552038198d,
								-5089.62769276751d, 24942.9548613986d, -81885.8282080333d, 182018.355331301d, -275322.350623966d,
								279413.005987244d, -182065.297593211d, 68803.0314713548d, -11456.4694989453d }));
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
	
	public static JColourtable getByName(String name) {
		if(name.endsWith("_r")) {
			JColourtable ct = pctables.get(name.substring(0, name.length()-2));
			if(ct!=null)
				return new ReverseJColourtable(ct);
		}
		return pctables.get(name);
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
	
	public static double oclamp(double vmin, double vmax, double v) {
		if(v<vmin) return vmin;
		if(v>vmax) return vmax;
		return v;
	}
	public static double opow(double v, int p) {
		if (p == 0)
			return 1d;
		double b = p < 0 ? 1d / v : v, res = 1d;
		int pc = p < 0 ? -p : p;
		while (pc > 0) {
			if (pc % 2 == 1)
				res *= b;
			b *= b;
			pc = (pc >> 1);
		}
		return res;
	}
	public static double polynom(double v, double[] c) {
		double p = 0d;
		for (int i=0; i<c.length; i++)
			p += opow(v,i) * c[i];
		return p;
	}
	public static int toColour(double r, double g, double b) {
		int ri = (int) (255.9999d * oclamp(0d,1d,r));
		int gi = (int) (255.9999d * oclamp(0d,1d,g));
		int bi = (int) (255.9999d * oclamp(0d,1d,b));
		return 0xff000000 | (ri<<16) | (gi<<8) | bi;
	}
}
