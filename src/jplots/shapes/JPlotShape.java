package jplots.shapes;

import java.io.PrintWriter;

import jplots.JPlot;
import jplots.JPlotConstants;
import processing.core.PGraphics;

public abstract class JPlotShape {

	public static boolean useFill = true;
	public static boolean useStroke = true;
	public static int fillColour = 0xffffffff;
	public static int strokeColour = 0xff000000;
	public static float strokeWeight = 1f;
	
	public final static float ROTATE_CLOCKWISE        = JPlotConstants.ROTATE_CLOCKWISE;
	public final static float ROTATE_COUNTERCLOCKWISE = JPlotConstants.ROTATE_COUNTERCLOCKWISE;

	public static void fill(int col) {
		fillColour = col;
		useFill = true;
	}

	public static void stroke(int col) {
		strokeColour = col;
		useStroke = true;
	}

	public static void strokeWeight(float w) {
		strokeWeight = w;
		useStroke = true;
	}

	public static void noFill() {
		useFill = false;
	}

	public static void noStroke() {
		useStroke = false;
	}

	public abstract void draw(JPlot plot, PGraphics g);
	
	public void printStack(PrintWriter pw, String off) {
		pw.println(off+this.getClass().getSimpleName());
	}
}
