package jplots.shapes;

import java.io.PrintWriter;

import jplots.JPlot;
import processing.core.PConstants;
import processing.core.PFont;
import processing.core.PGraphics;

public class JTextShape extends JPlotShape {

	private int halign, valign;
	private int txtColour;
	private float txtSize, x, y;
	private float rotation;
	private String txt;
	private PFont font;
	
	public JTextShape(String text) {
		this(text, 0f, 0f, 12f, PConstants.LEFT, PConstants.TOP, JPlotShape.fillColour, 0);
	}

	public JTextShape(String text, int colour) {
		this(text, 0f, 0f, 12f, PConstants.LEFT, PConstants.TOP, colour, 0);
	}

	public JTextShape(String text, float x_pos, float y_pos) {
		this(text, x_pos, y_pos, 12f, PConstants.LEFT, PConstants.TOP, JPlotShape.fillColour, 0);
	}

	public JTextShape(String text, float x_pos, float y_pos, int colour) {
		this(text, x_pos, y_pos, 12f, PConstants.LEFT, PConstants.TOP, colour, 0);
	}

	public JTextShape(String text, float x_pos, float y_pos, float size, int colour) {
		this(text, x_pos, y_pos, size, PConstants.LEFT, PConstants.TOP, colour, 0);
	}

	public JTextShape(String text, float x_pos, float y_pos, float size, int horizontal_alignement,
			int vertical_alignement, int colour, float rotate) {
		txt = text;
		x = x_pos;
		y = y_pos;
		halign = horizontal_alignement;
		valign = vertical_alignement;
		txtSize = size;
		txtColour = colour;
		rotation = rotate;
		font = null;
	}
	
	@Override
	public void draw(JPlot plot, PGraphics g) {
		g.pushMatrix();
		g.translate(x, y);
		g.rotate(rotation);
		g.fill(txtColour);
		g.noStroke();
		g.textAlign(halign, valign);
		if(font!=null)
			g.textFont(font);
		g.textSize(txtSize);
		switch (valign) {
			case PConstants.TOP:
				g.text(txt, 0f, 0.16666667f * txtSize);
				break;
			case PConstants.CENTER:
				g.text(txt, 0f, -0.10000000f * txtSize);
				break;
			case PConstants.BOTTOM:
				g.text(txt, 0f, -0.16666667f * txtSize);
				break;
		}
		g.popMatrix();
	}
	
	@Override
	public void printStack(PrintWriter pw, String off) {
		pw.println(off+this.getClass().getSimpleName()+" { text=\""+txt+"\" }");
	}
}
