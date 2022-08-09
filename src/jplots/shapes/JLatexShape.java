package jplots.shapes;

import java.io.PrintWriter;

import jplots.JPlot;
import latex.PTeX;
import processing.core.PConstants;
import processing.core.PFont;
import processing.core.PGraphics;
import processing.core.PShape;

public class JLatexShape extends JPlotShape {

	private int halign, valign;
	private int txtColour;
	private float txtSize, x, y;
	private float rotation;
	private String txt;
	private PFont font;
	
	public JLatexShape(String text) {
		this(text, 0f, 0f, 12f, PConstants.LEFT, PConstants.TOP, JPlotShape.fillColour, 0);
	}

	public JLatexShape(String text, int colour) {
		this(text, 0f, 0f, 12f, PConstants.LEFT, PConstants.TOP, colour, 0);
	}

	public JLatexShape(String text, float x_pos, float y_pos) {
		this(text, x_pos, y_pos, 12f, PConstants.LEFT, PConstants.TOP, JPlotShape.fillColour, 0);
	}

	public JLatexShape(String text, float x_pos, float y_pos, int colour) {
		this(text, x_pos, y_pos, 12f, PConstants.LEFT, PConstants.TOP, colour, 0);
	}

	public JLatexShape(String text, float x_pos, float y_pos, float size, int colour) {
		this(text, x_pos, y_pos, size, PConstants.LEFT, PConstants.TOP, colour, 0);
	}

	public JLatexShape(String text, float x_pos, float y_pos, float size, int horizontal_alignement,
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
		PShape sp = PTeX.toPShape(txt, (int) (txtSize+0.5d), 0x00999999, txtColour);
		System.out.println("\""+txt+"\" --> "+sp+"    {w/h: "+sp.width+"/"+sp.height+"}");
		g.pushMatrix();
		g.translate(x-0.5f*sp.width, y-0.5f*sp.height);
		g.rotate(rotation);
		float yoff = 0f;
		switch(valign) {
			case PConstants.TOP:    yoff =  0.16666667f * txtSize; break;
			case PConstants.CENTER: yoff = -0.10000000f * txtSize; break;
			case PConstants.BOTTOM: yoff = -0.16666667f * txtSize; break;
		}
		g.shape(sp, 0f, yoff);
		g.popMatrix();
	}
	
	@Override
	public void printStack(PrintWriter pw, String off) {
		pw.println(off+this.getClass().getSimpleName()+" { text=\""+txt+"\" }");
	}
}
