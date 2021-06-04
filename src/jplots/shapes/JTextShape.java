package jplots.shapes;

import jplots.JPlot;
import processing.core.PApplet;
import processing.core.PGraphics;

public class JTextShape extends JPlotShape {

	private int halign, valign;
	private int txtColour;
	private float txtSize, x, y;
	private float rotation;
	private String txt;

	public JTextShape(String text) {
		this(text,0f,0f,12f,PApplet.LEFT,PApplet.TOP,JPlotShape.fillColour,0); }
	public JTextShape(String text, int colour) {
		this(text,0f,0f,12f,PApplet.LEFT,PApplet.TOP,colour,0); }
	public JTextShape(String text, float x_pos, float y_pos) {
		this(text,x_pos,y_pos,12f,PApplet.LEFT,PApplet.TOP,JPlotShape.fillColour,0); }
	public JTextShape(String text, float x_pos, float y_pos, int colour) {
		this(text,x_pos,y_pos,12f,PApplet.LEFT,PApplet.TOP,colour,0); }
	public JTextShape(String text, float x_pos, float y_pos, float size, int colour) {
		this(text,x_pos,y_pos,size,PApplet.LEFT,PApplet.TOP,colour,0); }
	public JTextShape(String text, float x_pos, float y_pos, float size, int horizontal_alignement, int vertical_alignement, int colour, float rotate) {
		txt = text;
		x = x_pos;
		y = y_pos;
		halign = horizontal_alignement;
		valign = vertical_alignement;
		txtSize = size;
		txtColour = colour;
		rotation = rotate;
	}

	@Override
	public void draw(JPlot plot, PGraphics g) {
		g.fill(txtColour); g.noStroke();
		g.textAlign(halign,valign); g.textSize(txtSize);
		if(rotation==0f) {
			switch(valign) {
				case PApplet.TOP:    g.text(txt, x, y+0.16666667f*txtSize); break;
				case PApplet.CENTER: g.text(txt, x, y-0.10000000f*txtSize); break;
				case PApplet.BOTTOM: g.text(txt, x, y-0.16666667f*txtSize); break;
				default: g.text(txt, x, y); break;
			}
		} else {
			g.pushMatrix();
			g.translate(x, y);
			g.rotate(rotation);
			switch(valign) {
				case PApplet.TOP:    g.text(txt, 0f,  0.16666667f*txtSize); break;
				case PApplet.CENTER: g.text(txt, 0f, -0.10000000f*txtSize); break;
				case PApplet.BOTTOM: g.text(txt, 0f, -0.16666667f*txtSize); break;
			}
			g.popMatrix();
		}
	}
}
