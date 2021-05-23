package pplots.gfx;

import processing.core.PApplet;
import processing.core.PGraphics;

public class PTextShape extends PPlotShape {

	private int halign, valign;
	private int txtColour;
	private float txtSize, x, y;
	private String txt;

	public PTextShape(String text) {
		this(text,0f,0f,12f,PApplet.LEFT,PApplet.TOP,PPlotShape.fillColour); }
	public PTextShape(String text, int colour) {
		this(text,0f,0f,12f,PApplet.LEFT,PApplet.TOP,colour); }
	public PTextShape(String text, float x_pos, float y_pos) {
		this(text,x_pos,y_pos,12f,PApplet.LEFT,PApplet.TOP,PPlotShape.fillColour); }
	public PTextShape(String text, float x_pos, float y_pos, int colour) {
		this(text,x_pos,y_pos,12f,PApplet.LEFT,PApplet.TOP,colour); }
	public PTextShape(String text, float x_pos, float y_pos, float size, int colour) {
		this(text,x_pos,y_pos,size,PApplet.LEFT,PApplet.TOP,colour); }
	public PTextShape(String text, float x_pos, float y_pos, float size, int horizontal_alignement, int vertical_alignement, int colour) {
		txt = text;
		x = x_pos;
		y = y_pos;
		halign = horizontal_alignement;
		valign = vertical_alignement;
		txtSize = size;
		txtColour = colour;
	}

	@Override
	public void draw(PGraphics g) {
		g.fill(txtColour); g.noStroke();
		g.textAlign(halign,valign); g.textSize(txtSize);
		g.text(txt, x, y);
	}
}
