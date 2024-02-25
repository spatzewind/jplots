package jplots.shapes;

import java.awt.Color;
import java.awt.Dimension;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.lang.Character.UnicodeBlock;

import javax.swing.JLabel;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGeneratorContext;
import org.apache.batik.svggen.SVGGraphics2D;
import org.scilab.forge.jlatexmath.DefaultTeXFont;
import org.scilab.forge.jlatexmath.TeXConstants;
import org.scilab.forge.jlatexmath.TeXFormula;
import org.scilab.forge.jlatexmath.TeXIcon;
import org.scilab.forge.jlatexmath.cyrillic.CyrillicRegistration;
import org.scilab.forge.jlatexmath.greek.GreekRegistration;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;
import org.xml.sax.SAXException;

import jplots.JPlot;
import static jplots.JPlotConstants.*;
import processing.core.PGraphics;
import processing.core.PShape;
import processing.core.PShapeSVG;
import processing.data.XML;

public class JLatexShape extends JPlotShape {

	private static final DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
	private static final String svgNS = "http://www.w3.org/2000/svg";
	private static final Document document = domImpl.createDocument(svgNS, "svg", null);
	private static SVGGeneratorContext ctx = SVGGeneratorContext.createDefault(document);
	private static SVGGraphics2D svgGenerator = new SVGGraphics2D(ctx, true);
	
	static {
		DefaultTeXFont.registerAlphabet(new CyrillicRegistration());
		DefaultTeXFont.registerAlphabet(new GreekRegistration());
//		TeXFormula.registerExternalFont(UnicodeBlock.BASIC_LATIN, "Computer modern");
		TeXFormula.registerExternalFont(UnicodeBlock.BASIC_LATIN, "Computer modern");
	}

	private int halign, valign;
	private int txtColour;
	private float txtSize, x, y;
	private float rotation;
	private String txt, txtstyle;
	
	public JLatexShape(String text) {
		this(text, 0f, 0f, 12f, LEFT, TOP, JPlotShape.fillColour, 0, null);
	}

	public JLatexShape(String text, int colour) {
		this(text, 0f, 0f, 12f, LEFT, TOP, colour, 0, null);
	}

	public JLatexShape(String text, float x_pos, float y_pos) {
		this(text, x_pos, y_pos, 12f, LEFT, TOP, JPlotShape.fillColour, 0, null);
	}

	public JLatexShape(String text, float x_pos, float y_pos, int colour) {
		this(text, x_pos, y_pos, 12f, LEFT, TOP, colour, 0, null);
	}

	public JLatexShape(String text, float x_pos, float y_pos, float size, int colour) {
		this(text, x_pos, y_pos, size, LEFT, TOP, colour, 0, null);
	}

	public JLatexShape(String text, float x_pos, float y_pos, float size, int horizontal_alignement,
			int vertical_alignement, int colour, float rotate, String style) {
		txt = text;
		x = x_pos;
		y = y_pos;
		halign = horizontal_alignement;
		valign = vertical_alignement;
		txtSize = size;
		txtColour = colour;
		rotation = rotate;
		txtstyle = "\\text";
		if(style!=null)
			txtstyle = style;
	}
	
	@Override
	public void draw(JPlot plot, PGraphics g) {
		PShape sp = toPShape(txt, (int) (txtSize+0.5d), txtColour, txtstyle);
//		System.out.println("\""+txt+"\" --> "+sp+"    {w/h: "+sp.width+"/"+sp.height+"}");
		g.pushMatrix();
		g.translate(x, y);
		g.rotate(rotation);
		float xoff = 0f;
		switch(halign) {
			case LEFT:   xoff = -0.0f*sp.width; break;
			case CENTER: xoff = -0.5f*sp.width; break;
			case RIGHT:  xoff = -1.0f*sp.width; break;
		}
		float yoff = 0f;
		switch(valign) {
			case TOP:    yoff =  0.16666667f * txtSize - 0.0f*sp.height; break;
			case CENTER: yoff = -0.10000000f * txtSize - 0.5f*sp.height; break;
			case BOTTOM: yoff = -0.16666667f * txtSize - 1.0f*sp.height; break;
		}
		g.shape(sp, xoff, yoff);
		g.popMatrix();
	}
	
	@Override
	public void printStack(PrintWriter pw, String off) {
		pw.println(off+this.getClass().getSimpleName()+" { text=\""+txt+"\" }");
	}
	
	public static PShape toPShape(String content, int tsize, int tcolor, String txtstyle) {
		
		TeXFormula formula = new TeXFormula(txtstyle+"{"+content+"}");
		
		TeXIcon icon = formula.createTeXIcon(TeXConstants.STYLE_DISPLAY, tsize);
		svgGenerator.setSVGCanvasSize(new Dimension(icon.getIconWidth(), icon.getIconHeight()));
		
		JLabel jl = new JLabel();
		jl.setForeground(new Color(tcolor, true));
		icon.paintIcon(jl, svgGenerator, 0, 0);
		
		boolean useCSS = true;
		Writer out = new StringWriter();
		
		try {
			svgGenerator.stream(out, useCSS);
			return new PShapeSVG(XML.parse(out.toString()));
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (ParserConfigurationException e) {
			throw new RuntimeException(e);
		} catch (SAXException e) {
			throw new RuntimeException(e);
		}
	}
}
