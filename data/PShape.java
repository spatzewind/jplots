package processing.core;

import java.awt.Image;
import java.awt.image.BufferedImage;
import java.util.HashMap;
import java.util.Map;
import javax.swing.ImageIcon;
import javax.xml.bind.DatatypeConverter;

public class PShape implements PConstants {
	protected String name;
	protected Map<String, PShape> nameTable;
	public static final int PRIMITIVE = 101;
	public static final int PATH = 102;
	public static final int GEOMETRY = 103;
	protected int family;
	protected int kind;
	protected PMatrix matrix;
	protected int textureMode;
	protected PImage image;
	protected String imagePath = null;

	public static final String OUTSIDE_BEGIN_END_ERROR = "%1$s can only be called between beginShape() and endShape()";
	public static final String INSIDE_BEGIN_END_ERROR = "%1$s can only be called outside beginShape() and endShape()";
	public static final String NO_SUCH_VERTEX_ERROR = "%1$s vertex index does not exist";
	public static final String NO_VERTICES_ERROR = "getVertexCount() only works with PATH or GEOMETRY shapes";
	public static final String NOT_A_SIMPLE_VERTEX = "%1$s can not be called on quadratic or bezier vertices";
	public static final String PER_VERTEX_UNSUPPORTED = "This renderer does not support %1$s for individual vertices";

	public float width;
	public float height;
	public float depth;
	PGraphics g;
	protected boolean visible = true;
	protected boolean openShape = false;
	protected boolean openContour = false;
	protected boolean stroke;
	protected int strokeColor;
	protected float strokeWeight;
	protected int strokeCap;
	protected int strokeJoin;
	protected boolean fill;
	protected int fillColor;
	protected boolean tint;
	protected int tintColor;
	protected int ambientColor;
	protected boolean setAmbient;
	protected int specularColor;
	protected int emissiveColor;
	protected float shininess;
	protected int sphereDetailU;
	protected int sphereDetailV;
	protected int rectMode;
	protected int ellipseMode;
	protected boolean style = true;
	protected float[] params;
	protected int vertexCount;
	protected float[][] vertices;
	protected PShape parent;
	protected int childCount;
	protected PShape[] children;
	protected int vertexCodeCount;
	protected int[] vertexCodes;
	protected boolean close;
	protected float calcR;
	protected float calcG;
	protected float calcB;
	protected float calcA;
	protected int calcRi;
	protected int calcGi;
	protected int calcBi;
	protected int calcAi;
	protected int calcColor;
	protected boolean calcAlpha;
	public int colorMode;
	public float colorModeX;
	public float colorModeY;
	public float colorModeZ;
	public float colorModeA;
	boolean colorModeScale;
	boolean colorModeDefault;
	protected boolean is3D = false;
	protected boolean perVertexStyles = false;


	public PShape() { this.family = 0; }
	public PShape(int family) { this.family = family; }
	public PShape(PGraphics g, int family) {
		this.g = g;
		this.family = family;
		this.textureMode = g.textureMode;
		colorMode(g.colorMode, g.colorModeX, g.colorModeY, g.colorModeZ, g.colorModeA);
		this.fill = g.fill;
		this.fillColor = g.fillColor;
		this.stroke = g.stroke;
		this.strokeColor = g.strokeColor;
		this.strokeWeight = g.strokeWeight;
		this.strokeCap = g.strokeCap;
		this.strokeJoin = g.strokeJoin;
		this.tint = g.tint;
		this.tintColor = g.tintColor;
    
    this.setAmbient = g.setAmbient;
    this.ambientColor = g.ambientColor;
    this.specularColor = g.specularColor;
    this.emissiveColor = g.emissiveColor;
    this.shininess = g.shininess;
    
    this.sphereDetailU = g.sphereDetailU;
    this.sphereDetailV = g.sphereDetailV;




    
    this.rectMode = g.rectMode;
    this.ellipseMode = g.ellipseMode;
  }

  public PShape(PGraphics g, int kind, float... params) {
    this(g, 101);
    setKind(kind);
    setParams(params);
  }

  public void setFamily(int family) { this.family = family; }

  public void setKind(int kind) { this.kind = kind; }

  public void setName(String name) { this.name = name; }

  public String getName() { return this.name; }

  public boolean isVisible() { return this.visible; }

  public void setVisible(boolean visible) { this.visible = visible; }



















  
  public void disableStyle() {
    this.style = false;
    
    for (int i = 0; i < this.childCount; i++) {
      this.children[i].disableStyle();
    }
  }















  
  public void enableStyle() {
    this.style = true;
    
    for (int i = 0; i < this.childCount; i++) {
      this.children[i].enableStyle();
    }
  }


















  
  public float getWidth() { return this.width; }







  
  public float getHeight() { return this.height; }








  
  public float getDepth() { return this.depth; }




































  
  public boolean is2D() { return !this.is3D; }






  
  public boolean is3D() { return this.is3D; }



  
  public void set3D(boolean val) { this.is3D = val; }

















  
  public void textureMode(int mode) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "textureMode()" });
      
      return;
    } 
    this.textureMode = mode;
  }
  
  public void texture(PImage tex) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "texture()" });
      
      return;
    } 
    this.image = tex;
  }
  
  public void noTexture() {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "noTexture()" });
      
      return;
    } 
    this.image = null;
  }





  
  protected void solid(boolean solid) {}




  
  public void beginContour() {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "beginContour()" });
      
      return;
    } 
    if (this.family == 0) {
      PGraphics.showWarning("Cannot begin contour in GROUP shapes");
      
      return;
    } 
    if (this.openContour) {
      PGraphics.showWarning("Already called beginContour().");
      return;
    } 
    this.openContour = true;
    beginContourImpl();
  }

  
  protected void beginContourImpl() {
    if (this.vertexCodes == null) {
      this.vertexCodes = new int[10];
    } else if (this.vertexCodes.length == this.vertexCodeCount) {
      this.vertexCodes = PApplet.expand(this.vertexCodes);
    } 
    this.vertexCodes[this.vertexCodeCount++] = 4;
  }






  
  public void endContour() {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "endContour()" });
      
      return;
    } 
    if (this.family == 0) {
      PGraphics.showWarning("Cannot end contour in GROUP shapes");
      
      return;
    } 
    if (!this.openContour) {
      PGraphics.showWarning("Need to call beginContour() first.");
      return;
    } 
    endContourImpl();
    this.openContour = false;
  }


  
  protected void endContourImpl() {}

  
  public void vertex(float x, float y) {
    if (this.vertices == null) {
      this.vertices = new float[10][2];
    } else if (this.vertices.length == this.vertexCount) {
      this.vertices = (float[][])PApplet.expand(this.vertices);
    } 
    new float[2][0] = x; new float[2][1] = y; this.vertices[this.vertexCount++] = new float[2];
    
    if (this.vertexCodes == null) {
      this.vertexCodes = new int[10];
    } else if (this.vertexCodes.length == this.vertexCodeCount) {
      this.vertexCodes = PApplet.expand(this.vertexCodes);
    } 
    this.vertexCodes[this.vertexCodeCount++] = 0;
    
    if (x > this.width) {
      this.width = x;
    }
    if (y > this.height) {
      this.height = y;
    }
  }


  
  public void vertex(float x, float y, float u, float v) {}


  
  public void vertex(float x, float y, float z) { vertex(x, y); }



  
  public void vertex(float x, float y, float z, float u, float v) {}



  
  public void normal(float nx, float ny, float nz) {}



  
  public void attribPosition(String name, float x, float y, float z) {}



  
  public void attribNormal(String name, float nx, float ny, float nz) {}



  
  public void attribColor(String name, int color) {}


  
  public void attrib(String name, float... values) {}


  
  public void attrib(String name, int... values) {}


  
  public void attrib(String name, boolean... values) {}


  
  public void beginShape() { beginShape(20); }


  
  public void beginShape(int kind) {
    this.kind = kind;
    this.openShape = true;
  }






  
  public void endShape() { endShape(1); }


  
  public void endShape(int mode) {
    if (this.family == 0) {
      PGraphics.showWarning("Cannot end GROUP shape");
      
      return;
    } 
    if (!this.openShape) {
      PGraphics.showWarning("Need to call beginShape() first");
      
      return;
    } 
    this.close = (mode == 2);

    
    this.openShape = false;
  }






  
  public void strokeWeight(float weight) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "strokeWeight()" });
      
      return;
    } 
    this.strokeWeight = weight;
  }
  
  public void strokeJoin(int join) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "strokeJoin()" });
      
      return;
    } 
    this.strokeJoin = join;
  }
  
  public void strokeCap(int cap) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "strokeCap()" });
      
      return;
    } 
    this.strokeCap = cap;
  }






  
  public void noFill() {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "noFill()" });
      
      return;
    } 
    this.fill = false;
    this.fillColor = 0;
    
    if (!this.setAmbient) {
      this.ambientColor = this.fillColor;
    }
  }

  
  public void fill(int rgb) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "fill()" });
      
      return;
    } 
    this.fill = true;
    colorCalc(rgb);
    this.fillColor = this.calcColor;
    
    if (!this.setAmbient) {
      this.ambientColor = this.fillColor;
    }
  }

  
  public void fill(int rgb, float alpha) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "fill()" });
      
      return;
    } 
    this.fill = true;
    colorCalc(rgb, alpha);
    this.fillColor = this.calcColor;
    
    if (!this.setAmbient) {
      this.ambientColor = this.fillColor;
    }
  }

  
  public void fill(float gray) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "fill()" });
      
      return;
    } 
    this.fill = true;
    colorCalc(gray);
    this.fillColor = this.calcColor;
    
    if (!this.setAmbient) {
      this.ambientColor = this.fillColor;
    }
  }

  
  public void fill(float gray, float alpha) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "fill()" });
      
      return;
    } 
    this.fill = true;
    colorCalc(gray, alpha);
    this.fillColor = this.calcColor;
    
    if (!this.setAmbient) {
      ambient(this.fillColor);
      this.setAmbient = false;
    } 
    
    if (!this.setAmbient) {
      this.ambientColor = this.fillColor;
    }
  }

  
  public void fill(float x, float y, float z) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "fill()" });
      
      return;
    } 
    this.fill = true;
    colorCalc(x, y, z);
    this.fillColor = this.calcColor;
    
    if (!this.setAmbient) {
      this.ambientColor = this.fillColor;
    }
  }

  
  public void fill(float x, float y, float z, float a) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "fill()" });
      
      return;
    } 
    this.fill = true;
    colorCalc(x, y, z, a);
    this.fillColor = this.calcColor;
    
    if (!this.setAmbient) {
      this.ambientColor = this.fillColor;
    }
  }






  
  public void noStroke() {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "noStroke()" });
      
      return;
    } 
    this.stroke = false;
  }

  
  public void stroke(int rgb) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "stroke()" });
      
      return;
    } 
    this.stroke = true;
    colorCalc(rgb);
    this.strokeColor = this.calcColor;
  }

  
  public void stroke(int rgb, float alpha) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "stroke()" });
      
      return;
    } 
    this.stroke = true;
    colorCalc(rgb, alpha);
    this.strokeColor = this.calcColor;
  }

  
  public void stroke(float gray) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "stroke()" });
      
      return;
    } 
    this.stroke = true;
    colorCalc(gray);
    this.strokeColor = this.calcColor;
  }

  
  public void stroke(float gray, float alpha) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "stroke()" });
      
      return;
    } 
    this.stroke = true;
    colorCalc(gray, alpha);
    this.strokeColor = this.calcColor;
  }

  
  public void stroke(float x, float y, float z) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "stroke()" });
      
      return;
    } 
    this.stroke = true;
    colorCalc(x, y, z);
    this.strokeColor = this.calcColor;
  }

  
  public void stroke(float x, float y, float z, float alpha) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "stroke()" });
      
      return;
    } 
    this.stroke = true;
    colorCalc(x, y, z, alpha);
    this.strokeColor = this.calcColor;
  }






  
  public void noTint() {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "noTint()" });
      
      return;
    } 
    this.tint = false;
  }

  
  public void tint(int rgb) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "tint()" });
      
      return;
    } 
    this.tint = true;
    colorCalc(rgb);
    this.tintColor = this.calcColor;
  }

  
  public void tint(int rgb, float alpha) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "tint()" });
      
      return;
    } 
    this.tint = true;
    colorCalc(rgb, alpha);
    this.tintColor = this.calcColor;
  }

  
  public void tint(float gray) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "tint()" });
      
      return;
    } 
    this.tint = true;
    colorCalc(gray);
    this.tintColor = this.calcColor;
  }

  
  public void tint(float gray, float alpha) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "tint()" });
      
      return;
    } 
    this.tint = true;
    colorCalc(gray, alpha);
    this.tintColor = this.calcColor;
  }

  
  public void tint(float x, float y, float z) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "tint()" });
      
      return;
    } 
    this.tint = true;
    colorCalc(x, y, z);
    this.tintColor = this.calcColor;
  }

  
  public void tint(float x, float y, float z, float alpha) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "tint()" });
      
      return;
    } 
    this.tint = true;
    colorCalc(x, y, z, alpha);
    this.tintColor = this.calcColor;
  }





  
  public void ambient(int rgb) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "ambient()" });
      
      return;
    } 
    this.setAmbient = true;
    colorCalc(rgb);
    this.ambientColor = this.calcColor;
  }

  
  public void ambient(float gray) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "ambient()" });
      
      return;
    } 
    this.setAmbient = true;
    colorCalc(gray);
    this.ambientColor = this.calcColor;
  }

  
  public void ambient(float x, float y, float z) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "ambient()" });
      
      return;
    } 
    this.setAmbient = true;
    colorCalc(x, y, z);
    this.ambientColor = this.calcColor;
  }





  
  public void specular(int rgb) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "specular()" });
      
      return;
    } 
    colorCalc(rgb);
    this.specularColor = this.calcColor;
  }

  
  public void specular(float gray) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "specular()" });
      
      return;
    } 
    colorCalc(gray);
    this.specularColor = this.calcColor;
  }

  
  public void specular(float x, float y, float z) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "specular()" });
      
      return;
    } 
    colorCalc(x, y, z);
    this.specularColor = this.calcColor;
  }





  
  public void emissive(int rgb) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "emissive()" });
      
      return;
    } 
    colorCalc(rgb);
    this.emissiveColor = this.calcColor;
  }

  
  public void emissive(float gray) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "emissive()" });
      
      return;
    } 
    colorCalc(gray);
    this.emissiveColor = this.calcColor;
  }

  
  public void emissive(float x, float y, float z) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "emissive()" });
      
      return;
    } 
    colorCalc(x, y, z);
    this.emissiveColor = this.calcColor;
  }





  
  public void shininess(float shine) {
    if (!this.openShape) {
      PGraphics.showWarning("%1$s can only be called between beginShape() and endShape()", new Object[] { "shininess()" });
      
      return;
    } 
    this.shininess = shine;
  }






  
  public void bezierDetail(int detail) {}





  
  public void bezierVertex(float x2, float y2, float x3, float y3, float x4, float y4) {
    if (this.vertices == null) {
      this.vertices = new float[10][];
    } else if (this.vertexCount + 2 >= this.vertices.length) {
      this.vertices = (float[][])PApplet.expand(this.vertices);
    } 
    new float[2][0] = x2; new float[2][1] = y2; this.vertices[this.vertexCount++] = new float[2];
    new float[2][0] = x3; new float[2][1] = y3; this.vertices[this.vertexCount++] = new float[2];
    new float[2][0] = x4; new float[2][1] = y4; this.vertices[this.vertexCount++] = new float[2];

    
    if (this.vertexCodes.length == this.vertexCodeCount) {
      this.vertexCodes = PApplet.expand(this.vertexCodes);
    }
    this.vertexCodes[this.vertexCodeCount++] = 1;
    
    if (x4 > this.width) {
      this.width = x4;
    }
    if (y4 > this.height) {
      this.height = y4;
    }
  }



  
  public void bezierVertex(float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4) {}



  
  public void quadraticVertex(float cx, float cy, float x3, float y3) {
    if (this.vertices == null) {
      this.vertices = new float[10][];
    } else if (this.vertexCount + 1 >= this.vertices.length) {
      this.vertices = (float[][])PApplet.expand(this.vertices);
    } 
    new float[2][0] = cx; new float[2][1] = cy; this.vertices[this.vertexCount++] = new float[2];
    new float[2][0] = x3; new float[2][1] = y3; this.vertices[this.vertexCount++] = new float[2];

    
    if (this.vertexCodes.length == this.vertexCodeCount) {
      this.vertexCodes = PApplet.expand(this.vertexCodes);
    }
    this.vertexCodes[this.vertexCodeCount++] = 2;
    
    if (x3 > this.width) {
      this.width = x3;
    }
    if (y3 > this.height) {
      this.height = y3;
    }
  }






  
  public void quadraticVertex(float cx, float cy, float cz, float x3, float y3, float z3) {}






  
  public void curveDetail(int detail) {}






  
  public void curveTightness(float tightness) {}






  
  public void curveVertex(float x, float y) {}






  
  public void curveVertex(float x, float y, float z) {}





  
  protected void pre(PGraphics g) {
    if (this.matrix != null) {
      g.pushMatrix();
      g.applyMatrix(this.matrix);
    } 














    
    if (this.style) {
      g.pushStyle();
      styles(g);
    } 
  }




  
  protected void styles(PGraphics g) {
    if (this.stroke) {
      g.stroke(this.strokeColor);
      g.strokeWeight(this.strokeWeight);
      g.strokeCap(this.strokeCap);
      g.strokeJoin(this.strokeJoin);
    } else {
      g.noStroke();
    } 
    
    if (this.fill) {
      
      g.fill(this.fillColor);
    } else {
      g.noFill();
    } 
  }



















  
  protected void post(PGraphics g) {
    if (this.matrix != null) {
      g.popMatrix();
    }
    
    if (this.style) {
      g.popStyle();
    }
  }







  
  protected static PShape createShape(PApplet parent, PShape src) {
    PShape dest = null;
    if (src.family == 0) {
      dest = parent.createShape(0);
      copyGroup(parent, src, dest);
    } else if (src.family == 101) {
      dest = parent.createShape(src.kind, src.params);
      copyPrimitive(src, dest);
    } else if (src.family == 103) {
      dest = parent.createShape(src.kind);
      copyGeometry(src, dest);
    } else if (src.family == 102) {
      dest = parent.createShape(102);
      copyPath(src, dest);
    } 
    dest.setName(src.name);
    return dest;
  }


  
  protected static void copyGroup(PApplet parent, PShape src, PShape dest) {
    copyMatrix(src, dest);
    copyStyles(src, dest);
    copyImage(src, dest);
    for (int i = 0; i < src.childCount; i++) {
      PShape c = createShape(parent, src.children[i]);
      dest.addChild(c);
    } 
  }


  
  protected static void copyPrimitive(PShape src, PShape dest) {
    copyMatrix(src, dest);
    copyStyles(src, dest);
    copyImage(src, dest);
  }


  
  protected static void copyGeometry(PShape src, PShape dest) {
    dest.beginShape(src.getKind());
    
    copyMatrix(src, dest);
    copyStyles(src, dest);
    copyImage(src, dest);
    
    if (src.style) {
      for (int i = 0; i < src.vertexCount; i++) {
        float[] vert = src.vertices[i];
        
        dest.fill((int)(vert[6] * 255.0F) << 24 | 
            (int)(vert[3] * 255.0F) << 16 | 
            (int)(vert[4] * 255.0F) << 8 | 
            (int)(vert[5] * 255.0F));






        
        if (0.0F < PApplet.dist(vert[9], 
            vert[10], 
            vert[11], 0.0F, 0.0F, 0.0F)) {
          dest.normal(vert[9], 
              vert[10], 
              vert[11]);
        }
        dest.vertex(vert[0], vert[1], vert[2], 
            vert[7], 
            vert[8]);
      } 
    } else {
      for (int i = 0; i < src.vertexCount; i++) {
        float[] vert = src.vertices[i];
        if (vert[2] == 0.0F) {
          dest.vertex(vert[0], vert[1]);
        } else {
          dest.vertex(vert[0], vert[1], vert[2]);
        } 
      } 
    } 
    
    dest.endShape();
  }


  
  protected static void copyPath(PShape src, PShape dest) {
    copyMatrix(src, dest);
    copyStyles(src, dest);
    copyImage(src, dest);
    dest.close = src.close;
    dest.setPath(src.vertexCount, src.vertices, src.vertexCodeCount, src.vertexCodes);
  }


  
  protected static void copyMatrix(PShape src, PShape dest) {
    if (src.matrix != null) {
      dest.applyMatrix(src.matrix);
    }
  }


  
  protected static void copyStyles(PShape src, PShape dest) {
    dest.ellipseMode = src.ellipseMode;
    dest.rectMode = src.rectMode;
    
    if (src.stroke) {
      dest.stroke = true;
      dest.strokeColor = src.strokeColor;
      dest.strokeWeight = src.strokeWeight;
      dest.strokeCap = src.strokeCap;
      dest.strokeJoin = src.strokeJoin;
    } else {
      dest.stroke = false;
    } 
    
    if (src.fill) {
      dest.fill = true;
      dest.fillColor = src.fillColor;
    } else {
      dest.fill = false;
    } 
  }


  
  protected static void copyImage(PShape src, PShape dest) {
    if (src.image != null) {
      dest.texture(src.image);
    }
  }










  
  public void draw(PGraphics g) {
    if (this.visible) {
      pre(g);
      drawImpl(g);
      post(g);
    } 
  }




  
  protected void drawImpl(PGraphics g) {
    if (this.family == 0) {
      drawGroup(g);
    } else if (this.family == 101) {
      drawPrimitive(g);
    } else if (this.family == 103) {

      
      drawGeometry(g);
    } else if (this.family == 102) {
      drawPath(g);
    } 
  }

  
  protected void drawGroup(PGraphics g) {
    for (int i = 0; i < this.childCount; i++) {
      this.children[i].draw(g);
    }
  }

  
  protected void drawPrimitive(PGraphics g) {
    if (this.kind == 2) {
      g.point(this.params[0], this.params[1]);
    }
    else if (this.kind == 4) {
      if (this.params.length == 4) {
        g.line(this.params[0], this.params[1], 
            this.params[2], this.params[3]);
      } else {
        g.line(this.params[0], this.params[1], this.params[2], 
            this.params[3], this.params[4], this.params[5]);
      }
    
    } else if (this.kind == 8) {
      g.triangle(this.params[0], this.params[1], 
          this.params[2], this.params[3], 
          this.params[4], this.params[5]);
    }
    else if (this.kind == 16) {
      g.quad(this.params[0], this.params[1], 
          this.params[2], this.params[3], 
          this.params[4], this.params[5], 
          this.params[6], this.params[7]);
    }
    else if (this.kind == 30) {
      
      if (this.imagePath != null) {
        loadImage(g);
      }
      if (this.image != null) {
        int oldMode = g.imageMode;
        g.imageMode(0);
        g.image(this.image, this.params[0], this.params[1], this.params[2], this.params[3]);
        g.imageMode(oldMode);
      } else {
        int oldMode = g.rectMode;
        g.rectMode(this.rectMode);
        if (this.params.length == 4) {
          g.rect(this.params[0], this.params[1], 
              this.params[2], this.params[3]);
        } else if (this.params.length == 5) {
          g.rect(this.params[0], this.params[1], 
              this.params[2], this.params[3], 
              this.params[4]);
        } else if (this.params.length == 8) {
          g.rect(this.params[0], this.params[1], 
              this.params[2], this.params[3], 
              this.params[4], this.params[5], 
              this.params[6], this.params[7]);
        } 
        g.rectMode(oldMode);
      } 
    } else if (this.kind == 31) {
      int oldMode = g.ellipseMode;
      g.ellipseMode(this.ellipseMode);
      g.ellipse(this.params[0], this.params[1], 
          this.params[2], this.params[3]);
      g.ellipseMode(oldMode);
    }
    else if (this.kind == 32) {
      int oldMode = g.ellipseMode;
      g.ellipseMode(this.ellipseMode);
      if (this.params.length == 6) {
        g.arc(this.params[0], this.params[1], 
            this.params[2], this.params[3], 
            this.params[4], this.params[5]);
      } else if (this.params.length == 7) {
        g.arc(this.params[0], this.params[1], 
            this.params[2], this.params[3], 
            this.params[4], this.params[5], 
            (int)this.params[6]);
      } 
      g.ellipseMode(oldMode);
    }
    else if (this.kind == 41) {
      if (this.params.length == 1) {
        g.box(this.params[0]);
      } else {
        g.box(this.params[0], this.params[1], this.params[2]);
      }
    
    } else if (this.kind == 40) {
      g.sphere(this.params[0]);
    } 
  }


  
  protected void drawGeometry(PGraphics g) {
    g.beginShape(this.kind);
    if (this.style) {
      for (int i = 0; i < this.vertexCount; i++) {
        g.vertex(this.vertices[i]);
      }
    } else {
      for (int i = 0; i < this.vertexCount; i++) {
        float[] vert = this.vertices[i];
        if (vert[2] == 0.0F) {
          g.vertex(vert[0], vert[1]);
        } else {
          g.vertex(vert[0], vert[1], vert[2]);
        } 
      } 
    } 
    g.endShape(this.close ? 2 : 1);
  }

























































  
  protected void drawPath(PGraphics g) {
    if (this.vertices == null)
      return; 
    boolean insideContour = false;
    g.beginShape();
    
    if (this.vertexCodeCount == 0) {
      if (this.vertices[0].length == 2) {
        for (int i = 0; i < this.vertexCount; i++) {
          g.vertex(this.vertices[i][0], this.vertices[i][1]);
        }
      } else {
        for (int i = 0; i < this.vertexCount; i++) {
          g.vertex(this.vertices[i][0], this.vertices[i][1], this.vertices[i][2]);
        }
      } 
    } else {
      
      int index = 0;
      
      if (this.vertices[0].length == 2) {
        for (int j = 0; j < this.vertexCodeCount; j++) {
          switch (this.vertexCodes[j]) {
            
            case 0:
              g.vertex(this.vertices[index][0], this.vertices[index][1]);
              index++;
              break;
            
            case 2:
              g.quadraticVertex(this.vertices[index + 0][0], this.vertices[index + 0][1], 
                  this.vertices[index + 1][0], this.vertices[index + 1][1]);
              index += 2;
              break;
            
            case 1:
              g.bezierVertex(this.vertices[index + 0][0], this.vertices[index + 0][1], 
                  this.vertices[index + 1][0], this.vertices[index + 1][1], 
                  this.vertices[index + 2][0], this.vertices[index + 2][1]);
              index += 3;
              break;
            
            case 3:
              g.curveVertex(this.vertices[index][0], this.vertices[index][1]);
              index++;
              break;
            
            case 4:
              if (insideContour) {
                g.endContour();
              }
              g.beginContour();
              insideContour = true; break;
          } 
        } 
      } else {
        for (int j = 0; j < this.vertexCodeCount; j++) {
          switch (this.vertexCodes[j]) {
            
            case 0:
              g.vertex(this.vertices[index][0], this.vertices[index][1], this.vertices[index][2]);
              index++;
              break;
            
            case 2:
              g.quadraticVertex(this.vertices[index + 0][0], this.vertices[index + 0][1], this.vertices[index + 0][2], 
                  this.vertices[index + 1][0], this.vertices[index + 1][1], this.vertices[index + 0][2]);
              index += 2;
              break;

            
            case 1:
              g.bezierVertex(this.vertices[index + 0][0], this.vertices[index + 0][1], this.vertices[index + 0][2], 
                  this.vertices[index + 1][0], this.vertices[index + 1][1], this.vertices[index + 1][2], 
                  this.vertices[index + 2][0], this.vertices[index + 2][1], this.vertices[index + 2][2]);
              index += 3;
              break;
            
            case 3:
              g.curveVertex(this.vertices[index][0], this.vertices[index][1], this.vertices[index][2]);
              index++;
              break;
            
            case 4:
              if (insideContour) {
                g.endContour();
              }
              g.beginContour();
              insideContour = true; break;
          } 
        } 
      } 
    } 
    if (insideContour) {
      g.endContour();
    }
    g.endShape(this.close ? 2 : 1);
  }

  
  private void loadImage(PGraphics g) {
    if (this.imagePath.startsWith("data:image")) {
      loadBase64Image();
    }
    
    if (this.imagePath.startsWith("file://")) {
      loadFileSystemImage(g);
    }
    this.imagePath = null;
  }
  
  private void loadFileSystemImage(PGraphics g) {
    this.imagePath = this.imagePath.substring(7);
    PImage loadedImage = g.parent.loadImage(this.imagePath);
    if (loadedImage == null) {
      System.err.println("Error loading image file: " + this.imagePath);
    } else {
      setTexture(loadedImage);
    } 
  }
  
  private void loadBase64Image() {
    String[] parts = this.imagePath.split(";base64,");
    String extension = parts[0].substring(11);
    String encodedData = parts[1];
    
    byte[] decodedBytes = DatatypeConverter.parseBase64Binary(encodedData);
    
    if (decodedBytes == null) {
      System.err.println("Decode Error on image: " + this.imagePath.substring(0, 20));
      
      return;
    } 
    Image awtImage = (new ImageIcon(decodedBytes)).getImage();
    
    if (awtImage instanceof BufferedImage) {
      BufferedImage buffImage = (BufferedImage)awtImage;
      int space = buffImage.getColorModel().getColorSpace().getType();
      if (space == 9) {
        return;
      }
    } 
    
    PImage loadedImage = new PImage(awtImage);
    loadedImage.width;



    
    if (extension.equals("gif") || extension.equals("png") || 
      extension.equals("unknown")) {
      loadedImage.checkAlpha();
    }
    
    setTexture(loadedImage);
  }




  
  public PShape getParent() { return this.parent; }






  
  public int getChildCount() { return this.childCount; }




  
  protected void crop() {
    if (this.children.length != this.childCount) {
      this.children = (PShape[])PApplet.subset(this.children, 0, this.childCount);
    }
  }

  
  public PShape[] getChildren() {
    crop();
    return this.children;
  }














  
  public PShape getChild(int index) {
    crop();
    return this.children[index];
  }



  
  public PShape getChild(String target) {
    if (this.name != null && this.name.equals(target)) {
      return this;
    }
    if (this.nameTable != null) {
      PShape found = (PShape)this.nameTable.get(target);
      if (found != null) return found; 
    } 
    for (int i = 0; i < this.childCount; i++) {
      PShape found = this.children[i].getChild(target);
      if (found != null) return found; 
    } 
    return null;
  }





  
  public PShape findChild(String target) {
    if (this.parent == null) {
      return getChild(target);
    }
    
    return this.parent.findChild(target);
  }









  
  public void addChild(PShape who) {
    if (this.children == null) {
      this.children = new PShape[1];
    }
    if (this.childCount == this.children.length) {
      this.children = (PShape[])PApplet.expand(this.children);
    }
    this.children[this.childCount++] = who;
    who.parent = this;
    
    if (who.getName() != null) {
      addName(who.getName(), who);
    }
  }





  
  public void addChild(PShape who, int idx) {
    if (idx < this.childCount) {
      if (this.childCount == this.children.length) {
        this.children = (PShape[])PApplet.expand(this.children);
      }

      
      for (int i = this.childCount - 1; i >= idx; i--) {
        this.children[i + 1] = this.children[i];
      }
      this.childCount++;
      
      this.children[idx] = who;
      
      who.parent = this;
      
      if (who.getName() != null) {
        addName(who.getName(), who);
      }
    } 
  }




  
  public void removeChild(int idx) {
    if (idx < this.childCount) {
      PShape child = this.children[idx];

      
      for (int i = idx; i < this.childCount - 1; i++) {
        this.children[i] = this.children[i + 1];
      }
      this.childCount--;
      
      if (child.getName() != null && this.nameTable != null) {
        this.nameTable.remove(child.getName());
      }
    } 
  }




  
  public void addName(String nom, PShape shape) {
    if (this.parent != null) {
      this.parent.addName(nom, shape);
    } else {
      if (this.nameTable == null) {
        this.nameTable = new HashMap();
      }
      this.nameTable.put(nom, shape);
    } 
  }




  
  public int getChildIndex(PShape who) {
    for (int i = 0; i < this.childCount; i++) {
      if (this.children[i] == who) {
        return i;
      }
    } 
    return -1;
  }


  
  public PShape getTessellation() { return null; }







  
  public int getFamily() { return this.family; }



  
  public int getKind() { return this.kind; }



  
  public float[] getParams() { return getParams(null); }


  
  public float[] getParams(float[] target) {
    if (target == null || target.length != this.params.length) {
      target = new float[this.params.length];
    }
    PApplet.arrayCopy(this.params, target);
    return target;
  }


  
  public float getParam(int index) { return this.params[index]; }


  
  protected void setParams(float[] source) {
    if (this.params == null) {
      this.params = new float[source.length];
    }
    if (source.length != this.params.length) {
      PGraphics.showWarning("Wrong number of parameters");
      return;
    } 
    PApplet.arrayCopy(source, this.params);
  }


  
  public void setPath(int vcount, float[][] verts) { setPath(vcount, verts, 0, null); }


  
  protected void setPath(int vcount, float[][] verts, int ccount, int[] codes) {
    if (verts == null || verts.length < vcount)
      return;  if (ccount > 0 && (codes == null || codes.length < ccount))
      return; 
    int ndim = verts[0].length;
    this.vertexCount = vcount;
    this.vertices = new float[this.vertexCount][ndim];
    for (int i = 0; i < this.vertexCount; i++) {
      PApplet.arrayCopy(verts[i], this.vertices[i]);
    }
    
    this.vertexCodeCount = ccount;
    if (this.vertexCodeCount > 0) {
      this.vertexCodes = new int[this.vertexCodeCount];
      PApplet.arrayCopy(codes, this.vertexCodes, this.vertexCodeCount);
    } 
  }






  
  public int getVertexCount() {
    if (this.family == 0 || this.family == 101) {
      PGraphics.showWarning("getVertexCount() only works with PATH or GEOMETRY shapes");
    }
    return this.vertexCount;
  }









  
  public PVector getVertex(int index) { return getVertex(index, null); }





  
  public PVector getVertex(int index, PVector vec) {
    if (vec == null) {
      vec = new PVector();
    }
    float[] vert = this.vertices[index];
    vec.x = vert[0];
    vec.y = vert[1];
    if (vert.length > 2) {
      vec.z = vert[2];
    } else {
      vec.z = 0.0F;
    } 
    return vec;
  }


  
  public float getVertexX(int index) { return this.vertices[index][0]; }



  
  public float getVertexY(int index) { return this.vertices[index][1]; }



  
  public float getVertexZ(int index) { return this.vertices[index][2]; }











  
  public void setVertex(int index, float x, float y) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setVertex()" });
      
      return;
    } 
    this.vertices[index][0] = x;
    this.vertices[index][1] = y;
  }




  
  public void setVertex(int index, float x, float y, float z) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setVertex()" });
      
      return;
    } 
    this.vertices[index][0] = x;
    this.vertices[index][1] = y;
    this.vertices[index][2] = z;
  }




  
  public void setVertex(int index, PVector vec) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setVertex()" });
      
      return;
    } 
    this.vertices[index][0] = vec.x;
    this.vertices[index][1] = vec.y;
    
    if (this.vertices[index].length > 2) {
      this.vertices[index][2] = vec.z;
    } else if (vec.z != 0.0F && vec.z == vec.z) {
      throw new IllegalArgumentException("Cannot set a z-coordinate on a 2D shape");
    } 
  }


  
  public PVector getNormal(int index) { return getNormal(index, null); }


  
  public PVector getNormal(int index, PVector vec) {
    if (vec == null) {
      vec = new PVector();
    }
    vec.x = this.vertices[index][9];
    vec.y = this.vertices[index][10];
    vec.z = this.vertices[index][11];
    return vec;
  }


  
  public float getNormalX(int index) { return this.vertices[index][9]; }



  
  public float getNormalY(int index) { return this.vertices[index][10]; }



  
  public float getNormalZ(int index) { return this.vertices[index][11]; }


  
  public void setNormal(int index, float nx, float ny, float nz) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setNormal()" });
      
      return;
    } 
    this.vertices[index][9] = nx;
    this.vertices[index][10] = ny;
    this.vertices[index][11] = nz;
  }



  
  public void setAttrib(String name, int index, float... values) {}


  
  public void setAttrib(String name, int index, int... values) {}


  
  public void setAttrib(String name, int index, boolean... values) {}


  
  public float getTextureU(int index) { return this.vertices[index][7]; }



  
  public float getTextureV(int index) { return this.vertices[index][8]; }


  
  public void setTextureUV(int index, float u, float v) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setTextureUV()" });
      
      return;
    } 
    
    if (this.vertices == null || 
      index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "setTextureUV()" });
      
      return;
    } 
    
    this.vertices[index][7] = u;
    this.vertices[index][8] = v;
  }

  
  public void setTextureMode(int mode) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setTextureMode()" });
      
      return;
    } 
    this.textureMode = mode;
  }

  
  public void setTexture(PImage tex) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setTexture()" });
      
      return;
    } 
    this.image = tex;
  }


  
  public int getFill(int index) {
    if (this.vertices == null || 
      index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "getFill()" });
      return this.fillColor;
    } 
    
    if (this.image == null) {
      int a = (int)(this.vertices[index][6] * 255.0F);
      int r = (int)(this.vertices[index][3] * 255.0F);
      int g = (int)(this.vertices[index][4] * 255.0F);
      int b = (int)(this.vertices[index][5] * 255.0F);
      return a << 24 | r << 16 | g << 8 | b;
    } 
    return 0;
  }




  
  public void setFill(boolean fill) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setFill()" });
      
      return;
    } 
    this.fill = fill;
  }


















  
  public void setFill(int fill) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setFill()" });
      
      return;
    } 
    this.fillColor = fill;
    
    if (this.vertices != null && this.perVertexStyles) {
      for (int i = 0; i < this.vertexCount; i++) {
        setFill(i, fill);
      }
    }
  }



  
  public void setFill(int index, int fill) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setFill()" });
      
      return;
    } 
    if (!this.perVertexStyles) {
      PGraphics.showWarning("This renderer does not support %1$s for individual vertices", new Object[] { "setFill()" });
      
      return;
    } 
    
    if (this.vertices == null || index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "getFill()" });
      
      return;
    } 
    if (this.image == null) {
      this.vertices[index][6] = (fill >> 24 & 0xFF) / 255.0F;
      this.vertices[index][3] = (fill >> 16 & 0xFF) / 255.0F;
      this.vertices[index][4] = (fill >> 8 & 0xFF) / 255.0F;
      this.vertices[index][5] = (fill >> 0 & 0xFF) / 255.0F;
    } 
  }


  
  public int getTint(int index) {
    if (this.vertices == null || index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "getTint()" });
      return this.tintColor;
    } 
    
    if (this.image != null) {
      int a = (int)(this.vertices[index][6] * 255.0F);
      int r = (int)(this.vertices[index][3] * 255.0F);
      int g = (int)(this.vertices[index][4] * 255.0F);
      int b = (int)(this.vertices[index][5] * 255.0F);
      return a << 24 | r << 16 | g << 8 | b;
    } 
    return 0;
  }


  
  public void setTint(boolean tint) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setTint()" });
      
      return;
    } 
    this.tint = tint;
  }

  
  public void setTint(int fill) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setTint()" });
      
      return;
    } 
    this.tintColor = fill;
    
    if (this.vertices != null) {
      for (int i = 0; i < this.vertices.length; i++) {
        setFill(i, fill);
      }
    }
  }

  
  public void setTint(int index, int tint) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setTint()" });
      
      return;
    } 
    
    if (this.vertices == null || 
      index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "setTint()" });
      
      return;
    } 
    if (this.image != null) {
      this.vertices[index][6] = (tint >> 24 & 0xFF) / 255.0F;
      this.vertices[index][3] = (tint >> 16 & 0xFF) / 255.0F;
      this.vertices[index][4] = (tint >> 8 & 0xFF) / 255.0F;
      this.vertices[index][5] = (tint >> 0 & 0xFF) / 255.0F;
    } 
  }


  
  public int getStroke(int index) {
    if (this.vertices == null || 
      index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "getStroke()" });
      return this.strokeColor;
    } 
    
    int a = (int)(this.vertices[index][16] * 255.0F);
    int r = (int)(this.vertices[index][13] * 255.0F);
    int g = (int)(this.vertices[index][14] * 255.0F);
    int b = (int)(this.vertices[index][15] * 255.0F);
    return a << 24 | r << 16 | g << 8 | b;
  }



  
  public void setStroke(boolean stroke) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setStroke()" });
      
      return;
    } 
    this.stroke = stroke;
  }


















  
  public void setStroke(int stroke) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setStroke()" });
      
      return;
    } 
    this.strokeColor = stroke;
    
    if (this.vertices != null && this.perVertexStyles) {
      for (int i = 0; i < this.vertices.length; i++) {
        setStroke(i, stroke);
      }
    }
  }



  
  public void setStroke(int index, int stroke) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setStroke()" });
      
      return;
    } 
    if (!this.perVertexStyles) {
      PGraphics.showWarning("This renderer does not support %1$s for individual vertices", new Object[] { "setStroke()" });
      
      return;
    } 
    
    if (this.vertices == null || index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "setStroke()" });
      
      return;
    } 
    this.vertices[index][16] = (stroke >> 24 & 0xFF) / 255.0F;
    this.vertices[index][13] = (stroke >> 16 & 0xFF) / 255.0F;
    this.vertices[index][14] = (stroke >> 8 & 0xFF) / 255.0F;
    this.vertices[index][15] = (stroke >> 0 & 0xFF) / 255.0F;
  }


  
  public float getStrokeWeight(int index) {
    if (this.vertices == null || index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "getStrokeWeight()" });
      return this.strokeWeight;
    } 
    
    return this.vertices[index][17];
  }

  
  public void setStrokeWeight(float weight) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setStrokeWeight()" });
      
      return;
    } 
    this.strokeWeight = weight;
    
    if (this.vertices != null && this.perVertexStyles) {
      for (int i = 0; i < this.vertexCount; i++) {
        setStrokeWeight(i, weight);
      }
    }
  }

  
  public void setStrokeWeight(int index, float weight) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setStrokeWeight()" });
      
      return;
    } 
    if (!this.perVertexStyles) {
      PGraphics.showWarning("This renderer does not support %1$s for individual vertices", new Object[] { "setStrokeWeight()" });
      
      return;
    } 
    
    if (this.vertices == null || index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "setStrokeWeight()" });
      
      return;
    } 
    this.vertices[index][17] = weight;
  }

  
  public void setStrokeJoin(int join) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setStrokeJoin()" });
      
      return;
    } 
    this.strokeJoin = join;
  }

  
  public void setStrokeCap(int cap) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setStrokeCap()" });
      
      return;
    } 
    this.strokeCap = cap;
  }


  
  public int getAmbient(int index) {
    if (this.vertices == null || index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "getAmbient()" });
      return this.ambientColor;
    } 
    
    int r = (int)(this.vertices[index][25] * 255.0F);
    int g = (int)(this.vertices[index][26] * 255.0F);
    int b = (int)(this.vertices[index][27] * 255.0F);
    return 0xFF000000 | r << 16 | g << 8 | b;
  }

  
  public void setAmbient(int ambient) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setAmbient()" });
      
      return;
    } 
    this.ambientColor = ambient;
    
    if (this.vertices != null) {
      for (int i = 0; i < this.vertices.length; i++) {
        setAmbient(i, ambient);
      }
    }
  }

  
  public void setAmbient(int index, int ambient) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setAmbient()" });
      
      return;
    } 
    
    if (this.vertices == null || index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "setAmbient()" });
      
      return;
    } 
    this.vertices[index][25] = (ambient >> 16 & 0xFF) / 255.0F;
    this.vertices[index][26] = (ambient >> 8 & 0xFF) / 255.0F;
    this.vertices[index][27] = (ambient >> 0 & 0xFF) / 255.0F;
  }


  
  public int getSpecular(int index) {
    if (this.vertices == null || index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "getSpecular()" });
      return this.specularColor;
    } 
    
    int r = (int)(this.vertices[index][28] * 255.0F);
    int g = (int)(this.vertices[index][29] * 255.0F);
    int b = (int)(this.vertices[index][30] * 255.0F);
    return 0xFF000000 | r << 16 | g << 8 | b;
  }

  
  public void setSpecular(int specular) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setSpecular()" });
      
      return;
    } 
    this.specularColor = specular;
    
    if (this.vertices != null) {
      for (int i = 0; i < this.vertices.length; i++) {
        setSpecular(i, specular);
      }
    }
  }

  
  public void setSpecular(int index, int specular) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setSpecular()" });
      
      return;
    } 
    
    if (this.vertices == null || index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "setSpecular()" });
      
      return;
    } 
    this.vertices[index][28] = (specular >> 16 & 0xFF) / 255.0F;
    this.vertices[index][29] = (specular >> 8 & 0xFF) / 255.0F;
    this.vertices[index][30] = (specular >> 0 & 0xFF) / 255.0F;
  }


  
  public int getEmissive(int index) {
    if (this.vertices == null || index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "getEmissive()" });
      return this.emissiveColor;
    } 
    
    int r = (int)(this.vertices[index][32] * 255.0F);
    int g = (int)(this.vertices[index][33] * 255.0F);
    int b = (int)(this.vertices[index][34] * 255.0F);
    return 0xFF000000 | r << 16 | g << 8 | b;
  }

  
  public void setEmissive(int emissive) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setEmissive()" });
      
      return;
    } 
    this.emissiveColor = emissive;
    
    if (this.vertices != null) {
      for (int i = 0; i < this.vertices.length; i++) {
        setEmissive(i, emissive);
      }
    }
  }

  
  public void setEmissive(int index, int emissive) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setEmissive()" });
      
      return;
    } 
    
    if (this.vertices == null || 
      index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "setEmissive()" });
      
      return;
    } 
    this.vertices[index][32] = (emissive >> 16 & 0xFF) / 255.0F;
    this.vertices[index][33] = (emissive >> 8 & 0xFF) / 255.0F;
    this.vertices[index][34] = (emissive >> 0 & 0xFF) / 255.0F;
  }


  
  public float getShininess(int index) {
    if (this.vertices == null || 
      index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "getShininess()" });
      return this.shininess;
    } 
    
    return this.vertices[index][31];
  }

  
  public void setShininess(float shine) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setShininess()" });
      
      return;
    } 
    this.shininess = shine;
    
    if (this.vertices != null) {
      for (int i = 0; i < this.vertices.length; i++) {
        setShininess(i, shine);
      }
    }
  }

  
  public void setShininess(int index, float shine) {
    if (this.openShape) {
      PGraphics.showWarning("%1$s can only be called outside beginShape() and endShape()", new Object[] { "setShininess()" });
      
      return;
    } 
    
    if (this.vertices == null || 
      index >= this.vertices.length) {
      PGraphics.showWarning("%1$s vertex index does not exist (" + index + ")", new Object[] { "setShininess()" });
      
      return;
    } 
    
    this.vertices[index][31] = shine;
  }

  
  public int[] getVertexCodes() {
    if (this.vertexCodes == null) {
      return null;
    }
    if (this.vertexCodes.length != this.vertexCodeCount) {
      this.vertexCodes = PApplet.subset(this.vertexCodes, 0, this.vertexCodeCount);
    }
    return this.vertexCodes;
  }


  
  public int getVertexCodeCount() { return this.vertexCodeCount; }






  
  public int getVertexCode(int index) { return this.vertexCodes[index]; }



  
  public boolean isClosed() { return this.close; }









  
  public boolean contains(float x, float y) {
    if (this.family == 102) {
      
      PMatrix inverseCoords = this.matrix.get();
      inverseCoords.invert();
      inverseCoords.invert();
      PVector p = new PVector();
      inverseCoords.mult(new PVector(x, y), p);

      
      boolean c = false;
      for (int i = 0, j = this.vertexCount - 1; i < this.vertexCount; j = i++) {
        if (((this.vertices[i][1] > p.y) ? 1 : 0) != ((this.vertices[j][1] > p.y) ? 1 : 0))
          if (p.x < (
            this.vertices[j][0] - this.vertices[i][0]) * (
            y - this.vertices[i][1]) / (
            this.vertices[j][1] - this.vertices[i][1]) + 
            this.vertices[i][0]) {
            c = !c;
          } 
      } 
      return c;
    } 
    if (this.family == 0) {


      
      for (int i = 0; i < this.childCount; i++) {
        if (this.children[i].contains(x, y)) return true; 
      } 
      return false;
    } 

    
    throw new IllegalArgumentException("The contains() method is only implemented for paths.");
  }




































  
  public void translate(float x, float y) {
    checkMatrix(2);
    this.matrix.translate(x, y);
  }



  
  public void translate(float x, float y, float z) {
    checkMatrix(3);
    this.matrix.translate(x, y, z);
  }






























  
  public void rotateX(float angle) { rotate(angle, 1.0F, 0.0F, 0.0F); }
































  
  public void rotateY(float angle) { rotate(angle, 0.0F, 1.0F, 0.0F); }
































  
  public void rotateZ(float angle) { rotate(angle, 0.0F, 0.0F, 1.0F); }




























  
  public void rotate(float angle) {
    checkMatrix(2);
    this.matrix.rotate(angle);
  }



  
  public void rotate(float angle, float v0, float v1, float v2) {
    checkMatrix(3);
    float norm2 = v0 * v0 + v1 * v1 + v2 * v2;
    if (Math.abs(norm2 - 1.0F) > 1.0E-4F) {
      
      float norm = PApplet.sqrt(norm2);
      v0 /= norm;
      v1 /= norm;
      v2 /= norm;
    } 
    this.matrix.rotate(angle, v0, v1, v2);
  }



























  
  public void scale(float s) {
    checkMatrix(2);
    this.matrix.scale(s);
  }

  
  public void scale(float x, float y) {
    checkMatrix(2);
    this.matrix.scale(x, y);
  }





  
  public void scale(float x, float y, float z) {
    checkMatrix(3);
    this.matrix.scale(x, y, z);
  }

















  
  public void resetMatrix() {
    checkMatrix(2);
    this.matrix.reset();
  }

  
  public void applyMatrix(PMatrix source) {
    if (source instanceof PMatrix2D) {
      applyMatrix((PMatrix2D)source);
    } else if (source instanceof PMatrix3D) {
      applyMatrix((PMatrix3D)source);
    } 
  }

  
  public void applyMatrix(PMatrix2D source) {
    applyMatrix(source.m00, source.m01, 0.0F, source.m02, 
        source.m10, source.m11, 0.0F, source.m12, 
        0.0F, 0.0F, 1.0F, 0.0F, 
        0.0F, 0.0F, 0.0F, 1.0F);
  }


  
  public void applyMatrix(float n00, float n01, float n02, float n10, float n11, float n12) {
    checkMatrix(2);
    this.matrix.apply(n00, n01, n02, 
        n10, n11, n12);
  }

  
  public void applyMatrix(PMatrix3D source) {
    applyMatrix(source.m00, source.m01, source.m02, source.m03, 
        source.m10, source.m11, source.m12, source.m13, 
        source.m20, source.m21, source.m22, source.m23, 
        source.m30, source.m31, source.m32, source.m33);
  }




  
  public void applyMatrix(float n00, float n01, float n02, float n03, float n10, float n11, float n12, float n13, float n20, float n21, float n22, float n23, float n30, float n31, float n32, float n33) {
    checkMatrix(3);
    this.matrix.apply(n00, n01, n02, n03, 
        n10, n11, n12, n13, 
        n20, n21, n22, n23, 
        n30, n31, n32, n33);
  }








  
  protected void checkMatrix(int dimensions) {
    if (this.matrix == null) {
      if (dimensions == 2) {
        this.matrix = new PMatrix2D();
      } else {
        this.matrix = new PMatrix3D();
      } 
    } else if (dimensions == 3 && this.matrix instanceof PMatrix2D) {
      
      this.matrix = new PMatrix3D(this.matrix);
    } 
  }





























  
  public void colorMode(int mode) { colorMode(mode, this.colorModeX, this.colorModeY, this.colorModeZ, this.colorModeA); }





  
  public void colorMode(int mode, float max) { colorMode(mode, max, max, max, max); }








  
  public void colorMode(int mode, float maxX, float maxY, float maxZ) { colorMode(mode, maxX, maxY, maxZ, this.colorModeA); }





  
  public void colorMode(int mode, float maxX, float maxY, float maxZ, float maxA) {
    this.colorMode = mode;
    
    this.colorModeX = maxX;
    this.colorModeY = maxY;
    this.colorModeZ = maxZ;
    this.colorModeA = maxA;

    
    this.colorModeScale = !(
      maxA == 1.0F && maxX == maxY && maxY == maxZ && maxZ == maxA);


    
    this.colorModeDefault = (this.colorMode == 1 && 
      this.colorModeA == 255.0F && this.colorModeX == 255.0F && 
      this.colorModeY == 255.0F && this.colorModeZ == 255.0F);
  }

  
  protected void colorCalc(int rgb) {
    if ((rgb & 0xFF000000) == 0 && rgb <= this.colorModeX) {
      colorCalc(rgb);
    } else {
      
      colorCalcARGB(rgb, this.colorModeA);
    } 
  }

  
  protected void colorCalc(int rgb, float alpha) {
    if ((rgb & 0xFF000000) == 0 && rgb <= this.colorModeX) {
      colorCalc(rgb, alpha);
    } else {
      
      colorCalcARGB(rgb, alpha);
    } 
  }


  
  protected void colorCalc(float gray) { colorCalc(gray, this.colorModeA); }


  
  protected void colorCalc(float gray, float alpha) {
    if (gray > this.colorModeX) gray = this.colorModeX; 
    if (alpha > this.colorModeA) alpha = this.colorModeA;
    
    if (gray < 0.0F) gray = 0.0F; 
    if (alpha < 0.0F) alpha = 0.0F;
    
    this.calcR = this.colorModeScale ? (gray / this.colorModeX) : gray;
    this.calcG = this.calcR;
    this.calcB = this.calcR;
    this.calcA = this.colorModeScale ? (alpha / this.colorModeA) : alpha;
    
    this.calcRi = (int)(this.calcR * 255.0F); this.calcGi = (int)(this.calcG * 255.0F);
    this.calcBi = (int)(this.calcB * 255.0F); this.calcAi = (int)(this.calcA * 255.0F);
    this.calcColor = this.calcAi << 24 | this.calcRi << 16 | this.calcGi << 8 | this.calcBi;
    this.calcAlpha = (this.calcAi != 255);
  }


  
  protected void colorCalc(float x, float y, float z) { colorCalc(x, y, z, this.colorModeA); }

  
  protected void colorCalc(float x, float y, float z, float a) {
    float t, q, p, f, which;
    if (x > this.colorModeX) x = this.colorModeX; 
    if (y > this.colorModeY) y = this.colorModeY; 
    if (z > this.colorModeZ) z = this.colorModeZ; 
    if (a > this.colorModeA) a = this.colorModeA;
    
    if (x < 0.0F) x = 0.0F; 
    if (y < 0.0F) y = 0.0F; 
    if (z < 0.0F) z = 0.0F; 
    if (a < 0.0F) a = 0.0F;
    
    switch (this.colorMode) {
      case 1:
        if (this.colorModeScale) {
          this.calcR = x / this.colorModeX;
          this.calcG = y / this.colorModeY;
          this.calcB = z / this.colorModeZ;
          this.calcA = a / this.colorModeA; break;
        } 
        this.calcR = x; this.calcG = y; this.calcB = z; this.calcA = a;
        break;

      
      case 3:
        x /= this.colorModeX;
        y /= this.colorModeY;
        z /= this.colorModeZ;
        
        this.calcA = this.colorModeScale ? (a / this.colorModeA) : a;
        
        if (y == 0.0F) {
          this.calcR = this.calcG = this.calcB = z;
          break;
        } 
        which = (x - (int)x) * 6.0F;
        f = which - (int)which;
        p = z * (1.0F - y);
        q = z * (1.0F - y * f);
        t = z * (1.0F - y * (1.0F - f));
        
        switch ((int)which) { case 0:
            this.calcR = z; this.calcG = t; this.calcB = p; break;
          case 1: this.calcR = q; this.calcG = z; this.calcB = p; break;
          case 2: this.calcR = p; this.calcG = z; this.calcB = t; break;
          case 3: this.calcR = p; this.calcG = q; this.calcB = z; break;
          case 4: this.calcR = t; this.calcG = p; this.calcB = z; break;
          case 5: case 0: break; }  this.calcR = z; this.calcG = p; this.calcB = q;
        break;
    } 

    
    this.calcRi = (int)(255.0F * this.calcR); this.calcGi = (int)(255.0F * this.calcG);
    this.calcBi = (int)(255.0F * this.calcB); this.calcAi = (int)(255.0F * this.calcA);
    this.calcColor = this.calcAi << 24 | this.calcRi << 16 | this.calcGi << 8 | this.calcBi;
    this.calcAlpha = (this.calcAi != 255);
  }

  
  protected void colorCalcARGB(int argb, float alpha) {
    if (alpha == this.colorModeA) {
      this.calcAi = argb >> 24 & 0xFF;
      this.calcColor = argb;
    } else {
      this.calcAi = (int)((argb >> 24 & 0xFF) * alpha / this.colorModeA);
      this.calcColor = this.calcAi << 24 | argb & 0xFFFFFF;
    } 
    this.calcRi = argb >> 16 & 0xFF;
    this.calcGi = argb >> 8 & 0xFF;
    this.calcBi = argb & 0xFF;
    this.calcA = this.calcAi / 255.0F;
    this.calcR = this.calcRi / 255.0F;
    this.calcG = this.calcGi / 255.0F;
    this.calcB = this.calcBi / 255.0F;
    this.calcAlpha = (this.calcAi != 255);
  }
}
