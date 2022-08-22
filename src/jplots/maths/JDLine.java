package jplots.maths;

import java.util.ArrayList;
import java.util.List;

public class JDLine {
	
	private JDPoint[] points;
	private int[] colors;

	public JDLine(JDPoint... p) {
		points = p.clone();
		colors = new int[points.length*2];
		for(int i=0; i<colors.length; i++)
			colors[i] = 0x00999999;
	}
	private JDLine(JDPoint[] p, int[] c) {
		points = p.clone();
		colors = c.clone();
	}
	
	public void addVertexValues(int[] pc) {
		for(int i=0; i<points.length; i++) {
			colors[2*i] = pc[i];
			colors[2*i+1] = pc[i];
		}
	}
	public void addLineValues(int[] lc) {
		for(int i=0; i+1<points.length; i++) {
			colors[2*+1] = lc[i];
			colors[2*i+2] = lc[i];
		}
		colors[0] = colors[1];
		if(colors.length>1)
			colors[colors.length-1] = colors[colors.length-2];
	}
	
	public JDPoint[] getPoints() {
		return points;
	}
	
	public float[] getCoordsAsFloats() {
		return JPlotMath.toFloatArray1D(getCoords());
	}
	public double[] getCoords() {
		double[] coords = new double[points.length*2];
		for(int i=0; i<points.length; i++) {
			coords[ 2*i ] = points[i].x();
			coords[2*i+1] = points[i].y();
		}
		return coords;
	}
	public int[] allColors() {
		return colors;
	}
	public int[] getLineColors() {
		int[] res = new int[points.length-1];
		for(int i=0; i<res.length; i++)
			res[i] = colors[2*i+1];
		return res;
	}
	public int[] getVertexColors() {
		int[] res = new int[points.length-1];
		for(int i=0; i<res.length; i++)
			res[i] = colors[2*i];
		return res;
	}
	public List<JDEdge> toEdges() {
		List<JDEdge> edges = new ArrayList<JDEdge>();
		for(int i=1; i<points.length; i++)
			edges.add(new JDEdge(points[i-1].copy(),points[i].copy()));
		return edges;
	}
	
	public JDLine affine(double[][] transformationMatrix) {
		for(JDPoint p: points) p.affine(transformationMatrix);
		return this;
	}
	
	public boolean join(JDLine other) {
		return join(other, 0.0001d);
	}
	public boolean join(JDLine other, double delta) {
		JDPoint[] p1 = this.points.clone();
		int[] c1 = this.colors.clone();
		int len1 = p1.length;
		JDPoint[] p2 = other.points.clone(); //.getPoints();
		int[] c2 = other.colors.clone();
		int len2 = p2.length;
//		if(p1[p1.length-1].equals(p2[0], delta)) {
//			points = new JDPoint[len1+len2-1];
//			for(int i=0; i<len2; i++) points[i+len1-1] = p2[i];
//			for(int i=0; i<len1; i++) points[i]        = p1[i];
//			return true;
//		}
		if(p1[0].equals(p2[p2.length-1], delta)) {
			points = new JDPoint[len1+len2-1];
			colors = new int[2*points.length];
			for(int i=0; i<len2; i++) {
				points[i]        = p2[i];
				colors[2*i]      = c2[2*i];
				colors[2*i+1]    = c2[2*i+1];
			}
			for(int i=0; i<len1; i++) {
				points[i+len2-1] = p1[i];
				if(i>0) colors[2*(i+len2-1)] = c1[2*i];
				colors[2*(i+len2)-1] = c1[2*i+1];
			}
			return true;
		}
//		if(p1[0].equals(p2[0], delta)) {
//			points = new JDPoint[len1+len2-1];
//			for(int i=0; i<len2; i++) points[i]             = p2[i];
//			for(int i=0; i<len1; i++) points[len1+len2-i-2] = p1[i];
//			return true;
//		}
//		if(p1[p1.length-1].equals(p2[p2.length-1], delta)) {
//			points = new JDPoint[len1+len2-1];
//			for(int i=0; i<len2; i++) points[len1+len2-i-2] = p2[i];
//			for(int i=0; i<len1; i++) points[i]             = p1[i];
//			return true;
//		}
		return false;
	}
	
	public List<JDLine> intersectsCircle(JDPoint center, double radius) {
		List<JDLine> res = new ArrayList<>();
		double er2 = radius * radius;
		List<JDPoint> np = new ArrayList<>();
		List<int[]> nc = new ArrayList<>();
		double d1 = points[0].x*points[0].x + points[0].y*points[0].y - er2;
		if(d1<=0d) {
			np.add(points[0].copy());
			nc.add(new int[] {colors[0],colors[1]});
		}
		for(int i=1; i<points.length; i++) {
			double d0 = points[i].x*points[i].x + points[i].y*points[i].y - er2;
			if(d0<=0d) {
				if(d1>0d) {
					JDPoint c = center.circleCrossBetween(points[i], points[i-1], radius);
					if(c!=null) { 
						np.add(c);
						double frac = points[i].dist(c) / points[i].dist(points[i-1]);
						int col = JPlotMath.colorLerp(colors[2*i],colors[2*i-1],frac);
						nc.add(new int[] {col,col});
					}
//					np.add(points[i].fractionTowards(d0/(d0-d1), points[i-1]));
				}
				np.add(points[i].copy());
				nc.add(new int[] {colors[2*i],colors[2*i+1]});
				if(i+1==points.length) {
					res.add(new JDLine(np.toArray(new JDPoint[0]),toArray(nc)));
					np.clear();
					nc.clear();
				}
			} else
			if(d1<=0d) {
				JDPoint c = center.circleCrossBetween(points[i-1], points[i], radius);
				if(c!=null) {
					np.add(c);
					double frac = points[i-1].dist(c) / points[i-1].dist(points[i]);
					int col = JPlotMath.colorLerp(colors[2*i-1],colors[2*i],frac);
					nc.add(new int[] {col,col});
				}
//				np.add(points[i-1].fractionTowards(d1/(d1-d0), points[i]));
				res.add(new JDLine(np.toArray(new JDPoint[0]), toArray(nc)));
				np.clear();
				nc.clear();
			}
			d1 = d0;
		}
		if(res.size()>1)
			if(res.get(res.size()-1).join(res.get(0)))
				res.remove(0);
		return res;
	}
	public List<JDLine> intersectsAABB(double left, double right, double top, double bottom) {
		double xl = Math.min(left, right);
		double xr = Math.max(left, right);
		double yt = Math.min(top, bottom);
		double yb = Math.max(top, bottom);
		double dx = 1d / (xr-xl), dy = 1d / (yb-yt);
		List<JDPoint> np = new ArrayList<>();
		List<int[]>   nc = new ArrayList<>();
		List<JDPoint[]> sres = new ArrayList<>();
		List<int[]>     cres = new ArrayList<>();
		//check left & right
		double d0 = (points[0].x-xl)*dx;
		//d0 = 0d;
		if(0<=d0 && d0<=1d) {
			np.add(points[0].copy());
			nc.add(new int[] {colors[0],colors[1]});
		}
		for(int i=1; i<points.length; i++) {
			double di = (points[i].x-xl)*dx;
			double f0 = (0d-di)/(d0-di), f1 = (1d-di)/(d0-di);
			if(0<=di && di<=1d) {
				if(d0<0d) {
					np.add(points[i].fractionTowards(f0, points[i-1]));
					int col = JPlotMath.colorLerp(colors[2*i], colors[2*i-1], f0);
					nc.add(new int[] {col, col});
				}
				if(d0>1d) {
					np.add(points[i].fractionTowards(f1, points[i-1]));
					int col = JPlotMath.colorLerp(colors[2*i], colors[2*i-1], f1);
					nc.add(new int[] {col,col});
				}
				np.add(points[i].copy());
				nc.add(new int[] {colors[2*i],colors[2*i+1]});
				if(i+1==points.length) {
					sres.add(np.toArray(new JDPoint[0]));
					cres.add(toArray(nc));
					np.clear();
					nc.clear();
				}
			} else
			if(0d<=d0 && d0<=1d) {
				if(di<0d) {
					np.add(points[i].fractionTowards(f0, points[i-1]));
					int col = JPlotMath.colorLerp(colors[2*i], colors[2*i-1], f0);
					nc.add(new int[] {col, col});
				}
				if(di>1d) {
					np.add(points[i].fractionTowards(f1, points[i-1]));
					int col = JPlotMath.colorLerp(colors[2*i], colors[2*i-1], f1);
					nc.add(new int[] {col, col});
				}
				sres.add(np.toArray(new JDPoint[0]));
				cres.add(toArray(nc));
				np.clear();
				nc.clear();
			} else
			if(d0*di<0d) {
				np.add(points[i].fractionTowards(f0, points[i-1]));
				int col0 = JPlotMath.colorLerp(colors[2*i], colors[2*i-1], f0);
				np.add(points[i].fractionTowards(f1, points[i-1]));
				int col1 = JPlotMath.colorLerp(colors[2*i], colors[2*i-1], f1);
				sres.add(np.toArray(new JDPoint[0]));
				cres.add(new int[] {col0,col0,col1,col1});
				np.clear();
				nc.clear();
			}
			d0 = di;
		}
		if(sres.size()>1) {
			//check if line was a ring, then first and last line segment join together
			int idxF = 0, idxL = sres.size()-1;
			JDPoint[] first = sres.get(idxF), last = sres.get(idxL);
			int[] colF = cres.get(idxF), colL = cres.get(idxL);
			if(last[last.length-1].equals(first[0])) {
				JDPoint[] newline = new JDPoint[first.length+last.length-1];
				int[] newColors = new int[newline.length*2];
				for(int i=0; i<last.length; i++) {
					newline[i] = last[i];
					newColors[2*i] = colL[2*i];
					newColors[2*i+1] = colL[2*i+1];
				}
				for(int i=0; i<first.length; i++) {
					int j = i+last.length-1;
					newline[j] = first[i];
					newColors[2*j] = colF[2*i];
					newColors[2*j+1] = colF[2*i+1];
				}
				sres.add(newline);
				cres.add(newColors);
				sres.remove(idxL); sres.remove(idxF);
				cres.remove(idxL); cres.remove(idxF);
			}
		}
		
		//check top & bottom
		np.clear();
		nc.clear();
		List<JDLine> res = new ArrayList<>();
		for(int idx=0; idx<sres.size(); idx++) {
			JDPoint[] pnts = sres.get(idx);
			int[] cols = cres.get(idx);
			d0 = (pnts[0].y-yt)*dy;
			if(0<=d0 && d0<=1d) {
				np.add(pnts[0].copy());
				nc.add(new int[] {cols[0],cols[1]});
			}
			for(int i=1; i<pnts.length; i++) {
				double di = (pnts[i].y-yt)*dy;
				double f0 = (0d-di)/(d0-di), f1 = (1d-di)/(d0-di);
				if(0<=di && di<=1d) {
					if(d0<0d) {
						np.add(pnts[i].fractionTowards(f0, pnts[i-1]));
						int col = JPlotMath.colorLerp(cols[2*i], cols[2*i-1], f0);
						nc.add(new int[] {col,col});
					}
					if(d0>1d) {
						np.add(pnts[i].fractionTowards(f1, pnts[i-1]));
						int col = JPlotMath.colorLerp(cols[2*i], cols[2*i-1], f1);
						nc.add(new int[] {col,col});
					}
					np.add(pnts[i].copy());
					nc.add(new int[] {cols[2*i],cols[2*i+1]});
					if(i+1==pnts.length) {
						res.add(new JDLine(np.toArray(new JDPoint[0]), toArray(nc)));
						np.clear();
						nc.clear();
					}
				} else
				if(0d<=d0 && d0<=1d) {
					if(di<0d) {
						np.add(pnts[i].fractionTowards(f0, pnts[i-1]));
						int col = JPlotMath.colorLerp(cols[2*i], cols[2*i-1], f0);
						nc.add(new int[] {col,col});
					}
					if(di>1d) {
						np.add(pnts[i].fractionTowards(f1, pnts[i-1]));
						int col = JPlotMath.colorLerp(cols[2*i], cols[2*i-1], f1);
						nc.add(new int[] {col,col});
					}
					if(di<0d || 1d<di) {
						res.add(new JDLine(np.toArray(new JDPoint[0]), toArray(nc)));
						np.clear();
						nc.clear();
					}
				} else
				if(d0*di<0d) {
					np.add(pnts[i].fractionTowards(f0, pnts[i-1]));
					int col0 = JPlotMath.colorLerp(cols[2*i], cols[2*i-1], f0);
					np.add(pnts[i].fractionTowards(f1, pnts[i-1]));
					int col1 = JPlotMath.colorLerp(cols[2*i], cols[2*i-1], f1);
					res.add(new JDLine(np.toArray(new JDPoint[0]), new int[] {col0,col0,col1,col1}));
					np.clear();
					nc.clear();
				}
				d0 = di;
			}
		}
		if(res.size()>1)
			if(res.get(res.size()-1).join(res.get(0)))
				res.remove(0);
		return res;
	}
	
	public List<JDLine> cutByHalfPlane(double nx, double ny, double no) {
		List<JDLine> res = new ArrayList<>();
		List<JDPoint> np = new ArrayList<>();
		List<int[]> nc = new ArrayList<>();
		double d1 = points[0].x*nx + points[0].y*ny - no;
		if(d1>=0d) {
			np.add(points[0].copy());
			nc.add(new int[] {colors[0],colors[1]});
		}
		for(int i=1; i<points.length; i++) {
			double d0 = points[i].x*nx + points[i].y*ny - no;
			if(d0>=0d) {
				if(d1<0d) {
					np.add(points[i].fractionTowards(d0/(d0-d1), points[i-1]));
					int col = JPlotMath.colorLerp(colors[2*i], colors[2*i-1], d0/(d0-d1));
					nc.add(new int[] {col,col});
				}
				np.add(points[i].copy());
				nc.add(new int[] {colors[2*i],colors[2*i+1]});
				if(i+1==points.length) {
					res.add(new JDLine(np.toArray(new JDPoint[0]),toArray(nc)));
					np.clear();
					nc.clear();
				}
			} else
			if(d1>=0d) {
				np.add(points[i-1].fractionTowards(d1/(d1-d0), points[i]));
				int col = JPlotMath.colorLerp(colors[2*i-1], colors[2*i], d1/(d1-d0));
				nc.add(new int[] {col,col});
				res.add(new JDLine(np.toArray(new JDPoint[0]),toArray(nc)));
				np.clear();
				nc.clear();
			}
			d1 = d0;
		}
		return res;
	}
	
	private int[] toArray(List<int[]> pairs) {
		int[] res = new int[pairs.size()*2];
		for(int i=0; i<pairs.size(); i++) {
			int[] pair = pairs.get(i);
			res[2*i] = pair[0];
			res[2*i+1] = pair[1];
		}
		return res;
	}
}
