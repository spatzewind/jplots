package jplots.maths;

import java.util.Arrays;

public class JDTriangle {

	public int lev;
	public JDPoint a;
	public JDPoint b;
	public JDPoint c;
	public JDEdge ab;
	public JDEdge bc;
	public JDEdge ca;

	private Integer hash = null;

	public JDTriangle(JDPoint a, JDPoint b, JDPoint c) {
		JDPoint[] tmp = { a, b, c };
		Arrays.sort(tmp);
		this.a = tmp[0];
		this.b = tmp[1];
		this.c = tmp[2];
	}

	public void edges(JDEdge ab, JDEdge bc, JDEdge ca) {
		this.ab = ab.equals(a, b) ? ab : bc.equals(a, b) ? bc : ca;
		this.bc = ab.equals(b, c) ? ab : bc.equals(b, c) ? bc : ca;
		this.ca = ab.equals(c, a) ? ab : bc.equals(c, a) ? bc : ca;
	}
	
	public void contourIntervalLevel(double[] intervals) {
		JDPoint[] tmpP = { a, b, c };
		int[] tmpL = {-1, -1, -1};
		for(int c=0; c<3; c++) {
			if(Double.isNaN(tmpP[c].value))
				continue;
			tmpL[c] = 0;
			for(int l=0; l<intervals.length; l++)
				if(tmpP[c].value>intervals[l])
					tmpL[c]++;
		}
		lev = tmpL[0]+tmpL[1]+tmpL[2];
		if(lev>=0) {
			lev /= 3;
		} else {
			lev = -1;
		}
	}

	@Override
	public String toString() {
		return "t[" + a + " - " + b + " - " + c + "]";
	}

	@Override
	public int hashCode() {
		if (hash != null) {
			return hash;
		}
		return hash = hash(a, b, c);
	}

	public static int hash(JDPoint a, JDPoint b, JDPoint c) {
		final int prime = 31;
		int hash = 1;
		hash = prime * hash + a.hashCode();
		hash = prime * hash + b.hashCode();
		hash = prime * hash + c.hashCode();
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		JDTriangle A = this;
		JDTriangle B = (JDTriangle) obj;

		if (A.a.equals(B.a)) {
			return (A.b.equals(B.b) && A.c.equals(B.c)) || (A.b.equals(B.c) && A.c.equals(B.b));
		} else if (A.a.equals(B.b)) {
			return (A.b.equals(B.a) && A.c.equals(B.c)) || (A.b.equals(B.c) && A.c.equals(B.a));
		} else if (A.a.equals(B.c)) {
			return (A.b.equals(B.a) && A.c.equals(B.b)) || (A.b.equals(B.b) && A.c.equals(B.a));
		}

		return false;
	}
}
