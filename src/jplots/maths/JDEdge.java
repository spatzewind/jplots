package jplots.maths;

public class JDEdge {

	public JDPoint a;
	public JDPoint b;

	public JDEdge(JDPoint a, JDPoint b) {
		boolean swap = 0 < a.compareTo(b);
		this.a = swap ? b : a;
		this.b = swap ? a : b;
	}

	public JDTriangle A;
	public JDTriangle B;

	public JDTriangle[] getWing() {
		if (A != null && B != null) {
			return new JDTriangle[] { A, B };
		}
		return new JDTriangle[] { A };
	}

	public void wing(JDTriangle t) {
		if (this.A == null) {
			this.A = t;
		} else
		if (this.B == null) {
			this.B = t;
		} else {
			System.err.println("[ERR] error state in edge's wing triangle ...");
		}
	}
	
	public boolean has2triangles() {
		return (this.A!=null && this.B!=null);
	}

	@Override
	public String toString() {
		return "e[" + a + " - " + b + "]";
	}

	Integer hash = null;

	@Override
	public int hashCode() {
		if (hash != null) {
			return hash;
		}
		return hash = hash(a, b);
	}

	public static int hash(JDPoint a, JDPoint b) {
		int ahash = a.hashCode();
		int bhash = b.hashCode();
		return ahash * 31 + bhash;
	}

	public boolean equals(JDPoint a, JDPoint b) {
		if (this.a.equals(a) && this.b.equals(b)) {
			return true;
		}
		if (this.a.equals(b) && this.b.equals(a)) {
			return true;
		}
		return false;
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
		JDEdge A = this;
		JDEdge B = (JDEdge) obj;

		return (A.a.equals(B.a) && A.b.equals(B.b)) || (A.a.equals(B.b) && A.b.equals(B.a));
	}
}
