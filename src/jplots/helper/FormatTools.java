package jplots.helper;

import java.io.PrintStream;

public class FormatTools {
	
	public static void printMat(PrintStream stream, double[][] _mat) {
	    int la = _mat.length, lb = _mat[0].length;
	    String[][] mout = new String[la][lb];
	    //stream.println("    [PRINTMAT] "+la+"x"+lb+" matrix");
	    int ml = 0;
	    for(int a=0; a<la; a++) for(int b=0; b<lb; b++) {
	        mout[a][b] = " "+_mat[a][b]+" ";
	        ml = Math.max(ml,mout[a][b].length());
	    }
	    //stream.println("    [PRINTMAT] "+ml+"chars per entry");
	    for(int a=0; a<la; a++) for(int b=0; b<lb; b++) {
	        while(mout[a][b].length()<ml) { mout[a][b] = " "+mout[a][b]; }
	    }
	    for(int a=0; a<la; a++) {
	        String l = "";
	        for(int b=0; b<lb; b++) l += mout[a][b];
	        if(a==0) {
	        	if(a==la-1) {
	        		l = "("+l+")";
	        	} else {
	        		l = "/"+l+"\\";
	        	}
	        } else {
	        	if(a==la-1) {
	        		l = "\\"+l+"/";
	        	} else {
	        		l = "|"+l+"|";
	        	}
	        }
	        stream.println(l);
	    }
	}
	public static void printVec(PrintStream stream, double[] _vec, boolean transposed) {
	    int la = _vec.length, lb = 1;
	    String[] mout = new String[la];
	    //stream.println("    [PRINTMAT] "+la+"x"+lb+" matrix");
	    int ml = 0;
	    for(int a=0; a<la; a++) {
	        mout[a] = " "+_vec[a]+" ";
	        ml = Math.max(ml,mout[a].length());
	    }
	    //stream.println("    [PRINTMAT] "+ml+"chars per entry");
	    for(int a=0; a<la; a++) {
	        while(mout[a].length()<ml) { mout[a] = " "+mout[a]; }
	    }
	    if(transposed) { lb = la; la = 1; }
	    for(int a=0; a<la; a++) {
	        String l = "";
	        for(int b=0; b<lb; b++) { int i = a+b; l += mout[i]; }
	        if(a==0) {
	        	if(a==la-1) {
	        		l = "("+l+")"+(transposed?"T":"");
	        	} else {
	        		l = "/"+l+"\\"+(transposed?"T":"");
	        	}
	        } else {
	        	if(a==la-1) {
	        		l = "\\"+l+"/";
	        	} else {
	        		l = "|"+l+"|";
	        	}
	        }
	        stream.println(l);
	    }
	}
}
