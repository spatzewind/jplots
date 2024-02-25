package jplots;

import processing.core.PConstants;

public interface JPlotConstants extends PConstants {

	/** the version of the Library. */
	public final static String VERSION = "##library.prettyVersion##";
	
	public final static int NONE   =   0;
	
//	public final static int TOP    = 101; //PConstants.TOP;
//	public final static int LEFT   =  37; //PConstants.LEFT;
//	public final static int RIGHT  =  39; //PConstants.RIGHT;
//	public final static int BOTTOM = 102; //PConstants.BOTTOM;
	public final static int BOTH   = 139;
//	public final static int CENTER =   3; //PConstants.CENTER;
	
	public static final int VERTICAL   = 1;
	public static final int HORIZONTAL = 2;
	
	public final static float ROTATE_CLOCKWISE        =  1.5707964f; //PConstants.HALF_PI;
	public final static float ROTATE_COUNTERCLOCKWISE = -1.5707964f; //-PConstants.HALF_PI;
	
	public final static double DEG_TO_RAD = Math.PI / 180d;
	public final static double RAD_TO_DEG = 180d / Math.PI;

	public final static double EARTH_RADIUS_AEQU = 6378137.000d; // m
	public final static double EARTH_RADIUS_MEAN = 6371000.685d; // m
	public final static double EARTH_RADIUS_POL = 6356752.000d; // m
	public final static double EARTH_FLATTENING = 1d / 298.257222101d;
	
}
