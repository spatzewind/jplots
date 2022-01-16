import jplots.*;
import jplots.maths.JPlotMath.DateTime;

JPlot plt;
float[] temperature, minTemp, maxTemp;
float[] time;
String timeunit;

void setup() {
  size(800,400, P2D);
  smooth();
  
  //load data
  timeunit = "days since 1900-01-01";
  temperature_and_date();
  //PFont font = createFont("",40);
  
  //create plot
  plt = new JPlot(this);
  plt.debug(true);
  plt.figure(6.0,3.0);
  plt.plot(time,minTemp, color(0,0,255), 2f, "-",     "label", "minimum");
  plt.plot(time,maxTemp, color(255,0,0), 2f, "-",     "label", "maximum");
  plt.plot(time,temperature, color(0,255,0), 2f, "-", "label", "average");
  plt.setAsTimeAxis('x', timeunit, "gregorian", "m.yyyy");
  plt.setYRange(-25.0, 60.0);
  plt.setXTitle("Time");
  plt.setYTitle("Temperature");
  plt.setTitle("In Stuttgart Germany");
  //plt.setFont(font);
  plt.legend();
  plt.show();
  plt.redraw(true); //first plotting and a second draw call removes some issues
}

void draw() {
  background(255);
  image(plt.show(),0,0,width,height);
}

void temperature_and_date() {
  //load data from csv file
  String[] data = loadStrings("GM000002716.csv");
  //String[] data = loadStrings("https://www.ncei.noaa.gov/data/daily-summaries/access/GM000002716.csv");
  
  //first line is header gives name of each column
  time = new float[data.length-1];
  minTemp     = new float[data.length-1];
  temperature = new float[data.length-1];
  maxTemp     = new float[data.length-1];
  //but remember date is 
  for(int t=0; t<data.length; t++) {
    //so skip header line,
    if(t==0) continue;
    String entry = ""+data[t];
    //data file contains 16 columns
    String[] column = new String[16];
    int columnend = 0;
    for(int c=0; c<16; c++) {
      entry = entry.substring(columnend);
      if(entry.length()==0) {
        column[c] = "";
      } else
      if(entry.charAt(0)=='"') {
        columnend = entry.substring(1).indexOf("\"") + 1;
        column[c] = entry.substring(1,columnend).trim();
        columnend += 2;
      } else {
        columnend = entry.indexOf(",");
        if(columnend==0) {
          column[c] = "";
        } else {
          column[c] = entry.substring(0,columnend+1).trim();
        }
        columnend += 1;
      }
    }
    //date is 2nd column
    DateTime date = new DateTime(column[1], null);
    time[t-1] = (float) date.toDouble(timeunit, null);
    //minimal temperature is 11th column
    if(column[10].length()>0) minTemp[t-1] = 0.1*float(column[10].trim());
    else minTemp[t-1] = Float.NaN;
    //maximal temperature is 13th column
    if(column[12].length()>0) maxTemp[t-1] = 0.1*float(column[12].trim());
    else maxTemp[t-1] = Float.NaN;
    //temperature is 15th column
    if(column[14].length()>0) temperature[t-1] = 0.1*float(column[14].trim());
    else temperature[t-1] = Float.NaN;
  }
  //monthly average mean:
  
}
