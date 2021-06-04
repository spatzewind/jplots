import jplots.*;
import jplots.transform.*;

void setup() {
    size(1440,720,P2D);
    
    JPlot plt = new JPlot(this);
    plt.figure(4d,2d);
    plt.setRange(-PI,PI,-HALF_PI,HALF_PI);
    plt.predefImgShow("earth2");
    
    background(255);
    image(plt.show(),0,0,width,height);
    noLoop();
}
