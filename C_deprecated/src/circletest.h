#ifndef circletest_h
#define circletest_h

double sumsquare(double x, double y);

double distance(double x1, double x2);

double magnitude( double x0, double y0, double x1, double y1);

double determinant2(double e11, double e12, double e21, double e22);

double determinant3(double e11, double e12, double e13, 
		    double e21, double e22, double e23,
		    double e31, double e32, double e33);

int circle_test(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

#endif
