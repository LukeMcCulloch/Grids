#ifndef circletest_h
#define circletest_h

double sumsquare(double x, double y);

double distance(double x1, double x2);

double magnitude( double x0, double y0, double x1, double y1);

double determinant2(double e11, double e12, double e21, double e22);

double determinant3(double e11, double e12, double e13, 
		    double e21, double e22, double e23,
		    double e31, double e32, double e33);

bool ccw(double (&A)[2], double (&B)[2], double  (&C)[2]);

int circle_test_tlm(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

int circle_test(double (&a)[2], double (&b)[2], 
                      double (&c)[2], double (&d)[2]);
#endif
