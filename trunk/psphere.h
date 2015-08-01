#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_blas.h>
#include <time.h>

// #include <gsl/gsl_heapsort.h>

#define epsilon  0.00001
#define NEGATIVE(x)  ((x) < -epsilon)
#define ISZERO(x)    (((x) >= -epsilon) && ((x) <= epsilon))
#define POSITIVE(x)  ((x) > epsilon)
#define ISNOTZERO(x) (!ISZERO(x))

#define EPS1 0.001
#define EPS2 1.0E-8

#define epsilon2  0.01
#define NEGATIVE2(x)  ((x) < -epsilon2)
#define ISZERO2(x)  (((x) >= -epsilon2) && ((x) <= epsilon2))
#define POSITIVE2(x)  ((x) > epsilon2)

#define sqrt3 1.73205080759

/* Structure for (x,y,z) in 3-D Cartesian system */
struct point {
  double x;
  double y;
  double z;
};

struct Ppoint {
  struct point point;
  double param1;
};

double       solid_angle(struct point center,
                         struct point first, struct point second,struct point third);
struct point point_substract(struct point first, struct point second);
struct point point_add(struct point first, struct point second);
struct point point_multiply(struct point first, double a);

double       scalar_product(struct point first, struct point second);
struct point vector_product(struct point first, struct point second);
double       distance_point_to_plane(struct point top,
                                     struct point first, struct point second,struct point third);

struct point unit_vector(struct point first);
double       distance_to_origin(struct point first);
double       distance(struct point first, struct point second) ;
int          samepoint(struct point first, struct point second) ;
int          within_distance(struct point first, struct point second, double D) ;
int          calculate_fourth(struct point first, struct point second, struct point third,
                              double R1, double R2, double R3,
                              struct point * fourth0,struct point * fourth1,
                              double *r0, double *r1);
int          findmax(int mat[], int length);
int          findmin(int mat[], int length);
int          space_violation(struct point sphere, struct point ALLSHELL[],
                             int all_sphere_count, double d);
void         dump_centers(const char *filename, struct point ALLSHELL[],
                          int current_sphere_count);
int          precheck_triangle(struct point first, struct point second,
                               double distance);
int          isconvex(struct point center,
                      struct point first, struct point second,
                      struct point third, struct point fourth);
void         sort_ascending(double *a, int N);
void         sort_descending(double *a, int N);
void         swap(double *a, int i, int j) ;

void         sort_ascending_point(struct point *p, int N) ;
void         sort_ascending_Ppoint(struct Ppoint *p, int N) ;

double       angle(struct point first, struct point second, double r);
double       cos_angle(struct point first, struct point second);
int          compare_doubles(const double *a, const double *b);
int          compare_twopoints_ascending(const void *a, const void *b);
int          compare_twoPpoints_ascending(const void *a, const void *b);
int          compare_twopoints_descending(const void *a, const void *b);
int          isinside(struct point center,
                      struct point first, struct point second,struct point third, struct point fourth);
int          isinsidecircle(struct point center,
                            struct point first, struct point second,struct point third, struct point fourth);
double       volume(struct point top, struct point first, struct point second,struct point third) ;

int          blocked(struct point nshell[], int N, double r);

int          sbeta_point_sort(struct point *p, int Nc)  ;
int          random_index(int *c, int Nc, gsl_rng *s) ;

double       probks(double alam);
