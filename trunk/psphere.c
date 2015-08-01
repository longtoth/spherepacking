/*------------------------------------------------------------
 *  Source: 
 *          author =  {Long To and Zbigniew Stachurski},
 *          title =   {Random close packing of spheres in a round cell},
 *          journal = {Journal of Non-Crystalline Solids},
 *          year =    2004,
 *          volume =  333,
 *          pages =   {161-171}
 *-------------------------------------------------------------
 * DESCRIPTION:
 *  Dense packing of spheres inside a spherical cell
 * HISTORY: 
 * Last Update: 2006-07-11
 *-------------------------------------------------------------
 */


#include "psphere.h"


//--------------------------------------------------------------
//  solid_angle
// -------------------------------------------------------------
//  Desc: Calculate the solid angle form by three points on the
//        surface of a sphere
//
//  History:
//--------------------------------------------------------------
double solid_angle(struct point center,
                   struct point first, struct point second, struct point third)
{
  struct point a, b, c;
  double upper, lower, omega;
  a = unit_vector(point_substract(first,center))  ;
  b = unit_vector(point_substract(second,center)) ;
  c = unit_vector(point_substract(third,center))  ;
  
  if (ISZERO(distance_to_origin(vector_product(b,c))) ||
      ISZERO(distance_to_origin(vector_product(c,a))) ||
      ISZERO(distance_to_origin(vector_product(a,b))))
    omega = 0;
  else
    {
      upper = scalar_product(a,vector_product(b,c));
      lower = (1 + scalar_product(a,b) + scalar_product(b,c) + scalar_product(c,a));
      omega = 2*atan(fabs(upper)/lower);
      if (omega<0)
        omega = omega + 2*M_PI;
    }
  return omega;
}

//--------------------------------------------------------------
//  vector_product
// -------------------------------------------------------------
//  Desc: Calculate the vector product of two vectors to points
//        1 and 2
//
//  History:
//--------------------------------------------------------------
struct point vector_product(struct point first, struct point second)
{
  struct point out ;
  out.x = first.y*second.z - first.z*second.y;
  out.y = first.z*second.x - first.x*second.z;
  out.z = first.x*second.y - first.y*second.x;
  return out;
}

//--------------------------------------------------------------
//  scalar_product
// -------------------------------------------------------------
//  Desc: Calculate the scalar product of two vectors
//
//  History:
//--------------------------------------------------------------
double scalar_product(struct point first, struct point second)
{
  return (first.x*second.x + first.y*second.y + first.z*second.z);
}

//--------------------------------------------------------------
//  point_substract
// -------------------------------------------------------------
//  Desc: Subtract one vector from the other
//
//  History:
//--------------------------------------------------------------
struct point point_substract(struct point first, struct point second)
{
  struct point out ;
  out.x = first.x - second.x;
  out.y = first.y - second.y;
  out.z = first.z - second.z;
  return out;
}

//--------------------------------------------------------------
//   point_add
// -------------------------------------------------------------
//  Desc: Add two vectors
//
//  History:
//--------------------------------------------------------------
struct point point_add(struct point first, struct point second)
{
  struct point out ;
  out.x = first.x + second.x;
  out.y = first.y + second.y;
  out.z = first.z + second.z;
  return out;
}

//--------------------------------------------------------------
//   point_multiply
// -------------------------------------------------------------
//  Desc: Multiply a vector by a constant
//
//  History:
//--------------------------------------------------------------
struct point point_multiply(struct point first, double a)
{
  struct point out ;
  out.x = a*first.x;
  out.y = a*first.y;
  out.z = a*first.z;
  return out;
}

//--------------------------------------------------------------
//  distance_point_to_plane
// -------------------------------------------------------------
//  Desc: Calculate the distance of a point to a plane defined
//        by three other points
//  History:
//--------------------------------------------------------------
double distance_point_to_plane(struct point top,
                               struct point first, struct point second,struct point third)
{
  struct point A;
  A = vector_product(point_substract(second, first), point_substract(third, first));
  return fabs(scalar_product(A, point_substract(first, top))/distance_to_origin(A));
}
  
//--------------------------------------------------------------
//  isconvex
// -------------------------------------------------------------
//  Desc: Check if the polyhedron formed by 4 points contains the
//        center point
//  History:
//--------------------------------------------------------------
int isconvex(struct point center,
             struct point first, struct point second,
             struct point third, struct point fourth)
{
  if (ISZERO(distance_point_to_plane(center, first, second, third))  |
      ISZERO(distance_point_to_plane(center, first, second, fourth)) |
      ISZERO(distance_point_to_plane(center, first, third, fourth))  |
      ISZERO(distance_point_to_plane(center, second, third, fourth)) |
      (NEGATIVE((solid_angle(center, first, second, third) +
                 solid_angle(center, first, second, fourth) +
                 solid_angle(center, first, third, fourth) +
                 solid_angle(center, second, third, fourth)) - 4*M_PI)))
    return 0;
  else
    return 1;
}

//--------------------------------------------------------------
//  distance_to_origin
// -------------------------------------------------------------
//  Desc: Calculate the distance to the origin
//
//  History:
//--------------------------------------------------------------
double distance_to_origin(struct point first)
{
  return sqrt(first.x*first.x + first.y*first.y + first.z*first.z);
}

//--------------------------------------------------------------
//  unit_vector
// -------------------------------------------------------------
//  Desc: Calculate the coordinates of the unit vector along a
//        given direction
//  History:
//--------------------------------------------------------------
struct point unit_vector(struct point first)
{
  struct point out ;
  double d = distance_to_origin(first) ;
  out.x = first.x/d;
  out.y = first.y/d;
  out.z = first.z/d;
  return out;
}

//--------------------------------------------------------------
//  distance
// -------------------------------------------------------------
//  Desc: Calculate the distance between two points
//
//  History:
//--------------------------------------------------------------
double distance(struct point first, struct point second)
{
  double x1, y1, z1, x2, y2, z2;  
  x1 = first.x;
  y1 = first.y;  
  z1 = first.z;
  x2 = second.x;
  y2 = second.y;     
  z2 = second.z;
  return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}

//--------------------------------------------------------------
//  samepoint
// -------------------------------------------------------------
//  Desc: Check if two points co-incides
//
//  History:
//--------------------------------------------------------------
int samepoint(struct point first, struct point second)
{
  return (ISZERO(first.x - second.x) && ISZERO(first.y - second.y) &&
          ISZERO(first.z - second.z)) ;
}

//-------------------------------------------------------------
// within_distance
//-------------------------------------------------------------
// Desc: Check if two points are within the limit of each other
//
// History:
//------------------------------------------------------------
int within_distance(struct point first, struct point second, double D)
{
  if (!NEGATIVE2(D - distance(first, second)))
    return 1;
  else
    return 0;
}


//--------------------------------------------------------------
//  calculate_fourth
// -------------------------------------------------------------
//  Desc: Calculate the position of the fourth sphere which
//        touches three existing spheres
//
//  History:
//--------------------------------------------------------------
int calculate_fourth(struct point first, struct point second, struct point third,
                     double R1, double R2, double R3,
                     struct point * fourth0,struct point * fourth1,
                     double *r0, double *r1)
{
  int n = 0;
  double x[2], y[2], z[2];
  double a1, b1, c1, a2, b2, c2, d1, d2;
  double y_y1_squared, z_z1_squared;
  double x1, y1, z1, x2, y2, z2,x3, y3, z3;
  double A, B, C;
  double my, ny, mz, nz;
  int rn;
  double k;

  double a,b,c,p,Area,R_out;
  a = distance(second,third);
  b = distance(first,third);
  c = distance(first,second);
  p = 0.5*(a + b + c);
  Area = sqrt(p*(p-a)*(p-b)*(p-c));
  R_out = (a*b*c)/(4*Area);
  
  if (ISZERO(R1-R2) && ISZERO(R2-R3) && ISZERO(R1-R3))
    {
      if (NEGATIVE2(R1 - R_out))
        return 0;
      else if (R1 < R_out)
        {
          R1 = R2 = R3 = R_out+epsilon;
        }
    }
  
  x1 = first.x;
  y1 = first.y;  
  z1 = first.z;
  x2 = second.x;
  y2 = second.y;     
  z2 = second.z;
  x3 = third.x;
  y3 = third.y;  
  z3 = third.z;
  
  a1 = 2*(first.x - second.x);
  b1 = 2*(first.y - second.y);
  c1 = 2*(first.z - second.z);

  a2 = 2*(first.x - third.x);
  b2 = 2*(first.y - third.y);
  c2 = 2*(first.z - third.z);

  d1 = R2*R2 - R1*R1 + ((x1*x1 + y1*y1 + z1*z1) - (x2*x2 + y2*y2 + z2*z2));
  d2 = R3*R3 - R1*R1 + ((x1*x1 + y1*y1 + z1*z1) - (x3*x3 + y3*y3 + z3*z3));

  if ( ISZERO(b2) && ISZERO(c2) && ISZERO(a2) )
    n = 0;
  else if (ISZERO(b2) && ISZERO(c2) && ISNOTZERO(a2))
    {
      if ( ISZERO(b1) && ISZERO(c1) )
        n = 0;
      else if ( ISZERO(b1) && ISNOTZERO(c1) )
        {
          x[0] = (double) d2/a2;
          z[0] = (double) (d1 - a1*x[0])/c1;
          y_y1_squared = R1*R1 - (x[0]-x1)*(x[0]-x1) - (z[0]-z1)*(z[0]-z1);
          if (NEGATIVE(y_y1_squared))
            n = 0;
          else if (ISZERO(y_y1_squared))
            {
              n = 1;
              y[0] = y1; }
          else {
            n = 2;
            y[0] = y1 + sqrt(y_y1_squared);
            y[1] = y1 - sqrt(y_y1_squared);
            x[1] = x[0];
            z[1] = z[0];
          }
        }
      else {
        x[0] = (double) d2/a2;
        A = ((double) -c1/b1)*((double) -c1/b1)+1;
        B = (double) 2*(c1/b1)*((d2*a1/a2-d1)/b1+y1) - 2*z1;
        C = ((d2*a1/a2-d1)/b1+y1)*((d2*a1/a2-d1)/b1+y1) + z1*z1 - R1*R1 + (x[0]-x1)*(x[0]-x1);
        rn = gsl_poly_solve_quadratic(A,B,C,r0,r1);
        if (rn == 0)
          n = 0;
        else {
          n = 2;
          z[0] = *r0;
          z[1] = *r1;
          y[0] = (1/b1)*(d1 - d2*a1/a2 - c1*z[0]);
          y[1] = (1/b1)*(d1 - d2*a1/a2 - c1*z[1]);
          x[1] = x[0];
        }
      }
    }
  else if ( ISZERO(b2) && ISNOTZERO(c2) && ISZERO(b1) )
    {
      if ( ISZERO(c1) && ISZERO(a1) )
        n = 0;
      else if ( ISZERO(c1) && ISNOTZERO(a1) )
        {
          x[0] = d1/a1;
          z[0] = (d2-a2*x[0])/c2;
          y_y1_squared = R1*R1 - (x[0]-x1)*(x[0]-x1) - (z[0]-z1)*(z[0]-z1);
          if (NEGATIVE(y_y1_squared))
            n = 0;
          else if (ISZERO(y_y1_squared))
            {
              n = 1;
              y[0] = y1; }
          else {
            n = 2;
            y[0] = y1 + sqrt(y_y1_squared);
            y[1] = y1 - sqrt(y_y1_squared);
            x[1] = x[0];
            z[1] = z[0];        
          }
        }
      else {
        if (ISZERO (c2*a1 - c1*a2))
          n = 0;
        else {
          x[0] = (c2*d1 - c1*d2)/(c2*a1 - c1*a2);
          z[0] = (d2-a2*x[0])/c2;
          y_y1_squared = R1*R1 - (x[0]-x1)*(x[0]-x1) - (z[0]-z1)*(z[0]-z1);
          if (NEGATIVE(y_y1_squared))
            n = 0;
          else if (ISZERO(y_y1_squared))
            {
              n = 1;
              y[0] = y1; }
          else {
            n = 2;
            y[0] = y1 + sqrt(y_y1_squared);
            y[1] = y1 - sqrt(y_y1_squared);
            x[1] = x[0];
            z[1] = z[0];        
          }
        }
      }
    }
  else if ( ISNOTZERO(b2) && ISZERO(c2) && ISZERO(c1) )
    {
      if (ISZERO(b1) && ISZERO(a1) )
        n = 0;
      else if ( ISZERO(b1) && ISNOTZERO(a1) )
        {
          x[0] = d1/a1;
          y[0] = (d2-a2*x[0])/b2;
          z_z1_squared = R1*R1 - (x[0]-x1)*(x[0]-x1) - (y[0]-y1)*(y[0]-y1);
          if (NEGATIVE(z_z1_squared))
            n = 0;
          else if (ISZERO(z_z1_squared))
            {
              n = 1;
              z[0] = z1;
            }
          else {
            n = 2;
            z[0] = z1 + sqrt(z_z1_squared);
            z[1] = z1 - sqrt(z_z1_squared);
            x[1] = x[0];
            y[1] = y[0];        
          }
        }
      else {
        if (ISZERO(b2*a1 - b1*a2))
          n = 0;
        else {
          x[0] = (b2*d1-b1*d2)/(b2*a1-b1*a2) ;
          y[0] = (d2 - a2*x[0])/b2;
          z_z1_squared = R1*R1 - (x[0]-x1)*(x[0]-x1) - (y[0]-y1)*(y[0]-y1);
          if (NEGATIVE(z_z1_squared))
            n = 0;
          else if (ISZERO(z_z1_squared))
            {
              n = 1;
              z[0] = z1;
            }
          else {
            n = 2;
            z[0] = z1 + sqrt(z_z1_squared);
            z[1] = z1 - sqrt(z_z1_squared);
            x[1] = x[0];
            y[1] = y[0];    
          }
        }
      }
    }
  else if (ISZERO(b1*c2 - b2*c1))
    {
      k = b1/b2;
      if (ISZERO(k*a2 - a1))
        n = 0;
      else {
        x[0] = (k*d2-d1)/(k*a2-a1);
        A = (c2/b2)*(c2/b2) + 1;
        B = 2*(-c2/b2)*((d2-a2*x[0])/b2-y1) - 2*z1;
        C = ((d2-a2*x[0])/b2-y1)*((d2-a2*x[0])/b2-y1) + z1*z1 - R1*R1 + (x[0]-x1)*(x[0]-x1);
        rn = gsl_poly_solve_quadratic(A,B,C,r0,r1);
        if (rn == 0)
          n = 0;
        else {
          n = 2;
          x[1] = x[0];
          z[0] = *r0;
          z[1] = *r1;
          y[0] = (1/b2)*(d2 - a2*x[0] - c2*z[0]);
          y[1] = (1/b2)*(d2 - a2*x[1] - c2*z[1]);
        }
      }
    }
  else {
    my = (c2*d1-c1*d2)/(b1*c2-b2*c1);
    ny = (c1*a2-c2*a1)/(b1*c2-b2*c1);
    mz = (b2*d1-b1*d2)/(c1*b2-c2*b1);
    nz = (b1*a2-b2*a1)/(c1*b2-c2*b1);
    A = 1 + ny*ny + nz*nz;
    B = -2*x1 + 2*ny*(my-y1) + 2*nz*(mz-z1);
    C = x1*x1 + (my-y1)*(my-y1) + (mz-z1)*(mz-z1) - R1*R1;
    rn = gsl_poly_solve_quadratic(A,B,C,r0,r1);
    if (rn == 0)
      n = 0;
    else {
      n = 2;      
      x[0] = *r0;
      x[1] = *r1;
      y[0] = my + ny*x[0];
      y[1] = my + ny*x[1];
      z[0] = mz + nz*x[0];
      z[1] = mz + nz*x[1];
    }
  }
  if (n == 1)
    {
      fourth0->x = x[0];
      fourth0->y = y[0];
      fourth0->z = z[0];
    }
  else if (n == 2)
    {
      if ((x[0]*x[0] + y[0]*y[0] + z[0]*z[0]) > (x[1]*x[1] + y[1]*y[1] + z[1]*z[1]))
        {
          fourth0->x = x[0];
          fourth0->y = y[0];
          fourth0->z = z[0];
          fourth1->x = x[1];
          fourth1->y = y[1];
          fourth1->z = z[1];
        }
      else {
        fourth0->x = x[1];
        fourth0->y = y[1];
        fourth0->z = z[1];
        fourth1->x = x[0];
        fourth1->y = y[0];
        fourth1->z = z[0];
      }
    }
  else
    n = 0;
  return n;
}

//--------------------------------------------------------------
//  findmax
// -------------------------------------------------------------
//  Desc: Find the maximum element of a matrix
//
//  History:
//--------------------------------------------------------------
int findmax(int mat[], int length)
{
  int max = mat[0];
  int maxidx = 0;
  int i;
  for (i = 1; i<length; i++)
    {
      if (mat[i] > max)
        {
          max = mat[i];
          maxidx = i;
        }
    }
  return maxidx;
}

//--------------------------------------------------------------
//  findmin
// -------------------------------------------------------------
//  Desc: Find the minimum element of a matrix
//
//  History:
//--------------------------------------------------------------
int findmin(int mat[], int length)
{
  int min = mat[0];
  int minidx = 0;
  int i;
  for (i = 1; i<length; i++)
    {
      if (mat[i] < min)
        {
          min = mat[i];
          minidx = i;
        }
    }
  return minidx;
}

//--------------------------------------------------------------
//  space_violation
// -------------------------------------------------------------
//  Desc: Find out if a sphere collides with the rest of the
//        structure
//
//  History:
//--------------------------------------------------------------
int space_violation(struct point sphere, struct point *ALLSHELL,
                    int all_sphere_count, double d)
{
  int i = 0;

  for (i = all_sphere_count-1; i >=0; i--)
    {
      if (precheck_triangle(sphere,ALLSHELL[i], d) &&
          NEGATIVE2(distance(sphere,ALLSHELL[i]) - d))
        {
          return 1;
        }
    }
  return 0 ;
}

//--------------------------------------------------------------
//  precheck_triangle
// -------------------------------------------------------------
//  Desc: Check (preliminary) if two spheres are close enough
//        to have a third sphere touching both of them
//  History:
//--------------------------------------------------------------
int precheck_triangle(struct point first, struct point second, double distance)
{
  if (POSITIVE2(fabs(first.x - second.x) - distance) ||
      POSITIVE2(fabs(first.y - second.y) - distance) ||
      POSITIVE2(fabs(first.z - second.z) - distance) ||
      POSITIVE2(fabs(first.x - second.x) + fabs(first.y - second.y) + fabs(first.z - second.z) - distance*sqrt3))
    return 0;
  else
    return 1;
}

//--------------------------------------------------------------
//  dump_centers
// -------------------------------------------------------------
//  Desc: Print coordinates of the spheres
//
//  History:
//--------------------------------------------------------------
void dump_centers(const char *filename, struct point ALLSHELL[], int current_sphere_count)
{
  int i;
  FILE *f;

  if ((f = fopen(filename, "w")) == (FILE *)NULL)
    abort();
  for (i = 0; i < current_sphere_count; ++i)
    fprintf(f, "%g %g %g\n", ALLSHELL[i].x, ALLSHELL[i].y, ALLSHELL[i].z);
  fclose(f);
}


//--------------------------------------------------------------
//  sort_ascending
// -------------------------------------------------------------
//  Desc: Sort an array in ascending order
//
//  History:
//--------------------------------------------------------------
void sort_ascending(double *a, int N)
{
  int i, j;
  double temp;
  
  for (i = 1; i<N; i++)
    {
      if (a[i] < a[0])
        {
          temp = a[i];
          memcpy(a+1,a,i*sizeof(double));
          a[0] = temp;
        }
      else if (a[i] < a[i-1])
        {
          j = 0;
          while (a[i] >= a[j])
            j++;

          temp = a[i];
          memcpy(a+j+1,a+j,(i-j)*sizeof(double));
          a[j] = temp;
        }
    }
}

//--------------------------------------------------------------
//  sort_ascending_point
// -------------------------------------------------------------
//  Desc: Sort an array of points according to distance to the
//        origin (ascending)
//
//  History:
//--------------------------------------------------------------
void sort_ascending_point(struct point *p, int N)
{
  qsort(p,N,sizeof(struct point),compare_twopoints_ascending);
}

//--------------------------------------------------------------
//  sort_ascending_Ppoint
// -------------------------------------------------------------
//  Desc: Sort an array of Ppoints according to the fourth parameter
//
//  History:
//--------------------------------------------------------------
void sort_ascending_Ppoint(struct Ppoint *p, int N)
{
  qsort(p,N,sizeof(struct Ppoint),compare_twoPpoints_ascending);
}

//--------------------------------------------------------------
//  sort_descending_point
// -------------------------------------------------------------
//  Desc: Sort an array of points according to distance to the
//        origin (decending)
//
//  History:
//--------------------------------------------------------------
void sort_descending_point(struct point *p, int N)
{
  qsort(p,N,sizeof(struct point),compare_twopoints_descending);
}

//--------------------------------------------------------------
//  swap
// -------------------------------------------------------------
//  Desc: Swap two elements of an array
//
//  History:
//--------------------------------------------------------------
void swap(double *a, int i, int j)
{
  double temp;
  if (i == j)
    return;
  else {
    temp = a[i];
    a[i] = a[j];
    a[j] = temp;
  }
}

//--------------------------------------------------------------
//   angle
// -------------------------------------------------------------
//  Desc: Calculate the angle form by two points on a sphere
//
//  History:
//--------------------------------------------------------------
double angle(struct point first, struct point second, double r)
{
  return acos((first.x*second.x + first.y*second.y + first.z*second.z)/(r*r));
}

//--------------------------------------------------------------
//  cos_angle
// -------------------------------------------------------------
//  Desc: Calculate the cosine of an angle between two points
//
//  History:
//--------------------------------------------------------------
double cos_angle(struct point first, struct point second)
{
  return scalar_product(first, second)/(distance_to_origin(first)*distance_to_origin(second));
}


//--------------------------------------------------------------
//  compare_doubles
// -------------------------------------------------------------
//  Desc: Compare two doubles
//
//  History:
//--------------------------------------------------------------
int compare_doubles(const double *a, const double *b)
{
  if (*a > *b)
    return 1;
  else if (*a < *b)
    return -1;
  else
    return 0;
}

//--------------------------------------------------------------
//  compare_twopoints_ascending
// -------------------------------------------------------------
//  Desc: Compare two points - for use in sort functions
//
//  History:
//--------------------------------------------------------------
int compare_twopoints_ascending(const void *a, const void *b)
{
  double x = distance_to_origin(*(struct point *)a);
  double y = distance_to_origin(*(struct point *)b);  
  if (x > y)
    return 1;
  else if (x < y)
    return -1;
  else
    return 0;
}

//--------------------------------------------------------------
//  compare_twoPpoints_ascending
// -------------------------------------------------------------
//  Desc: Compare two Ppoints - for use in sort functions
//
//  History:
//--------------------------------------------------------------
int compare_twoPpoints_ascending(const void *a, const void *b)
{
  double x = ((struct Ppoint *)a)->param1;
  double y = ((struct Ppoint *)b)->param1;
  if (x > y)
    return 1;
  else if (x < y)
    return -1;
  else
    return 0;
}


//--------------------------------------------------------------
//  compare_twopoints_descending
// -------------------------------------------------------------
//  Desc: Compare two points - for use in sort functions
//
//  History:
//--------------------------------------------------------------
int compare_twopoints_descending(const void *a, const void *b)
{
  return -compare_twopoints_ascending((struct point *)a, (struct point *)b) ;
}

//--------------------------------------------------------------
//  isinside
// -------------------------------------------------------------
//  Desc: Another implementation of isconvex
//
//  History:
//--------------------------------------------------------------
int isinside(struct point center,
             struct point first, struct point second,struct point third, struct point fourth)
{

  struct point a, b, c, d;
  a = point_substract(first,center);
  b = point_substract(second,center);
  c = point_substract(third,center);
  d = point_substract(fourth,center);

  if (ISZERO(solid_angle(center, first, second, third)))
    return 1;
  else if (ISZERO(distance_to_origin(vector_product(d,a))) ||
           ISZERO(distance_to_origin(vector_product(d,b))) ||
           ISZERO(distance_to_origin(vector_product(d,c))))
    return 0;
  else if (ISZERO(solid_angle(center, first, second, third) - solid_angle(center, first, second, fourth) -
                  solid_angle(center, first, third, fourth) - solid_angle(center, second, third, fourth)))
    return 1;
  else
    return 0;
}

//--------------------------------------------------------------
//  spherical_harmonic
// -------------------------------------------------------------
//  Desc: Calculate the spherical harmonic
//
//  History:
//--------------------------------------------------------------
void spherical_harmonic(int l, int m, double theta, double phi, double *Yc, double *Ys)
{
  if (m >= 0)
    {
      *Yc = gsl_sf_legendre_sphPlm(l,m,cos(theta))*cos(m*phi) ;
      *Ys = gsl_sf_legendre_sphPlm(l,m,cos(theta))*sin(m*phi) ;
    }
  else // spherical harmonic for negative m
    {
      m = -m;
      *Yc = gsl_sf_legendre_sphPlm(l,m,cos(theta))*cos(m*phi) ;
      *Ys = gsl_sf_legendre_sphPlm(l,m,cos(theta))*sin(m*phi) ;
      *Yc = pow(-1,m)*(*Yc);
      *Ys = -pow(-1,m)*(*Ys);
    }
}

//--------------------------------------------------------------
//  volume
// -------------------------------------------------------------
//  Desc: Calculate the volume of a polyhedron form by 4 points
//
//  History:
//--------------------------------------------------------------
double volume(struct point top, struct point first, struct point second,struct point third)
{
  double a,b,c,h,p,Area;
  a = distance(second,third);
  b = distance(first,third);
  c = distance(first,second);
  p = 0.5*(a + b + c);
  Area = sqrt(p*(p-a)*(p-b)*(p-c));
  h = distance_point_to_plane(top, first, second, third);

  return (h*Area);
}

//--------------------------------------------------------------
//  spherical_coordinate
// -------------------------------------------------------------
//  Desc: Convert from Cartesian to spherical coordinates
//
//  History:
//--------------------------------------------------------------
void spherical_coordinate(struct point p, double *theta, double *phi)
{
  double r;
  r = distance_to_origin(p);
  *theta = acos(p.z/r);      // polar coordinate

  if (p.x==0)
    if (p.y==0)
      *phi = 0;
    else if (p.y > 0)
      *phi = M_PI/2;
    else
      *phi = 3*M_PI/2;
  else
    *phi   = atan(p.y/p.x);   // azimuthal coordinate

  if ((p.x<0) && (p.y<=0))
    *phi = *phi + M_PI ;
  else if ((p.x<0) && (p.y>0))
    *phi = *phi + M_PI ;  
  else if ((p.x>0) && (p.y<0))
    *phi = *phi + 2*M_PI ;  
}


//--------------------------------------------------------------
//  cartesian_to_spherical_coordinate
// -------------------------------------------------------------
//  Desc: Convert from spherical coordinates to Cartesian
//
//  History:
//--------------------------------------------------------------
void cartesian_to_spherical_coordinate(double theta, double phi, struct point * p)
{
  p->z = cos(theta);
  p->x = sin(theta)*cos(phi);
  p->y = sin(theta)*sin(phi);
}


//--------------------------------------------------------------
//   MAIN PROGRAM
// -------------------------------------------------------------
//  Desc: Pack a large number of spheres into a round cell.
//        The maximum number of spheres is set by N_spheres_limit
//
//  History:
//--------------------------------------------------------------
int main(int argc, char *argv[])
{
  char *allcenter =   "allcenter" ;
  char *alltouch  =   "alltouch" ;
  char *otherinfo =   "otherinfo";
  char *oldcenter =   "seed" ;
  char *buffered_center =  "buffered_center" ;

  FILE *file_center, *file_other, *file_touch;
  FILE *file_center_old;
  FILE *file_center_buffer;

  double r = 0.5 ;
  int N_shells =  100;            // number of shells
  int N_old_spheres = 200000;
  int N_spheres_limit = 120000;
  int N_spheres = (N_spheres_limit*7)/5;
  int TMPIDX = N_spheres_limit/3;
  int N_shells_limit = 100;

  struct point zero = {0, 0, 0};
  struct point first, second, third;
  struct point * fourth0 = (struct point *)malloc(sizeof(struct point));
  struct point * fourth1 = (struct point *)malloc(sizeof(struct point));

  int i, j, rn, cnt = 0, j_n, exist_flag;
  double t;
  int flip[6];
  
  int S, maxidx;
  int current_sphere_count = 0, N, M;
  int last_sphere_count ;
  
  int check_range;
  int *tooshort,  *mark ;
  struct point *SHELL_temp;
  struct Ppoint *SHELL_Ptemp;
  
  struct point *OLDSHELL;
  struct point *ALLSHELL;  // coordinates of all spheres in the cells, i.e all shells
  int **ALLTOUCH;          // touching spheres

  double gap ;
  
  int SPCOUNT[N_shells] ;      // number of spheres in each shells
  int SPCOUNT_pass1[N_shells];
  int SPCOUNT_pass2[N_shells];
  int SPCOUNT_pass3[N_shells];
  int SPCOUNT_pass1a[N_shells];
  int SPCOUNT_pass2a[N_shells];
  int SPCOUNT_pass3a[N_shells];
  int L = 0;
  
  int N_next_points = 0, removed_pts = 0;
  double max_radius, min_radius;
  int min_idx,max_idx;

  double *r0 = (double *)malloc(sizeof(double));
  double *r1 = (double *)malloc(sizeof(double));
  
  int m, n , k;

  double D = 2*r;
  double D_as = D*sqrt((double)2/3)*pow((0.74/0.64),((double)1/3));
  double shell_radius[N_shells] ;

  double misc_shellinfo[N_shells][8];

  gsl_rng *s ;
  const gsl_rng_type * Ts;
  int dt[6][3];
  double Dt[3];
  int passno;

  double radius[N_old_spheres], core_radius[N_old_spheres];
  double shift_x, shift_y, shift_z,  shell_widthdiff, core_sample_radius, shift_range;
  float x, y, z;

  double d1, d2;
  
  if ((file_center_old = fopen(oldcenter, "r")) == (FILE *)NULL)
    {
      printf("Initial center file missing!\n");
      abort();
    }

  // pointer initialization
  ALLSHELL = (struct point *)malloc(N_spheres * sizeof(struct point)) ;
  OLDSHELL = (struct point *)malloc(N_old_spheres * sizeof(struct point)) ;
  ALLTOUCH = (int **)malloc(N_spheres * sizeof(int *));
  for(i = 0; i < N_spheres; i++)
    ALLTOUCH[i] = (int *)malloc(13 * sizeof(int));
  tooshort = (int *)malloc(TMPIDX * sizeof(int));
  mark = (int *)malloc(TMPIDX * sizeof(int));

  // interpreting input argurments
  if (argc != 8)
    {
       printf("Not enough parameters!\n");
       printf("Sample command line: psphere 000 000 000 000 000 000 111111\n");
       return 1;
    }
       
  for (i=0; i<18; i++)
    {
      dt[i/3][i%3] = (atoi(argv[(i/3)+1])/pow(10,2-i%3));
      dt[i/3][i%3] = (dt[i/3][i%3] % 10);
      printf("%2d %2d %2d \n", i/3, i%3, dt[i/3][i%3]);
    }

  for (i=0; i<6; i++)
    {
      flip[i] = (atoi(argv[7])/pow(10,5-i));
      flip[i] = flip[i] % 10;
      printf("%2d ", flip[i]);
    }
  printf("\n");


  if ((file_center = fopen(allcenter, "w")) == (FILE *)NULL)
    abort();
  if ((file_center_buffer = fopen(buffered_center, "w")) == (FILE *)NULL)
    abort();
  if ((file_other = fopen(otherinfo, "w")) == (FILE *)NULL)
    abort();
  if ((file_touch = fopen(alltouch, "w")) == (FILE *)NULL)
    abort();

  core_sample_radius = 5;
  max_radius = 0;
  i = 0;
  while (fscanf(file_center_old, "%f %f %f", &x, &y, &z) != EOF)
    {
      OLDSHELL[i].x = x;
      OLDSHELL[i].y = y;
      OLDSHELL[i].z = z;
      radius[i] = distance_to_origin(OLDSHELL[i]);
      if (radius[i] > max_radius)
        max_radius = radius[i];
      i++ ;
    }
  fclose(file_center_old);
  N = i;
  printf("Number of spheres: %d\n",N);
  
  shell_widthdiff = 2;

  shift_range = max_radius - shell_widthdiff - core_sample_radius;
  
  gsl_rng_env_setup();
  Ts = gsl_rng_default;
  s = gsl_rng_alloc(Ts);
  for (i = 0; i < 100; i++)
    t = gsl_rng_uniform_pos(s) ;
  shift_x = shift_range*gsl_rng_uniform_pos(s) ;
  shift_y = shift_range*gsl_rng_uniform_pos(s) ;
  shift_z = shift_range*gsl_rng_uniform_pos(s) ;
  printf("Shifts: %g %g %g \n", shift_x, shift_y, shift_z);
  
  j = 0;
  for (i = 0; i < N; i++)
    {
      OLDSHELL[i].x = OLDSHELL[i].x - shift_x;
      OLDSHELL[i].y = OLDSHELL[i].y - shift_y;
      OLDSHELL[i].z = OLDSHELL[i].z -   shift_z;
      radius[i] = distance_to_origin(OLDSHELL[i]);
      if (radius[i] < core_sample_radius)
        {
          ALLSHELL[j] = OLDSHELL[i] ;
          core_radius[j] = radius[i];
          j++;
        }
    }
  free(OLDSHELL);
  
  M = j;
  last_sphere_count = 0 ;
  current_sphere_count = M;
  SPCOUNT[0] = 0 ;
  SPCOUNT[1] = M ;

  printf("Spheres in new core: %d \n", M);

  sort_ascending_point(ALLSHELL, M);
  
  for (i = 0; i < N_spheres; i++)
    {
      ALLTOUCH[i][13] =  0;
      for (j = 0; j < 12; j++)
        ALLTOUCH[i][j] = -1;
    }
    
  // Calculate the shell distances
  shell_radius[0] = 0.5;
  shell_radius[1] = 1.5;
  for (i = 2; i<N_shells; i++)
    {
      shell_radius[i] = shell_radius[i-1] + D_as + (D - D_as)/i;
    }

  // gsl_ieee_env_setup();  // read GSL_IEEE_MODE
  gsl_rng_env_setup();
  Ts = gsl_rng_default;
  s = gsl_rng_alloc (Ts);

  L = 2;

  while ((current_sphere_count < N_spheres_limit) && (L < N_shells_limit))
    {
      printf("\n------>> Procecessing shell %2d \n", L);
    
      // * * * * * * * * * * * * * * * * * * * * *
      //              FIRST PASS
      // * * * * * * * * * * * * * * * * * * * * *

      // Given shell (L-1), build a maximum set of tri-contact balls in shell L
      cnt = 0;
      SHELL_temp  = (struct point  *)malloc(TMPIDX * sizeof(struct point));
      SHELL_Ptemp = (struct Ppoint *)malloc(TMPIDX * sizeof(struct Ppoint));

      if (flip[0] == 0)
        {
          SPCOUNT_pass1[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      else {
        for (m = 0; m < SPCOUNT[L-1]-2; m++)
          {
            first = ALLSHELL[m + last_sphere_count];
            for (n = m+1; n < SPCOUNT[L-1]-1; n++)
              {
                second = ALLSHELL[n + last_sphere_count];
                if (precheck_triangle(first,second, 4*r) &&
                    !POSITIVE(distance(first,second) - 4*r))
                  {
                    for (k = n+1; k < SPCOUNT[L-1]; k++)
                      {
                        third = ALLSHELL[k + last_sphere_count];
                        if (precheck_triangle(first,third, 4*r) && precheck_triangle(second,third, 4*r) &&
                            !POSITIVE(distance(first,third)-4*r) && !POSITIVE(distance(second,third)-4*r))
                          {
                
                            passno = 0;
                            for (i = 0; i<3; i++)
                              {
                                Dt[i] = D;
                                if (dt[passno][i])
                                  Dt[i] = D + gsl_rng_uniform(s);
                              }
                            rn = calculate_fourth(first, second, third, Dt[0], Dt[1], Dt[2], fourth0, fourth1, r0, r1);
                
                            if (L == 2)
                              check_range = SPCOUNT[L-2]+SPCOUNT[L-1];
                            else
                              check_range = SPCOUNT[L-3]+SPCOUNT[L-2]+SPCOUNT[L-1];

                            if ((rn == 2) &&
                                !space_violation(*fourth0, ALLSHELL+current_sphere_count-check_range, check_range, D))
                              if (cnt==0)
                                {
                                  SHELL_temp[cnt] = *fourth0;               
                                  SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                  d1 = distance_to_origin(SHELL_temp[cnt]);
                                  d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                  SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;

                                  cnt++;
                                }
                              else {
                                exist_flag = 0;
                                for (j_n = 0; j_n<cnt; j_n++)
                                  {
                                    if (samepoint(*fourth0,SHELL_temp[j_n]))
                                      {
                                        exist_flag = 1;
                                        j_n = cnt;
                                      }
                                  }
                                if (exist_flag == 0)
                                  {
                                    SHELL_temp[cnt] = *fourth0;
                                    SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                    d1 = distance_to_origin(SHELL_temp[cnt]);
                                    d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                    SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                    cnt++;
                                  }
                              }
                          }
                      }
                  }
              }
          }

        N_next_points = cnt;
    
        if (N_next_points > 0)
          {
            if (flip[0] == 1)
              {
                sort_ascending_point(SHELL_temp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_temp[0];
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass1[L] = N_next_points - removed_pts;
              }
        
            else if (flip[0] == 4)
              {
                sort_ascending_Ppoint(SHELL_Ptemp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_Ptemp[0].point;
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_Ptemp[i].point, ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_Ptemp[i].point;
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass1[L] = N_next_points - removed_pts;
              }
        
            else if (flip[0] == 3)
              {
                sort_descending_point(SHELL_temp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_temp[0];
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass1[L] = N_next_points - removed_pts;
              }
            else if (flip[0] == 2)
              {
                for (i = 0; i < N_next_points; i++)
                  {
                    tooshort[i] = 0; mark[i] = 1;
                  }
                // Find out those pair , btw whom the distance is less than 2R
                // Count the frequency of apperance of each point in those pairs
                // S is the total number of times a point appear in such a pair -
                // hence twice the number of pair
                S = 0;
                for (i = 0; i < N_next_points - 1; i++)
                  {
                    for (j = i+1; j < N_next_points; j++)
                      {
                        if (within_distance(SHELL_temp[i], SHELL_temp[j], D))
                          {
                            tooshort[i]++;
                            tooshort[j]++;
                            S = S + 2;
                          }
                      }
                  }
                // Eliminate the most frequent-appearing points
                removed_pts = 0;
                while (S > 0)
                  {
                    maxidx = findmax(tooshort,N_next_points); // find the maximum appearance
                    tooshort[maxidx] = 0 ;                    // reset to 0
                    mark[maxidx] = 0 ;                        // mark the points to be removed
                    removed_pts++ ;                           // keep count of how many points will be removed
                    for (j = 0; j < N_next_points; j++)
                      {
                        if (mark[j] && (within_distance(SHELL_temp[j], SHELL_temp[maxidx], D)))
                          {
                            tooshort[j]-- ;
                            S = S - 2;
                          }
                      }
                  }
                //  Record the points selected after the first pass
                SPCOUNT_pass1[L] = N_next_points - removed_pts;
                for (i = 0; i < N_next_points; i++)
                  {
                    if (mark[i] == 1)
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++;
                      }
                  }
              }
            else {
              SPCOUNT_pass1[L] = 0;
              removed_pts = N_next_points;
            }
          }
        else {
          SPCOUNT_pass1[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      }
      free(SHELL_temp);
      printf("All Pts: %5d   Remd Pts: %5d       Pass 1: %5d    Total: %5d \n",
             N_next_points, removed_pts, SPCOUNT_pass1[L], current_sphere_count);
      // * * * * * * * * * * * * * * * * * * * * *
      //           END OF FIRST PASS
      // * * * * * * * * * * * * * * * * * * * * *

      // * * * * * * * * * * * * * * * * * * * * *
      //             SECOND PASS
      // * * * * * * * * * * * * * * * * * * * * *
      cnt = 0;
      SHELL_temp  = (struct point  *)malloc(TMPIDX * sizeof(struct point));
      SHELL_Ptemp = (struct Ppoint *)malloc(TMPIDX * sizeof(struct Ppoint));

      if (flip[1] == 0)
        {
          SPCOUNT_pass2[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      else {
        for (m = 0; m < SPCOUNT[L-1]-1; m++)
          {
            first =  ALLSHELL[m + last_sphere_count] ;
            for (n = m+1; n < SPCOUNT[L-1]; n++)
              {
                second = ALLSHELL[n + last_sphere_count] ;
                if (precheck_triangle(first,second, 4*r) &&
                    !POSITIVE(distance(first,second)-4*r))
                  {
                    for (k = 0; k < SPCOUNT_pass1[L]; k++)
                      {
                        third = ALLSHELL[k + current_sphere_count - SPCOUNT_pass1[L]] ;
                        if (precheck_triangle(first,third, 4*r) && precheck_triangle(second,third, 4*r) &&
                            !POSITIVE(distance(first,third)-4*r) && !POSITIVE(distance(second,third)-4*r))
                          {
                            // Eliminate duplication in the listsample_codes/

                            passno = 1;
                            for (i = 0; i<3; i++)
                              {
                                Dt[i] = D;
                                if (dt[passno][i])
                                  Dt[i] = D+gsl_rng_uniform(s);
                              }
                            rn = calculate_fourth(first, second, third, Dt[0], Dt[1], Dt[2], fourth0, fourth1, r0, r1);

                            if (L==2)
                              check_range = SPCOUNT[L-2]+SPCOUNT[L-1]+SPCOUNT_pass1[L];
                            else
                              check_range = SPCOUNT[L-3]+SPCOUNT[L-2]+SPCOUNT[L-1]+SPCOUNT_pass1[L];
                
                            if ((rn == 2) &&
                                !space_violation(*fourth0, ALLSHELL+current_sphere_count-check_range, check_range, D))
                              if (cnt==0)
                                {
                                  SHELL_temp[cnt] = *fourth0;
                                  SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                  d1 = distance_to_origin(SHELL_temp[cnt]);
                                  d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                  SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                  cnt++;
                                }
                              else {
                                exist_flag = 0;
                                for (j_n = 0; j_n<cnt; j_n++)
                                  {
                                    if (samepoint(*fourth0,SHELL_temp[j_n]))
                                      {
                                        exist_flag = 1;
                                        j_n = cnt;
                                      }
                                  }
                                SHELL_temp[cnt] = *fourth0;
                                SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                d1 = distance_to_origin(SHELL_temp[cnt]);
                                d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                cnt++;
                              }
                          }
                      }
                  }
              }
          }

        N_next_points = cnt;

        if (N_next_points > 0)
          {

            if (flip[1] ==1)
              {
                sort_ascending_point(SHELL_temp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_temp[0];
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                //  Record the points selected after the second pass
                SPCOUNT_pass2[L] = N_next_points - removed_pts;
              }

            else if (flip[1] == 4)
              {
                sort_ascending_Ppoint(SHELL_Ptemp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_Ptemp[0].point;
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_Ptemp[i].point, ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_Ptemp[i].point;
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass2[L] = N_next_points - removed_pts;
              }
        
            else if (flip[1] ==3)
              {
                sort_descending_point(SHELL_temp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_temp[0];
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                //  Record the points selected after the second pass
                SPCOUNT_pass2[L] = N_next_points - removed_pts;
              }
            else if (flip[1] ==2)
              {
                for (i = 0; i < N_next_points; i++)
                  {
                    tooshort[i] = 0; mark[i] = 1;
                  }
                // Find out those pair , btw whom the distance is less than 2R
                // Count the frequency of apperance of each point in those pairs
                // S is the total number of times a point appear in such a pair -
                // hence twice the number of pair
                S = 0;
                for (i = 0; i < N_next_points-1; i++)
                  {
                    for (j = i+1; j < N_next_points; j++)
                      {
                        if (within_distance(SHELL_temp[i], SHELL_temp[j], D))
                          {
                            tooshort[i]++;
                            tooshort[j]++;
                            S = S + 2;
                          }
                      }
                  }
                // Eliminate the most frequent-appearing points
                removed_pts = 0;
                while (S > 0)
                  {
                    maxidx = findmax(tooshort,N_next_points); // find the maximum appearance
                    tooshort[maxidx] = 0 ;                    // reset to 0
                    mark[maxidx] = 0 ;                        // mark the points to be removed
                    removed_pts++ ;                           // keep count of how many points will be removed
                    for (j = 0; j < N_next_points; j++)
                      {
                        if (mark[j] && (within_distance(SHELL_temp[j], SHELL_temp[maxidx], D)))
                          {
                            tooshort[j]-- ;
                            S = S - 2;
                          }
                      }
                  }
                //  Record the points selected after the second pass
                SPCOUNT_pass2[L] = N_next_points - removed_pts;
                for (i = 0; i < N_next_points; i++)
                  {
                    if (mark[i] == 1)
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i] ;
                        current_sphere_count ++ ;
                      }
                  }
              }
            else {
              SPCOUNT_pass2[L] = 0;
              removed_pts = N_next_points;
            }
          }
        else {
          SPCOUNT_pass2[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      }
      free(SHELL_temp);
      printf("All Pts: %5d   Remd Pts: %5d       Pass 2: %5d    Total: %5d \n",
             N_next_points, removed_pts, SPCOUNT_pass2[L], current_sphere_count);
      // * * * * * * * * * * * * * * * * * * * * *
      //           END OF SECOND PASS
      // * * * * * * * * * * * * * * * * * * * * *

      // * * * * * * * * * * * * * * * * * * * * *
      //               THIRD PASS
      // * * * * * * * * * * * * * * * * * * * * *
      cnt = 0;
      SHELL_temp  = (struct point  *)malloc(TMPIDX * sizeof(struct point));
      SHELL_Ptemp = (struct Ppoint *)malloc(TMPIDX * sizeof(struct Ppoint));

      if (flip[2] == 0)
        {
          SPCOUNT_pass3[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      else
        {
          // Eliminate duplication in the list
          for (m = 0; m < SPCOUNT[L-1]; m++)
            {
              first = ALLSHELL[m + last_sphere_count];
              for (n = 0; n < (SPCOUNT_pass1[L] + SPCOUNT_pass2[L]-1); n++)
                {
                  second = ALLSHELL[n + current_sphere_count - (SPCOUNT_pass1[L] + SPCOUNT_pass2[L])];
                  if (precheck_triangle(first,second, 4*r) &&
                      !POSITIVE(distance(first,second)-4*r))
                    {
                      for (k = n+1; k < (SPCOUNT_pass1[L] + SPCOUNT_pass2[L]); k++)
                        {
                          third = ALLSHELL[k + current_sphere_count - (SPCOUNT_pass1[L] + SPCOUNT_pass2[L])];
                          if (precheck_triangle(first,third, 4*r) && precheck_triangle(second,third, 4*r) &&
                              !POSITIVE(distance(first,third)-4*r) && !POSITIVE(distance(second,third)-4*r))
                            {
                              passno = 2;
                              for (i = 0; i<3; i++)
                                {
                                  Dt[i] = D;
                                  if (dt[passno][i])
                                    Dt[i] = D+gsl_rng_uniform(s);
                                }
                              rn = calculate_fourth(first, second, third, Dt[0], Dt[1], Dt[2], fourth0, fourth1, r0, r1);
    
                              if (L==2)
                                check_range = SPCOUNT[L-2]+SPCOUNT[L-1]+SPCOUNT_pass1[L]+SPCOUNT_pass2[L];
                              else
                                check_range = SPCOUNT[L-3]+SPCOUNT[L-2]+SPCOUNT[L-1]+SPCOUNT_pass1[L]+SPCOUNT_pass2[L];
                
                              if ((rn == 2) &&
                                  !space_violation(*fourth0, ALLSHELL+current_sphere_count-check_range, check_range, D))
                                if (cnt==0)
                                  {
                                    SHELL_temp[cnt] = *fourth0;
                                    SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                    d1 = distance_to_origin(SHELL_temp[cnt]);
                                    d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                    SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                    cnt++;
                                  }
                                else
                                  {
                                    exist_flag = 0;
                                    for (j_n = 0; j_n<cnt; j_n++)
                                      {
                                        if (samepoint(*fourth0,SHELL_temp[j_n]))
                                          {
                                            exist_flag = 1;
                                            j_n = cnt;
                                          }
                                      }
                                    SHELL_temp[cnt] = *fourth0;
                                    SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                    d1 = distance_to_origin(SHELL_temp[cnt]);
                                    d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                    SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                    cnt++;
                                  }
                            }
                        }
                    }
                }
            }

          N_next_points = cnt;
  
          if (N_next_points > 0)
            {

              if (flip[2]==1)
                {
                  sort_ascending_point(SHELL_temp,N_next_points);
                  removed_pts = 0;
                  ALLSHELL[current_sphere_count] = SHELL_temp[0];
                  current_sphere_count++;
                  j = 1;
                  for (i = 1; i < N_next_points; i++)
                    {
                      if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                        {
                          ALLSHELL[current_sphere_count] = SHELL_temp[i];
                          current_sphere_count++ ;
                          j++ ;
                        }
                      else
                        removed_pts++ ;
                    }
                  SPCOUNT_pass3[L] = N_next_points - removed_pts;
                }
              else if (flip[2] == 4)
                {
                  sort_ascending_Ppoint(SHELL_Ptemp,N_next_points);
                  removed_pts = 0;
                  ALLSHELL[current_sphere_count] = SHELL_Ptemp[0].point;
                  current_sphere_count++;
                  j = 1;
                  for (i = 1; i < N_next_points; i++)
                    {
                      if (!space_violation(SHELL_Ptemp[i].point, ALLSHELL + current_sphere_count - j, j, D))
                        {
                          ALLSHELL[current_sphere_count] = SHELL_Ptemp[i].point;
                          current_sphere_count++ ;
                          j++ ;
                        }
                      else
                        removed_pts++ ;
                    }
                  SPCOUNT_pass3[L] = N_next_points - removed_pts;
                }
              else if (flip[2]==3)
                {
                  sort_descending_point(SHELL_temp,N_next_points);
                  removed_pts = 0;
                  ALLSHELL[current_sphere_count] = SHELL_temp[0];
                  current_sphere_count++;
                  j = 1;
                  for (i = 1; i < N_next_points; i++)
                    {
                      if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                        {
                          ALLSHELL[current_sphere_count] = SHELL_temp[i];
                          current_sphere_count++ ;
                          j++ ;
                        }
                      else
                        removed_pts++ ;
                    }
                  SPCOUNT_pass3[L] = N_next_points - removed_pts;
                }
              else if (flip[2]==2)
                {
                  for (i = 0; i < N_next_points; i++)
                    {
                      tooshort[i] = 0; mark[i] = 1;
                    }
                  // Find out those pair, between whom the distance is less than 2R
                  // Count the frequency of apperance of each point in those pairs
                  // S is the total number of times a point appear in such a pair -
                  // hence twice the number of pair
                  S = 0;
                  for (i = 0; i < N_next_points-1; i++)
                    {
                      for (j = i+1; j < N_next_points; j++)
                        {
                          if (within_distance(SHELL_temp[i], SHELL_temp[j], D))
                            {
                              tooshort[i]++;
                              tooshort[j]++;
                              S = S + 2;
                            }
                        }
                    }
                  // Eliminate the most frequent-appearing points
                  removed_pts = 0;
                  while (S > 0)
                    {
                      maxidx = findmax(tooshort,N_next_points); // find the maximum appearance
                      tooshort[maxidx] = 0 ;                    // reset to 0
                      mark[maxidx] = 0 ;                        // mark the points to be removed
                      removed_pts++ ;                           // keep count of how many points will be removed
                      for (j = 0; j < N_next_points; j++)
                        {
                          if (mark[j] && (within_distance(SHELL_temp[j], SHELL_temp[maxidx], D)))
                            {
                              tooshort[j]-- ;
                              S = S - 2;
                            }
                        }
                    }
                  //  Record the points selected after the third pass
                  SPCOUNT_pass3[L] = N_next_points - removed_pts;
                  for (i = 0; i < N_next_points; i++)
                    {
                      if (mark[i] == 1)
                        {
                          ALLSHELL[current_sphere_count] = SHELL_temp[i];
                          current_sphere_count++ ;
                        }
                    }
                }
              else {
                SPCOUNT_pass3[L] = 0;
                removed_pts = N_next_points;
              }
            }
          else
            {
              SPCOUNT_pass3[L] = 0;
              N_next_points = 0;
              removed_pts = N_next_points;
            }
        }
      free(SHELL_temp);
      printf("All Pts: %5d   Remd Pts: %5d       Pass 3: %5d    Total: %5d \n",
             N_next_points, removed_pts, SPCOUNT_pass3[L], current_sphere_count);
      // * * * * * * * * * * * * * * * * * * * * *
      //            END OF THIRD PASS
      // * * * * * * * * * * * * * * * * * * * * *

      // Update shell information
      SPCOUNT[L] = SPCOUNT_pass1[L] + SPCOUNT_pass2[L]+ SPCOUNT_pass3[L];

      last_sphere_count = last_sphere_count+SPCOUNT[L-1];
      max_radius = 0;
      max_idx = 1000;
      min_radius = 1000;
      min_idx = 1000;
      for (i = 0; i<SPCOUNT[L]; i++)
        {
          if (POSITIVE(distance(ALLSHELL[i+last_sphere_count],zero) - max_radius))
            {
              max_radius = distance(ALLSHELL[i+last_sphere_count],zero);
              max_idx = i;
            }
          if (NEGATIVE(distance(ALLSHELL[i+last_sphere_count],zero) - min_radius))
            {
              min_radius = distance(ALLSHELL[i+last_sphere_count],zero);
              min_idx = i;
            }
        }

      max_radius = max_radius + r;
      min_radius = min_radius + r;

      printf("In shell: %2d; in cell: %2d, packing = %6.3f\n",
             SPCOUNT[L], current_sphere_count,
             current_sphere_count*pow((r/max_radius),3));
      printf("Shell radius: Max %6.3f, Min %6.3f; Calculated %6.3f - %6.3f\n",
             max_radius, min_radius, shell_radius[L-1],shell_radius[L]);
  
      // * * * * * * * * * * * * * * * * * * * * *
      //              FOURTH PASS
      // * * * * * * * * * * * * * * * * * * * * *
      cnt = 0;
      SHELL_temp  = (struct point  *)malloc(TMPIDX * sizeof(struct point));
      SHELL_Ptemp = (struct Ppoint *)malloc(TMPIDX * sizeof(struct Ppoint));

      if (flip[3] == 0)
        {
          SPCOUNT_pass1a[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      else {
        for (m = 0; m < SPCOUNT[L]-2; m++)
          {
            first = ALLSHELL[m + last_sphere_count];
            for (n = m+1; n < SPCOUNT[L]-1; n++)
              {
                second = ALLSHELL[n + last_sphere_count];
                if (precheck_triangle(first,second, 4*r) &&
                    !POSITIVE(distance(first,second) - 4*r))
                  {
                    for (k = n+1; k < SPCOUNT[L]; k++)
                      {
                        third = ALLSHELL[k + last_sphere_count];
                        if (precheck_triangle(first,third, 4*r) && precheck_triangle(second,third, 4*r) &&
                            !POSITIVE(distance(first,third)-4*r) && !POSITIVE(distance(second,third)-4*r))
                          {
                            passno = 3;
                            for (i = 0; i<3; i++)
                              {
                                Dt[i] = D;
                                if (dt[passno][i])
                                  Dt[i] = D+gsl_rng_uniform(s);
                              }
                            rn = calculate_fourth(first, second, third, Dt[0], Dt[1], Dt[2], fourth0, fourth1, r0, r1);
                
                            check_range = SPCOUNT[L-2]+SPCOUNT[L-1]+SPCOUNT[L] ;
                
                            if ((rn == 2) && NEGATIVE(distance(*fourth0, zero) - max_radius + r) &&
                                !space_violation(*fourth0, ALLSHELL+current_sphere_count-check_range, check_range, D))
                              if (cnt==0)
                                {
                                  SHELL_temp[cnt] = *fourth0;
                                  SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                  d1 = distance_to_origin(SHELL_temp[cnt]);
                                  d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                  SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                  cnt++;
                                }
                              else
                                {
                                  exist_flag = 0;
                                  for (j_n = 0; j_n<cnt; j_n++)
                                    {
                                      if (samepoint(*fourth0,SHELL_temp[j_n]))
                                        {
                                          exist_flag = 1;
                                          j_n = cnt;
                                        }
                                    }
                                  if (exist_flag == 0)
                                    {
                                      SHELL_temp[cnt] = *fourth0;
                                      SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                      d1 = distance_to_origin(SHELL_temp[cnt]);
                                      d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                      SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                      cnt++;
                                    }
                                }
                          }
                      }
                  }
              }
          }
        N_next_points = cnt;

        if (N_next_points > 0)
          {
            if (flip[3] == 1)
              {
                sort_ascending_point(SHELL_temp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_temp[0];
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass1a[L] = N_next_points - removed_pts;
              }

            else if (flip[3] == 4)
              {
                sort_ascending_Ppoint(SHELL_Ptemp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_Ptemp[0].point;
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_Ptemp[i].point, ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_Ptemp[i].point;
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass1a[L] = N_next_points - removed_pts;
              }
        
            else if (flip[3] == 3)
              {
                sort_descending_point(SHELL_temp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_temp[0];
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass1a[L] = N_next_points - removed_pts;
              }
            else if (flip[3]==2)
              {
                for (i = 0; i < N_next_points; i++)
                  {
                    tooshort[i] = 0; mark[i] = 1;
                  }
                S = 0;
                for (i = 0; i < N_next_points-1; i++)
                  {
                    for (j = i+1; j < N_next_points; j++)
                      {
                        if (within_distance(SHELL_temp[i], SHELL_temp[j], D))
                          {
                            tooshort[i]++;
                            tooshort[j]++;
                            S = S + 2;
                          }
                      }
                  }
                // Eliminate the most frequent-appearing points
                removed_pts = 0;
                while (S > 0)
                  {
                    maxidx = findmax(tooshort,N_next_points); // find the maximum appearance
                    tooshort[maxidx] = 0 ;                    // reset to 0
                    mark[maxidx] = 0 ;                        // mark the points to be removed
                    removed_pts++ ;                           // keep count of how many points will be removed
                    for (j = 0; j < N_next_points; j++)
                      {
                        if (mark[j] && (within_distance(SHELL_temp[j], SHELL_temp[maxidx], D)))
                          {
                            tooshort[j]-- ;
                            S = S - 2;
                          }
                      }
                  }
                SPCOUNT_pass1a[L] = N_next_points - removed_pts;
                for (i = 0; i < N_next_points; i++)
                  {
                    if (mark[i] == 1)
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                      }
                  }
              }
            else {
              SPCOUNT_pass1a[L] = 0;
              removed_pts = N_next_points;
            }
          }
        else {
          SPCOUNT_pass1a[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      }
      free(SHELL_temp);
      printf("All Pts: %5d   Remd Pts: %5d       Pass 1: %5d    Total: %5d \n",
             N_next_points, removed_pts, SPCOUNT_pass1a[L], current_sphere_count);
      // * * * * * * * * * * * * * * * * * * * * *
      //           END OF FOURTH PASS
      // * * * * * * * * * * * * * * * * * * * * *


      // * * * * * * * * * * * * * * * * * * * * *
      //             FIFTH PASS
      // * * * * * * * * * * * * * * * * * * * * *
      cnt = 0;
      SHELL_temp  = (struct point  *)malloc(TMPIDX * sizeof(struct point));
      SHELL_Ptemp = (struct Ppoint *)malloc(TMPIDX * sizeof(struct Ppoint));

      if (flip[4] == 0)
        {
          SPCOUNT_pass2a[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      else {
        for (m = 0; m < SPCOUNT[L]-1; m++)
          {
            first =  ALLSHELL[m + last_sphere_count] ;
            for (n = m+1; n < SPCOUNT[L]; n++)
              {
                second = ALLSHELL[n + last_sphere_count] ;
                if (precheck_triangle(first,second, 4*r) &&
                    !POSITIVE(distance(first,second)-4*r))
                  {
                    for (k = 0; k < SPCOUNT_pass1a[L]; k++)
                      {
                        third = ALLSHELL[k + current_sphere_count - SPCOUNT_pass1a[L]] ;
                        if (precheck_triangle(first,third, 4*r) && precheck_triangle(second,third, 4*r) &&
                            !POSITIVE(distance(first,third)-4*r) && !POSITIVE(distance(second,third)-4*r))
                          {
                            // Eliminate duplication in the list

                            passno = 4;
                            for (i = 0; i<3; i++)
                              {
                                Dt[i] = D;
                                if (dt[passno][i])
                                  Dt[i] = D+gsl_rng_uniform(s);
                              }
                            rn = calculate_fourth(first, second, third, Dt[0], Dt[1], Dt[2], fourth0, fourth1, r0, r1);

                            check_range = SPCOUNT[L-2]+SPCOUNT[L-1]+SPCOUNT[L]+SPCOUNT_pass1a[L];
                            if ((rn == 2) &&
                                NEGATIVE(distance(*fourth0, zero) - max_radius + r) &&
                                !space_violation(*fourth0, ALLSHELL+current_sphere_count-check_range, check_range, D))
                  
                              if (cnt==0)
                                {
                                  SHELL_temp[cnt] = *fourth0;
                                  SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                  d1 = distance_to_origin(SHELL_temp[cnt]);
                                  d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                  SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                  cnt++;
                                }
                              else
                                {
                                  exist_flag = 0;
                                  for (j_n = 0; j_n<cnt; j_n++)
                                    {
                                      if (samepoint(*fourth0,SHELL_temp[j_n]))
                                        {
                                          exist_flag = 1;
                                          j_n = cnt;
                                        }
                                    }
                                  SHELL_temp[cnt] = *fourth0;
                                  SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                  d1 = distance_to_origin(SHELL_temp[cnt]);
                                  d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                  SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                  cnt++;
                                }
                          }
                      }
                  }
              }
          }

        N_next_points = cnt;

        if (N_next_points > 0)
          {
            if (flip[4] ==1)
              {
                sort_ascending_point(SHELL_temp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_temp[0];
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                //  Record the points selected after the second pass
                SPCOUNT_pass2a[L] = N_next_points - removed_pts;
              }

            else if (flip[4] == 4)
              {
                sort_ascending_Ppoint(SHELL_Ptemp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_Ptemp[0].point;
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_Ptemp[i].point, ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_Ptemp[i].point;
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass2a[L] = N_next_points - removed_pts;
              }
        
            else if (flip[4] == 3)
              {
                sort_descending_point(SHELL_temp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_temp[0];
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                //  Record the points selected after the second pass
                SPCOUNT_pass2a[L] = N_next_points - removed_pts;
              }
            else if (flip[4] == 2)
              {
                for (i = 0; i < N_next_points; i++)
                  {
                    tooshort[i] = 0; mark[i] = 1;
                  }
                S = 0;
                for (i = 0; i < N_next_points-1; i++)
                  {
                    for (j = i+1; j < N_next_points; j++)
                      {
                        if (within_distance(SHELL_temp[i], SHELL_temp[j], D))
                          {
                            tooshort[i]++;
                            tooshort[j]++;
                            S = S + 2;
                          }
                      }
                  }
                removed_pts = 0;
                while (S > 0)
                  {
                    maxidx = findmax(tooshort,N_next_points); // find the maximum appearance
                    tooshort[maxidx] = 0 ;                    // reset to 0
                    mark[maxidx] = 0 ;                        // mark the points to be removed
                    removed_pts++ ;                           // keep count of how many points will be removed
                    for (j = 0; j < N_next_points; j++)
                      {
                        if (mark[j] && (within_distance(SHELL_temp[j], SHELL_temp[maxidx], D)))
                          {
                            tooshort[j]-- ;
                            S = S - 2;
                          }
                      }
                  }
                SPCOUNT_pass2a[L] = N_next_points - removed_pts;
                for (i = 0; i < N_next_points; i++)
                  {
                    if (mark[i] == 1)
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                      }
                  }
              }
            else {
              SPCOUNT_pass2a[L] = 0;
              removed_pts = N_next_points;
            }
          }
        else {
          SPCOUNT_pass2a[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      }
      free(SHELL_temp);
      printf("All Pts: %5d   Remd Pts: %5d       Pass 2: %5d    Total: %5d \n",
             N_next_points, removed_pts, SPCOUNT_pass2a[L], current_sphere_count);
      // * * * * * * * * * * * * * * * * * * * * *
      //           END OF FIFTH PASS
      // * * * * * * * * * * * * * * * * * * * * *

      // * * * * * * * * * * * * * * * * * * * * *
      //               SIXTH PASS
      // * * * * * * * * * * * * * * * * * * * * *
      cnt = 0;
      SHELL_temp  = (struct point  *)malloc(TMPIDX * sizeof(struct point));
      SHELL_Ptemp = (struct Ppoint *)malloc(TMPIDX * sizeof(struct Ppoint));

      // Eliminate duplication in the list
      if (flip[5] == 0)
        {
          SPCOUNT_pass3a[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      else {
        for (m = 0; m < SPCOUNT[L]; m++)
          {
            first = ALLSHELL[m + last_sphere_count];
            for (n = 0; n < (SPCOUNT_pass1a[L] + SPCOUNT_pass2a[L]-1); n++)
              {
                second = ALLSHELL[n + current_sphere_count - (SPCOUNT_pass1a[L] + SPCOUNT_pass2a[L])];
                if (precheck_triangle(first,second, 4*r) &&
                    !POSITIVE(distance(first,second)-4*r))
                  {
                    for (k = n+1; k < (SPCOUNT_pass1a[L] + SPCOUNT_pass2a[L]); k++)
                      {
                        third = ALLSHELL[k + current_sphere_count - (SPCOUNT_pass1a[L] + SPCOUNT_pass2a[L])];
                        if (precheck_triangle(first,third, 4*r) && precheck_triangle(second,third, 4*r) &&
                            !POSITIVE(distance(first,third)-4*r) && !POSITIVE(distance(second,third)-4*r))
                          {

                            passno = 5;
                            for (i = 0; i<3; i++)
                              {
                                Dt[i] = D;
                                if (dt[passno][i])
                                  Dt[i] = D+gsl_rng_uniform(s);
                              }
                            rn = calculate_fourth(first, second, third, Dt[0], Dt[1], Dt[2], fourth0, fourth1, r0, r1);
                            check_range = SPCOUNT[L-2]+SPCOUNT[L-1]+SPCOUNT[L]+SPCOUNT_pass1a[L]+SPCOUNT_pass2a[L];
                            if ((rn == 2) &&
                                NEGATIVE(distance(*fourth0, zero) - max_radius + r) &&
                                !space_violation(*fourth0, ALLSHELL+current_sphere_count-check_range, check_range, D))
                              if (cnt==0)
                                {
                                  SHELL_temp[cnt] = *fourth0;
                                  SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                  d1 = distance_to_origin(SHELL_temp[cnt]);
                                  d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                  SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                  cnt++;
                                }
                              else
                                {
                                  exist_flag = 0;
                                  for (j_n = 0; j_n<cnt; j_n++)
                                    {
                                      if (samepoint(*fourth0,SHELL_temp[j_n]))
                                        {
                                          exist_flag = 1;
                                          j_n = cnt;
                                        }
                                    }
                                  SHELL_temp[cnt] = *fourth0;
                                  SHELL_Ptemp[cnt].point = SHELL_temp[cnt];
                                  d1 = distance_to_origin(SHELL_temp[cnt]);
                                  d2 = distance_point_to_plane(SHELL_temp[cnt],first,second,third);
                                  SHELL_Ptemp[cnt].param1 = d1*d1+d2*d2;
                                  cnt++;
                                }
                          }
                      }
                  }
              }
          }

        N_next_points = cnt;

        if (N_next_points > 0)
          {
            if (flip[5]==1)
              {
                sort_ascending_point(SHELL_temp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_temp[0];
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass3a[L] = N_next_points - removed_pts;
              }

            else if (flip[5] == 4)
              {
                sort_ascending_Ppoint(SHELL_Ptemp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_Ptemp[0].point;
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_Ptemp[i].point, ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_Ptemp[i].point;
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass3a[L] = N_next_points - removed_pts;
              }

            else if (flip[5]==3)
              {
                sort_descending_point(SHELL_temp,N_next_points);
                removed_pts = 0;
                ALLSHELL[current_sphere_count] = SHELL_temp[0];
                current_sphere_count++;
                j = 1;
                for (i = 1; i < N_next_points; i++)
                  {
                    if (!space_violation(SHELL_temp[i], ALLSHELL + current_sphere_count - j, j, D))
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                        j++ ;
                      }
                    else
                      removed_pts++ ;
                  }
                SPCOUNT_pass3a[L] = N_next_points - removed_pts;
              }
            else if (flip[5] == 2)
              {
                for (i = 0; i < N_next_points; i++)
                  {
                    tooshort[i] = 0; mark[i] = 1;
                  }
                S = 0;
                for (i = 0; i < N_next_points-1; i++)
                  {
                    for (j = i+1; j < N_next_points; j++)
                      {
                        if (within_distance(SHELL_temp[i], SHELL_temp[j], D))
                          {
                            tooshort[i]++;
                            tooshort[j]++;
                            S = S + 2;
                          }
                      }
                  }
                removed_pts = 0;
                while (S > 0)
                  {
                    maxidx = findmax(tooshort,N_next_points); // find the maximum appearance
                    tooshort[maxidx] = 0 ;                    // reset to 0
                    mark[maxidx] = 0 ;                        // mark the points to be removed
                    removed_pts++ ;                           // keep count of how many points will be removed
                    for (j = 0; j < N_next_points; j++)
                      {
                        if (mark[j] && (within_distance(SHELL_temp[j], SHELL_temp[maxidx], D)))
                          {
                            tooshort[j]-- ;
                            S = S - 2;
                          }
                      }
                  }
                SPCOUNT_pass3a[L] = N_next_points - removed_pts;
                for (i = 0; i < N_next_points; i++)
                  {
                    if (mark[i] == 1)
                      {
                        ALLSHELL[current_sphere_count] = SHELL_temp[i];
                        current_sphere_count++ ;
                      }
                  }
              }
            else {
              SPCOUNT_pass3a[L] = 0;
              removed_pts = N_next_points;
            }
          }
        else {
          SPCOUNT_pass3a[L] = 0;
          N_next_points = 0;
          removed_pts = N_next_points;
        }
      }
      free(SHELL_temp);
      printf("All Pts: %5d   Remd Pts: %5d       Pass 3: %5d    Total: %5d \n",
             N_next_points, removed_pts, SPCOUNT_pass3a[L], current_sphere_count);
      // * * * * * * * * * * * * * * * * * * * * *
      //            END OF SIXTH PASS
      // * * * * * * * * * * * * * * * * * * * * *

      SPCOUNT[L] = SPCOUNT[L] + SPCOUNT_pass1a[L] + SPCOUNT_pass2a[L]+ SPCOUNT_pass3a[L];

      max_radius = 0;
      max_idx = 1000;
      min_radius = 1000;
      min_idx = 1000;
      for (i = SPCOUNT[L] - (SPCOUNT_pass1a[L] + SPCOUNT_pass2a[L]+ SPCOUNT_pass3a[L]); i<SPCOUNT[L]; i++)
        {
          if (POSITIVE(distance(ALLSHELL[i+last_sphere_count],zero) - max_radius))
            {
              max_radius = distance(ALLSHELL[i+last_sphere_count],zero);
              max_idx = i;
            }
          if (NEGATIVE(distance(ALLSHELL[i+last_sphere_count],zero) - min_radius))
            {
              min_radius = distance(ALLSHELL[i+last_sphere_count],zero);
              min_idx = i;
            }
        }
      max_radius = max_radius + r;
      min_radius = min_radius + r;

      printf("In shell: %2d; in cell: %2d, packing = %6.3f\n",
             SPCOUNT[L], current_sphere_count,
             current_sphere_count*pow((r/max_radius),3));
      printf("Shell radius: Max %6.3f, Min %6.3f\n",
             max_radius, min_radius);
    
      misc_shellinfo[L][0] = SPCOUNT[L] ;
      misc_shellinfo[L][1] = SPCOUNT_pass1[L] ;
      misc_shellinfo[L][2] = SPCOUNT_pass2[L] ;
      misc_shellinfo[L][3] = SPCOUNT_pass3[L] ;
      misc_shellinfo[L][4] = current_sphere_count;
      misc_shellinfo[L][5] = max_radius ;
      misc_shellinfo[L][6] = min_radius ;
      misc_shellinfo[L][7] = current_sphere_count*pow((r/max_radius),3);

      for (i = 0; i < SPCOUNT[L]; i++)
        {
          fprintf(file_center_buffer, "%g %g %g\n",
                  ALLSHELL[i+last_sphere_count].x, ALLSHELL[i+last_sphere_count].y,
                  ALLSHELL[i+last_sphere_count].z);
        }
      L++;
    }
  fclose(file_center_buffer);

  sort_ascending_point(ALLSHELL,current_sphere_count);
  for (i = 0; i < current_sphere_count; i++)
    {
      fprintf(file_center, "%g %g %g\n",
              ALLSHELL[i].x, ALLSHELL[i].y, ALLSHELL[i].z);
    }
  fclose(file_center);
  
  // Check for validity of ball placements
  printf("Checking for invalid entries in the data...");
  for (i = 0; i < current_sphere_count-1; i++)
    {
      for (j = i+1; j < current_sphere_count; j++)
        {
          if (precheck_triangle(ALLSHELL[i],ALLSHELL[j], 2*r))
            {
              gap = distance(ALLSHELL[i],ALLSHELL[j]) - D ;
              if (NEGATIVE2(gap) || ISZERO2(gap)){
                ALLTOUCH[i][ALLTOUCH[i][12]] = j;
                ALLTOUCH[j][ALLTOUCH[j][12]] = i;
                ALLTOUCH[i][12]++;
                ALLTOUCH[j][12]++;
                if (NEGATIVE2(gap))
                  {
                    printf("!!! ERROR !!! %2d %2d %6.4f\n",i,j, distance(ALLSHELL[i],ALLSHELL[j]));
                    printf("%6.4f %6.4f %6.4f x\n", ALLSHELL[i].x,ALLSHELL[i].y,ALLSHELL[i].z);
                    printf("%6.4f %6.4f %6.4f x\n", ALLSHELL[j].x,ALLSHELL[j].y,ALLSHELL[j].z);
                  }
              }
              if (POSITIVE2(distance_to_origin(ALLSHELL[j]) - distance_to_origin(ALLSHELL[i]) - D))
                break;
            }
        }
    }

  // write the touch file
  for (i = 0; i < current_sphere_count; i++)
    {
      for (j = 0; j <13; j++)
        {
          fprintf(file_touch, "%6d ",  ALLTOUCH[i][j]);
        }
      fprintf(file_touch, "\n");
    }
  fclose(file_touch);


  // dump other information
  for (i = 2; i < N_shells; ++i)
    {
      for (j = 0; j <8; j++)
        {
          fprintf(file_other, "%g ",  misc_shellinfo[i][j]); }
      fprintf(file_other, "\n");
    }
  fclose(file_other);

  return 0;
}
