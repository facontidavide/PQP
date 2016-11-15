/*************************************************************************\

  Copyright 1999 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             S. Gottschalk
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:               geom@cs.unc.edu


\**************************************************************************/

#ifndef PQP_MATVEC_H
#define PQP_MATVEC_H

#include <math.h>
#include <stdio.h>
#include "PQP_Compile.h"

#ifndef M_PI
const PQP_REAL M_PI = (PQP_REAL)3.14159265359;
#endif

#ifdef gnu
#include "zzzz.h"

#ifdef hppa
#define myfabs(x) \
 ({double __value, __arg = (x); \
   asm("fabs,dbl %1, %0": "=f" (__value): "f" (__arg)); \
   __value; \
});
#endif

#ifdef mips
#define myfabs(x) \
 ({double __value, __arg = (x); \
   asm("abs.d %0, %1": "=f" (__value): "f" (__arg)); \
   __value; \
});
#endif

#else  

#define myfabs(x) ((x < 0) ? -x : x)

#endif


inline
void
Mprintg(const Matrix& M)
{
  printf("%g %g %g\n%g %g %g\n%g %g %g\n",
         M(0,0), M(0,1), M(0,2),
         M(1,0), M(1,1), M(1,2),
         M(2,0), M(2,1), M(2,2));
}


inline
void
Mfprint(FILE *f, const Matrix& M)
{
  fprintf(f, "%g %g %g\n%g %g %g\n%g %g %g\n",
         M(0,0), M(0,1), M(0,2),
         M(1,0), M(1,1), M(1,2),
         M(2,0), M(2,1), M(2,2));
}

inline
void
Mprint(const Matrix& M)
{
  printf("%g %g %g\n%g %g %g\n%g %g %g\n",
         M(0,0), M(0,1), M(0,2),
         M(1,0), M(1,1), M(1,2),
         M(2,0), M(2,1), M(2,2));
}

inline
void
Vprintg(const Vector& V)
{
  printf("%g %g %g\n", V[0], V[1], V[2]);
}

inline
void
Vfprint(FILE *f, const Vector& V)
{
  fprintf(f, "%g %g %g\n", V[0], V[1], V[2]);
}

inline
void
Vprint(const Vector& V)
{
  printf("%g %g %g\n", V[0], V[1], V[2]);
}

inline
void
Midentity(Matrix& M)
{
  M(0,0) = M(1,1) = M(2,2) = 1.0;
  M(0,1) = M(1,2) = M(2,0) = 0.0;
  M(0,2) = M(1,0) = M(2,1) = 0.0;
}

inline
void
Videntity(Vector& T)
{
  T[0] = T[1] = T[2] = 0.0;
}

inline
void
McM(Matrix& Mr, const Matrix& M)
{
  Mr = M;
}

inline
void
MTcM(Matrix& Mr, const Matrix& M)
{
  Mr = M.transpose();
}

inline
void
VcV(Vector& Vr, const Vector& V)
{
  Vr = V;
}

inline
void
McolcV(Vector& Vr, const Matrix& M, int c)
{
  Vr = M.col(c);
}

inline
void
McolcMcol(Matrix& Mr, int cr, const Matrix& M, int c)
{
  Mr.col(cr) = M.col(c);
}

inline
void
MxMpV(Matrix& Mr, const Matrix& M1, const Matrix& M2, const Vector& T)
{
  Mr(0,0) = (M1(0,0) * M2(0,0) +
              M1(0,1) * M2(1,0) +
              M1(0,2) * M2(2,0) +
	      T[0]);
  Mr(1,0) = (M1(1,0) * M2(0,0) +
              M1(1,1) * M2(1,0) +
              M1(1,2) * M2(2,0) +
	      T[1]);
  Mr(2,0) = (M1(2,0) * M2(0,0) +
              M1(2,1) * M2(1,0) +
              M1(2,2) * M2(2,0) +
	      T[2]);
  Mr(0,1) = (M1(0,0) * M2(0,1) +
              M1(0,1) * M2(1,1) +
              M1(0,2) * M2(2,1) +
	      T[0]);
  Mr(1,1) = (M1(1,0) * M2(0,1) +
              M1(1,1) * M2(1,1) +
              M1(1,2) * M2(2,1) +
	      T[1]);
  Mr(2,1) = (M1(2,0) * M2(0,1) +
              M1(2,1) * M2(1,1) +
              M1(2,2) * M2(2,1) +
	      T[2]);
  Mr(0,2) = (M1(0,0) * M2(0,2) +
              M1(0,1) * M2(1,2) +
              M1(0,2) * M2(2,2) +
	      T[0]);
  Mr(1,2) = (M1(1,0) * M2(0,2) +
              M1(1,1) * M2(1,2) +
              M1(1,2) * M2(2,2) +
	      T[1]);
  Mr(2,2) = (M1(2,0) * M2(0,2) +
              M1(2,1) * M2(1,2) +
              M1(2,2) * M2(2,2) +
	      T[2]);
}

inline
void
MxM(Matrix& Mr, const Matrix& M1, const Matrix& M2)
{
  Mr = M1*M2;
}


inline
void
MxMT(Matrix& Mr, const Matrix& M1, const Matrix& M2)
{
  Mr = M1*M2.transpose();
}

inline
void
MTxM(Matrix& Mr, const Matrix& M1, const Matrix& M2)
{
  Mr = M1.transpose()*M2;
}

inline
void
MxV(Vector& Vr, const Matrix& M1, const Vector& V1)
{
  Vr = M1*V1;
}


inline
void
MxVpV(Vector& Vr, const Matrix& M1, const Vector& V1, const Vector& V2)
{
    Vr = M1*V1 + V2;
}


inline
void
sMxVpV(Vector& Vr, PQP_REAL s1, const Matrix& M1, const Vector& V1, const Vector& V2)
{
  Vr = s1 * (M1 * V1);
}

inline
void
MTxV(Vector& Vr, const Matrix& M1, const Vector& V1)
{
    Vr = M1.transpose()*V1;
}

inline
void
sMTxV(Vector& Vr, PQP_REAL s1, const Matrix& M1, const Vector& V1)
{
    Vr = M1.transpose()*V1;
}

inline
void
sMxV(Vector& Vr, PQP_REAL s1, const Matrix& M1, const Vector& V1)
{
    Vr = s1* ( M1*V1);
}


inline
void
VmV(Vector& Vr, const Vector& V1, const Vector& V2)
{
  Vr = V1 - V2;
}

inline
void
VpV(Vector& Vr, const Vector& V1, const Vector& V2)
{
  Vr = V1 + V2;
}

inline
void
VpVxS(Vector& Vr, const Vector& V1, const Vector& V2, PQP_REAL s)
{
  Vr = V1 + V2 * s;
}

inline 
void
MskewV(Matrix& M, const Vector& v)
{
  M(0,0) = M(1,1) = M(2,2) = 0.0;
  M(1,0) = v[2];
  M(0,1) = -v[2];
  M(0,2) = v[1];
  M(2,0) = -v[1];
  M(1,2) = -v[0];
  M(2,1) = v[0];
}


inline
void
VcrossV(Vector& Vr, const Vector& V1, const Vector& V2)
{
  Vr[0] = V1[1]*V2[2] - V1[2]*V2[1];
  Vr[1] = V1[2]*V2[0] - V1[0]*V2[2];
  Vr[2] = V1[0]*V2[1] - V1[1]*V2[0];
}

inline
PQP_REAL
Vlength(Vector& V)
{
  return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
}

inline
void
Vnormalize(Vector& V)
{
  PQP_REAL d = (PQP_REAL)1.0 / sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
  V *= d;
}

inline
PQP_REAL
VdotV(const Vector& V1, const Vector& V2)
{
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}

inline
PQP_REAL
VdistV2(const Vector& V1, const Vector& V2)
{
  return ( (V1[0]-V2[0]) * (V1[0]-V2[0]) + 
	   (V1[1]-V2[1]) * (V1[1]-V2[1]) + 
	   (V1[2]-V2[2]) * (V1[2]-V2[2]));
}

inline
void
VxS(Vector& Vr, const Vector& V, PQP_REAL s)
{
  Vr = V * s;
}

inline
void
MRotZ(Matrix& Mr, PQP_REAL t)
{
  Mr(0,0) = cos(t);
  Mr(1,0) = sin(t);
  Mr(0,1) = -Mr(1,0);
  Mr(1,1) = Mr(0,0);
  Mr(2,0) = Mr(2,1) = 0.0;
  Mr(0,2) = Mr(1,2) = 0.0;
  Mr(2,2) = 1.0;
}

inline
void
MRotX(Matrix& Mr, PQP_REAL t)
{
  Mr(1,1) = cos(t);
  Mr(2,1) = sin(t);
  Mr(1,2) = -Mr(2,1);
  Mr(2,2) = Mr(1,1);
  Mr(0,1) = Mr(0,2) = 0.0;
  Mr(1,0) = Mr(2,0) = 0.0;
  Mr(0,0) = 1.0;
}

inline
void
MRotY(Matrix& Mr, PQP_REAL t)
{
  Mr(2,2) = cos(t);
  Mr(0,2) = sin(t);
  Mr(2,0) = -Mr(0,2);
  Mr(0,0) = Mr(2,2);
  Mr(1,2) = Mr(1,0) = 0.0;
  Mr(2,1) = Mr(0,1) = 0.0;
  Mr(1,1) = 1.0;
}

inline
void
MVtoOGL(double oglm[16], const Matrix& R, const Vector& T)
{
  oglm[0] = (double)R(0,0);
  oglm[1] = (double)R(1,0);
  oglm[2] = (double)R(2,0);
  oglm[3] = 0.0;
  oglm[4] = (double)R(0,1);
  oglm[5] = (double)R(1,1);
  oglm[6] = (double)R(2,1);
  oglm[7] = 0.0;
  oglm[8] = (double)R(0,2);
  oglm[9] = (double)R(1,2);
  oglm[10] = (double)R(2,2);
  oglm[11] = 0.0;
  oglm[12] = (double)T[0];
  oglm[13] = (double)T[1];
  oglm[14] = (double)T[2];
  oglm[15] = 1.0;
}

inline 
void
OGLtoMV(Matrix& R, Vector& T, const double oglm[16])
{
  R(0,0) = (PQP_REAL)oglm[0];
  R(1,0) = (PQP_REAL)oglm[1];
  R(2,0) = (PQP_REAL)oglm[2];

  R(0,1) = (PQP_REAL)oglm[4];
  R(1,1) = (PQP_REAL)oglm[5];
  R(2,1) = (PQP_REAL)oglm[6];

  R(0,2) = (PQP_REAL)oglm[8];
  R(1,2) = (PQP_REAL)oglm[9];
  R(2,2) = (PQP_REAL)oglm[10];

  T[0] = (PQP_REAL)oglm[12];
  T[1] = (PQP_REAL)oglm[13];
  T[2] = (PQP_REAL)oglm[14];
}

// taken from quatlib, written by Richard Holloway
const int QX = 0;
const int QY = 1;
const int QZ = 2;
const int QW = 3;

inline
void 
MRotQ(Matrix& destMatrix, Eigen::Vector4d srcQuat)
{
  PQP_REAL  s;
  PQP_REAL  xs, ys, zs,
    	    wx, wy, wz,
	        xx, xy, xz,
	        yy, yz, zz;

  /* 
   * For unit srcQuat, just set s = 2.0; or set xs = srcQuat[QX] + 
   *   srcQuat[QX], etc. 
   */

  s = (PQP_REAL)2.0 / (srcQuat(QX)*srcQuat(QX) + srcQuat(QY)*srcQuat(QY) +
             srcQuat(QZ)*srcQuat(QZ) + srcQuat(QW)*srcQuat(QW));

  xs = srcQuat(QX) * s;   ys = srcQuat(QY) * s;   zs = srcQuat(QZ) * s;
  wx = srcQuat(QW) * xs;  wy = srcQuat(QW) * ys;  wz = srcQuat(QW) * zs;
  xx = srcQuat(QX) * xs;  xy = srcQuat(QX) * ys;  xz = srcQuat(QX) * zs;
  yy = srcQuat(QY) * ys;  yz = srcQuat(QY) * zs;  zz = srcQuat(QZ) * zs;

  destMatrix(QX,QX) = (PQP_REAL)1.0 - (yy + zz);
  destMatrix(QX,QY) = xy + wz;
  destMatrix(QX,QZ) = xz - wy;

  destMatrix(QY,QX) = xy - wz;
  destMatrix(QY,QY) = (PQP_REAL)1.0 - (xx + zz);
  destMatrix(QY,QZ) = yz + wx;

  destMatrix(QZ,QX) = xz + wy;
  destMatrix(QZ,QY) = yz - wx;
  destMatrix(QZ,QZ) = (PQP_REAL)1.0 - (xx + yy);
} 

inline
void
Mqinverse(Matrix& Mr, Matrix& m)
{
  int i,j;

  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
    {
      int i1 = (i+1)%3;
      int i2 = (i+2)%3;
      int j1 = (j+1)%3;
      int j2 = (j+2)%3;
      Mr(i,j) = (m(j1,i1)*m(j2,i2) - m(j1,i2)*m(j2,i1));
    }
}

// Meigen from Numerical Recipes in C

#if 0

#define rfabs(x) ((x < 0) ? -x : x)

#define ROT(a,i,j,k,l) g=a[i][j]; h=a[k][l]; a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);

int
inline
Meigen(Matrix& vout, PQP_REAL dout[3], Matrix& a)
{
  int i;
  PQP_REAL tresh,theta,tau,t,sm,s,h,g,c;
  int nrot;
  PQP_REAL b[3];
  PQP_REAL z[3];
  Vector& V[3];
  PQP_REAL d[3];

  v(0,0) = v(1,1) = v(2,2) = 1.0;
  v(0,1) = v(1,2) = v(2,0) = 0.0;
  v(0,2) = v(1,0) = v(2,1) = 0.0;
  
  b[0] = a(0,0); d[0] = a(0,0); z[0] = 0.0;
  b[1] = a(1,1); d[1] = a(1,1); z[1] = 0.0;
  b[2] = a(2,2); d[2] = a(2,2); z[2] = 0.0;

  nrot = 0;

  
  for(i=0; i<50; i++)
    {

      printf("2\n");

      sm=0.0; sm+=fabs(a(0,1)); sm+=fabs(a(0,2)); sm+=fabs(a(1,2));
      if (sm == 0.0) { McM(vout,v); VcV(dout,d); return i; }
      
      if (i < 3) tresh=0.2*sm/(3*3); else tresh=0.0;
      
      {
        g = 100.0*rfabs(a(0,1));
	if (i>3 && rfabs(d[0])+g==rfabs(d[0]) && rfabs(d[1])+g==rfabs(d[1]))
          a(0,1)=0.0;
        else if (rfabs(a(0,1))>tresh)
	  {
	    h = d[1]-d[0];
            if (rfabs(h)+g == rfabs(h)) t=(a(0,1))/h;
	    else
	      {
                theta=0.5*h/(a(0,1));
		t=1.0/(rfabs(theta)+sqrt(1.0+theta*theta));
		if (theta < 0.0) t = -t;
	      }
            c=1.0/sqrt(1+t*t); s=t*c; tau=s/(1.0+c); h=t*a(0,1);
	    z[0] -= h; z[1] += h; d[0] -= h; d[1] += h;
            a(0,1)=0.0;
	    ROT(a,0,2,1,2); ROT(v,0,0,0,1); ROT(v,1,0,1,1); ROT(v,2,0,2,1); 
	    nrot++;
	  }
      }

      {
        g = 100.0*rfabs(a(0,2));
	if (i>3 && rfabs(d[0])+g==rfabs(d[0]) && rfabs(d[2])+g==rfabs(d[2]))
          a(0,2)=0.0;
        else if (rfabs(a(0,2))>tresh)
	  {
	    h = d[2]-d[0];
            if (rfabs(h)+g == rfabs(h)) t=(a(0,2))/h;
	    else
	      {
                theta=0.5*h/(a(0,2));
		t=1.0/(rfabs(theta)+sqrt(1.0+theta*theta));
		if (theta < 0.0) t = -t;
	      }
            c=1.0/sqrt(1+t*t); s=t*c; tau=s/(1.0+c); h=t*a(0,2);
	    z[0] -= h; z[2] += h; d[0] -= h; d[2] += h;
            a(0,2)=0.0;
	    ROT(a,0,1,1,2); ROT(v,0,0,0,2); ROT(v,1,0,1,2); ROT(v,2,0,2,2); 
	    nrot++;
	  }
      }


      {
        g = 100.0*rfabs(a(1,2));
	if (i>3 && rfabs(d[1])+g==rfabs(d[1]) && rfabs(d[2])+g==rfabs(d[2]))
          a(1,2)=0.0;
        else if (rfabs(a(1,2))>tresh)
	  {
	    h = d[2]-d[1];
            if (rfabs(h)+g == rfabs(h)) t=(a(1,2))/h;
	    else
	      {
                theta=0.5*h/(a(1,2));
		t=1.0/(rfabs(theta)+sqrt(1.0+theta*theta));
		if (theta < 0.0) t = -t;
	      }
            c=1.0/sqrt(1+t*t); s=t*c; tau=s/(1.0+c); h=t*a(1,2);
	    z[1] -= h; z[2] += h; d[1] -= h; d[2] += h;
            a(1,2)=0.0;
	    ROT(a,0,1,0,2); ROT(v,0,1,0,2); ROT(v,1,1,1,2); ROT(v,2,1,2,2); 
	    nrot++;
	  }
      }

      b[0] += z[0]; d[0] = b[0]; z[0] = 0.0;
      b[1] += z[1]; d[1] = b[1]; z[1] = 0.0;
      b[2] += z[2]; d[2] = b[2]; z[2] = 0.0;
      
    }

  fprintf(stderr, "eigen: too many iterations in Jacobi transform (%d).\n", i);

  return i;
}

#else



#define ROTATE(a,i,j,k,l) g=a(i,j); h=a(k,l); a(i,j)=g-s*(h+g*tau); a(k,l)=h+s*(g-h*tau);

void
inline
Meigen(Matrix& vout, Vector& dout, Matrix& a)
{
  int n = 3;
  int j,iq,ip,i;
  PQP_REAL tresh,theta,tau,t,sm,s,h,g,c;
  int nrot;
  Vector b;
  Vector z;
  Matrix v;
  Vector d;
  
  Midentity(v);
  for(ip=0; ip<n; ip++) 
    {
      b[ip] = a(ip,ip);
      d[ip] = a(ip,ip);
      z[ip] = 0.0;
    }
  
  nrot = 0;
  
  for(i=0; i<50; i++)
    {

      sm=0.0;
      for(ip=0;ip<n;ip++) for(iq=ip+1;iq<n;iq++) sm+=fabs(a(ip,iq));
      if (sm == 0.0)
	{
	  McM(vout, v);
	  VcV(dout, d);
	  return;
	}
      
      
      if (i < 3) tresh=(PQP_REAL)0.2*sm/(n*n);
      else tresh=0.0;
      
      for(ip=0; ip<n; ip++) for(iq=ip+1; iq<n; iq++)
	{
          g = (PQP_REAL)100.0*fabs(a(ip,iq));
	  if (i>3 && 
	      fabs(d[ip])+g==fabs(d[ip]) && 
	      fabs(d[iq])+g==fabs(d[iq]))
            a(ip,iq)=0.0;
          else if (fabs(a(ip,iq))>tresh)
	    {
	      h = d[iq]-d[ip];
              if (fabs(h)+g == fabs(h)) t=(a(ip,iq))/h;
	      else
		{
                  theta=(PQP_REAL)0.5*h/(a(ip,iq));
		  t=(PQP_REAL)(1.0/(fabs(theta)+sqrt(1.0+theta*theta)));
		  if (theta < 0.0) t = -t;
		}
	      c=(PQP_REAL)1.0/sqrt(1+t*t);
	      s=t*c;
	      tau=s/((PQP_REAL)1.0+c);
              h=t*a(ip,iq);
	      z[ip] -= h;
	      z[iq] += h;
	      d[ip] -= h;
	      d[iq] += h;
              a(ip,iq)=0.0;
	      for(j=0;j<ip;j++) { ROTATE(a,j,ip,j,iq); } 
	      for(j=ip+1;j<iq;j++) { ROTATE(a,ip,j,j,iq); } 
	      for(j=iq+1;j<n;j++) { ROTATE(a,ip,j,iq,j); } 
	      for(j=0;j<n;j++) { ROTATE(v,j,ip,j,iq); } 
	      nrot++;
	    }
	}
      for(ip=0;ip<n;ip++)
	{
	  b[ip] += z[ip];
	  d[ip] = b[ip];
	  z[ip] = 0.0;
	}
    }

  fprintf(stderr, "eigen: too many iterations in Jacobi transform.\n");

  return;
}


#endif

#endif
// MATVEC_H
