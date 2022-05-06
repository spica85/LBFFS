//------------------------------------------------------------------------------
//
// kernel:  updatePUVW
//
// Purpose: Compute the elementwise sum c = a+b
//
// input: a and b float vectors of length count
//
// output: c float vector of length count holding the sum a + b
//

__kernel void updateP
(
   __global float* p, __global float* u, __global float* v, __global float* w,
   __global float* pTmp,
   __global unsigned* upID,
   const unsigned elements,
   const float Ap
)
{
   int i = get_global_id(0);

   float wt[19] = {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
   float cx[19] = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
   float cy[19] = {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
   float cz[19] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0};

   float pUp[19];
   float uUp[19];
   float vUp[19];
   float wUp[19];

   pTmp[i] = 0.f;
   for(int q = 0; q < 19; q++)
   {
      int qid = q*elements+i;
      int uID = upID[qid];
      pUp[q] = p[uID];
      uUp[q] = u[uID];
      vUp[q] = v[uID];
      wUp[q] = w[uID];

      float UUpDotC = cx[q]*uUp[q] + cy[q]*vUp[q] + cz[q]*wUp[q];
      float uSqr = uUp[q]*uUp[q] + vUp[q]*vUp[q] + wUp[q]*wUp[q];

      float f_eq_in_Utmp = wt[q]*(3.0f*UUpDotC + 4.5f*UUpDotC*UUpDotC - 1.5*uSqr);
      float f_eq_in = 3.0f*wt[q]*pUp[q] +f_eq_in_Utmp;
      // f_eq_in_U[qid] = f_eq_in_Utmp;
      // f_eq_in_U[qid] = wt[q]*(3.0f*UUpDotC + 4.5f*UUpDotC*UUpDotC - 1.5*uSqr);
      // float f_eq_in = 3.0f*wt[q]*pUp[q] +f_eq_in_U[qid];

      float UDotC = cx[q]*u[i] + cy[q]*v[i] + cz[q]*w[i];

      float f_Delta_Utmp = wt[q]*(UDotC -UUpDotC);
      float f = f_eq_in + 3.0f*Ap*f_Delta_Utmp;
      // f_Delta_U[qid] = f_Delta_Utmp;
      // f_Delta_U[qid] = wt[q]*(UDotC -UUpDotC);
      // float f = f_eq_in + 3.0f*Ap*f_Delta_U[qid];

      pTmp[i] += f;
      // printf("i=%d, pUp=%f", i, pUp[q]);
   }
   pTmp[i] /= 3.0f;
   // printf("i=%d, p=%f", i, pTmp[i]);
}

__kernel void BC_P
(
   __global float* p,
   __global unsigned* innerID, __global unsigned* wallID
)
{
   int i = get_global_id(0);
   p[wallID[i]] = p[innerID[i]];
}

__kernel void updateU
(
   __global float* pTmp, __global float* u, __global float* v, __global float* w,
   __global float* uTmp, __global float* vTmp, __global float* wTmp,
   __global unsigned* upID,
   const unsigned elements,
   const float Au
)
{
   int i = get_global_id(0);

   float wt[19] = {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
   float cx[19] = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
   float cy[19] = {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
   float cz[19] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0};

   uTmp[i] = 0.f;
   vTmp[i] = 0.f;
   wTmp[i] = 0.f;

   float pUp[19];
   float uUp[19];
   float vUp[19];
   float wUp[19];

   for(int q = 0; q < 19; q++)
   {
      int qid = q*elements+i;
      int uID = upID[qid];
      pUp[q] = pTmp[uID];
      uUp[q] = u[uID];
      vUp[q] = v[uID];
      wUp[q] = w[uID];

      float UUpDotC = cx[q]*uUp[q] + cy[q]*vUp[q] + cz[q]*wUp[q];
      float uSqr = uUp[q]*uUp[q] + vUp[q]*vUp[q] + wUp[q]*wUp[q];

      float f_eq_in_Utmp = wt[q]*(3.0f*UUpDotC + 4.5f*UUpDotC*UUpDotC - 1.5*uSqr);
      float f_eq_in = 3.0f*wt[q]*pUp[q] +f_eq_in_Utmp;

      float UDotC = cx[q]*u[i] + cy[q]*v[i] + cz[q]*w[i];

      float f_Delta_Utmp = wt[q]*(UDotC -UUpDotC);

      float f = f_eq_in + 3.0f*Au*f_Delta_Utmp;
      uTmp[i] += cx[q]*f;
      vTmp[i] += cy[q]*f;
      wTmp[i] += cz[q]*f;
   }
}

__kernel void BC_U_fixedWall
(
   __global float* u, __global float* v, __global float* w,
   __global unsigned* fixedWallID
)
{
   int i = get_global_id(0);
   u[fixedWallID[i]] = 0.f;
   v[fixedWallID[i]] = 0.f;
   w[fixedWallID[i]] = 0.f;
}

__kernel void BC_U_movingWall
(
   __global float* u, __global float* v, __global float* w,
   __global unsigned* movingWallID,
   const float u0
)
{
   int i = get_global_id(0);
   u[movingWallID[i]] = u0;
   v[movingWallID[i]] = 0.f;
   w[movingWallID[i]] = 0.f;
}