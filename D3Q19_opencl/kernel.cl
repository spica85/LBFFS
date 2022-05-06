__kernel void k_collision
(
   __global float* f, __global float* fTmp,
   const unsigned elements,
   const float omega
)
{
    int i = get_global_id(0);

    float wt[19] = {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    float cx[19] = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
    float cy[19] = {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
    float cz[19] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0};

    float rho = 0.0;
    float u = 0.0;
    float v = 0.0;
    float w = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int qic = q*elements +i;

        rho += f[qic];

        u += f[qic]*cx[q];
        v += f[qic]*cy[q];
        w += f[qic]*cz[q];
    }
    u /= rho;
    v /= rho;
    w /= rho;

    for(int q = 0; q < 19; q++)
    {
        float uSqr =u*u+v*v+w*w;
        float uDotC = u*cx[q]+v*cy[q]+w*cz[q];
        float feq = (1.0+3.0*uDotC +4.5*uDotC*uDotC -1.5*uSqr)*wt[q]*rho;

        int qic = q*elements +i;

        fTmp[qic] = f[qic] - omega *(f[qic] -feq);
        // printf("feq: %.1f", feq);
    }
}

__kernel void k_externalForce
(
   __global float* f, __global float* fTmp,
   const unsigned elements,
   const float dpdx
)
{
    int i = get_global_id(0);

    float wt[19] = {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    float cx[19] = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};

    float rho = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int qic = q*elements +i;

        rho += f[qic];
    }

    for(int q = 0; q < 19; q++)
    {
        int qic = q*elements +i;

        fTmp[qic] += rho*wt[q]*3.0f*dpdx*cx[q];
    }
}

__kernel void k_streaming
(
   __global float* f, __global float* fTmp,
   __global unsigned* upQID
)
{
    int qi = get_global_id(0);
    f[qi] = fTmp[upQID[qi]];
}

__kernel void k_bounceBack
(
   __global float* f,
   __global int* BBID, __global int* BBQID,
   const unsigned elements, const unsigned nBB
)
{
    int i = get_global_id(0);
    for(int q = 0; q < 19; q++)
    {
        int qic = q*elements +BBID[i];
        int qid = q*nBB +i;
        if(BBQID[qid] >= 0)
        {
            f[qic] = f[BBQID[qid]];
        }
        else
        {
            f[qic] = 1.0f/36.0f;
        }
    }   
}

__kernel void k_bounceBackMovingWall
(
    __global float* f,
    __global int* BBmovingWallID,
    __global float* u0, __global float* v0, __global float* w0,
    __global int* normal,
    const unsigned elements
)
{
    int i = get_global_id(0);
    int ic = BBmovingWallID[i];
    float u = u0[ic];
    float v = v0[ic];
    float w = w0[ic];
    int nVec = normal[ic];  

    float* Vc    = &f[0*elements +ic];//(0)
    float* Vin   = &f[1*elements +ic];//(+x)
    float* Vip   = &f[2*elements +ic];//(-x)
    float* Vjn   = &f[3*elements +ic];//(+y)
    float* Vjp   = &f[4*elements +ic];//(-y)
    float* Vkn   = &f[5*elements +ic];//(+z)
    float* Vkp   = &f[6*elements +ic];//(-z)
    float* Vinjn = &f[7*elements +ic];//(+x,+y)
    float* Vipjp = &f[8*elements +ic];//(-x,-y)
    float* Vinjp = &f[9*elements +ic];//(+x,-y)
    float* Vipjn = &f[10*elements +ic];//(-x,+y)
    float* Vinkn = &f[11*elements +ic];//(+x,+z)
    float* Vipkp = &f[12*elements +ic];//(-x,-z)
    float* Vinkp = &f[13*elements +ic];//(+x,-z)
    float* Vipkn = &f[14*elements +ic];//(-x,+z)
    float* Vjnkn = &f[15*elements +ic];//(+y,+z)
    float* Vjpkp = &f[16*elements +ic];//(-y,-z)
    float* Vjnkp = &f[17*elements +ic];//(+y,-z)
    float* Vjpkn = &f[18*elements +ic];//(-y,+z)

    if(nVec == -1)
    {
        float rho =( (*Vc)+(*Vjn)+(*Vjp)+(*Vkn)+(*Vkp)+(*Vjnkn)+(*Vjpkp)+(*Vjnkp)+(*Vjpkn) + 2.0f*((*Vip)+(*Vipjn)+(*Vipjp)+(*Vipkn)+(*Vipkp)))/( 1.0f-u);
        (*Vin) = (*Vip) +rho*u/ 3.0f;
        float Nyx = -rho*v/ 3.0f + 0.5f*((*Vjn)+(*Vjnkn)+(*Vjnkp)-((*Vjp)+(*Vjpkn)+(*Vjpkp)));
        float Nzx = -rho*w/ 3.0f + 0.5f*((*Vkn)+(*Vjnkn)+(*Vjpkn)-((*Vkp)+(*Vjnkp)+(*Vjpkp)));
        (*Vinjn) = (*Vipjp) +rho*(u+v)/ 6.0f -Nyx;
        (*Vinjp) = (*Vipjn) +rho*(u-v)/ 6.0f +Nyx;
        (*Vinkn) = (*Vipkp) +rho*(u+w)/ 6.0f -Nzx;
        (*Vinkp) = (*Vipkn) +rho*(u-w)/ 6.0f +Nzx;
    }
    else if(nVec == 1)
    {
        float rho =( (*Vc)+(*Vjn)+(*Vjp)+(*Vkn)+(*Vkp)+(*Vjnkn)+(*Vjpkp)+(*Vjnkp)+(*Vjpkn) + 2.0f*((*Vin)+(*Vinjn)+(*Vinjp)+(*Vinkn)+(*Vinkp)))/( 1.0f+u);
        (*Vip) = (*Vin) -rho*u/ 3.0f;
        float Nyx = -rho*v/ 3.0f + 0.5f*((*Vjn)+(*Vjnkn)+(*Vjnkp)-((*Vjp)+(*Vjpkn)+(*Vjpkp)));
        float Nzx = -rho*w/ 3.0f + 0.5f*((*Vkn)+(*Vjnkn)+(*Vjpkn)-((*Vkp)+(*Vjnkp)+(*Vjpkp)));
        (*Vipjp) = (*Vinjn) -rho*(u+v)/ 6.0f +Nyx;
        (*Vipjn) = (*Vinjp) -rho*(u-v)/ 6.0f -Nyx;
        (*Vipkp) = (*Vinkn) -rho*(u+w)/ 6.0f +Nzx;
        (*Vipkn) = (*Vinkp) -rho*(u-w)/ 6.0f -Nzx;
    }
    else if(nVec == -2)
    {
        float rho = ( (*Vc)+(*Vin)+(*Vip)+(*Vkn)+(*Vkp)+(*Vinkn)+(*Vipkp)+(*Vipkn)+(*Vinkp) + 2.0f*((*Vjp)+(*Vinjp)+(*Vipjp)+(*Vjpkn)+(*Vjpkp)))/( 1.0f-v);
        (*Vjn) = (*Vjp) +rho*v/ 3.0f;
        float Nxy = -rho*u/ 3.0f + 0.5f*((*Vin)+(*Vinkn)+(*Vinkp)-((*Vip)+(*Vipkn)+(*Vipkp)));
        float Nzy = -rho*w/ 3.0f + 0.5f*((*Vkn)+(*Vinkn)+(*Vipkn)-((*Vkp)+(*Vinkp)+(*Vipkp)));
        (*Vinjn) = (*Vipjp) +rho*(v+u)/ 6.0f -Nxy;
        (*Vipjn) = (*Vinjp) +rho*(v-u)/ 6.0f +Nxy;
        (*Vjnkn) = (*Vjpkp) +rho*(v+w)/ 6.0f -Nzy;
        (*Vjnkp) = (*Vjpkn) +rho*(v-w)/ 6.0f +Nzy;
    }
    else if(nVec == 2)
    {
        float rho =( (*Vc)+(*Vin)+(*Vip)+(*Vkn)+(*Vkp)+(*Vinkn)+(*Vipkp)+(*Vipkn)+(*Vinkp) + 2.0f*((*Vjn)+(*Vinjn)+(*Vipjn)+(*Vjnkn)+(*Vjnkp)))/( 1.0f+v);
        (*Vjp) = (*Vjn) -rho*v/ 3.0f;
        float Nxy = -rho*u/ 3.0f + 0.5f*((*Vin)+(*Vinkn)+(*Vinkp)-((*Vip)+(*Vipkn)+(*Vipkp)));
        float Nzy = -rho*w/ 3.0f + 0.5f*((*Vkn)+(*Vinkn)+(*Vipkn)-((*Vkp)+(*Vinkp)+(*Vipkp)));
        (*Vipjp) = (*Vinjn)  -rho*(v+u)/ 6.0f +Nxy;
        (*Vinjp) = (*Vipjn)  -rho*(v-u)/ 6.0f -Nxy;
        (*Vjpkp) = (*Vjnkn)  -rho*(v+w)/ 6.0f +Nzy;
        (*Vjpkn) = (*Vjnkp)  -rho*(v-w)/ 6.0f -Nzy;
    }
    else if(nVec == -3)
    {
        float rho = ( (*Vc)+(*Vin)+(*Vip)+(*Vjn)+(*Vjp)+(*Vinjn)+(*Vipjp)+(*Vinjp)+(*Vipjn)+ 2.0f*((*Vkp)+(*Vipkp)+(*Vinkp)+(*Vjpkp)+(*Vjnkp)))/( 1.0f-w);
        (*Vkn) = (*Vkp) +rho*w/ 3.0f;
        float Nxz = -rho*u/ 3.0f + 0.5f*((*Vin)+(*Vinjn)+(*Vinjp)-((*Vip)+(*Vipjn)+(*Vipjp)));
        float Nyz = -rho*v/ 3.0f + 0.5f*((*Vjn)+(*Vinjn)+(*Vipjp)-((*Vjp)+(*Vinjp)+(*Vipjp)));
        (*Vinkn) = (*Vipkp) +rho*(w+u)/ 6.0f -Nxz;
        (*Vipkn) = (*Vinkp) +rho*(w-u)/ 6.0f +Nxz;
        (*Vjnkn) = (*Vjpkp) +rho*(w+v)/ 6.0f -Nyz;
        (*Vjpkn) = (*Vjnkn) +rho*(w-v)/ 6.0f +Nyz;
    }
    else if(nVec == 3)
    {
        float rho = ( (*Vc)+(*Vin)+(*Vip)+(*Vjn)+(*Vjp)+(*Vinjn)+(*Vipjp)+(*Vinjp)+(*Vipjn)+ 2.0f*((*Vkn)+(*Vipkn)+(*Vinkn)+(*Vjpkn)+(*Vjnkn)))/( 1.0f+w);
        (*Vkp) = (*Vkn) -rho*w/ 3.0f;
        float Nxz = -rho*u/ 3.0f + 0.5f*((*Vin)+(*Vinjn)+(*Vinjp)-((*Vip)+(*Vipjn)+(*Vipjp)));
        float Nyz = -rho*v/ 3.0f + 0.5f*((*Vjn)+(*Vinjn)+(*Vipjp)-((*Vjp)+(*Vinjp)+(*Vipjp)));
        (*Vipkp) = (*Vinkn) -rho*(w+u)/ 6.0f +Nxz;
        (*Vinkp) = (*Vipkn) -rho*(w-u)/ 6.0f -Nxz;
        (*Vjpkp) = (*Vjnkn) -rho*(w+v)/ 6.0f +Nyz;
        (*Vjnkp) = (*Vjpkn) -rho*(w-v)/ 6.0f -Nyz;
    }
}