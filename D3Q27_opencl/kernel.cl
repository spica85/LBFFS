// -- D3Q27 functions
inline int ic2i(int ic, int nx, int ny)
{
    return (ic%(nx*ny))%nx;
}

inline int ic2j(int ic, int nx, int ny)
{
    return ic%(nx*ny)/nx;
}

inline int ic2k(int ic, int nx, int ny)
{
    return ic/(nx*ny);
}

inline int index1d(int i, int j, int k, int nx, int ny)
{
    return nx*ny*k+nx*j+i;
}

inline int index1df(int q, int i, int j, int k, int nx, int ny, int nz)
{
    return q*nx*ny*nz+nx*ny*k+nx*j+i;
}

inline int idf(int q, int i, int nx, int ny, int nz)
{
    return q*nx*ny*nz+i;
}

inline int c2q(const float cxi, const float cyi, const float czi)
{
    float cx[27] = {0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f};
    float cy[27] = {0.f, 0.f,  0.f, 1.f, -1.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f};
    float cz[27] = {0.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f};

    #pragma unroll
    for(int q = 0; q < 27; q++)
    {
        if(cxi == cx[q] && cyi == cy[q] && czi == cz[q])
        {
            return q;
        }
    }
}

inline int reflectQ(const int q)
{
    if(q == 0)
    {
        return 0;
    }
    else if(q == 1)
    {
        return 2;
    }
    else if(q == 2)
    {
        return 1;
    }
    else if(q == 3)
    {
        return 4;
    }
    else if(q == 4)
    {
        return 3;
    }
    else if(q == 5)
    {
        return 6;
    }
    else if(q == 6)
    {
        return 5;
    }
    else if(q == 7)
    {
        return 8;
    }
    else if(q == 8)
    {
        return 7;
    }
    else if(q == 9)
    {
        return 10;
    }
    else if(q == 10)
    {
        return 9;
    }
    else if(q == 11)
    {
        return 12;
    }
    else if(q == 12)
    {
        return 11;
    }
    else if(q == 13)
    {
        return 14;
    }
    else if(q == 14)
    {
        return 13;
    }
    else if(q == 15)
    {
        return 16;
    }
    else if(q == 16)
    {
        return 15;
    }
    else if(q == 17)
    {
        return 18;
    }
    else if(q == 18)
    {
        return 17;
    }
    else if(q == 19)
    {
        return 20;
    }
    else if(q == 20)
    {
        return 19;
    }
    else if(q == 21)
    {
        return 22;
    }
    else if(q == 22)
    {
        return 21;
    }
    else if(q == 23)
    {
        return 24;
    }
    else if(q == 24)
    {
        return 23;
    }
    else if(q == 25)
    {
        return 26;
    }
    else if(q == 26)
    {
        return 25;
    }
    // std::cerr << "Error q = [0:18], q = " << q << std::endl;
    // exit(EXIT_FAILURE);    
}

inline int reflectQs(const int q) //Simple but slow
{
    float cx[27] = {0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f};
    float cy[27] = {0.f, 0.f,  0.f, 1.f, -1.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f};
    float cz[27] = {0.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f};

    float cxi = -cx[q];
    float cyi = -cy[q];
    float czi = -cz[q];

    return c2q(cxi,cyi,czi);
}

inline int mirrorQxs(const int q) //Simple but slow
{
    float cx[27] = {0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f};
    float cy[27] = {0.f, 0.f,  0.f, 1.f, -1.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f};
    float cz[27] = {0.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f};

    float cxi = -cx[q];
    float cyi = cy[q];
    float czi = cz[q];
    return c2q(cxi, cyi, czi);
}

inline int mirrorQx(const int q)
{
    if(q == 0)
    {
        return 0;
    }
    else if(q == 1)
    {
        return 2;
    }
    else if(q == 2)
    {
        return 1;
    }
    else if(q == 3)
    {
        return 4;
    }
    else if(q == 4)
    {
        return 3;
    }
    else if(q == 5)
    {
        return 6;
    }
    else if(q == 6)
    {
        return 5;
    }
    else if(q == 7)
    {
        return 10;
    }
    else if(q == 8)
    {
        return 9;
    }
    else if(q == 9)
    {
        return 8;
    }
    else if(q == 10)
    {
        return 7;
    }
    else if(q == 11)
    {
        return 14;
    }
    else if(q == 12)
    {
        return 13;
    }
    else if(q == 13)
    {
        return 12;
    }
    else if(q == 14)
    {
        return 11;
    }
    else if(q == 15)
    {
        return 16;
    }
    else if(q == 16)
    {
        return 15;
    }
    else if(q == 17)
    {
        return 18;
    }
    else if(q == 18)
    {
        return 17;
    }
    else if(q == 19)
    {
        return 21;
    }
    else if(q == 20)
    {
        return 22;
    }
    else if(q == 21)
    {
        return 19;
    }
    else if(q == 22)
    {
        return 20;
    }
    else if(q == 23)
    {
        return 26;
    }
    else if(q == 24)
    {
        return 25;
    }
    else if(q == 25)
    {
        return 24;
    }
    else if(q == 26)
    {
        return 23;
    }
    // std::cerr << "Error q = [0:18], q = " << q << std::endl;
    // exit(EXIT_FAILURE);    
}

inline int mirrorQy(const int q)
{
    if(q == 0)
    {
        return 0;
    }
    else if(q == 1)
    {
        return 2;
    }
    else if(q == 2)
    {
        return 1;
    }
    else if(q == 3)
    {
        return 4;
    }
    else if(q == 4)
    {
        return 3;
    }
    else if(q == 5)
    {
        return 6;
    }
    else if(q == 6)
    {
        return 5;
    }
    else if(q == 7)
    {
        return 9;
    }
    else if(q == 8)
    {
        return 10;
    }
    else if(q == 9)
    {
        return 7;
    }
    else if(q == 10)
    {
        return 8;
    }
    else if(q == 11)
    {
        return 12;
    }
    else if(q == 12)
    {
        return 11;
    }
    else if(q == 13)
    {
        return 14;
    }
    else if(q == 14)
    {
        return 13;
    }
    else if(q == 15)
    {
        return 18;
    }
    else if(q == 16)
    {
        return 17;
    }
    else if(q == 17)
    {
        return 16;
    }
    else if(q == 18)
    {
        return 15;
    }
    else if(q == 19)
    {
        return 23;
    }
    else if(q == 20)
    {
        return 24;
    }
    else if(q == 21)
    {
        return 26;
    }
    else if(q == 22)
    {
        return 25;
    }
    else if(q == 23)
    {
        return 19;
    }
    else if(q == 24)
    {
        return 20;
    }
    else if(q == 25)
    {
        return 22;
    }
    else if(q == 26)
    {
        return 21;
    }
    // std::cerr << "Error q = [0:18], q = " << q << std::endl;
    // exit(EXIT_FAILURE);    
}

inline int mirrorQz(const int q)
{
    if(q == 0)
    {
        return 0;
    }
    else if(q == 1)
    {
        return 2;
    }
    else if(q == 2)
    {
        return 1;
    }
    else if(q == 3)
    {
        return 4;
    }
    else if(q == 4)
    {
        return 3;
    }
    else if(q == 5)
    {
        return 6;
    }
    else if(q == 6)
    {
        return 5;
    }
    else if(q == 7)
    {
        return 8;
    }
    else if(q == 8)
    {
        return 7;
    }
    else if(q == 9)
    {
        return 10;
    }
    else if(q == 10)
    {
        return 9;
    }
    else if(q == 11)
    {
        return 13;
    }
    else if(q == 12)
    {
        return 14;
    }
    else if(q == 13)
    {
        return 11;
    }
    else if(q == 14)
    {
        return 12;
    }
    else if(q == 15)
    {
        return 17;
    }
    else if(q == 16)
    {
        return 18;
    }
    else if(q == 17)
    {
        return 15;
    }
    else if(q == 18)
    {
        return 16;
    }
    else if(q == 19)
    {
        return 25;
    }
    else if(q == 20)
    {
        return 26;
    }
    else if(q == 21)
    {
        return 24;
    }
    else if(q == 22)
    {
        return 23;
    }
    else if(q == 23)
    {
        return 22;
    }
    else if(q == 24)
    {
        return 21;
    }
    else if(q == 25)
    {
        return 19;
    }
    else if(q == 26)
    {
        return 20;
    }
    // std::cerr << "Error q = [0:18], q = " << q << std::endl;
    // exit(EXIT_FAILURE);    
}


int upwindID_B(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz, const int boundary1, const int boundary2, const int boundary3)
{
    bool isnotWall1 = (boundary1 == 1 || boundary1 == 4 || boundary1 == 5 || boundary1 == 6 || boundary1 == 7) ? false : true;
    bool isnotWall2 = (boundary2 == 1 || boundary2 == 4 || boundary2 == 5 || boundary2 == 6 || boundary2 == 7) ? false : true;
    bool isnotWall3 = (boundary3 == 1 || boundary3 == 4 || boundary3 == 5 || boundary3 == 6 || boundary3 == 7) ? false : true;

    float cx[27] = {0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f};
    float cy[27] = {0.f, 0.f,  0.f, 1.f, -1.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f};
    float cz[27] = {0.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f};

    int upi = i -cx[q];
    int upj = j -cy[q];
    int upk = k -cz[q];

    upi = (upi == -1) ? (isnotWall1 ? nx-1 : -1) : ((upi == nx) ? (isnotWall1 ? 0 : -1) : upi);
    upj = (upj == -1) ? (isnotWall2 ? ny-1 : -1) : ((upj == ny) ? (isnotWall2 ? 0 : -1) : upj);
    upk = (upk == -1) ? (isnotWall3 ? nz-1 : -1) : ((upk == nz) ? (isnotWall3 ? 0 : -1) : upk);

    return (upi == -1 || upj == -1 || upk == -1) ? -1 : index1d(upi,upj,upk,nx,ny);
}

inline int upwindID(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz)
{
    int ic = index1d(i, j, k, nx, ny);

    if(q == 0)
    {
        return index1d(i,j,k,nx,ny);
    }
    else if(q == 1)
    {
        return i != 0 ? index1d(i-1,j,k,nx,ny) : index1d(nx-1,j,k,nx,ny);
    }
    else if(q == 2)
    {
        return i != nx-1 ? index1d(i+1,j,k,nx,ny) : index1d(0,j,k,nx,ny);
    }
    else if(q == 3)
    {
        if(j != 0)
        {
            return index1d(i,j-1,k,nx,ny);
        }
        else
        {
            return index1d(i,ny-1,k,nx,ny);
        }
    }
    else if(q == 4)
    {
        if(j != ny-1)
        {
            return index1d(i,j+1,k,nx,ny);
        }
        else
        {
            return index1d(i,0,k,nx,ny);
        }
    }
    else if(q == 5)
    {
        if(k != 0)
        {
            return index1d(i,j,k-1,nx,ny);
        }
        else
        {
            return index1d(i,j,nz-1,nx,ny);
        }
    }
    else if(q == 6)
    {
        if(k != nz-1)
        {
            return index1d(i,j,k+1,nx,ny);
        }
        else
        {
            return index1d(i,j,0,nx,ny);
        }
    }
    else if(q == 7)
    {
        return (i != 0 && j != 0) ? index1d(i-1, j-1, k, nx, ny) :
               (i == 0 && j != 0) ? index1d(nx-1, j-1, k, nx, ny) :
               (i != 0 && j == 0) ? index1d(i-1, ny-1, k, nx, ny) :
               index1d(nx-1, ny-1, k, nx, ny);
    }
    else if(q == 8) 
    {
        return (i != nx-1 && j != ny-1) ? index1d(i+1, j+1, k, nx, ny) :
               (i == nx-1 && j != ny-1) ? index1d(0, j+1, k, nx, ny) :
               (i != nx-1 && j == ny-1) ? index1d(i+1, 0, k, nx, ny) :
               index1d(0, 0, k, nx, ny);
    }
    else if(q == 9)
    {
        return (i != 0 && j != ny-1) ? index1d(i-1, j+1, k, nx, ny) :
               (i == 0 && j != ny-1) ? index1d(nx-1, j+1, k, nx, ny) :
               (i != 0 && j == ny-1) ? index1d(i-1, 0, k, nx, ny) :
               index1d(nx-1, 0, k, nx, ny);
    }
    else if(q == 10)
    {
        return (i != nx-1 && j != 0) ? index1d(i+1, j-1, k, nx, ny) :
               (i == nx-1 && j != 0) ? index1d(0, j-1, k, nx, ny) :
               (i != nx-1 && j == 0) ? index1d(i+1, ny-1, k, nx, ny) :
               index1d(0, ny-1, k, nx, ny);
    }
    else if(q == 11)
    {
        if(i != 0 && k != 0)
        {
            return index1d(i-1, j, k-1, nx, ny);
        }
        else if(i == 0 && k != 0)
        {
            return index1d(nx-1, j, k-1, nx, ny);
        }
        else if(i != 0 && k == 0)
        {
            return index1d(i-1, j, nz-1, nx, ny);
        }
        else
        {
            return index1d(nx-1, j, nz-1, nx, ny);
        }
    }
    else if(q == 12)
    {
        if(i != nx-1 && k != nz-1)
        {
            return index1d(i+1, j, k+1, nx, ny);
        }
        else if(i == nx-1 && k != nz-1)
        {
            return index1d(0, j, k+1, nx, ny);
        }
        else if(i != nx-1 && k == nz-1)
        {
            return index1d(i+1, j, 0, nx, ny);
        }
        else
        {
            return index1d(0, j, 0, nx, ny);
        }
    }
    else if(q == 13)
    {
        if(i != 0 && k != nz-1)
        {
            return index1d(i-1, j, k+1, nx, ny);
        }
        else if(i == 0 && k != nz-1)
        {
            return index1d(nx-1, j, k+1, nx, ny);
        }
        else if(i != 0 && k == nz-1)
        {
            return index1d(i-1, j, 0, nx, ny);
        }
        else
        {
            return index1d(nx-1, j, 0, nx, ny);
        }
    }
    else if(q == 14)
    {
        if(i != nx-1 && k != 0)
        {
            return index1d(i+1, j, k-1, nx, ny);
        }
        else if(i == nx-1 && k != 0)
        {
            return index1d(0, j, k-1, nx, ny);
        }
        else if(i != nx-1 && k == 0)
        {
            return index1d(i+1, j, nz-1, nx, ny);
        }
        else
        {
            return index1d(0, j, nz-1, nx, ny);
        }
    }
    else if(q == 15)
    {
        if(j != 0 && k != 0)
        {
            return index1d(i, j-1, k-1, nx, ny);
        }
        else if(j == 0 && k != 0)
        {
            return index1d(i, ny-1, k-1, nx, ny);
        }
        else if(j != 0 && k == 0)
        {
            return index1d(i, j-1, nz-1, nx, ny);
        }
        else
        {
            return index1d(i, ny-1, nz-1, nx, ny);
        }
    }
    else if(q == 16)
    {
        if(j != ny-1 && k != nz-1)
        {
            return index1d(i, j+1, k+1, nx, ny);
        }
        else if(j == ny-1 && k != nz-1)
        {
            return index1d(i, 0, k+1, nx, ny);
        }
        else if(j != ny-1 && k == nz-1)
        {
            return index1d(i, j+1, 0, nx, ny);
        }
        else
        {
            return index1d(i, 0, 0, nx, ny);
        }
    }
    else if(q == 17)
    {
        if(j != 0 && k != nz-1)
        {
            return index1d(i, j-1, k+1, nx, ny);
        }
        else if(j == 0 && k != nz-1)
        {
            return index1d(i, ny-1, k+1, nx, ny);
        }
        else if(j != 0 && k == nz-1)
        {
            return index1d(i, j-1, 0, nx, ny);
        }
        else
        {
            return index1d(i, ny-1, 0, nx, ny);
        }
    }
    else if(q == 18)
    {
        if(j != ny-1 && k != 0)
        {
            return index1d(i, j+1, k-1, nx, ny);
        }
        else if(j == ny-1 && k != 0)
        {
            return index1d(i, 0, k-1, nx, ny);
        }
        else if(j != ny-1 && k == 0)
        {
            return index1d(i, j+1, nz-1, nx, ny);
        }
        else
        {
            return index1d(i, 0, nz-1, nx, ny);
        }
    }
    else
    {
        return 0;
    }
}

inline int upwindID_internal(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz)
{
    float cx[27] = {0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f};
    float cy[27] = {0.f, 0.f,  0.f, 1.f, -1.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f};
    float cz[27] = {0.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f};

    int upi = i -cx[q];
    int upj = j -cy[q];
    int upk = k -cz[q];

    return index1d(upi,upj,upk,nx,ny);
}

inline int upwindID_2D(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz)
{
    float cx[27] = {0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f};
    float cy[27] = {0.f, 0.f,  0.f, 1.f, -1.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f};
    float cz[27] = {0.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f};

    int upi = i -cx[q];
    int upj = j -cy[q];

    return index1d(upi,upj,k,nx,ny);
}

int icNear(int ic, int iNear, int nx, int ny, int nz)
{
    int i = ic2i(ic,nx,ny);
    int j = ic2j(ic,nx,ny);
    int k = ic2k(ic,nx,ny);

    // if(i == 0 || i == nx-1 || j == 0 || j == ny-1 || k == 0 || k == nz-1)
    // {
    //     printf("Error: Missing near points!\n");
    // }

    if(iNear == 0)
    {
        return ic;
    }
    else if(iNear == 1)
    {
        return i+1 <= nx ? index1d(i+1, j, k, nx, ny) : -1;
    }
    else if(iNear == 2)
    {
        return i-1 >= 0 ? index1d(i-1, j, k, nx, ny) : -1;
    }
    else if(iNear == 3)
    {
        return j+1 <= ny ? index1d(i, j+1, k, nx, ny) : -1;
    }
    else if(iNear == 4)
    {
        return j-1 >= 0 ? index1d(i, j-1, k, nx, ny) : -1;
    }
    else if(iNear == 5)
    {
        return k+1 <= nz ? index1d(i, j, k+1, nx, ny) : -1;
    }
    else if(iNear == 6)
    {
        return k-1 >= 0 ? index1d(i, j, k-1, nx, ny) : -1;
    }
    else if(iNear == 7)
    {
        return (i+1 <= nx && j+1 <= ny) ? index1d(i+1, j+1, k, nx, ny) : -1;
    }
    else if(iNear == 8)
    {
        return (i-1 >= 0 && j-1 >= 0) ? index1d(i-1, j-1, k, nx, ny) : -1;
    }
    else if(iNear == 9)
    {
        return (i+1 <= nx && j-1 >= 0) ? index1d(i+1, j-1, k, nx, ny) : -1;
    }
    else if(iNear == 10)
    {
        return (i-1 >= 0 && j+1 <= ny) ? index1d(i-1, j+1, k, nx, ny) : -1;
    }
    else if(iNear == 11)
    {
        return (i+1 <= nx && k+1 <= nz) ? index1d(i+1, j, k+1, nx, ny) : -1;
    }
    else if(iNear == 12)
    {
        return (i-1 >= 0 && k-1 >= 0) ? index1d(i-1, j, k-1, nx, ny) : -1;
    }
    else if(iNear == 13)
    {
        return (i+1 <= nx && k-1 >= 0) ? index1d(i+1, j, k-1, nx, ny) : -1;
    }
    else if(iNear == 14)
    {
        return (i-1 >= 0 && k+1 <= nz) ? index1d(i-1, j, k+1, nx, ny) : -1;
    }
    else if(iNear == 15)
    {
        return (j+1 <= ny && k+1 <= nz) ? index1d(i, j+1, k+1, nx, ny) : -1;
    }
    else if(iNear == 16)
    {
        return (j-1 >= 0 && k-1 >= 0) ? index1d(i, j-1, k-1, nx, ny) : -1;
    }
    else if(iNear == 17)
    {
        return (j+1 <= ny && k-1 >= 0) ? index1d(i, j+1, k-1, nx, ny) : -1;
    }
    else if(iNear == 18)
    {
        return (j-1 >= 0 && k+1 <= nz) ? index1d(i, j-1, k+1, nx, ny) : -1;
    }
    else if(iNear == 19)
    {
        return (i+1 <= nx && j+1 <= ny && k+1 <= nz) ? index1d(i+1, j+1, k+1, nx, ny) : -1;
    }
    else if(iNear == 20)
    {
        return (i-1 >= 0 && j-1 >= 0 && k-1 >= 0) ? index1d(i-1, j-1, k-1, nx, ny) : -1;
    }
    else if(iNear == 21)
    {
        return (i+1 <= nx && j-1 >= 0 && k+1 <= nz) ? index1d(i+1, j-1, k+1, nx, ny) : -1;
    }
    else if(iNear == 22)
    {
        return (i-1 >= 0 && j+1 <= ny && k-1 >= 0) ? index1d(i-1, j+1, k-1, nx, ny) : -1;
    }
    else if(iNear == 23)
    {
        return (i+1 <= nx && j+1 <= ny && k-1 >= 0) ? index1d(i+1, j+1, k-1, nx, ny) : -1;
    }
    else if(iNear == 24)
    {
        return (i-1 >= 0 && j-1 >= 0 && k+1 <= nz) ? index1d(i-1, j-1, k+1, nx, ny) : -1;
    }
    else if(iNear == 25)
    {
        return (i+1 <= nx && j-1 >= 0 && k-1 >= 0) ? index1d(i+1, j-1, k-1, nx, ny) : -1;
    }
    else if(iNear == 26)
    {
        return (i-1 >= 0 && j+1 <= ny && k+1 <= nz) ? index1d(i-1, j+1, k+1, nx, ny) : -1;
    }
}

int icBox(int ic, int iBox, int nx, int ny, int nz)
{
    int i = ic2i(ic,nx,ny);
    int j = ic2j(ic,nx,ny);
    int k = ic2k(ic,nx,ny);

    // if(i == 0 || i == nx-1 || j == 0 || j == ny-1 || k == 0 || k == nz-1)
    // {
    //     printf("Error: Missing near points!\n");
    // }

    if(iBox == 0)
    {
        return ic;
    }
    else if(iBox == 1)
    {
        return i+1 < nx ? index1d(i+1, j, k, nx, ny) : -1;
    }
    else if(iBox == 2)
    {
        return j+1 < ny ? index1d(i, j+1, k, nx, ny) : -1;
    }
    else if(iBox == 3)
    {
        return k+1 < nz ? index1d(i, j, k+1, nx, ny) : -1;
    }
    else if(iBox == 4)
    {
        return (i+1 < nx && j+1 < ny) ? index1d(i+1, j+1, k, nx, ny) : -1;
    }
    else if(iBox == 5)
    {
        return (i+1 < nx && k+1 < nz) ? index1d(i+1, j, k+1, nx, ny) : -1;
    }
    else if(iBox == 6)
    {
        return (j+1 < ny && k+1 < nz) ? index1d(i, j+1, k+1, nx, ny) : -1;
    }
    else if(iBox == 7)
    {
        return (i+1 < nx && j+1 < ny && k+1 < nz) ? index1d(i+1, j+1, k+1, nx, ny) : -1;
    }
}

void cal_rhoUVW(const float* f, float* rho, float* u, float* v, float* w)
{
    float cx[27] = {0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f};
    float cy[27] = {0.f, 0.f,  0.f, 1.f, -1.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f};
    float cz[27] = {0.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, 0.f, 0.f, 0.f, 0.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f};

    *rho = 0.0f;
    *u = 0.0f;
    *v = 0.0f;
    *w = 0.0f;

    #pragma unroll
    for(int q = 0; q < 27; q++)
    {
        *rho += f[q];
        *u += f[q]*cx[q];
        *v += f[q]*cy[q];
        *w += f[q]*cz[q];
    }
    *u /= *rho;
    *v /= *rho;
    *w /= *rho;
}

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
float atom_add_float(__global float* const address, const float value)
{
  uint oldval, newval, readback;
  
  *(float*)&oldval = *address;
  *(float*)&newval = (*(float*)&oldval + value);
  while ((readback = atom_cmpxchg((__global uint*)address, oldval, newval)) != oldval) {
    oldval = readback;
    *(float*)&newval = (*(float*)&oldval + value);
  }
  return *(float*)&oldval;
}

void Urot
(
    const float wallX, const float wallY, const float wallZ, 
    const float rotX, const float rotY, const float rotZ, 
    const float rotAxisX, const float rotAxisY, const float rotAxisZ, 
    const float rotOmega, 
    float* uRot, float* vRot, float* wRot
)
{
    const float magRotAxis = sqrt(rotAxisX*rotAxisX +rotAxisY*rotAxisY +rotAxisZ*rotAxisZ);
    const float e0_x = rotAxisX;
    const float e0_y = rotAxisY;
    const float e0_z = rotAxisZ;

    const float OX_x = wallX - rotX;
    const float OX_y = wallY - rotY;
    const float OX_z = wallZ - rotZ;

    const float OXdotE0 = (OX_x*e0_x +OX_y*e0_y +OX_x*e0_z);
    const float OOp_x = OXdotE0*e0_x;
    const float OOp_y = OXdotE0*e0_y;
    const float OOp_z = OXdotE0*e0_z;

    const float OpX_x = OX_x -OOp_x;
    const float OpX_y = OX_y -OOp_y;
    const float OpX_z = OX_z -OOp_z;
    const float magOpX = sqrt(OpX_x*OpX_x +OpX_y*OpX_y +OpX_z*OpX_z);

    const float e1_x = OpX_x/magOpX;
    const float e1_y = OpX_y/magOpX;
    const float e1_z = OpX_z/magOpX;

    const float e2_x = e0_y*e1_z -e0_z*e1_y;
    const float e2_y = e0_z*e1_x -e0_x*e1_z;
    const float e2_z = e0_x*e1_y -e0_y*e1_x;

    
    const float r = sqrt((OX_x -OOp_x)*(OX_x -OOp_x) +(OX_y -OOp_y)*(OX_y -OOp_y) +(OX_z -OOp_z)*(OX_z -OOp_z));
    *uRot = r*rotOmega*e2_x;
    *vRot = r*rotOmega*e2_y;
    *wRot = r*rotOmega*e2_z;
}

void Xrot
(
    const float wallX, const float wallY, const float wallZ, 
    const float rotX, const float rotY, const float rotZ, 
    const float rotAxisX, const float rotAxisY, const float rotAxisZ, 
    const float rotOmega, 
    float* Xrot_x, float* Xrot_y, float* Xrot_z
)
{
    const float magRotAxis = sqrt(rotAxisX*rotAxisX +rotAxisY*rotAxisY +rotAxisZ*rotAxisZ);
    const float e0_x = rotAxisX;
    const float e0_y = rotAxisY;
    const float e0_z = rotAxisZ;

    const float OX_x = wallX - rotX;
    const float OX_y = wallY - rotY;
    const float OX_z = wallZ - rotZ;

    const float OXdotE0 = (OX_x*e0_x +OX_y*e0_y +OX_x*e0_z);
    const float OOp_x = OXdotE0*e0_x;
    const float OOp_y = OXdotE0*e0_y;
    const float OOp_z = OXdotE0*e0_z;

    const float OpX_x = OX_x -OOp_x;
    const float OpX_y = OX_y -OOp_y;
    const float OpX_z = OX_z -OOp_z;
    const float magOpX = sqrt(OpX_x*OpX_x +OpX_y*OpX_y +OpX_z*OpX_z);

    const float e1_x = OpX_x/magOpX;
    const float e1_y = OpX_y/magOpX;
    const float e1_z = OpX_z/magOpX;

    const float e2_x = e0_y*e1_z -e0_z*e1_y;
    const float e2_y = e0_z*e1_x -e0_x*e1_z;
    const float e2_z = e0_x*e1_y -e0_y*e1_x;

    
    const float r = sqrt((OX_x -OOp_x)*(OX_x -OOp_x) +(OX_y -OOp_y)*(OX_y -OOp_y) +(OX_z -OOp_z)*(OX_z -OOp_z));
    *Xrot_x = rotX +OOp_x +r*(cos(rotOmega)*e1_x +sin(rotOmega)*e2_x);
    *Xrot_y = rotY +OOp_y +r*(cos(rotOmega)*e1_y +sin(rotOmega)*e2_y);
    *Xrot_z = rotZ +OOp_z +r*(cos(rotOmega)*e1_z +sin(rotOmega)*e2_z);
}

__attribute__((always_inline))
void streaming(float* ft, const __global float* f, int* upID, const int boundary1, const int boundary2, const int boundary3, const int ic, const int i, const int j, const int k, const int nx, const int ny, const int nz, const int elements, __global float* Fwx, __global float* Fwy, __global float* Fwz, float* cx, float* cy, float* cz)
{
    Fwx[ic] = 0.f;
    Fwy[ic] = 0.f;
    Fwz[ic] = 0.f;

    ft[0] = f[ic];
    upID[0] = ic;

    #pragma unroll
    for(int q = 1; q < 27; q++)
    {
        upID[q] = upwindID_B(q,i,j,k,nx,ny,nz,boundary1,boundary2,boundary3);

        if(upID[q] != -1) // Streaming
        {   
            int upQID = idf(q, upID[q], nx, ny, nz);
            ft[q] = f[upQID];
        }
        else // Bounce-Back or Symmetry for boundary wall
        {
            int qbb = reflectQ(q);
            if(boundary1 == 4)
            {
                qbb = mirrorQx(q);
            }
            if(boundary2 == 4)
            {
                qbb = mirrorQy(q);
            }
            if(boundary3 == 4)
            {
                qbb = mirrorQz(q);
            }
            int bbQID = idf(qbb, ic, nx, ny, nz);
            ft[q] = f[bbQID];
            
            Fwx[ic] += -2.f*ft[q]*cx[q];
            Fwy[ic] += -2.f*ft[q]*cy[q];
            Fwz[ic] += -2.f*ft[q]*cz[q];
            // printf("q: %d\n",q);
        }
    }
}

__attribute__((always_inline))
void streamingInternal(float* ft, const __global float* f, int* upID, const int ic, const int i, const int j, const int k, const int nx, const int ny, const int nz, const int elements)
{
    ft[0] = f[ic];
    upID[0] = ic;

    #pragma unroll
    for(int q = 1; q < 27; q++)
    {
        upID[q] = upwindID_internal(q,i,j,k,nx,ny,nz);

        int upQID = idf(q, upID[q], nx, ny, nz);
        ft[q] = f[upQID];
    }
}

__attribute__((always_inline))
void streaming2D(float* ft, const __global float* f, int* upID, const int ic, const int i, const int j, const int k, const int nx, const int ny, const int nz, const int elements)
{
    ft[0] = f[ic];
    upID[0] = ic;

    #pragma unroll
    for(int q = 1; q < 27; q++)
    {
        upID[q] = upwindID_2D(q,i,j,k,nx,ny,nz);

        int upQID = idf(q, upID[q], nx, ny, nz);
        ft[q] = f[upQID];
    }
}

int cornerFlag(const int boundary3, const int i, const int j, const int k, const int nx, const int ny, const int nz)
{
    if(nz-1 == 0 || boundary3 == 0)
    {
        if
        (
            i == 0 && j == 0
            ||
            i == 0 && j == ny-1
            ||
            i == nx-1 && j == 0
            ||
            i == nx-1 && j == ny-1
        )
        {
            return 1;
        }
    }
    else
    {
        if
        (
            i == 0 && j == 0
            ||
            i == 0 && j == ny-1
            ||
            i == nx-1 && j == 0
            ||
            i == nx-1 && j == ny-1
            ||
            i == 0 && k == 0
            ||
            i == 0 && k == nz-1
            ||
            i == nx-1 && k == 0
            ||
            i == nx-1 && k == nz-1
            ||
            j == 0 && k == 0
            ||
            j == 0 && k == nz-1
            ||
            j == ny-1 && k == 0
            ||
            j == ny-1 && k == nz-1
        )
        {
            return 1;   
        }
    }
    return 0;
}

void fixedVelocityBC(float* ft, const __global float* rhoList, const float u0, const float v0, const float w0, const int boundary1, const int boundary2, const int boundary3, const int ic, const int i, const int j, const int k, const int nx, const int ny, const int nz, const int corner)
{
    if(boundary1 == 5 || boundary2 == 5 || boundary3 == 5)
    {
        if(corner == 0)
        {
            if(i == 0 && boundary1 == 5)
            {
                float rhow = (ft[0]+ft[3]+ft[4]+ft[5]+ft[6]+ft[15]+ft[16]+ft[17]+ft[18] +2.f*(ft[2]+ft[8]+ft[10]+ft[12]+ft[14]+ft[20]+ft[21]+ft[24]+ft[26]))/(1.f-u0);
                if(u0 < 0.f)
                {
                    int innerID = index1d(1,j,k,nx,ny);
                    rhow = rhoList[innerID];
                }
                ft[1] += 4.f*rhow*u0/9.f;
                // float Nyx = -rhow*v0/3.f +0.5f*(ft[3]+ft[15]+ft[17]-(ft[4]+ft[18]+ft[16]));
                // float Nzx = -rhow*w0/3.f +0.5f*(ft[5]+ft[15]+ft[18]-(ft[6]+ft[17]+ft[16]));
                float Nyx = 0.f;
                float Nzx = 0.f;
                ft[7]  += rhow*(u0+v0)/9.f -Nyx;
                ft[9]  += rhow*(u0-v0)/9.f +Nyx;
                ft[11] += rhow*(u0+w0)/9.f -Nzx;
                ft[13] += rhow*(u0-w0)/9.f +Nzx;
                ft[19] += rhow*(u0+v0+w0)/36.f -Nyx -Nzx;
                ft[22] += rhow*(u0-v0-w0)/36.f +Nyx +Nzx;
                ft[23] += rhow*(u0-v0+w0)/36.f +Nyx -Nzx;
                ft[25] += rhow*(u0+v0-w0)/36.f -Nyx +Nzx;
            }
            else if(i == nx-1 && boundary1 == 5)
            {
                float rhow = (ft[0]+ft[3]+ft[4]+ft[5]+ft[6]+ft[15]+ft[16]+ft[17]+ft[18] +2.f*(ft[1]+ft[7]+ft[9]+ft[11]+ft[13]+ft[19]+ft[22]+ft[23]+ft[25]))/(1.f+u0);
                if(u0 > 0.f)
                {
                    int innerID = index1d(nx-2,j,k,nx,ny);
                    rhow = rhoList[ic];
                }
                ft[2] += -4.f*rhow*u0/9.f;
                // float Nyx = -rhow*v0/3.f +0.5f*(ft[3]+ft[15]+ft[17]-(ft[4]+ft[18]+ft[16]));
                // float Nzx = -rhow*w0/3.f +0.5f*(ft[5]+ft[15]+ft[18]-(ft[6]+ft[17]+ft[16]));
                float Nyx = 0.f;
                float Nzx = 0.f;
                ft[8]  += -rhow*(u0+v0)/9.f +Nyx;
                ft[10] += -rhow*(u0-v0)/9.f -Nyx;
                ft[12] += -rhow*(u0+w0)/9.f +Nzx;
                ft[14] += -rhow*(u0-w0)/9.f -Nzx;
                ft[20] += -rhow*(u0+v0+w0)/36.f -Nyx -Nzx;
                ft[21] += -rhow*(u0-v0-w0)/36.f +Nyx +Nzx;
                ft[24] += -rhow*(u0-v0+w0)/36.f -Nyx +Nzx;
                ft[26] += -rhow*(u0+v0-w0)/36.f +Nyx -Nzx;
            }
            if(j == 0 && boundary2 == 5)
            {
                float rhow = (ft[0]+ft[1]+ft[2]+ft[5]+ft[6]+ft[11]+ft[12]+ft[14]+ft[13] +2.f*(ft[4]+ft[8]+ft[9]+ft[16]+ft[18]+ft[20]+ft[22]+ft[23]+ft[26]))/(1.f-v0);
                if(v0 < 0.f)
                {
                    int innerID = index1d(i,1,k,nx,ny);
                    rhow = rhoList[innerID];
                }
                ft[3] += 4.f*rhow*v0/9.f;
                // float Nxy = -rhow*u0/3.f +0.5f*(ft[1]+ft[11]+ft[13]-(ft[2]+ft[14]+ft[12]));
                // float Nzy = -rhow*w0/3.f +0.5f*(ft[5]+ft[11]+ft[14]-(ft[6]+ft[13]+ft[12]));
                float Nxy = 0.f;
                float Nzy = 0.f;
                ft[7]  += rhow*(v0+u0)/9.f -Nxy;
                ft[10] += rhow*(v0-u0)/9.f +Nxy;
                ft[15] += rhow*(v0+w0)/9.f -Nzy;
                ft[17] += rhow*(v0-w0)/9.f +Nzy;
                ft[19] += rhow*(u0+v0+w0)/36.f -Nxy -Nzy;
                ft[21] += rhow*(-u0+v0+w0)/36.f +Nxy -Nzy;
                ft[24] += rhow*(-u0+v0-w0)/36.f +Nxy +Nzy;
                ft[25] += rhow*(u0+v0-w0)/36.f -Nxy +Nzy;
            }
            else if(j == ny-1 && boundary2 == 5)
            {
                float rhow = (ft[0]+ft[1]+ft[2]+ft[5]+ft[6]+ft[11]+ft[12]+ft[14]+ft[13] +2.f*(ft[3]+ft[7]+ft[10]+ft[15]+ft[17]+ft[19]+ft[21]+ft[24]+ft[25]))/(1.f+v0);
                if(v0 > 0.f)
                {
                    int innerID = index1d(i,ny-2,k,nx,ny);
                    rhow = rhoList[innerID];
                }
                ft[4] += -4.f*rhow*v0/9.f;
                // float Nxy = -rhow*u0/3.f +0.5f*(ft[1]+ft[11]+ft[13]-(ft[2]+ft[14]+ft[12]));
                // float Nzy = -rhow*w0/3.f +0.5f*(ft[5]+ft[11]+ft[14]-(ft[6]+ft[13]+ft[12]));
                float Nxy = 0.f;
                float Nzy = 0.f;
                ft[8]  += -rhow*(v0+u0)/9.f +Nxy;
                ft[9]  += -rhow*(v0-u0)/9.f -Nxy;
                ft[16] += -rhow*(v0+w0)/9.f +Nzy;
                ft[18] += -rhow*(v0-w0)/9.f -Nzy;
                ft[20] += -rhow*(u0+v0+w0)/36.f +Nxy +Nzy;
                ft[22] += -rhow*(-u0+v0+w0)/36.f -Nxy +Nzy;
                ft[23] += -rhow*(-u0+v0-w0)/36.f -Nxy -Nzy;
                ft[26] += -rhow*(u0+v0-w0)/36.f +Nxy -Nzy;
            }
            if(k == 0 && boundary3 == 5)
            {
                float rhow =(ft[0]+ft[1]+ft[2]+ft[3]+ft[4]+ft[7]+ft[8]+ft[9]+ft[10]+2.f*(ft[6]+ft[12]+ft[13]+ft[16]+ft[17]+ft[20]+ft[22]+ft[24]+ft[25]))/(1.f-w0);
                if(w0 < 0.f)
                {
                    int innerID = index1d(i,j,1,nx,ny);
                    rhow = rhoList[innerID];
                }
                ft[5] += 4.f*rhow*w0/9.f;
                // float Nxz = -rhow*u0/3.f +0.5f*(ft[1]+ft[7]+ft[9]-(ft[2]+ft[10]+ft[8]));
                // float Nyz = -rhow*v0/3.f +0.5f*(ft[3]+ft[7]+ft[8]-(ft[4]+ft[9]+ft[8]));
                float Nxz = 0.f;
                float Nyz = 0.f;
                ft[11] += rhow*(w0+u0)/9.f -Nxz;
                ft[14] += rhow*(w0-u0)/9.f +Nxz;
                ft[15] += rhow*(w0+v0)/9.f -Nyz;
                ft[18] += rhow*(w0-v0)/9.f +Nyz;
                ft[19] += rhow*(u0+v0+w0)/36.f -Nxz -Nyz;
                ft[21] += rhow*(-u0+v0+w0)/36.f +Nxz -Nyz;
                ft[23] += rhow*(u0-v0+w0)/36.f -Nxz +Nyz;
                ft[26] += rhow*(-u0-v0+w0)/36.f +Nxz +Nyz;
            }
            else if(k == nz-1 && boundary3 == 5)
            {
                float rhow = (ft[0]+ft[1]+ft[2]+ft[3]+ft[4]+ft[7]+ft[8]+ft[9]+ft[10]+2.f*(ft[5]+ft[11]+ft[14]+ft[15]+ft[18]+ft[19]+ft[21]+ft[23]+ft[26]))/(1.f+w0);
                if(w0 > 0.f)
                {
                    int innerID = index1d(i,j,nz-2,nx,ny);
                    rhow = rhoList[innerID];
                }
                ft[6] += -4.f*rhow*w0/9.f;
                // float Nxz = -rhow*u0/3.f +0.5f*(ft[1]+ft[7]+ft[9]-(ft[2]+ft[10]+ft[8]));
                // float Nyz = -rhow*v0/3.f +0.5f*(ft[3]+ft[7]+ft[8]-(ft[4]+ft[9]+ft[8]));
                float Nxz = 0.f;
                float Nyz = 0.f;
                ft[12] += -rhow*(w0+u0)/9.f +Nxz;
                ft[13] += -rhow*(w0-u0)/9.f -Nxz;
                ft[16] += -rhow*(w0+v0)/9.f +Nyz;
                ft[17] += -rhow*(w0-v0)/9.f -Nyz;
                ft[20] += -rhow*(u0+v0+w0)/9.f +Nxz +Nyz;
                ft[22] += -rhow*(-u0+v0+w0)/9.f -Nxz +Nyz;
                ft[24] += -rhow*(u0-v0+w0)/9.f +Nxz -Nyz;
                ft[25] += -rhow*(-u0-v0+w0)/9.f -Nxz -Nyz;
            }
        }    
    }
}

void fixedDensityBC(float* ft, const float rhow, const __global float* uList, const __global float* vList, const __global float* wList, const float* cx, const float* cy, const float* cz, const float* wt, const int boundary1, const int boundary2, const int boundary3, const int ic, const int i, const int j, const int k, const int nx, const int ny, const int nz, const int corner)
{
    if(boundary1 == 6 || boundary2 == 6 || boundary3 == 6)
    {
        if(corner == 0)
        {
            if(boundary1 == 6)
            {
                int iIn = (i == 0) ? 1 : nx-2;
                int innerID = index1d(iIn,j,k,nx,ny);
                float u0 = uList[ic] +0.5f*(uList[ic] -uList[innerID]);
                float v0 = vList[ic] +0.5f*(vList[ic] -vList[innerID]);
                float w0 = wList[ic] +0.5f*(wList[ic] -wList[innerID]);
                float uSqr =u0*u0+v0*v0+w0*w0;
                
                if(i == 0)
                {
                    int qList[9] = {1,7,9,11,13,19,22,23,25};
                    
                    for(int qid = 0; qid < 9; qid++)
                    {
                        int q = qList[qid];
                        float uDotC = u0*cx[q]+v0*cy[q]+w0*cz[q];
                        ft[q] = -ft[q] +2.f*wt[q]*rhow*(1.f +4.5f*uDotC*uDotC -1.5f*uSqr);
                    }
                }
                if(i == nx-1)
                {
                    int qList[9] = {2,8,10,12,14,20,21,24,26};
                    
                    for(int qid = 0; qid < 9; qid++)
                    {
                        int q = qList[qid];
                        float uDotC = u0*cx[q]+v0*cy[q]+w0*cz[q];
                        ft[q] = -ft[q] +2.f*wt[q]*rhow*(1.f +4.5f*uDotC*uDotC -1.5f*uSqr);
                    }
                }
                
            }
            if(boundary2 == 6)
            {
                int jIn = (j == 0) ? 1 : ny-2;
                int innerID = index1d(i,jIn,k,nx,ny);
                float u0 = uList[ic] +0.5f*(uList[ic] -uList[innerID]);
                float v0 = vList[ic] +0.5f*(vList[ic] -vList[innerID]);
                float w0 = wList[ic] +0.5f*(wList[ic] -wList[innerID]);
                float uSqr =u0*u0+v0*v0+w0*w0;
                if(j == 0)
                {
                    int qList[9] = {3,7,10,15,17,19,21,24,25};                    
                    for(int qid = 0; qid < 9; qid++)
                    {
                        int q = qList[qid];
                        float uDotC = u0*cx[q]+v0*cy[q]+w0*cz[q];
                        ft[q] = -ft[q] +2.f*wt[q]*rhow*(1.f +4.5f*uDotC*uDotC -1.5f*uSqr);
                    }
                }
                if(j == ny-1)
                {
                    int qList[9] = {4,8,9,16,18,20,22,23,26};    
                    for(int qid = 0; qid < 9; qid++)
                    {
                        int q = qList[qid];
                        float uDotC = u0*cx[q]+v0*cy[q]+w0*cz[q];
                        ft[q] = -ft[q] +2.f*wt[q]*rhow*(1.f +4.5f*uDotC*uDotC -1.5f*uSqr);
                    }
                }
            }

            if(boundary3 == 6)
            {
                int kIn = (k == 0) ? 1: nz-2;
                int innerID = index1d(i,j,kIn,nx,ny);
                float u0 = uList[ic] +0.5f*(uList[ic] -uList[innerID]);
                float v0 = vList[ic] +0.5f*(vList[ic] -vList[innerID]);
                float w0 = wList[ic] +0.5f*(wList[ic] -wList[innerID]);
                float uSqr =u0*u0+v0*v0+w0*w0;
                if(k == 0 && boundary3 == 6)
                {
                    int qList[9] = {5,11,14,15,18,19,21,23,26};
                    for(int qid = 0; qid < 9; qid++)
                    {
                        int q = qList[qid];
                        float uDotC = u0*cx[q]+v0*cy[q]+w0*cz[q];
                        ft[q] = -ft[q] +2.f*wt[q]*rhow*(1.f +4.5f*uDotC*uDotC -1.5f*uSqr);
                    }
                }
                if(k == nz-1 && boundary3 == 6)
                {
                    int qList[9] = {6,12,13,16,17,20,22,24,25};
                    for(int qid = 0; qid < 9; qid++)
                    {
                        int q = qList[qid];
                        float uDotC = u0*cx[q]+v0*cy[q]+w0*cz[q];
                        ft[q] = -ft[q] +2.f*wt[q]*rhow*(1.f +4.5f*uDotC*uDotC -1.5f*uSqr);
                    }
                }
            }
        }
    }
}

float calcTauw(const float rhow, const float magUp, const float nu, const float nuEff)
{
    const float y = 1.5f;
    const float kappa = 0.42f;
    const float E = 9.1f;

    float ut = sqrt(nuEff*magUp/y);   

    int iter = 0;
    float err = 1.f;

    if (ut > 1e-18f)
    {
        // printf("nu: %f, nuEff: %f, magUp: %f, ut: %f\n", nu,nuEff,magUp,ut);

        do
        {
            float kUu = min(kappa*magUp/ut,50.f);
            float fkUu = exp(kUu) -1.f -kUu*(1.f +0.5f*kUu);

            float func = -ut*y/nu +magUp/ut +(fkUu -kUu*kUu*kUu/6.f)/E;
            float dfunc = y/nu +magUp/(ut*ut) +(kUu*fkUu/ut)/E;

            float utNew = ut +func/dfunc;
            err = fabs((ut-utNew)/ut);

            ut = utNew;
            // printf("func: %f, dfunc: %f, ut: %f\n", func,dfunc,ut);
        } while (ut > 1e-18f && err > 0.01f && ++iter < 10);
    }
    ut = max(ut, 0.f);
    // float yPlus = y*ut/nu;
    float tauw = rhow*ut*ut;
    // printf("iter: %d, tauw: %f, y+: %f\n", iter,tauw,yPlus);
    return tauw;
}

void wallFunctionBC(float* ft, const __global float* rhoList, const __global float* uList, const __global float* vList, const __global float* wList, const int boundary1, const int boundary2, const int boundary3, const int ic, const int i, const int j, const int k, const int nx, const int ny, const int nz, const int corner, const float nu, const float nuEff, const float Fwx, const float Fwy, const float Fwz) // Stationary wall
{
    if(boundary1 == 7 || boundary2 == 7 || boundary3 == 7)
    {
        if(corner == 0)
        {
            if(i == 0 && boundary1 == 7)
            {
                const float u0 = 0.f;
                float rhow = (ft[0]+ft[3]+ft[4]+ft[5]+ft[6]+ft[15]+ft[16]+ft[17]+ft[18] +2.f*(ft[2]+ft[8]+ft[10]+ft[12]+ft[14]+ft[20]+ft[21]+ft[24]+ft[26]))/(1.f-u0);

                float v0 = 0.f;
                float w0 = 0.f;
                float u = uList[index1d(i+1,j,k,nx,ny)];
                float v = vList[index1d(i+1,j,k,nx,ny)];
                float w = wList[index1d(i+1,j,k,nx,ny)];
                float magUp = sqrt(v*v+w*w);
                
                if(magUp != 0.f)
                {
                    float tauw = calcTauw(rhow,magUp,nu,nuEff);
                    u = uList[index1d(i,j,k,nx,ny)];
                    v = vList[index1d(i,j,k,nx,ny)];
                    w = wList[index1d(i,j,k,nx,ny)];
                    magUp = sqrt(v*v+w*w);
                    if(magUp != 0.f)
                    {
                        float ex = 0.f/magUp;
                        float ey = v/magUp;
                        float ez = w/magUp;

                        float Fuw_x = tauw*ex -Fwx;
                        float Fuw_y = tauw*ey -Fwy;
                        float Fuw_z = tauw*ez -Fwz;

                        v0 = -3.f*Fuw_y/rhow;
                        w0 = -3.f*Fuw_z/rhow;
                    }
                }
                            
                ft[1] += 4.f*rhow*u0/9.f;
                // float Nyx = -rhow*v0/3.f +0.5f*(ft[3]+ft[15]+ft[17]-(ft[4]+ft[18]+ft[16]));
                // float Nzx = -rhow*w0/3.f +0.5f*(ft[5]+ft[15]+ft[18]-(ft[6]+ft[17]+ft[16]));
                float Nyx = 0.f;
                float Nzx = 0.f;
                ft[7]  += rhow*(u0+v0)/9.f -Nyx;
                ft[9]  += rhow*(u0-v0)/9.f +Nyx;
                ft[11] += rhow*(u0+w0)/9.f -Nzx;
                ft[13] += rhow*(u0-w0)/9.f +Nzx;
                ft[19] += rhow*(u0+v0+w0)/36.f -Nyx -Nzx;
                ft[22] += rhow*(u0-v0-w0)/36.f +Nyx +Nzx;
                ft[23] += rhow*(u0-v0+w0)/36.f +Nyx -Nzx;
                ft[25] += rhow*(u0+v0-w0)/36.f -Nyx +Nzx;
            }
            else if(i == nx-1 && boundary1 == 7)
            {
                const float u0 = 0.f;
                float rhow = (ft[0]+ft[3]+ft[4]+ft[5]+ft[6]+ft[15]+ft[16]+ft[17]+ft[18] +2.f*(ft[1]+ft[7]+ft[9]+ft[11]+ft[13]+ft[19]+ft[22]+ft[23]+ft[25]))/(1.f+u0);

                float v0 = 0.f;
                float w0 = 0.f;
                float u = uList[index1d(i-1,j,k,nx,ny)];
                float v = vList[index1d(i-1,j,k,nx,ny)];
                float w = wList[index1d(i-1,j,k,nx,ny)];
                float magUp = sqrt(v*v+w*w);
                
                if(magUp != 0.f)
                {
                    float tauw = calcTauw(rhow,magUp,nu,nuEff);
                    u = uList[index1d(i,j,k,nx,ny)];
                    v = vList[index1d(i,j,k,nx,ny)];
                    w = wList[index1d(i,j,k,nx,ny)];
                    magUp = sqrt(v*v+w*w);
                    if(magUp != 0.f)
                    {
                        float ex = 0.f/magUp;
                        float ey = v/magUp;
                        float ez = w/magUp;

                        float Fuw_x = tauw*ex -Fwx;
                        float Fuw_y = tauw*ey -Fwy;
                        float Fuw_z = tauw*ez -Fwz;

                        v0 = -3.f*Fuw_y/rhow;
                        w0 = -3.f*Fuw_z/rhow;
                    }
                }
                
                ft[2] += -4.f*rhow*u0/9.f;
                // float Nyx = -rhow*v0/3.f +0.5f*(ft[3]+ft[15]+ft[17]-(ft[4]+ft[18]+ft[16]));
                // float Nzx = -rhow*w0/3.f +0.5f*(ft[5]+ft[15]+ft[18]-(ft[6]+ft[17]+ft[16]));
                float Nyx = 0.f;
                float Nzx = 0.f;
                ft[8]  += -rhow*(u0+v0)/9.f +Nyx;
                ft[10] += -rhow*(u0-v0)/9.f -Nyx;
                ft[12] += -rhow*(u0+w0)/9.f +Nzx;
                ft[14] += -rhow*(u0-w0)/9.f -Nzx;
                ft[20] += -rhow*(u0+v0+w0)/36.f -Nyx -Nzx;
                ft[21] += -rhow*(u0-v0-w0)/36.f +Nyx +Nzx;
                ft[24] += -rhow*(u0-v0+w0)/36.f -Nyx +Nzx;
                ft[26] += -rhow*(u0+v0-w0)/36.f +Nyx -Nzx;
            }
            if(j == 0 && boundary2 == 7)
            {
                const float v0 = 0.f;
                float rhow = (ft[0]+ft[1]+ft[2]+ft[5]+ft[6]+ft[11]+ft[12]+ft[14]+ft[13] +2.f*(ft[4]+ft[8]+ft[9]+ft[16]+ft[18]+ft[20]+ft[22]+ft[23]+ft[26]))/(1.f-v0);

                float u0 = 0.f;
                float w0 = 0.f;
                float u = uList[index1d(i,j+1,k,nx,ny)];
                float v = vList[index1d(i,j+1,k,nx,ny)];
                float w = wList[index1d(i,j+1,k,nx,ny)];
                float magUp = sqrt(u*u+w*w);

                if(magUp != 0.f)
                {
                    float tauw = calcTauw(rhow,magUp,nu,nuEff);
                    u = uList[index1d(i,j,k,nx,ny)];
                    v = vList[index1d(i,j,k,nx,ny)];
                    w = wList[index1d(i,j,k,nx,ny)];
                    magUp = sqrt(u*u+w*w);
                    // printf("magUp1: %.12f, u1: %.12f,w1: %.12f\n", magUp, u, w);
                    if(magUp != 0.f)
                    {
                        float ex = u/magUp;
                        float ey = 0.f/magUp;
                        float ez = w/magUp;

                        float Fuw_x = tauw*ex -Fwx;
                        float Fuw_y = tauw*ey -Fwy;
                        float Fuw_z = tauw*ez -Fwz;

                        u0 = -3.f*Fuw_x/rhow;
                        w0 = -3.f*Fuw_z/rhow;
                    }
                }
                
                ft[3] += 4.f*rhow*v0/9.f;
                // float Nxy = -rhow*u0/3.f +0.5f*(ft[1]+ft[11]+ft[13]-(ft[2]+ft[14]+ft[12]));
                // float Nzy = -rhow*w0/3.f +0.5f*(ft[5]+ft[11]+ft[14]-(ft[6]+ft[13]+ft[12]));
                float Nxy = 0.f;
                float Nzy = 0.f;
                ft[7]  += rhow*(v0+u0)/9.f -Nxy;
                ft[10] += rhow*(v0-u0)/9.f +Nxy;
                ft[15] += rhow*(v0+w0)/9.f -Nzy;
                ft[17] += rhow*(v0-w0)/9.f +Nzy;
                ft[19] += rhow*(u0+v0+w0)/36.f -Nxy -Nzy;
                ft[21] += rhow*(-u0+v0+w0)/36.f +Nxy -Nzy;
                ft[24] += rhow*(-u0+v0-w0)/36.f +Nxy +Nzy;
                ft[25] += rhow*(u0+v0-w0)/36.f -Nxy +Nzy;
            }
            else if(j == ny-1 && boundary2 == 7)
            {
                const float v0 = 0.f;
                float rhow = (ft[0]+ft[1]+ft[2]+ft[5]+ft[6]+ft[11]+ft[12]+ft[14]+ft[13] +2.f*(ft[3]+ft[7]+ft[10]+ft[15]+ft[17]+ft[19]+ft[21]+ft[24]+ft[25]))/(1.f+v0);

                float u0 = 0.f;
                float w0 = 0.f;
                float u = uList[index1d(i,j-1,k,nx,ny)];
                float v = vList[index1d(i,j-1,k,nx,ny)];
                float w = wList[index1d(i,j-1,k,nx,ny)];
                float magUp = sqrt(u*u+w*w);

                if(magUp != 0.f)
                {
                    float tauw = calcTauw(rhow,magUp,nu,nuEff);
                    u = uList[index1d(i,j,k,nx,ny)];
                    v = vList[index1d(i,j,k,nx,ny)];
                    w = wList[index1d(i,j,k,nx,ny)];
                    magUp = sqrt(u*u+w*w);
                    if(magUp != 0.f)
                    {
                        float ex = u/magUp;
                        float ey = 0.f/magUp;
                        float ez = w/magUp;

                        float Fuw_x = tauw*ex -Fwx;
                        float Fuw_y = tauw*ey -Fwy;
                        float Fuw_z = tauw*ez -Fwz;

                        u0 = -3.f*Fuw_x/rhow;
                        w0 = -3.f*Fuw_z/rhow;
                    }
                }
                
                ft[4] += -4.f*rhow*v0/9.f;
                // float Nxy = -rhow*u0/3.f +0.5f*(ft[1]+ft[11]+ft[13]-(ft[2]+ft[14]+ft[12]));
                // float Nzy = -rhow*w0/3.f +0.5f*(ft[5]+ft[11]+ft[14]-(ft[6]+ft[13]+ft[12]));
                float Nxy = 0.f;
                float Nzy = 0.f;
                ft[8]  += -rhow*(v0+u0)/9.f +Nxy;
                ft[9]  += -rhow*(v0-u0)/9.f -Nxy;
                ft[16] += -rhow*(v0+w0)/9.f +Nzy;
                ft[18] += -rhow*(v0-w0)/9.f -Nzy;
                ft[20] += -rhow*(u0+v0+w0)/36.f +Nxy +Nzy;
                ft[22] += -rhow*(-u0+v0+w0)/36.f -Nxy +Nzy;
                ft[23] += -rhow*(-u0+v0-w0)/36.f -Nxy -Nzy;
                ft[26] += -rhow*(u0+v0-w0)/36.f +Nxy -Nzy;
            }
            if(k == 0 && boundary3 == 7)
            {
                const float w0 = 0.f;
                float rhow =(ft[0]+ft[1]+ft[2]+ft[3]+ft[4]+ft[7]+ft[8]+ft[9]+ft[10]+2.f*(ft[6]+ft[12]+ft[13]+ft[16]+ft[17]+ft[20]+ft[22]+ft[24]+ft[25]))/(1.f-w0);

                float u0 = 0.f;
                float v0 = 0.f;
                float u = nz-1 != 0 ? uList[index1d(i,j,k+1,nx,ny)] : uList[ic];
                float v = nz-1 != 0 ? vList[index1d(i,j,k+1,nx,ny)] : vList[ic];
                float w = nz-1 != 0 ? wList[index1d(i,j,k+1,nx,ny)] : wList[ic];
                float magUp = sqrt(u*u+v*v);
                
                if(magUp != 0.f)
                {
                    float tauw = calcTauw(rhow,magUp,nu,nuEff);
                    u = uList[index1d(i,j,k,nx,ny)];
                    v = vList[index1d(i,j,k,nx,ny)];
                    w = wList[index1d(i,j,k,nx,ny)];
                    magUp = sqrt(u*u+v*v);
                    if(magUp != 0.f)
                    {
                        float ex = u/magUp;
                        float ey = v/magUp;
                        float ez = 0.f/magUp;

                        float Fuw_x = tauw*ex -Fwx;
                        float Fuw_y = tauw*ey -Fwy;
                        float Fuw_z = tauw*ez -Fwz;

                        u0 = -3.f*Fuw_x/rhow;
                        v0 = -3.f*Fuw_y/rhow;
                    }
                }
                
                ft[5] += 4.f*rhow*w0/9.f;
                // float Nxz = -rhow*u0/3.f +0.5f*(ft[1]+ft[7]+ft[9]-(ft[2]+ft[10]+ft[8]));
                // float Nyz = -rhow*v0/3.f +0.5f*(ft[3]+ft[7]+ft[8]-(ft[4]+ft[9]+ft[8]));
                float Nxz = 0.f;
                float Nyz = 0.f;
                ft[11] += rhow*(w0+u0)/9.f -Nxz;
                ft[14] += rhow*(w0-u0)/9.f +Nxz;
                ft[15] += rhow*(w0+v0)/9.f -Nyz;
                ft[18] += rhow*(w0-v0)/9.f +Nyz;
                ft[19] += rhow*(u0+v0+w0)/36.f -Nxz -Nyz;
                ft[21] += rhow*(-u0+v0+w0)/36.f +Nxz -Nyz;
                ft[23] += rhow*(u0-v0+w0)/36.f -Nxz +Nyz;
                ft[26] += rhow*(-u0-v0+w0)/36.f +Nxz +Nyz;
            }
            else if(k == nz-1 && boundary3 == 7)
            {
                const float w0 = 0.f;
                float rhow = (ft[0]+ft[1]+ft[2]+ft[3]+ft[4]+ft[7]+ft[8]+ft[9]+ft[10]+2.f*(ft[5]+ft[11]+ft[14]+ft[15]+ft[18]+ft[19]+ft[21]+ft[23]+ft[26]))/(1.f+w0);

                float u0 = 0.f;
                float v0 = 0.f;
                float u = nz-1 != 0 ? uList[index1d(i,j,k+1,nx,ny)] : uList[ic];
                float v = nz-1 != 0 ? vList[index1d(i,j,k+1,nx,ny)] : vList[ic];
                float w = nz-1 != 0 ? wList[index1d(i,j,k+1,nx,ny)] : wList[ic];
                float magUp = sqrt(u*u+v*v);
                
                if(magUp != 0.f)
                {
                    float tauw = calcTauw(rhow,magUp,nu,nuEff);
                    u = uList[index1d(i,j,k,nx,ny)];
                    v = vList[index1d(i,j,k,nx,ny)];
                    w = wList[index1d(i,j,k,nx,ny)];
                    magUp = sqrt(u*u+v*v);
                    if(magUp != 0.f)
                    {
                        float ex = u/magUp;
                        float ey = v/magUp;
                        float ez = 0.f/magUp;

                        float Fuw_x = tauw*ex -Fwx;
                        float Fuw_y = tauw*ey -Fwy;
                        float Fuw_z = tauw*ez -Fwz;

                        u0 = -3.f*Fuw_x/rhow;
                        v0 = -3.f*Fuw_y/rhow;
                    }
                }
                
                ft[6] += -4.f*rhow*w0/9.f;
                // float Nxz = -rhow*u0/3.f +0.5f*(ft[1]+ft[7]+ft[9]-(ft[2]+ft[10]+ft[8]));
                // float Nyz = -rhow*v0/3.f +0.5f*(ft[3]+ft[7]+ft[8]-(ft[4]+ft[9]+ft[8]));
                float Nxz = 0.f;
                float Nyz = 0.f;
                ft[12] += -rhow*(w0+u0)/9.f +Nxz;
                ft[13] += -rhow*(w0-u0)/9.f -Nxz;
                ft[16] += -rhow*(w0+v0)/9.f +Nyz;
                ft[17] += -rhow*(w0-v0)/9.f -Nyz;
                ft[20] += -rhow*(u0+v0+w0)/9.f +Nxz +Nyz;
                ft[22] += -rhow*(-u0+v0+w0)/9.f -Nxz +Nyz;
                ft[24] += -rhow*(u0-v0+w0)/9.f +Nxz -Nyz;
                ft[25] += -rhow*(-u0-v0+w0)/9.f -Nxz -Nyz;
            }
        }
    }
}

void internalWallBC(float* ft, const __global float* f, __global float* Fwx, __global float* Fwy, __global float* Fwz, const __global unsigned char* solidList, const __global unsigned char* neiSolidList, const __global float* sdfList, const float sdf, const int* upID, const float omega, const float tauSGS, const __global float* rhoList, const __global float* uList, const __global float* vList, const __global float* wList, const float* cx, const float* cy, const float* cz, const float* wt, const int ic, const int i, const int j, const int k, const int nx, const int ny, const int nz, const int elements)
{
    float u = uList[ic];
    float v = vList[ic];
    float w = wList[ic];
    float rho = rhoList[ic];

    Fwx[ic] = 0.f;
    Fwy[ic] = 0.f;
    Fwz[ic] = 0.f;

    #pragma unroll
    for(int q = 1; q < 27; q++)
    {
        if(neiSolidList[ic] == 1)
        {
            if(solidList[upID[q]] == 1)
            {
                int qbb = reflectQ(q);
                int bbQID = idf(qbb, ic, nx, ny, nz);
                
                // ft[q] = f[bbQID]; // Simple Bounce-Back               

                //- Interpolated Bounce Back (Unstable for hi Reynolds number flows)
                const float sdf0 = sdf;
                const float sdf1 = sdfList[upID[q]];
                const float qf = fabs(sdf0)/(fabs(sdf0)+fabs(sdf1));

                float tau = 1.f/omega;
                float omegaEff = 1.f/(tau +tauSGS);

                int upQID = idf(q, upID[qbb], nx, ny, nz);
                int upQBBID = idf(qbb, upID[qbb], nx, ny, nz);

                float uSqr =u*u+v*v+w*w;

                if(qf <= 0.5f)
                {
                    // ft[q] = (1.f -2.f*qf)*ft[qbb] +(qf*f[bbQID])*2.f; // Bouzidi et al.'s Interpolated Bounce-Back
                    // ft[q] = (1.f -2.f*qf)*f[upQBBID] +(qf*f[bbQID])*2.f; // Bouzidi et al.'s Interpolated Bounce-Back (local)

                    //- Filippova & Hanel's Interpolated Bounce-Back (FHIBB) (physically local)
                    // float uDotC = -u*cx[q]-v*cy[q]-w*cz[q];
                    // float feq = (rho+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q];
                    // float chi = omegaEff*(2.f*qf -1.f)/(1.f-omegaEff);

                    //- Mei, Luo and Shyy's modification for FHIBB
                    float uDotC = -u*cx[q]-v*cy[q]-w*cz[q];
                    float ufDotC = -uList[upID[qbb]]*cx[q]-vList[upID[qbb]]*cy[q]-wList[upID[qbb]]*cz[q];
                    float feq = (rho+3.0f*ufDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q];
                    float chi = omegaEff*(2.f*qf -1.f)/(1.f-2.f*omegaEff);

                    ft[q] = (1.f -chi)*f[bbQID] +chi*feq;
                }
                else
                {
                    // ft[q] = (1.f -0.5f/qf)*ft[upQID] +(0.5f/qf)*f[bbQID]; // Bouzidi et al.'s Interpolated Bounce-Back
                    // ft[q] = (1.f -0.5f/qf)*f[q] +(0.5f/qf)*f[bbQID]; // Bouzidi et al.'s Interpolated Bounce-Back (local)

                    // - Filippova & Hanel's Interpolated Bounce-Back (physically local)
                    float uDotC = -u*cx[q]-v*cy[q]-w*cz[q];
                    float feq = (rho+3.0f*(1.f -1.f/qf)*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q];
                    float chi = omegaEff*(2.f*qf -1.f);

                    ft[q] = (1.f -chi)*f[bbQID] +chi*feq; 
                }
                //--

                Fwx[ic] += -(f[bbQID] + ft[q])*cx[q];
                Fwy[ic] += -(f[bbQID] + ft[q])*cy[q];
                Fwz[ic] += -(f[bbQID] + ft[q])*cz[q];
            }
        }
    }
}

// Geier et al., Comput. Math. Appl. (2015), Appendix F
void outflowBC(float* ft, const __global float* f, const float u, const float v, const float w, const int boundary1, const int boundary2, const int boundary3, const int ic, const int i, const int j, const int k, const int nx, const int ny, const int nz)
{
    if(boundary1 == 3 || boundary2 == 3 || boundary3 == 3)
    {
        #pragma unroll
        for(int q = 1; q < 27; q++)
        {
            int i = ic2i(ic,nx,ny);
            int j = ic2j(ic,nx,ny);
            int k = ic2k(ic,nx,ny);

            if(i == 0)
            {
                if(boundary1 == 3)
                {
                    if(q == 1 || q == 7 || q == 9 || q == 11 || q == 13 || q == 19 || q == 22 || q == 23 || q == 25)
                    {
                        int innerID = index1d(1,j,k,nx,ny);
                        int qinic = idf(q, innerID, nx, ny, nz);
                        int qic = idf(q, ic, nx, ny, nz);
                        ft[q] = (1.f/sqrt(3.f) -u)*f[qinic] +(1.f -(1.f/sqrt(3.f) -u))*f[qic];
                    }
                }
            }
            if(i == nx-1)
            {
                if(boundary1 == 3)
                {
                    if(q == 2 || q == 8 || q == 10 || q == 12 || q == 14 || q == 20 || q == 21 || q == 24 || q == 26)
                    {
                        int innerID = index1d(nx-2,j,k,nx,ny);
                        int qinic = idf(q, innerID, nx, ny, nz);
                        int qic = idf(q, ic, nx, ny, nz);
                        ft[q] = (1.f/sqrt(3.f) -u)*f[qinic] +(1.f -(1.f/sqrt(3.f) -u))*f[qic];
                    }
                }
            }
            if(j == 0)
            {
                if(boundary2 == 3)
                {
                    if(q == 3 || q == 7 || q == 10 || q == 15 || q == 17 || q == 19 || q == 21 || q == 24 || q == 25)
                    {
                        int innerID = index1d(i,1,k,nx,ny);
                        int qinic = idf(q, innerID, nx, ny, nz);
                        int qic = idf(q, ic, nx, ny, nz);
                        ft[q] = (1.f/sqrt(3.f) -v)*f[qinic] +(1.f -(1.f/sqrt(3.f) -v))*f[qic];
                    }
                }
            }
            if(j == ny-1)
            {
                if(boundary2 == 3)
                {
                    if(q == 4 || q == 8 || q == 9 || q == 16 || q == 18 || q == 20|| q == 22 || q == 23 || q == 26)
                    {
                        int innerID = index1d(i,ny-2,k,nx,ny);
                        int qinic = idf(q, innerID, nx, ny, nz);
                        int qic = idf(q, ic, nx, ny, nz);
                        ft[q] = (1.f/sqrt(3.f) -v)*f[qinic] +(1.f -(1.f/sqrt(3.f) -v))*f[qic];
                    }
                }
            }
            if(k == 0)
            {
                if(boundary3 == 3)
                {
                    if(q == 5 || q == 11 || q == 14 || q == 15 || q == 18 || q == 19 || q == 21 || q == 23 || q == 26)
                    {
                        int innerID = index1d(i,j,1,nx,ny);
                        int qinic = idf(q, innerID, nx, ny, nz);
                        int qic = idf(q, ic, nx, ny, nz);
                        ft[q] = (1.f/sqrt(3.f) -w)*f[qinic] +(1.f -(1.f/sqrt(3.f) -w))*f[qic];
                    }
                }
            }
            if(k == nz-1)
            {
                if(boundary3 == 3)
                {
                    if(q == 6 || q == 12 || q == 13 || q == 16 || q == 17 || q == 20 || q == 22 || q == 24 || q == 25)
                    {
                        int innerID = index1d(i,j,nz-2,nx,ny);
                        int qinic = idf(q, innerID, nx, ny, nz);
                        int qic = idf(q, ic, nx, ny, nz);
                        ft[q] = (1.f/sqrt(3.f) -w)*f[qinic] +(1.f -(1.f/sqrt(3.f) -w))*f[qic];
                    }
                }
            }
        }
    }    
}

void gradU(float* dudx, float* dudy, float* dudz, float* dvdx, float* dvdy, float* dvdz, float* dwdx, float* dwdy, float* dwdz, const int i, const int j, const int k, const int nx, const int ny, const int nz, const __global float* uList, const __global float* vList, const __global float* wList)
{
    *dudx = (i > 0 && i < nx-1) ? 0.5f*(uList[index1d(i+1,j,k,nx,ny)] -uList[index1d(i-1,j,k,nx,ny)]) : 
                    (i == 0) ? uList[index1d(1,j,k,nx,ny)]-uList[index1d(0,j,k,nx,ny)] :
                    uList[index1d(nx-1,j,k,nx,ny)]-uList[index1d(nx-2,j,k,nx,ny)];
    *dudy = (j > 0 && j < ny-1) ? 0.5f*(uList[index1d(i,j+1,k,nx,ny)] -uList[index1d(i,j-1,k,nx,ny)]) :
                    (j == 0) ? uList[index1d(i,1,k,nx,ny)]-uList[index1d(i,0,k,nx,ny)] :
                    uList[index1d(i,ny-1,k,nx,ny)]-uList[index1d(i,ny-2,k,nx,ny)];
    *dudz = (k > 0 && k < nz-1) ? 0.5f*(uList[index1d(i,j,k+1,nx,ny)] -uList[index1d(i,j,k-1,nx,ny)]) : 
                    (nz-1 == 0) ? 0.f:
                    (k == 0) ? uList[index1d(i,j,1,nx,ny)]-uList[index1d(i,j,1,nx,ny)] :
                    uList[index1d(i,j,nz-1,nx,ny)]-uList[index1d(i,j,nz-2,nx,ny)];
    *dvdx = (i > 0 && i < nx-1) ? 0.5f*(vList[index1d(i+1,j,k,nx,ny)] -vList[index1d(i-1,j,k,nx,ny)]) : 
                    (i == 0) ? vList[index1d(1,j,k,nx,ny)]-vList[index1d(0,j,k,nx,ny)] :
                    vList[index1d(nx-1,j,k,nx,ny)]-vList[index1d(nx-2,j,k,nx,ny)];
    *dvdy = (j > 0 && j < ny-1) ? 0.5f*(vList[index1d(i,j+1,k,nx,ny)] -vList[index1d(i,j-1,k,nx,ny)]) : 
                    (j == 0) ? vList[index1d(i,1,k,nx,ny)]-vList[index1d(i,0,k,nx,ny)] :
                    vList[index1d(i,ny-1,k,nx,ny)]-vList[index1d(i,ny-2,k,nx,ny)];;
    *dvdz = (k > 0 && k < nz-1) ? 0.5f*(vList[index1d(i,j,k+1,nx,ny)] -vList[index1d(i,j,k-1,nx,ny)]) : 
                    (nz-1 == 0) ? 0.f:
                    (k == 0) ? vList[index1d(i,j,1,nx,ny)]-vList[index1d(i,j,1,nx,ny)] :
                    vList[index1d(i,j,nz-1,nx,ny)]-vList[index1d(i,j,nz-2,nx,ny)];
    *dwdx = (i > 0 && i < nx-1) ? 0.5f*(wList[index1d(i+1,j,k,nx,ny)] -wList[index1d(i-1,j,k,nx,ny)]) : 
                    (i == 0) ? wList[index1d(1,j,k,nx,ny)]-wList[index1d(0,j,k,nx,ny)] :
                    wList[index1d(nx-1,j,k,nx,ny)]-wList[index1d(nx-2,j,k,nx,ny)];
    *dwdy = (j > 0 && j < ny-1) ? 0.5f*(wList[index1d(i,j+1,k,nx,ny)] -wList[index1d(i,j-1,k,nx,ny)]) : 
                    (j == 0) ? wList[index1d(i,1,k,nx,ny)]-wList[index1d(i,0,k,nx,ny)] :
                    wList[index1d(i,ny-1,k,nx,ny)]-wList[index1d(i,ny-2,k,nx,ny)];
    *dwdz = (k > 0 && k < nz-1) ? 0.5f*(wList[index1d(i,j,k+1,nx,ny)] -wList[index1d(i,j,k-1,nx,ny)]) : 
                    (nz-1 == 0) ? 0.f:
                    (k == 0) ? wList[index1d(i,j,1,nx,ny)]-wList[index1d(i,j,1,nx,ny)] :
                    wList[index1d(i,j,nz-1,nx,ny)]-wList[index1d(i,j,nz-2,nx,ny)];
}

float smagorinskyTauSGS(const float* ft, const float rho, const float u, const float v, const float w, const float omega, const float LES, const float sdf, const float* cx, const float* cy, const float* cz, const float* wt, const float Cs, const int damping)
{
    // float Cs = 0.1f;// 0.1--0.2
    // float Cs = 0.2f;// 0.1--0.2
    // float Cs = 0.33f;// 0.1--0.2

    float tau = 1.f/omega;
            
    float PIxx = 0.f;
    float PIxy = 0.f;
    float PIxz = 0.f;
    float PIyy = 0.f;
    float PIyz = 0.f;
    float PIzz = 0.f;

    #pragma unroll
    for(int q = 0; q < 27; q++)
    {
        float uSqr =u*u+v*v+w*w;
        float uDotC = u*cx[q]+v*cy[q]+w*cz[q];
        // float feq = (1.0f+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rho;
        float feq = (rho+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q];
        
        PIxx += cx[q]*cx[q]*(ft[q] -feq);
        PIxy += cx[q]*cy[q]*(ft[q] -feq);
        PIxz += cx[q]*cz[q]*(ft[q] -feq);
        PIyy += cy[q]*cy[q]*(ft[q] -feq);
        PIyz += cy[q]*cz[q]*(ft[q] -feq);
        PIzz += cz[q]*cz[q]*(ft[q] -feq);
    }
    float sqrtPIPI = sqrt(PIxx*PIxx+PIyy*PIyy+PIzz*PIzz+2.f*(PIxy*PIxy+PIxz*PIxz+PIyz*PIyz));
    
    float tauSGS = LES*0.5f*(-tau +sqrt(tau*tau +18.f*sqrt(2.f)*Cs*Cs*sqrtPIPI/rho));
    // tauSGS[ic] = LES*0.5f*(-tau +sqrt(tau*tau +18.f*sqrt(2.f)*Cs*Cs*sqrtPIPI/rho));
    // tauSGS[ic] = 3.f*(Cs*Cs)*sqrt(2.f)*sqrtPIPI*0.5f/rho*3.0f/tau;

    //-- Damping of nuSGS (tauSGS) near wall
    if(damping == 1)
    {
        const float y = sdf;
        if(y != 10000.f && y > 0.f)
        {
            const float kappa = 0.41f;
            tauSGS *= pow(min(1.f, kappa*y/Cs),2.f); // Treatment in Fluent

            // //-- van Driest damping function
            // const float Aplus = 26.f;
            // const float Cdelta = 0.158f;
            // const float nu = (tau -0.5f)/3.f;
            // const float E = 9.8f;

            // float yPlus = 1.f;
            // float yPlusN = yPlus;
            // const float epsilon = 0.001f;
            // int iNew = 0;
            // do
            // {
            //     float U = sqrt(u*u+v*v+w*w);
            //     yPlusN = yPlus;
            //     yPlus = ((kappa*U*y/nu) +yPlus)/(1.f +log(E*yPlus));
            //     iNew++;
            // }while(fabs(yPlus -yPlusN) > epsilon && iNew < 10);
            
            // tauSGS[ic] *= pow(min(1.f, (kappa/Cs)*((1.f +1e-10f)-exp(-yPlus/Aplus))*y),2.f);
            // //--
        }
    }
    return tauSGS;
}

float vremanTauSGS(const int i, const int j, const int k, const int nx, const int ny, const int nz, const __global float* uList, const __global float* vList, const __global float* wList, const float LES)
{   
    float Cs = 0.1f;// 0.1--0.2
    // float Cs = 0.2f;// 0.1--0.2
    // float Cs = 0.33f;// 0.1--0.2

    float Cv = 2.5f*Cs*Cs;

    float dudx;
    float dudy;
    float dudz;
    float dvdx;
    float dvdy;
    float dvdz;
    float dwdx;
    float dwdy;
    float dwdz;
    gradU(&dudx,&dudy,&dudz,&dvdx,&dvdy,&dvdz,&dwdx,&dwdy,&dwdz, i, j, k, nx, ny, nz, uList, vList, wList);

    float b11 = dudx*dudx+dvdx*dvdx+dwdx*dwdx;
    float b12 = dudx*dudy+dvdx*dvdy+dwdx*dwdy;
    float b13 = dudx*dudz+dvdx*dvdz+dwdx*dwdz;
    float b22 = dudy*dudy+dvdy*dvdy+dwdy*dwdy;
    float b23 = dudy*dudz+dvdy*dvdz+dwdy*dwdz;
    float b33 = dudz*dudz+dvdz*dvdz+dwdz*dwdz;

    float abeta = dudx*dudx+dvdx*dvdx+dwdx*dwdx 
                    +dudy*dudy+dvdy*dvdy+dwdy*dwdy
                    +dudz*dudz+dvdz*dvdz+dwdz*dwdz;
    
    float bbeta = b11*b22 -b12*b12 +b11*b33 -b13*b13 +b22*b33 -b23*b23;

    float tauSGS = (bbeta > 0.f && abeta > 0.f) ? 3.f*LES*(Cv*sqrt(bbeta/abeta))*3.f : 0.f;

    return tauSGS;
}

float waleTauSGS(const int i, const int j, const int k, const int nx, const int ny, const int nz, const __global float* uList, const __global float* vList, const __global float* wList, const float LES)
{   
    // float Cs = 0.1f;// 0.1--0.2
    // float Cs = 0.2f;// 0.1--0.2
    // float Cs = 0.33f;// 0.1--0.2

    float Cw = 0.325f;

    float dudx;
    float dudy;
    float dudz;
    float dvdx;
    float dvdy;
    float dvdz;
    float dwdx;
    float dwdy;
    float dwdz;

    gradU(&dudx,&dudy,&dudz,&dvdx,&dvdy,&dvdz,&dwdx,&dwdy,&dwdz, i, j, k, nx, ny, nz, uList, vList, wList);

    float S11 = dudx;
    float S12 = 0.5f*(dudy+dvdx);
    float S13 = 0.5f*(dudz+dwdx);
    float S21 = S12;
    float S22 = dvdy;
    float S23 = 0.5f*(dvdz+dwdy);
    float S31 = S13;
    float S32 = S23;
    float S33 = dwdz;

    float O11 = 0.f;
    float O12 = 0.5f*(dudy-dvdx);
    float O13 = 0.5f*(dudz-dwdx);
    float O21 = O12;
    float O22 = 0.f;
    float O23 = 0.5f*(dvdz-dwdy);
    float O31 = O13;
    float O32 = O23;
    float O33 = 0.f;
    
    float SS = S11*S11+S12*S12+S13*S13
                +S21*S21+S22*S22+S23*S23
                +S31*S31+S32*S32+S33*S33;

    float OO = O11*O11+O12*O12+O13*O13
                +O21*O21+O22*O22+O23*O23
                +O31*O31+O32*O32+O33*O33;

    float gradUdotGradU11 = dudx*dudx+dvdx*dudy+dwdx*dudz;
    float gradUdotGradU12 = dudx*dvdx+dvdx*dvdy+dwdx*dvdz;
    float gradUdotGradU13 = dudx*dwdx+dvdx*dwdy+dwdx*dwdz;
    float gradUdotGradU21 = dudy*dudx+dvdy*dudy+dwdy*dudz;
    float gradUdotGradU22 = dudy*dvdx+dvdy*dvdy+dwdy*dvdz;
    float gradUdotGradU23 = dudy*dwdx+dvdy*dwdy+dwdy*dwdz;
    float gradUdotGradU31 = dudz*dudx+dvdz*dudy+dwdz*dudz;
    float gradUdotGradU32 = dudz*dvdx+dvdz*dvdy+dwdz*dvdz;
    float gradUdotGradU33 = dudz*dwdx+dvdz*dwdy+dwdz*dwdz;
    float sqrGradU = dudx*dudx+dudy*dvdx+dudz*dwdx
                        +dvdx*dudy+dvdy*dvdy+dvdz*dwdy
                        +dwdx*dudz+dwdy*dvdz+dwdz*dwdz;

    float Sd11 = gradUdotGradU11 -sqrGradU/3.f;
    float Sd12 = 0.5f*(gradUdotGradU12 +gradUdotGradU21);
    float Sd13 = 0.5f*(gradUdotGradU13 +gradUdotGradU31);
    float Sd21 = Sd12;
    float Sd22 = gradUdotGradU22 -sqrGradU/3.f;
    float Sd23 = 0.5f*(gradUdotGradU23 +gradUdotGradU32);
    float Sd31 = Sd13;
    float Sd32 = Sd23;
    float Sd33 = gradUdotGradU33 -sqrGradU/3.f;

    float SdSd = Sd11*Sd11+Sd12*Sd12+Sd13*Sd13
                +Sd21*Sd21+Sd22*Sd22+Sd23*Sd23
                +Sd31*Sd31+Sd32*Sd32+Sd33*Sd33;

    
    float tauSGS = (SS > 1e-24f && SdSd > 1e-24f ) ? 3.f*(Cw*Cw*pow(SdSd,1.5f)/(pow(SS,2.5f)+pow(SdSd,1.25f)))*3.f : 0.f;

    return tauSGS;
}

float CSsmagorinskyTauSGS(const float* ft, const float rho, const float u, const float v, const float w, const float omega, const float LES, const float sdf, const float* cx, const float* cy, const float* cz, const float* wt, const int i, const int j, const int k, const int nx, const int ny, const int nz, const __global float* uList, const __global float* vList, const __global float* wList)
{   
    // float Cs = 0.1f;// 0.1--0.2
    // float Cs = 0.2f;// 0.1--0.2
    // float Cs = 0.33f;// 0.1--0.2

    float C1 = 1.f/22.f;

    float dudx;
    float dudy;
    float dudz;
    float dvdx;
    float dvdy;
    float dvdz;
    float dwdx;
    float dwdy;
    float dwdz;

    gradU(&dudx,&dudy,&dudz,&dvdx,&dvdy,&dvdz,&dwdx,&dwdy,&dwdz, i, j, k, nx, ny, nz, uList, vList, wList);

    float S11 = dudx;
    float S12 = 0.5f*(dudy+dvdx);
    float S13 = 0.5f*(dudz+dwdx);
    float S21 = S12;
    float S22 = dvdy;
    float S23 = 0.5f*(dvdz+dwdy);
    float S31 = S13;
    float S32 = S23;
    float S33 = dwdz;

    float O11 = 0.f;
    float O12 = 0.5f*(dudy-dvdx);
    float O13 = 0.5f*(dudz-dwdx);
    float O21 = O12;
    float O22 = 0.f;
    float O23 = 0.5f*(dvdz-dwdy);
    float O31 = O13;
    float O32 = O23;
    float O33 = 0.f;
    
    float SS = S11*S11+S12*S12+S13*S13
                +S21*S21+S22*S22+S23*S23
                +S31*S31+S32*S32+S33*S33;

    float OO = O11*O11+O12*O12+O13*O13
                +O21*O21+O22*O22+O23*O23
                +O31*O31+O32*O32+O33*O33;

    float E = 0.5f*(OO+SS);
    float Q = 0.5f*(OO-SS);
    float Fcs = E > 0.f ? Q/E : 0.f;
    float Cs = sqrt(C1*pow(fabs(Fcs),1.5f)*(1.f -Fcs));
    
    float tauSGS = smagorinskyTauSGS(ft, rho, u, v, w, omega, LES, sdf, cx, cy, cz, wt, Cs, 0);

    return tauSGS;
}

float tauSpongeZone(const float spzWidth, const __global int* boundary1List, const __global int* boundary2List, const __global int* boundary3List, const float omega, const float tauSGS, const int i, const int j, const int k, const int nx, const int ny, const int nz)
{
    float tau = 1.f/omega;
    int icX0 = index1d(0,ny/2,nz/2,nx,ny);
    int icXE = index1d(nx-1,ny/2,nz/2,nx,ny);
    int icY0 = index1d(nx/2,0,nz/2,nx,ny);
    int icYE = index1d(nx/2,ny-1,nz/2,nx,ny);
    int icZ0 = index1d(nx/2,ny/2,0,nx,ny);
    int icZE = index1d(nx/2,ny/2,nz-1,nx,ny);
    
    if(boundary1List[icX0] == 3 || boundary1List[icX0] == 6)
    {
        if(i < nx*spzWidth)
        {
            const float fx = i/(nx*spzWidth);
            tau = (1.f -fx)*(1.f -tauSGS) +fx*tau;
        }
    }
    if(boundary1List[icXE] == 3 || boundary1List[icXE] == 6)
    {
        if(i > nx*(1.f -spzWidth))
        {
            const float fx = (nx-i)/(nx*spzWidth);
            tau = (1.f -fx)*(1.f -tauSGS) +fx*tau;
        }
    }
    if(boundary2List[icY0] == 3 || boundary2List[icY0] == 6)
    {
        if(j < ny*spzWidth)
        {
            const float fx = j/(ny*spzWidth);
            tau = (1.f -fx)*(1.f -tauSGS) +fx*tau;
        }
    }
    if(boundary2List[icYE] == 3 || boundary2List[icYE] == 6)
    {
        if(j > ny*(1.f -spzWidth))
        {
            const float fx = (ny-j)/(ny*spzWidth);
            tau = (1.f -fx)*(1.f -tauSGS) +fx*tau;
        }
    }
    if(boundary3List[icZ0] == 3 || boundary3List[icZ0] == 6)
    {
        if(k < nz*spzWidth)
        {
            const float fx = k/(nz*spzWidth);
            tau = (1.f -fx)*(1.f -tauSGS) +fx*tau;
        }
    }
    if(boundary3List[icZE] == 3 || boundary3List[icZE] == 6)
    {
        if(k > nz*(1.f -spzWidth))
        {
            const float fx = (nz-k)/(nz*spzWidth);
            tau = (1.f -fx)*(1.f -tauSGS) +fx*tau;
        }
    }
    return tau;    
}

float collisionBGK(float* fTmp, const float* ft, const float rho, const float u, const float v, const float w, const float tau, const float tauSGS, const float dpdx, const float* cx, const float* cy, const float* cz, const float* wt, const int ic, const int elements)
{
    const float omegaEff = 1.f/(tau +tauSGS);
    #pragma unroll
    for(int q = 0; q < 27; q++)
    {
        float uSqr =u*u+v*v+w*w;
        float uDotC = u*cx[q]+v*cy[q]+w*cz[q];
        float feq = (1.0f+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rho;

        int qic = q*elements +ic;

        float S = rho*wt[q]*(1.f -0.5f*omegaEff)*(3.0f*cx[q] +9.0f*cx[q]*uDotC -3.0f*u)*dpdx;
        fTmp[qic] = (1.0f -omegaEff)*ft[q] + omegaEff *feq +S;
    }
}


float collisionCumulantGeier(__global float* fTmp, const float* ft, const float rho, float u, float v, float w, const float tau, const float tauSGS, const float omegaB, const float clim, const float dpdx, const float* cx, const float* cy, const float* cz, const float* wt, const int ic, const int elements)
{
    const float sqrCs = 1.f/3.f;
    const float quadCs = 1.f/9.f;
    const float omegaEff = 1.f/(tau +tauSGS);

    float K000 = rho;
    float K100 = 0.f;
    float K010 = 0.f;
    float K001 = 0.f;

    float K200 = 0.f;
    float K020 = 0.f;
    float K002 = 0.f;
    float K110 = 0.f;
    float K101 = 0.f;
    float K011 = 0.f;
    
    float K210 = 0.f;
    float K201 = 0.f;
    float K021 = 0.f;
    float K120 = 0.f;
    float K102 = 0.f;
    float K012 = 0.f;
    float K111 = 0.f;

    float K220 = 0.f;
    float K202 = 0.f;
    float K022 = 0.f;
    float K211 = 0.f;  
    float K121 = 0.f;
    float K112 = 0.f;

    float K221 = 0.f;
    float K212 = 0.f;
    float K122 = 0.f;

    float K222 = 0.f;


    #pragma unroll
    for(int q = 0; q < 27; q++)
    {
        float cxq = cx[q] -u;
        float cyq = cy[q] -v;
        float czq = cz[q] -w;

        K100 += cxq*ft[q];
        K010 += cyq*ft[q];
        K001 += czq*ft[q];

        K200 += cxq*cxq*ft[q];
        K020 += cyq*cyq*ft[q];
        K002 += czq*czq*ft[q];
        K110 += cxq*cyq*ft[q];
        K101 += cxq*czq*ft[q];
        K011 += cyq*czq*ft[q];
        
        K210 += cxq*cxq*cyq*ft[q];
        K201 += cxq*cxq*czq*ft[q];
        K021 += cyq*cyq*czq*ft[q];
        K120 += cxq*cyq*cyq*ft[q];
        K102 += cxq*czq*czq*ft[q];
        K012 += cyq*czq*czq*ft[q];
        K111 += cxq*cyq*czq*ft[q];

        K220 += cxq*cxq*cyq*cyq*ft[q];
        K202 += cxq*cxq*czq*czq*ft[q];
        K022 += cyq*cyq*czq*czq*ft[q];

        K211 += cxq*cxq*cyq*czq*ft[q];
        K121 += cxq*cyq*cyq*czq*ft[q];
        K112 += cxq*cyq*czq*czq*ft[q];

        K221 += cxq*cxq*cyq*cyq*czq*ft[q];
        K212 += cxq*cxq*cyq*czq*czq*ft[q];
        K122 += cxq*cyq*cyq*czq*czq*ft[q];

        K222 += cxq*cxq*cyq*cyq*czq*czq*ft[q];
    }   

    float invRho = 1.f/rho;

    float C200 = K200;
    float C020 = K020;
    float C002 = K002;
    float C110 = K110;
    float C101 = K101;
    float C011 = K011;
    float C210 = K210;
    float C201 = K201;
    float C021 = K021;
    float C120 = K120;
    float C102 = K102;
    float C012 = K012;
    float C111 = K111;
    float C211 = K211 -(K200*K011 +2.f*K110*K101)*invRho;
    float C121 = K121 -(K020*K101 +2.f*K110*K011)*invRho;
    float C112 = K112 -(K002*K110 +2.f*K101*K011)*invRho;
    float C220 = K220 -(K200*K020 +2.f*K110*K110)*invRho;
    float C202 = K202 -(K002*K200 +2.f*K101*K101)*invRho;
    float C022 = K022 -(K020*K002 +2.f*K011*K011)*invRho;
    float C122 = K122 -(K002*K120 +K020*K102 +4.f*K011*K111 +2.f*(K101*K021 +K110*K012))*invRho;
    float C212 = K212 -(K200*K012 +K002*K210 +4.f*K101*K111 +2.f*(K110*K102 +K011*K201))*invRho;
    float C221 = K221 -(K020*K201 +K200*K021 +4.f*K110*K111 +2.f*(K011*K210 +K101*K120))*invRho;
    float C222 = K222 -(4.f*K111*K111 +K200*K022 +K020*K202 +K002*K220 +4.f*(K011*K211 +K101*K121 +K110*K112) +2.f*(K120*K102 +K210*K012 +K201*K021))*invRho +(16.f*K110*K101*K011 +4.f*(K101*K101*K020 +K011*K011*K200 +K110*K110*K002) +2.f*K200*K020*K002)*invRho*invRho;


    float omega1 = omegaEff;
    float omega2 = omegaB;

    float omega3 = 8.f*(omega1 -2.f)*(omega2*(3.f*omega1 -1.f) -5.f*omega1)/(8.f*(5.f -2.f*omega1)*omega1 +omega2*(8.f +omega1*(9.f*omega1 -26.f)));
    float omega4 = 8.f*(omega1 -2.f)*(omega1 +omega2*(3.f*omega1 -7.f))/(omega2*(56.f -42.f*omega1 +9.f*omega1*omega1) -8.f*omega1);
    float omega5 = 24.f*(omega1 -2.f)*(4.f*omega1*omega1 +omega1*omega2*(18.f -13.f*omega1) +omega2*omega2*(2.f +omega1*(6.f*omega1 -11.f)))/(16.f*omega1*omega1*(omega1 -6.f) -2.f*omega1*omega2*(216.f +5.f*omega1*(9.f*omega1 -46.0)) +omega2*omega2*(omega1*(3.f*omega1 -10.f)*(15.f*omega1 -28.f) -48.f));
    float omega3alim = fabs(C120 +C102) > 0.f ? omega3 +(1.f -omega3)*fabs(C120 +C102)/(fabs(C120 +C102) +clim*rho) : 1.f;
    float omega3blim = fabs(C210 +C012) > 0.f ? omega3 +(1.f -omega3)*fabs(C210 +C012)/(fabs(C210 +C012) +clim*rho) : 1.f;
    float omega3clim = fabs(C201 +C021) > 0.f ? omega3 +(1.f -omega3)*fabs(C201 +C021)/(fabs(C201 +C021) +clim*rho) : 1.f;
    float omega4alim = fabs(C120 -C102) > 0.f ? omega4 +(1.f -omega4)*fabs(C120 -C102)/(fabs(C120 -C102) +clim*rho) : 1.f;
    float omega4blim = fabs(C210 -C012) > 0.f ? omega4 +(1.f -omega4)*fabs(C210 -C012)/(fabs(C210 -C012) +clim*rho) : 1.f;
    float omega4clim = fabs(C201 -C021) > 0.f ? omega4 +(1.f -omega4)*fabs(C201 -C021)/(fabs(C201 -C021) +clim*rho) : 1.f;
    float omega5lim = fabs(C111) > 0.f ? omega5 +(1.f -omega5)*fabs(C111)/(fabs(C111) +clim*rho) : 1.f;

    float omega6 = 1.f;
    float omega7 = 1.f;
    float omega8 = 1.f;
    float omega9 = 1.f;
    float omega10 = 1.f;

    float A = (4.f*omega1*omega1 +2.f*omega1*omega2*(omega1 -6.f)+omega2*omega2*(omega1*(10.f -3.f*omega1) -4.f))/((omega1 -omega2)*(omega2*(2.f +3.f*omega1) -8.f*omega1));
    float B = (4.f*omega1*omega2*(9.f*omega1 -16.f) -4.f*omega1*omega1 -2.f*omega2*omega2*(2.f +9.f*omega1*(omega1 -2.f)))/(3.f*(omega1 -omega2)*(omega2*(2.f +3.f*omega1)- 8.f*omega1));

    A = max(A,0.f);
    B = max(B,0.f);

    float Cs110 = (1.f -omega1)*C110;
    float Cs101 = (1.f -omega1)*C101;
    float Cs011 = (1.f -omega1)*C011;

    float Dxu = -0.5f*(omega1*(2.f*C200 -C020 -C002) +omega2*(C200 +C020 +C002 -K000))*invRho;
    float Dyv = Dxu +1.5f*omega1*(C200 -C020)*invRho;
    float Dzw = Dxu +1.5f*omega1*(C200 -C002)*invRho;
    float Dxv_Dyu = -3.f*omega1*C110*invRho;
    float Dxw_Dzu = -3.f*omega1*C101*invRho;
    float Dyw_Dzv = -3.f*omega1*C011*invRho;

    float X1 = (1.f -omega1)*(C200 -C020) -3.f*rho*(1.f -0.5f*omega1)*(u*u*Dxu -v*v*Dyv);
    float X2 = (1.f -omega1)*(C200 -C002) -3.f*rho*(1.f -0.5f*omega1)*(u*u*Dxu -w*w*Dzw);
    float X3 =  K000*omega2 +(1.f -omega2)*(C200 +C020 +C002) -3.f*rho*(1.f -0.5f*omega2)*(u*u*Dxu +v*v*Dyv +w*w*Dzw);
    float Cs200 = (X1 +X2 +X3)/3.f;
    float Cs020 = Cs200 -X1;
    float Cs002 = Cs200 -X2;

    X1 = (1.f -omega3alim)*(C120 +C102);
    X2 = (1.f -omega4alim)*(C120 -C102);
    float Cs120 = 0.5f*(X1 +X2);
    float Cs102 = 0.5f*(X1 -X2);

    X1 = (1.f -omega3blim)*(C210 +C012);
    X2 = (1.f -omega4blim)*(C210 -C012);
    float Cs210 = 0.5f*(X1 +X2);
    float Cs012 = 0.5f*(X1 -X2);

    X1 = (1.f -omega3clim)*(C201 +C021);
    X2 = (1.f -omega4clim)*(C201 -C021);
    float Cs201 = 0.5f*(X1 +X2);
    float Cs021 = 0.5f*(X1 -X2);

    float Cs111 = (1.f -omega5lim)*C111;

    X1 = (1.f -omega6)*(C220 -2.f*C202 +C022) +2.f*(1.f/omega1 -0.5f)*omega6*A*rho*(Dxu -2.f*Dyv +Dzw)/3.f;
    X2 = (1.f -omega6)*(C220 +C202 -2.f*C022) +2.f*(1.f/omega1 -0.5f)*omega6*A*rho*(Dxu +Dyv -2.f*Dzw)/3.f;
    X3 = (1.f -omega7)*(C220 +C202 +C022) -4.f*(1.f/omega1 -0.5f)*omega7*A*rho*(Dxu +Dyv +Dzw)/3.f;
    float Cs220 = (X1 +X2 +X3)/3.f;
    float Cs202 = -(X1 -X3)/3.f;
    float Cs022 = -(X2 -X3)/3.f;

    float Cs211 = (1.f -omega8)*C211 -(1.f/omega1 -0.5f)*omega8*B*rho*Dyw_Dzv/3.f;
    float Cs121 = (1.f -omega8)*C121 -(1.f/omega1 -0.5f)*omega8*B*rho*Dxw_Dzu/3.f;
    float Cs112 = (1.f -omega8)*C112 -(1.f/omega1 -0.5f)*omega8*B*rho*Dxv_Dyu/3.f;

    float Cs221 = (1.f -omega9)*C221;
    float Cs212 = (1.f -omega9)*C212;
    float Cs122 = (1.f -omega9)*C122;

    float Cs222 = (1.f -omega10)*C222;


    float Ks000 = K000; // ?
    float Ks100 = -K100;
    float Ks010 = -K010;
    float Ks001 = -K001;

    float Ks200 = Cs200;
    float Ks020 = Cs020;
    float Ks002 = Cs002;
    float Ks110 = Cs110;
    float Ks101 = Cs101;
    float Ks011 = Cs011;
    float Ks111 = Cs111;
    float Ks210 = Cs210;
    float Ks201 = Cs201;
    float Ks021 = Cs021;
    float Ks120 = Cs120;
    float Ks102 = Cs102;
    float Ks012 = Cs012;

    float Ks211 = Cs211 +(Ks200*Ks011 +2.f*Ks110*Ks101)*invRho;
    float Ks121 = Cs121 +(Ks020*Ks101 +2.f*Ks011*Ks110)*invRho;
    float Ks112 = Cs112 +(Ks002*Ks110 +2.f*Ks101*Ks011)*invRho;

    float Ks220 = Cs220 +(Ks200*Ks020 +2.f*Ks110*Ks110)*invRho;
    float Ks202 = Cs202 +(Ks002*Ks200 +2.f*Ks101*Ks101)*invRho;
    float Ks022 = Cs022 +(Ks020*Ks002 +2.f*Ks011*Ks011)*invRho;

    float Ks122 = Cs122 +(Ks022*Ks120 +Ks020*Ks102 +4.f*Ks011*Ks111 +2.f*(Ks101*Ks021 +Ks110*Ks012))*invRho;
    float Ks212 = Cs212 +(Ks202*Ks012 +Ks002*Ks210 +4.f*Ks101*Ks111 +2.f*(Ks110*Ks102 +Ks011*Ks201))*invRho;
    float Ks221 = Cs221 +(Ks220*Ks201 +Ks200*Ks021 +4.f*Ks110*Ks111 +2.f*(Ks011*Ks210 +Ks101*Ks120))*invRho;

    float Ks222 = Cs222 +(4.f*Ks111*Ks111 +Ks200*Ks022 +Ks020*Ks202 +Ks002*Ks220 +4.f*(Ks011*Ks211 +Ks101*Ks121 +Ks110*Ks112) +2.f*(Ks120*Ks102 +Ks210*Ks012 +Ks201*Ks021))*invRho -(16.f*Ks110*Ks101*Ks011 +4.f*(Ks101*Ks101*Ks020 +Ks011*Ks011*Ks200 +Ks110*Ks110*Ks002) +2.f*Ks200*Ks020*Ks002)*invRho*invRho;

    float Ks0_00 = Ks000*(1.f -u*u) -2.f*u*Ks100 -Ks200;
    float Ks0_10 = Ks010*(1.f -u*u) -2.f*u*Ks110 -Ks210;
    float Ks0_01 = Ks001*(1.f -u*u) -2.f*u*Ks101 -Ks201;
    float Ks0_11 = Ks011*(1.f -u*u) -2.f*u*Ks111 -Ks211;
    float Ks0_20 = Ks020*(1.f -u*u) -2.f*u*Ks120 -Ks220;
    float Ks0_02 = Ks002*(1.f -u*u) -2.f*u*Ks102 -Ks202;
    float Ks0_21 = Ks021*(1.f -u*u) -2.f*u*Ks121 -Ks221;
    float Ks0_12 = Ks012*(1.f -u*u) -2.f*u*Ks112 -Ks212;
    float Ks0_22 = Ks022*(1.f -u*u) -2.f*u*Ks122 -Ks222;

    float Ksm1_00 = 0.5f*(Ks000*(u*u -u) +Ks100*(2.f*u -1.f) +Ks200);
    float Ksm1_10 = 0.5f*(Ks010*(u*u -u) +Ks110*(2.f*u -1.f) +Ks210);
    float Ksm1_01 = 0.5f*(Ks001*(u*u -u) +Ks101*(2.f*u -1.f) +Ks201);
    float Ksm1_11 = 0.5f*(Ks011*(u*u -u) +Ks111*(2.f*u -1.f) +Ks211);
    float Ksm1_20 = 0.5f*(Ks020*(u*u -u) +Ks120*(2.f*u -1.f) +Ks220);
    float Ksm1_02 = 0.5f*(Ks002*(u*u -u) +Ks102*(2.f*u -1.f) +Ks202);
    float Ksm1_21 = 0.5f*(Ks021*(u*u -u) +Ks121*(2.f*u -1.f) +Ks221);
    float Ksm1_12 = 0.5f*(Ks012*(u*u -u) +Ks112*(2.f*u -1.f) +Ks212);
    float Ksm1_22 = 0.5f*(Ks022*(u*u -u) +Ks122*(2.f*u -1.f) +Ks222);

    float Ks1_00 = 0.5f*(Ks000*(u*u +u) +Ks100*(2.f*u +1.f) +Ks200);
    float Ks1_10 = 0.5f*(Ks010*(u*u +u) +Ks110*(2.f*u +1.f) +Ks210);
    float Ks1_01 = 0.5f*(Ks001*(u*u +u) +Ks101*(2.f*u +1.f) +Ks201);
    float Ks1_11 = 0.5f*(Ks011*(u*u +u) +Ks111*(2.f*u +1.f) +Ks211);
    float Ks1_20 = 0.5f*(Ks020*(u*u +u) +Ks120*(2.f*u +1.f) +Ks220);
    float Ks1_02 = 0.5f*(Ks002*(u*u +u) +Ks102*(2.f*u +1.f) +Ks202);
    float Ks1_21 = 0.5f*(Ks021*(u*u +u) +Ks121*(2.f*u +1.f) +Ks221);
    float Ks1_12 = 0.5f*(Ks012*(u*u +u) +Ks112*(2.f*u +1.f) +Ks212);
    float Ks1_22 = 0.5f*(Ks022*(u*u +u) +Ks122*(2.f*u +1.f) +Ks222);

    float Ks00_0 = Ks0_00*(1.f -v*v) -2.f*v*Ks0_10 -Ks0_20;
    float Ks00_1 = Ks0_01*(1.f -v*v) -2.f*v*Ks0_11 -Ks0_21;
    float Ks00_2 = Ks0_02*(1.f -v*v) -2.f*v*Ks0_12 -Ks0_22;
    float Ks10_0 = Ks1_00*(1.f -v*v) -2.f*v*Ks1_10 -Ks1_20;
    float Ks10_1 = Ks1_01*(1.f -v*v) -2.f*v*Ks1_11 -Ks1_21;
    float Ks10_2 = Ks1_02*(1.f -v*v) -2.f*v*Ks1_12 -Ks1_22;
    float Ksm10_0 = Ksm1_00*(1.f -v*v) -2.f*v*Ksm1_10 -Ksm1_20;
    float Ksm10_1 = Ksm1_01*(1.f -v*v) -2.f*v*Ksm1_11 -Ksm1_21;
    float Ksm10_2 = Ksm1_02*(1.f -v*v) -2.f*v*Ksm1_12 -Ksm1_22;

    float Ks0m1_0 = 0.5f*(Ks0_00*(v*v -v) +Ks0_10*(2.f*v -1.f) +Ks0_20);
    float Ks0m1_1 = 0.5f*(Ks0_01*(v*v -v) +Ks0_11*(2.f*v -1.f) +Ks0_21);
    float Ks0m1_2 = 0.5f*(Ks0_02*(v*v -v) +Ks0_12*(2.f*v -1.f) +Ks0_22);
    float Ks1m1_0 = 0.5f*(Ks1_00*(v*v -v) +Ks1_10*(2.f*v -1.f) +Ks1_20);
    float Ks1m1_1 = 0.5f*(Ks1_01*(v*v -v) +Ks1_11*(2.f*v -1.f) +Ks1_21);
    float Ks1m1_2 = 0.5f*(Ks1_02*(v*v -v) +Ks1_12*(2.f*v -1.f) +Ks1_22);
    float Ksm1m1_0 = 0.5f*(Ksm1_00*(v*v -v) +Ksm1_10*(2.f*v -1.f) +Ksm1_20);
    float Ksm1m1_1 = 0.5f*(Ksm1_01*(v*v -v) +Ksm1_11*(2.f*v -1.f) +Ksm1_21);
    float Ksm1m1_2 = 0.5f*(Ksm1_02*(v*v -v) +Ksm1_12*(2.f*v -1.f) +Ksm1_22);

    float Ks01_0 = 0.5f*(Ks0_00*(v*v +v) +Ks0_10*(2.f*v +1.f) +Ks0_20);
    float Ks01_1 = 0.5f*(Ks0_01*(v*v +v) +Ks0_11*(2.f*v +1.f) +Ks0_21);
    float Ks01_2 = 0.5f*(Ks0_02*(v*v +v) +Ks0_12*(2.f*v +1.f) +Ks0_22);
    float Ks11_0 = 0.5f*(Ks1_00*(v*v +v) +Ks1_10*(2.f*v +1.f) +Ks1_20);
    float Ks11_1 = 0.5f*(Ks1_01*(v*v +v) +Ks1_11*(2.f*v +1.f) +Ks1_21);
    float Ks11_2 = 0.5f*(Ks1_02*(v*v +v) +Ks1_12*(2.f*v +1.f) +Ks1_22);
    float Ksm11_0 = 0.5f*(Ksm1_00*(v*v +v) +Ksm1_10*(2.f*v +1.f) +Ksm1_20);
    float Ksm11_1 = 0.5f*(Ksm1_01*(v*v +v) +Ksm1_11*(2.f*v +1.f) +Ksm1_21);
    float Ksm11_2 = 0.5f*(Ksm1_02*(v*v +v) +Ksm1_12*(2.f*v +1.f) +Ksm1_22);


    float f000 = Ks00_0*(1.f -w*w) -2.f*w*Ks00_1 -Ks00_2;
    float f100 = Ks10_0*(1.f -w*w) -2.f*w*Ks10_1 -Ks10_2;
    float fm100 = Ksm10_0*(1.f -w*w) -2.f*w*Ksm10_1 -Ksm10_2;
    float f010 = Ks01_0*(1.f -w*w) -2.f*w*Ks01_1 -Ks01_2;
    float f0m10 = Ks0m1_0*(1.f -w*w) -2.f*w*Ks0m1_1 -Ks0m1_2;
    float f110 = Ks11_0*(1.f -w*w) -2.f*w*Ks11_1 -Ks11_2;
    float fm110 = Ksm11_0*(1.f -w*w) -2.f*w*Ksm11_1 -Ksm11_2;
    float f1m10 = Ks1m1_0*(1.f -w*w) -2.f*w*Ks1m1_1 -Ks1m1_2;
    float fm1m10 = Ksm1m1_0*(1.f -w*w) -2.f*w*Ksm1m1_1 -Ksm1m1_2;


    float f00m1 = 0.5f*(Ks00_0*(w*w -w) +Ks00_1*(2.f*w -1.f) +Ks00_2);
    float f10m1 = 0.5f*(Ks10_0*(w*w -w) +Ks10_1*(2.f*w -1.f) +Ks10_2);
    float fm10m1 = 0.5f*(Ksm10_0*(w*w -w) +Ksm10_1*(2.f*w -1.f) +Ksm10_2);
    float f01m1 = 0.5f*(Ks01_0*(w*w -w) +Ks01_1*(2.f*w -1.f) +Ks01_2);
    float f0m1m1 = 0.5f*(Ks0m1_0*(w*w -w) +Ks0m1_1*(2.f*w -1.f) +Ks0m1_2);
    float f11m1 = 0.5f*(Ks11_0*(w*w -w) +Ks11_1*(2.f*w -1.f) +Ks11_2);
    float f1m1m1 = 0.5f*(Ks1m1_0*(w*w -w) +Ks1m1_1*(2.f*w -1.f) +Ks1m1_2);
    float fm11m1 = 0.5f*(Ksm11_0*(w*w -w) +Ksm11_1*(2.f*w -1.f) +Ksm11_2);
    float fm1m1m1 = 0.5f*(Ksm1m1_0*(w*w -w) +Ksm1m1_1*(2.f*w -1.f) +Ksm1m1_2);
    
    float f001 = 0.5f*(Ks00_0*(w*w +w) +Ks00_1*(2.f*w +1.f) +Ks00_2);
    float f101 = 0.5f*(Ks10_0*(w*w +w) +Ks10_1*(2.f*w +1.f) +Ks10_2);
    float fm101 = 0.5f*(Ksm10_0*(w*w +w) +Ksm10_1*(2.f*w +1.f) +Ksm10_2);
    float f011 = 0.5f*(Ks01_0*(w*w +w) +Ks01_1*(2.f*w +1.f) +Ks01_2);
    float f0m11 = 0.5f*(Ks0m1_0*(w*w +w) +Ks0m1_1*(2.f*w +1.f) +Ks0m1_2);
    float f111 = 0.5f*(Ks11_0*(w*w +w) +Ks11_1*(2.f*w +1.f) +Ks11_2);
    float fm111 = 0.5f*(Ksm11_0*(w*w +w) +Ksm11_1*(2.f*w +1.f) +Ksm11_2);
    float f1m11 = 0.5f*(Ks1m1_0*(w*w +w) +Ks1m1_1*(2.f*w +1.f) +Ks1m1_2);
    float fm1m11 = 0.5f*(Ksm1m1_0*(w*w +w) +Ksm1m1_1*(2.f*w +1.f) +Ksm1m1_2);

    fTmp[ 0*elements +ic] = f000;
    fTmp[ 1*elements +ic] = f100;
    fTmp[ 2*elements +ic] = fm100;
    fTmp[ 3*elements +ic] = f010;
    fTmp[ 4*elements +ic] = f0m10;
    fTmp[ 5*elements +ic] = f001;
    fTmp[ 6*elements +ic] = f00m1;
    fTmp[ 7*elements +ic] = f110;
    fTmp[ 8*elements +ic] = fm1m10;
    fTmp[ 9*elements +ic] = f1m10;
    fTmp[10*elements +ic] = fm110;
    fTmp[11*elements +ic] = f101;
    fTmp[12*elements +ic] = fm10m1;
    fTmp[13*elements +ic] = f10m1;
    fTmp[14*elements +ic] = fm101;
    fTmp[15*elements +ic] = f011;
    fTmp[16*elements +ic] = f0m1m1;
    fTmp[17*elements +ic] = f01m1;
    fTmp[18*elements +ic] = f0m11;
    fTmp[19*elements +ic] = f111;
    fTmp[20*elements +ic] = fm1m1m1;
    fTmp[21*elements +ic] = fm111;
    fTmp[22*elements +ic] = f1m1m1;
    fTmp[23*elements +ic] = f1m11;
    fTmp[24*elements +ic] = fm11m1;
    fTmp[25*elements +ic] = f11m1;
    fTmp[26*elements +ic] = fm1m11;
}

void updateRhoUVW(const float* ft, __global float* rhoList, __global float* uList, __global float* vList, __global float* wList, const int ic, const float dpdx)
{
    float rho = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    float w = 0.0f;

    cal_rhoUVW(ft, &rho, &u, &v, &w);

    rhoList[ic] = rho;
    uList[ic] = u;
    vList[ic] = v;
    wList[ic] = w;

    uList[ic] += 0.5f*dpdx/rho;
}

__kernel void k_streamingCollision // Pull
(
   __global float* f, __global float* fTmp,
   __global int* boundaryList,
   __global int* boundary1List, __global int* boundary2List, __global int* boundary3List,
   __global float* sdfList, __global unsigned char* solidList, __global unsigned char* neiSolidList, 
   __global float* u0List, __global float* v0List, __global float* w0List,
   __global float* Fwx, __global float* Fwy, __global float* Fwz,
   __global float* tauSGSList,
   __global float* rhoList,
   __global float* uList, __global float* vList, __global float* wList,
   __global float* omegaList,
   __global float* GxIBM, __global float* GyIBM, __global float* GzIBM,
   const unsigned elements,
   const float dpdx,
   const float rho_av,
   const int nx, const int ny, const int nz,
   const float LES,
   const int isReadMovingWalls,
   const float omegaB,
   const float spzWidth,
   float clim
)
{
    
    int ic = get_global_id(0);
    if(solidList[ic] == 0)
    {
        int i = ic2i(ic,nx,ny);
        int j = ic2j(ic,nx,ny);
        int k = ic2k(ic,nx,ny);

        float wt[27] = {8.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0};

        //                 0    1     2    3    4     5     6    7     8     9    10    11    12    13    14   15    16    17    18   19    20    21    22    23    24    25   26
        float cx[27] = {0.f, 1.f, -1.f, 0.f,  0.f,  0.f,  0.f, 1.f, -1.f,  1.f, -1.f,  1.f, -1.f,  1.f, -1.f, 0.f,  0.f,  0.f,  0.f, 1.f, -1.f, -1.f,  1.f,  1.f, -1.f,  1.f, -1.f};
        float cy[27] = {0.f, 0.f,  0.f, 1.f, -1.f,  0.f,  0.f, 1.f, -1.f, -1.f,  1.f,  0.f,  0.f,  0.f,  0.f, 1.f, -1.f,  1.f, -1.f, 1.f, -1.f,  1.f, -1.f, -1.f,  1.f,  1.f, -1.f};
        float cz[27] = {0.f, 0.f,  0.f, 0.f,  0.f,  1.f, -1.f, 0.f,  0.f,  0.f,  0.f,  1.f, -1.f, -1.f,  1.f, 1.f, -1.f, -1.f,  1.f, 1.f, -1.f,  1.f, -1.f,  1.f, -1.f, -1.f,  1.f};

        float ft[27] = {0.f};
        int upID[27];

        float u = uList[ic];
        float v = vList[ic];
        float w = wList[ic];
        const int boundary = boundaryList[ic];
        const int boundary1 = boundary1List[ic];
        const int boundary2 = boundary2List[ic];
        const int boundary3 = boundary3List[ic];
        const float sdf = sdfList[ic];
        const float u0 = u0List[ic];
        const float v0 = v0List[ic];
        const float w0 = w0List[ic];
        const float omega = omegaList[ic];
        float tauSGS = tauSGSList[ic];
        const float nu = (1.f/omega -0.5f)/3.f;
        const float nuEff = nu+tauSGS/3.f;

        if(boundary == 1)
        {
            streaming(ft, f, upID, boundary1, boundary2, boundary3, ic, i, j, k, nx, ny, nz, elements, Fwx, Fwy, Fwz, cx, cy, cz);
        
            const int corner = cornerFlag(boundary3, i, j, k, nx, ny, nz);
            fixedVelocityBC(ft, rhoList, u0, v0, w0, boundary1, boundary2, boundary3, ic, i, j, k, nx, ny, nz, corner);

            fixedDensityBC(ft, 1.f, uList, vList, wList, cx, cy, cz, wt, boundary1, boundary2, boundary3, ic, i, j, k, nx, ny, nz, corner);        

            wallFunctionBC(ft, rhoList, uList, vList, wList, boundary1, boundary2, boundary3, ic, i, j, k, nx, ny, nz, corner, nu, nuEff, Fwx[ic], Fwy[ic], Fwz[ic]);

            outflowBC(ft, f, u, v, w, boundary1, boundary2, boundary3, ic, i, j, k, nx, ny, nz);
        }
        else
        {
            if(nz-1 == 0)
            {
                streaming2D(ft, f, upID, ic, i, j, k, nx, ny, nz, elements);
            }
            else
            {
                streamingInternal(ft, f, upID, ic, i, j, k, nx, ny, nz, elements);
            }
        }

        internalWallBC(ft, f, Fwx, Fwy, Fwz, solidList, neiSolidList, sdfList, sdf, upID, omega, tauSGS, rhoList, uList, vList, wList, cx, cy, cz, wt, ic, i, j, k, nx, ny, nz, elements);
        updateRhoUVW(ft, rhoList, uList, vList, wList, ic, dpdx);
        float rho = rhoList[ic];
        u = uList[ic];
        v = vList[ic];
        w = wList[ic];

        GxIBM[ic] = 0.f;
        GyIBM[ic] = 0.f;
        GzIBM[ic] = 0.f;            

        float dudx;
        float dudy;
        float dudz;
        float dvdx;
        float dvdy;
        float dvdz;
        float dwdx;
        float dwdy;
        float dwdz;

        // tauSGSList[ic] = smagorinskyTauSGS(ft, rho, u, v, w, omega, LES, sdf, cx, cy, cz, wt, 0.1f, 1); //base
        // tauSGSList[ic] = vremanTauSGS(i,j,k,nx,ny,nz,uList,vList,wList,LES); //Performance -3.7%
        tauSGSList[ic] = waleTauSGS(i,j,k,nx,ny,nz,uList,vList,wList,LES); //Performance -8.6%
        // tauSGSList[ic] = CSsmagorinskyTauSGS(ft, rho, u, v, w, omega, LES, sdf, cx, cy, cz, wt,i,j,k,nx,ny,nz,uList,vList,wList); // Performance -9.12%

        tauSGS = tauSGSList[ic];

        float tau = tauSpongeZone(spzWidth, boundary1List, boundary2List, boundary3List, omega, tauSGS, i, j, k, nx, ny, nz);
        
        // collisionBGK(fTmp, ft, rho, u, v, w, tau, tauSGS, dpdx, cx, cy, cz, wt, ic, elements);
        collisionCumulantGeier(fTmp, ft, rho, u, v, w, tau, tauSGS, omegaB, clim, dpdx, cx, cy, cz, wt, ic, elements);
    }
}

// __kernel void k_Gwall
// (
//     __global float* rho,
//     __global float* u, __global float* v, __global float* w,
//     __global float* movingSTLcList,
//     __global float* GxMovingWall, __global float* GyMovingWall, __global float* GzMovingWall,
//     const int nMovingSTL,
//     const int nx, const int ny, const int nz,
//     const float uMovingTrans, const float vMovingTrans, const float wMovingTrans,
//     const float rotX, const float rotY, const float rotZ,
//     const float rotAxisX, const float rotAxisY, const float rotAxisZ,
//     const float rotOmega
// )
// {
//     int iMSTL = get_global_id(0);

//     float u0 = uMovingTrans;
//     float v0 = vMovingTrans;
//     float w0 = wMovingTrans;

//     float uMovingWall = 0.f;
//     float vMovingWall = 0.f;
//     float wMovingWall = 0.f;

//     float wallX = movingSTLcList[3*iMSTL];
//     float wallY = movingSTLcList[3*iMSTL+1];
//     float wallZ = movingSTLcList[3*iMSTL+2];

//     int i = (int)(wallX);
//     int j = (int)(wallY);
//     int k = (int)(wallZ);

//     float uMovingRot = 0.f;
//     float vMovingRot = 0.f;
//     float wMovingRot = 0.f;
//     Urot(wallX, wallY, wallZ, rotX, rotY, rotZ, rotAxisX, rotAxisY, rotAxisZ, rotOmega, &uMovingRot, &vMovingRot, &wMovingRot);

//     u0 += uMovingRot;
//     v0 += vMovingRot;
//     w0 += wMovingRot;

//     if(i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
//     {
//         int icM = index1d(i,j,k,nx,ny);

//         for(int iBox = 0; iBox < 8; iBox++)
//         {
//             int icBoxPoint = icBox(icM, iBox, nx, ny, nz);

//             if(icBoxPoint != -1)
//             {
//                 int iBox = ic2i(icBoxPoint,nx,ny);
//                 int jBox = ic2j(icBoxPoint,nx,ny); 
//                 int kBox = ic2k(icBoxPoint,nx,ny); 

//                 float delta =   (1.f -fabs(iBox -wallX))
//                                 *(1.f -fabs(jBox -wallY))
//                                 *(1.f -fabs(kBox -wallZ));
                
//                 uMovingWall += u[icBoxPoint]*delta;
//                 vMovingWall += v[icBoxPoint]*delta;
//                 wMovingWall += w[icBoxPoint]*delta;
//             }
//         }
//     }
//     GxMovingWall[iMSTL] = u0 -uMovingWall;
//     GyMovingWall[iMSTL] = v0 -vMovingWall;
//     GzMovingWall[iMSTL] = w0 -wMovingWall;
// }

// __kernel void k_Gibm
// (
//     __global float* rho,
//     __global float* GxIBM, __global float* GyIBM, __global float* GzIBM,
//     __global float* movingSTLcList,
//     __global float* GxMovingWall, __global float* GyMovingWall, __global float* GzMovingWall,
//     const int nMovingSTL,
//     const int nx, const int ny, const int nz,
//     const float uMovingTrans, const float vMovingTrans, const float wMovingTrans,
//     const float rotX, const float rotY, const float rotZ,
//     const float rotAxisX, const float rotAxisY, const float rotAxisZ,
//     const float rotOmega
// )
// {
//     int iMSTL = get_global_id(0);

//     int i = (int)(movingSTLcList[3*iMSTL]);
//     int j = (int)(movingSTLcList[3*iMSTL+1]);
//     int k = (int)(movingSTLcList[3*iMSTL+2]);

//     float u0 = uMovingTrans;
//     float v0 = vMovingTrans;
//     float w0 = wMovingTrans;

//     if(i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
//     {
//         int icM = index1d(i,j,k,nx,ny);

//         for(int iBox = 0; iBox < 8; iBox++)
//         {
//             int icBoxPoint = icBox(icM, iBox, nx, ny, nz);

//             if(icBoxPoint != -1)
//             {
//                 int iBox = ic2i(icBoxPoint,nx,ny);
//                 int jBox = ic2j(icBoxPoint,nx,ny); 
//                 int kBox = ic2k(icBoxPoint,nx,ny); 

//                 float delta =   (1.f -fabs(iBox -movingSTLcList[3*iMSTL]))
//                                 *(1.f -fabs(jBox -movingSTLcList[3*iMSTL+1]))
//                                 *(1.f -fabs(kBox -movingSTLcList[3*iMSTL+2]));
                
//                 GxIBM[icBoxPoint] = atom_add_float(&GxIBM[icBoxPoint],GxMovingWall[iMSTL]*delta);
//                 GyIBM[icBoxPoint] = atom_add_float(&GyIBM[icBoxPoint],GyMovingWall[iMSTL]*delta);
//                 GzIBM[icBoxPoint] = atom_add_float(&GzIBM[icBoxPoint],GzMovingWall[iMSTL]*delta);
//             }
//         }
//     }
    
//     movingSTLcList[3*iMSTL] += u0;
//     movingSTLcList[3*iMSTL+1] += v0;
//     movingSTLcList[3*iMSTL+2] += w0;

//     float wallX = movingSTLcList[3*iMSTL];
//     float wallY = movingSTLcList[3*iMSTL+1];
//     float wallZ = movingSTLcList[3*iMSTL+2];

//     Xrot(wallX, wallY, wallZ, rotX, rotY, rotZ, rotAxisX, rotAxisY, rotAxisZ, rotOmega, &wallX, &wallY, &wallZ);
//     movingSTLcList[3*iMSTL] = wallX;
//     movingSTLcList[3*iMSTL+1] = wallY;
//     movingSTLcList[3*iMSTL+2] = wallZ;
// }

// __kernel void k_Force
// (
//     __global float* fTmp,
//     __global float* solid,
//     __global float* rho,
//     __global float* GxIBM, __global float* GyIBM, __global float* GzIBM,
//     const unsigned elements
// )
// {
//     int ic = get_global_id(0);

//     float wt[19] = {1.0f/3.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f};

//     //                 0     1      2     3      4     5      6     7      8      9     10    11     12     13     14    15     16     17     18
//     float cx[19] = {0.0f, 1.0f, -1.0f, 0.0f,  0.0f, 0.0f,  0.0f, 1.0f, -1.0f,  1.0f, -1.0f, 1.0f, -1.0f,  1.0f, -1.0f, 0.0f,  0.0f,  0.0f,  0.0f};
//     float cy[19] = {0.0f, 0.0f,  0.0f, 1.0f, -1.0f, 0.0f,  0.0f, 1.0f, -1.0f, -1.0f,  1.0f, 0.0f,  0.0f,  0.0f,  0.0f, 1.0f, -1.0f,  1.0f, -1.0f};
//     float cz[19] = {0.0f, 0.0f,  0.0f, 0.0f,  0.0f, 1.0f, -1.0f, 0.0f,  0.0f,  0.0f,  0.0f, 1.0f, -1.0f, -1.0f,  1.0f, 1.0f, -1.0f, -1.0f,  1.0f};

//     if(solid[ic] == 0)
//     {
//         for(int q = 0; q < 19; q++)
//         {
//             int qic = q*elements +ic;
//             fTmp[qic] += rho[ic]*wt[q]*3.0f*(GxIBM[ic]*cx[q] +GyIBM[ic]*cy[q] +GzIBM[ic]*cz[q]);
//         }
//     }
// }