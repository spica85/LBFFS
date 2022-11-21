
#ifndef base_H
#define base_H

inline int ic2i(int ic, int nx, int ny)
{
    return int((ic%(nx*ny))%nx);
}

inline int ic2j(int ic, int nx, int ny)
{
    return int(ic%(nx*ny)/nx);
}

inline int ic2k(int ic, int nx, int ny)
{
    return int(ic/(nx*ny));
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

char* asciiToBinary(char* str, const float x)
{
    str[0] = ((char*)&x)[3 + 0];
    str[1] = ((char*)&x)[2 + 0];
    str[2] = ((char*)&x)[1 + 0];
    str[3] = ((char*)&x)[3 + 0];

    return str;    
}

#endif
