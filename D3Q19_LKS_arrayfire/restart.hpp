std::ifstream restartFile;

std::cout << "Reading from ./out/data.dat" <<std::endl;
restartFile.open("./out/data.dat", std::ios::in | std::ios::binary);

restartFile.read((char *) &startTimeStep,sizeof(int));
nextOutTime = startTimeStep+outInterval;
std::cout << "nextOutTime: " << nextOutTime << std::endl;

restartFile.read((char *) &nx,sizeof(int));
restartFile.read((char *) &ny,sizeof(int));
restartFile.read((char *) &nz,sizeof(int));

// read u
for(int k = 0; k < nz; k++)
{
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);
            restartFile.read((char *) &u[ic],sizeof(double));
        }
    }
}

// read v
for(int k = 0; k < nz; k++)
{
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);
            restartFile.read((char *) &v[ic],sizeof(double));
        }
    }
}

// read w
for(int k = 0; k < nz; k++)
{
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);
            restartFile.read((char *) &w[ic],sizeof(double));
        }
    }
}

// read p
for(int k = 0; k < nz; k++)
{
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);
            restartFile.read((char *) &p[ic],sizeof(double));
        }
    }
}
