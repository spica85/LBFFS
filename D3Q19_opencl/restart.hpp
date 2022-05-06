std::ifstream restartFile;

std::cout << "Reading from ./out/data.dat" <<std::endl;
restartFile.open("./out/data.dat", std::ios::in | std::ios::binary);

restartFile.read((char *) &startTimeStep,sizeof(int));
nextOutTime = startTimeStep;

restartFile.read((char *) &nx,sizeof(int));
restartFile.read((char *) &ny,sizeof(int));
restartFile.read((char *) &nz,sizeof(int));

for(int q = 0; q < 19; q++)
{
    for(int k = 0; k < nz; k++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int i = 0; i < nx; i++)
            {
                int ic = index1d(i,j,k,nx,ny);
                int qic = idf(q,ic,nx,ny,nz);
                restartFile.read((char *) &f[qic],sizeof(float));
                //std::cout << "f[" << qic << "] = " << f[qic] << std::endl; // For debug
            }
        }
    }
}
