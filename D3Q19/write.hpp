if(Fwrite && nextOutTime < nt +1)
{
    std::ostringstream ss;
    ss << std::setw(8) << std::setfill('0') << std::to_string(nt);
    std::string filename = ss.str();
    std::ofstream writeFile;
    if(writeBinary)
    {
        writeFile.open("./out/res"+filename+".vtk", std::ios::out | std::ios::binary);        
    }
    else
    {
        writeFile.open("./out/res"+filename+".vtk", std::ios::out);
    }
    char str[300];// for binary

    std::cout << "Time = " << nt << std::endl;

    writeFile << "# vtk DataFile Version 3.0\n";
    writeFile << "vtk output\n";
    if(writeBinary)
    {
        writeFile << "BINARY\n";
    }
    else
    {
        writeFile << "ASCII\n";        
    }
    writeFile << "DATASET STRUCTURED_POINTS\n";
    writeFile << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";

    writeFile << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << "\n";
    writeFile << "SPACING " << 1.0 << " " << 1.0 << " " << 1.0 << "\n\n";


    writeFile << "POINT_DATA " << nx*ny*nz << "\n";

    // u output
    {
        writeFile << "SCALARS u float\n";
        writeFile << "LOOKUP_TABLE default\n";
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {                                
                    int ic = index1d(i,j,k,nx,ny);
                    if(writeBinary)
                    {
                        asciiToBinary(str,u[ic]);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << u[ic] << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }

    // v output
    {
        writeFile << "SCALARS v float\n";
        writeFile << "LOOKUP_TABLE default\n";
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {                                
                    int ic = index1d(i,j,k,nx,ny);
                    if(writeBinary)
                    {
                        asciiToBinary(str,v[ic]);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << v[ic] << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }

    // w output
    {
        writeFile << "SCALARS w float\n";
        writeFile << "LOOKUP_TABLE default\n";
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {                                
                    int ic = index1d(i,j,k,nx,ny);
                    if(writeBinary)
                    {
                        asciiToBinary(str,w[ic]);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << w[ic] << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }                    

        
    // rho output
    {
        writeFile << "SCALARS rho float\n";
        writeFile << "LOOKUP_TABLE default\n";

        for(int k = 0; k < nz; k++)                            
        {                        
            for(int j = 0; j < ny; j++)                            
            {
                for(int i = 0; i < nx; i++)
                {
                    int ic = index1d(i,j,k,nx,ny);
                    if(writeBinary)
                    {
                        asciiToBinary(str,rho[ic]);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << rho[ic] << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }

    // obst output
    {
        writeFile << "SCALARS obst float\n";
        writeFile << "LOOKUP_TABLE default\n";
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {                                
                    int ic = index1d(i,j,k,nx,ny);
                    if(writeBinary)
                    {
                        asciiToBinary(str,obst[ic].boundary);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << obst[ic].boundary << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }                    

    writeFile.close();

    std::string filename_data = std::to_string(nt);
    std::ofstream writeFile_data;
    writeFile_data.open("./out/data"+filename_data+".dat", std::ios::out | std::ios::binary);

    //nx ny nz output
    {
        writeFile_data.write((char *) &nx,sizeof(int));
        writeFile_data.write((char *) &ny,sizeof(int));
        writeFile_data.write((char *) &nz,sizeof(int));
    }

    //f output
    {
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
                        writeFile_data.write((char *) &f[qic],sizeof(float));
                    }
                }
            }
        }

    }

    nextOutTime += outInterval;
}
