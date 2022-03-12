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
    char str[300];// For binary

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
                        asciiToBinary(str,(float)u[ic]);
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
                        asciiToBinary(str,(float)v[ic]);
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
                        asciiToBinary(str,(float)w[ic]);
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

        
    // p output
    {
        writeFile << "SCALARS p float\n";
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
                        asciiToBinary(str,(float)p[ic]);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << p[ic] << "\n";
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
                        asciiToBinary(str,(float)obst[ic].boundary);
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

    // Writing for restart
    std::string filename_data = std::to_string(nt);
    std::ofstream writeFile_data;
    writeFile_data.open("./out/data.dat", std::ios::out | std::ios::binary);

    // timeStep
    {
        writeFile_data.write((char *) &nt,sizeof(int));
    }

    //nx ny nz output
    {
        writeFile_data.write((char *) &nx,sizeof(int));
        writeFile_data.write((char *) &ny,sizeof(int));
        writeFile_data.write((char *) &nz,sizeof(int));
    }

    //u output
    {
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {
                    int ic = index1d(i,j,k,nx,ny);
                    writeFile_data.write((char *) &u[ic],sizeof(double));
                }
            }
        }
    }

    //v output
    {
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {
                    int ic = index1d(i,j,k,nx,ny);
                    writeFile_data.write((char *) &v[ic],sizeof(double));
                }
            }
        }
    }

    //w output
    {
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {
                    int ic = index1d(i,j,k,nx,ny);
                    writeFile_data.write((char *) &w[ic],sizeof(double));
                }
            }
        }
    }

    //p output
    {
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {
                    int ic = index1d(i,j,k,nx,ny);
                    writeFile_data.write((char *) &p[ic],sizeof(double));
                }
            }
        }
    }
    writeFile_data.close();

    nextOutTime += outInterval;
}
