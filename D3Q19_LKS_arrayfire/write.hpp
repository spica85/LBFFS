if(Fwrite && nextOutTime < nt +1)
{   
    float* pHost = p_af.host<float>();
    p = std::vector<float>(pHost,&pHost[nx*ny*nz]);
    float* uHost = u_af.host<float>();
    u = std::vector<float>(uHost,&uHost[nx*ny*nz]);
    float* vHost = v_af.host<float>();
    v = std::vector<float>(vHost,&vHost[nx*ny*nz]);
    float* wHost = w_af.host<float>();
    w = std::vector<float>(wHost,&wHost[nx*ny*nz]);

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

    af::sync();
    auto stop_clock = af::timer::stop();
    float time = static_cast<float>(stop_clock);
    std::cout << "Execution time: " << time << " (s)" << std::endl;
    std::cout << "Speed: " << float(outInterval)*float(nx*ny*nz)/time << " (LUPS)" << std::endl;

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
                    writeFile_data.write((char *) &u[ic],sizeof(float));
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
                    writeFile_data.write((char *) &v[ic],sizeof(float));
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
                    writeFile_data.write((char *) &w[ic],sizeof(float));
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
                    writeFile_data.write((char *) &p[ic],sizeof(float));
                }
            }
        }
    }
    writeFile_data.close();

    nextOutTime += outInterval;

    af::timer::start();
}
