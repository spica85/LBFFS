if(Fwrite && nextOutTime < nt +1)
{
    cl::copy(queue, f_d, f.begin(), f.end());
    cl::copy(queue, Fwx_d, Fwx.begin(), Fwx.end());
    cl::copy(queue, Fwy_d, Fwy.begin(), Fwy.end());
    cl::copy(queue, Fwz_d, Fwz.begin(), Fwz.end());
    cl::copy(queue, tauSGS_d, tauSGS.begin(), tauSGS.end());
    rho_av = 0.0;
    float SumFwx = 0.f;
    float SumFwy = 0.f;
    float SumFwz = 0.f;
    #pragma omp parallel for
    for(int ic =0; ic <nx*ny*nz; ic++)
    {
        cal_rhoUVW(ic, nx, ny, nz, f, cx, cy, cz, rho[ic], u[ic], v[ic], w[ic]);
        rho_av += rho[ic];
        SumFwx += Fwx[ic];
        SumFwy += Fwy[ic];
        SumFwz += Fwz[ic];
    }
    rho_av /= float(elements);
    float rhoMax = *std::max_element(rho.begin(), rho.end());
    float rhoMin = *std::min_element(rho.begin(), rho.end());

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

    std::cout << "Time = " << nt*deltaT << " (s), " << nt << " steps" << std::endl;

    end = std::chrono::system_clock::now();
    float time = static_cast<float>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() *1e-6);
    std::cout << "Execution time: " << time << " (s)" << std::endl;
    std::cout << "Speed: " << float(nt)*float(nx*ny*nz)/time*1e-6 << " (MLUPS)" << std::endl;

    std::cout << "rhoMax: " << rhoMax << ", rhoMin: " << rhoMin << std::endl;

    //- For Flow around cylinder
    // std::cout << "Cx: " << SumFwx*c*c*L*L/(0.5*rho_av*uMax*uMax*d*L) << ", Cy: " << SumFwy*c*c*L*L/(0.5*rho_av*uMax*uMax*d*L) << ", Cz: " << SumFwz*c*c*L*L/(0.5*rho_av*uMax*uMax*d*L) << std::endl;
    //--

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

    // writeFile << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << "\n";
    writeFile << "ORIGIN " << L/2 << " " << L/2 << " " << L/2 << "\n";
    writeFile << "SPACING " << L << " " << L << " " << L << "\n\n";


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
                        asciiToBinary(str,(float)(u[ic]*c));
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << u[ic]*c << "\n";
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
                        asciiToBinary(str,(float)(v[ic]*c));
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << v[ic]*c << "\n";
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
                        asciiToBinary(str,(float)(w[ic]*c));
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << w[ic]*c << "\n";
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
                        asciiToBinary(str,(float)rho[ic]);
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

    // tauSGS output
    {
        writeFile << "SCALARS tauSGS float\n";
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
                        asciiToBinary(str,(float)tauSGS[ic]);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << tauSGS[ic] << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }

    // obst output
    {
        writeFile << "SCALARS boundary1 float\n";
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
                        asciiToBinary(str,(float)obst[ic].boundary1);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << obst[ic].boundary1 << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }                    

    // obst output
    {
        writeFile << "SCALARS boundary2 float\n";
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
                        asciiToBinary(str,(float)obst[ic].boundary2);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << obst[ic].boundary2 << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }                    

    // obst output
    {
        writeFile << "SCALARS boundary3 float\n";
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
                        asciiToBinary(str,(float)obst[ic].boundary3);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << obst[ic].boundary3 << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }                    

    // sdf output
    {
        writeFile << "SCALARS sdf float\n";
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
                        asciiToBinary(str,(float)(sdf[ic]*L));
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << sdf[ic]*L << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }

    // solid output
    {
        writeFile << "SCALARS solid float\n";
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
                        asciiToBinary(str,(float)(solid[ic]));
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << (float)(solid[ic]) << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }


    // neiSolid output
    {
        writeFile << "SCALARS neiSolid float\n";
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
                        asciiToBinary(str,(float)(neiSolid[ic]));
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << (float)(neiSolid[ic]) << "\n";
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
