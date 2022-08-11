if(Fwrite && nextOutTime < nt +1)
{
    cl::copy(queue, f_d, f.begin(), f.end());
    cl::copy(queue, Fwx_d, Fwx.begin(), Fwx.end());
    cl::copy(queue, Fwy_d, Fwy.begin(), Fwy.end());
    cl::copy(queue, Fwz_d, Fwz.begin(), Fwz.end());
    cl::copy(queue, GxIBM_d, GxIBM.begin(), GxIBM.end());
    cl::copy(queue, GyIBM_d, GyIBM.begin(), GyIBM.end());
    cl::copy(queue, GzIBM_d, GzIBM.begin(), GzIBM.end());
    cl::copy(queue, tauSGS_d, tauSGS.begin(), tauSGS.end());
    rho_av = 0.f;
    int eleFluid = 0;
    float SumFwx = 0.f;
    float SumFwy = 0.f;
    float SumFwz = 0.f;
    float SumGxIBM = 0.f;
    float SumGyIBM = 0.f;
    float SumGzIBM = 0.f;

    for(int ic =0; ic < elements; ic++)
    {
        cal_rhoUVW(ic, nx, ny, nz, f, cx, cy, cz, rho[ic], u[ic], v[ic], w[ic]);
        if(solid[ic] == 0)
        {
            rho_av += rho[ic] -1.f;
            SumFwx += Fwx[ic];
            SumFwy += Fwy[ic];
            SumFwz += Fwz[ic];
            SumGxIBM += GxIBM[ic];
            SumGyIBM += GyIBM[ic];
            SumGzIBM += GzIBM[ic];
            eleFluid++;
            // std::cout << "eleFluid: " << eleFluid
            //           << ", rho_av: " << rho_av
            //           << ", rho: " << rho[ic]
            //           << std::endl;
        }
    }
    rho_av = rho_av/float(eleFluid) +1.f;
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

    std::cout << "rhoMax: " << rhoMax << ", rhoMin: " << rhoMin << ", rhoAve: " << rho_av << std::endl;

    //- For Flow around cylinder
    if(forceCoeffs)
    {
        std::cout << "Cx: " << SumFwx*c*c*L*L/(0.5*rho_av*uMax*uMax*Dref*L) << ", Cy: " << SumFwy*c*c*L*L/(0.5*rho_av*uMax*uMax*Dref*L) << ", Cz: " << SumFwz*c*c*L*L/(0.5*rho_av*uMax*uMax*Dref*L) << std::endl;
        std::cout << "CxIBM: " << -SumGxIBM*c*c*L*L/(0.5*rho_av*uMax*uMax*Dref*L) << ", CyIBM: " << -SumGyIBM*c*c*L*L/(0.5*rho_av*uMax*uMax*Dref*L) << ", CzIBM: " << -SumGzIBM*c*c*L*L/(0.5*rho_av*uMax*uMax*Dref*L) << std::endl;
    }
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

    // boundary1 conditions output
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
                        asciiToBinary(str,(float)boundary1[ic]);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << boundary1[ic] << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }                    

    // boundary2 conditions output
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
                        asciiToBinary(str,(float)boundary2[ic]);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << boundary2[ic] << "\n";
                    }
                }
            }
        }
        writeFile << "\n";
    }                    

    // boundary3 conditions output
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
                        asciiToBinary(str,(float)boundary3[ic]);
                        writeFile.write(str,sizeof(char)*4);
                    }
                    else
                    {
                        writeFile << boundary3[ic] << "\n";
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
