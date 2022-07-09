//-- For channel flow
float Retau = 10;
float Retau = 1;
float utau = 0.005;
float nu = utau*0.5*ny/Retau;
float dpdx = utau*utau/(0.5*ny);
float h = 1.f;
const float L = h/float(ny);
