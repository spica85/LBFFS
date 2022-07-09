//- For Poiseuille flow
float h = 1.f;
float nu = uMax*h/Re;
const float L = h/float(ny);
float dpdx = 8.0*nu*uMax/(h*h);
