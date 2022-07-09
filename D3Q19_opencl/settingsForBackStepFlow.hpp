//- For backward facing step flow
float h = 1.f;
float Ly = 3.f*h; // 5.f*h// Flow around Digital Science
float Lx = 40.f*h;
float nu = uMax*h/Re;
const float L = Ly/float(ny);
float dpdx = 0.f;
