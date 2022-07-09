//- For cavity flow
const float a = 1.0; //Dimensional length of system (m)
const float L = a/float(nx); //Representative length (-)  
float nu = uMax*a/Re; //Dimensional kinematic viscosity
float dpdx = 0.0;//Dimensionless external force
//--