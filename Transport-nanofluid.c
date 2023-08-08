// This UDF simulates the transport and deposition of nanoparticles in a nanofluid
// A nanofluid is a fluid that contains suspended nanoparticles that can enhance its thermal and electrical properties
// The UDF uses various effects such as buoyancy, thermophoresis, Brownian motion, Dufour effect, and Soret effect
// The UDF also demonstrates how to use macros, aliases, and user-defined functions in Fluent UDFs

#include "udf.h"

// User inputs and constants
#define s0 0 // Initial salinity of the mixture (%)
#define t0 300 // Initial temperature of the mixture (K)
#define beta_t 0.0002 // Coefficient of thermal expansion of the base fluid (1/K)
#define beta_s 0.0008 // Coefficient of salinity of the base fluid (1/%)
#define beta_t_p 0.0000008 // Coefficient of thermal expansion of the nanoparticle (1/K)
#define ks 0.0000000014	// Mass diffusivity of the nanoparticle (m^2/s)
#define rho_p 3950 // Density of the nanoparticle (kg/m^3)
#define k_bf 0.5974 // Thermal conductivity of the base fluid (W/m.K)
#define C_P_p 765 // Specific heat capacity of the nanoparticle (J/kg.K)
#define mu_bf 0.001 // Dynamic viscosity of the base fluid (kg/m.s)
#define d_p 0.00000008 // Diameter of the nanoparticle (m)
#define density_nf 1029.5 // Density of the nanofluid (kg/m^3)
#define spheat_nf 4050.896 // Specific heat capacity of the nanofluid (J/kg.K)
#define thercond_nf 0.58058 // Thermal conductivity of the nanofluid (W/m.K)
#define dyvisco_nf 0.00102544 // Dynamic viscosity of the nanofluid (kg/m.s)
#define beta_t_nf 0.2315 // Coefficient of thermal expansion of the nanofluid (1/K)
#define D_T 3.7568e-14 // Thermal diffusivity of the nanofluid (m^2/s)

/* Defined constants */


#define k_b 1.38064852e-23 /*Boltz Const UNIT:m^2kgs^-2K^-1*/
#define pi 3.14159 /* pi */
#define g 9.81
#define phi(c,t) C_UDSI(c,t,0)

DEFINE_PROPERTY(nf_density, c, t)
{ 
       real rho_nf;
       rho_nf = density_nf;
       return rho_nf;
}

DEFINE_SPECIFIC_HEAT(nf_spheat, T, Tref, h, yi)
{
       real cp=4050.896 ;
       *h = cp*(T-Tref);
       return cp;
}

DEFINE_PROPERTY(nf_th_conductivity, c, t)
{
       real k_nf;
       k_nf = thercond_nf;
       return k_nf;
}

DEFINE_PROPERTY(nf_viscosity, c, t)
{
       real mu_nf;
       mu_nf = dyvisco_nf ;
       return mu_nf;
}

DEFINE_DIFFUSIVITY(nf_diffusivity, c, t)
{
       real gamma;
       gamma = ks; /*MassDiffusivity Unit:m^2/s*/
       return gamma;
}

DEFINE_DIFFUSIVITY(UDS_diffusivity, c, t, i)
{
       real gamma_UDS, temp, mu_nf, rho_nf, d_b, phi;	
       temp = C_T(c, t);
       phi = C_UDSI(c,t,0);
       d_b = (k_b * temp)/(3*pi*dyvisco_nf*d_p);
       /*d_b Brownian Diffusion Coefficient UNIT:m^2/s*/
       gamma_UDS = density_nf * d_b;
       return gamma_UDS;
}

DEFINE_SOURCE(ymom_source, c, t, dS, eqn)
{
       real x[ND_ND], momt_source, delta_t, delta_s;
       real conc_salt, temp;
       temp = C_T(c, t);
       conc_salt = C_YI(c, t, 0);
       delta_t = temp - t0; /* temperature difference */
       delta_s = conc_salt - s0;
       /* C_CENTROID(x,c,t);*/
       momt_source = - density_nf * (1 - beta_t_nf * delta_t + beta_s * delta_s)*g;
       dS[eqn] = 0;
       return momt_source;
}

DEFINE_SOURCE(UDS_source, c, t, dS, eqn)
{
       real x[ND_ND], vof_source, temp, phi, D_ttemp;
       temp = C_T(c, t);
       phi = C_UDSI(c,t,0);								
       D_ttemp = phi/temp;
       C_UDSI(c,t,1) = D_ttemp * C_T_G(c,t)[0];
       C_UDSI(c,t,2) = D_ttemp * C_T_G(c,t)[1]; 
       vof_source = density_nf * D_T *( C_UDSI_G(c,t,1)[0] + C_UDSI_G(c,t,1)[1] );
       C_CENTROID(x,c,t);
       dS[eqn] = 0;
       return vof_source;
}

DEFINE_SOURCE(my_energy_source, c, t, dS, eqn)
{
       real x[ND_ND], energy_source, temp, d_b;
       real term1, term11, term12, term2, term3, term4;
       real phi, D_ttemp;
       real d_dufour = 0.00000001;
       temp = C_T(c, t);
       phi = C_UDSI(c,t,0);
       d_b = (k_b * temp)/(3*pi*dyvisco_nf*d_p); 
       term11 = (C_UDSI_G(c,t,0)[0] * C_T_G(c,t)[0]);
       term12 = (C_UDSI_G(c,t,0)[1] * C_T_G(c,t)[1]);
       term1 = rho_p * C_P_p * d_b * (term11 + term12);
       D_ttemp = phi/temp;
       term2 = rho_p * C_P_p * (D_ttemp) * (pow (C_T_G(c,t)[0],2) + pow (C_T_G(c,t)[0],2));
       C_UDSI(c,t,3) = C_YI_G(c,t,0)[0]; 
       C_UDSI(c,t,4) = C_YI_G(c,t,0)[1];
       term3 = density_nf * spheat_nf * d_dufour * (C_UDSI_G(c,t,3)[0]);
       term4 = density_nf * spheat_nf * d_dufour * (C_UDSI_G(c,t,4)[1]);
       C_CENTROID(x,c,t);
       energy_source = term1 + term2 + term3 + term4;
       dS[eqn] = 0;
       return energy_source;
}

DEFINE_SOURCE(my_concentration_source, c, t, dS, eqn)
{
       real x[ND_ND],concentration_source, phi;
       real d_soret = 0.000000001;
       phi = C_UDSI(c,t,0);
       C_UDSI(c,t,5) = C_T_G(c,t)[0];
       C_UDSI(c,t,6) = C_T_G(c,t)[1];
       concentration_source = density_nf * d_soret * ((C_UDSI_G(c,t,5)[0]) + (C_UDSI_G(c,t,6)[1]));
       C_CENTROID(x,c,t);
       dS[eqn] = 0;
       return concentration_source;
}
