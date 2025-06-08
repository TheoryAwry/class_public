double phi_init;
double S_init;
double lambda;
double V_lambda;
double m_S;
double beta;
double alpha_s;
double m0_s;


enum background_quantities {
  ...
  phi_D_fld,
  phi_D_prime_fld,
  S_fld,
  S_prime_fld,
  ...
};

y[phi_D_fld] = pba->phi_init;
y[phi_D_prime_fld] = 0.0;
y[S_fld] = pba->S_init;
y[S_prime_fld] = 0.0;

// Virtualization/dark field potentials
double V_phi = pba->V_lambda * exp(-pba->lambda * y[phi_D_fld]);
double dV_dphi = -pba->lambda * V_phi;

double V_S = 0.5 * pba->m_S * pba->m_S * y[S_fld] * y[S_fld];
double dV_dS = pba->m_S * pba->m_S * y[S_fld];

// Equations of motion
dy[phi_D_fld] = y[phi_D_prime_fld];
dy[phi_D_prime_fld] = -3.0 * aH * y[phi_D_prime_fld] - dV_dphi + pba->beta * y[S_fld];

dy[S_fld] = y[S_prime_fld];
dy[S_prime_fld] = -3.0 * aH * y[S_prime_fld] - dV_dS + pba->beta * y[phi_D_fld];

double rho_phi = 0.5 * pow(y[phi_D_prime_fld], 2) + V_phi;
double p_phi = 0.5 * pow(y[phi_D_prime_fld], 2) - V_phi;

double rho_S = 0.5 * pow(y[S_prime_fld], 2) + V_S;
double p_S = 0.5 * pow(y[S_prime_fld], 2) - V_S;

*pba->rho_tot += rho_phi + rho_S;
*pba->p_tot += p_phi + p_S;

double m_s = pba->m0_s + pba->alpha_s * y[S_fld];

class_read_double("phi_init", pba->phi_init);
class_read_double("S_init", pba->S_init);
class_read_double("lambda", pba->lambda);
class_read_double("V_lambda", pba->V_lambda);
class_read_double("m_S", pba->m_S);
class_read_double("beta", pba->beta);
class_read_double("alpha_s", pba->alpha_s);
class_read_double("m0_s", pba->m0_s);

make clean && make
./class explanatory/your_model.ini

