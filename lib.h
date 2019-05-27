double square (double x );
double f_funcx(double cons,double x,double y);
double f_funcy(double cons,double x,double y);
double v_func(double v);
void RK (double& cons,double& x,double& y,double& vx,double& vy, double& h);
void leapfrog (double &cons,double& x,double& y,double& vx,double& vy, double& h);
void Euler (double &cons,double& x,double& y,double& vx,double& vy, double& h);
void Mid_point_RK2 (double &cons,double& x,double& y,double& vx,double& vy, double& h);
