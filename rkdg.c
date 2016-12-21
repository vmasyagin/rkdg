#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define M_PI		3.14159265358979323846
#define R_GAS       8.314472

#define N_BF 3
#define NX 80
#define NY 800

#define N_GP_EDGE 2
#define N_GP_CELL 4

#define XMIN  0.0
#define XMAX  0.04
#define YMIN -0.1
#define YMAX  0.3

#define TAU   5.0e-8
#define TMAX  2.5e-3

#define SAVE_STEP 100

#define LIM_ALPHA 2.0

#define _SIGN_(X) (fabs(X)/(X))
#define _MIN_(X,Y) ((X)<(Y) ? (X) : (Y))
#define _MAX_(X,Y) ((X)>(Y) ? (X) : (Y))

typedef struct {
    union {
        double flds[4][N_BF];
        struct {
            double ro[N_BF];
            double ru[N_BF];
            double rv[N_BF];
            double re[N_BF];
        };
    };
} data_t;

typedef struct {
    double r, p, u, v, e, e_tot, t, cz;
} prim_t;

typedef struct {
    double x,y;
} point_t;


double HX, HY;

data_t data[NX][NY], r_int[NX][NY];

point_t gp_edg_x[NX+1][NY][N_GP_EDGE];
point_t gp_edg_y[NX][NY+1][N_GP_EDGE];
double gw_edg_x[NX+1][NY][N_GP_EDGE];
double gw_edg_y[NX][NY+1][N_GP_EDGE];
double gj_edg_x[NX+1][NY];
double gj_edg_y[NX][NY+1];

point_t gp_cell[NX][NY][N_GP_CELL];
double gw_cell[NX][NY][N_GP_CELL];
double gj_cell[NX][NY];
point_t c_cell[NX][NY];

double matr_a[NX][NY][N_BF][N_BF];

void init();
void calc_flx();
void calc_int();
void calc_new();
void calc_lim();

void zero_r();

void save_vtk(int step);

double bf    (int i_func, int i, int j, double x, double y);
double bf_dx (int i_func, int i, int j, double x, double y);
double bf_dy (int i_func, int i, int j, double x, double y);

double get_fld(int i_fld, int i, int j, double x, double y);

void cons_to_prim(int i, int j, double x, double y, prim_t *prim);



int main(int argc, char** argv)
{
    double t = 0.0;
    int step = 0;
    init();
    save_vtk(step);
    while (t < TMAX) {
        t += TAU;
        step++;
        zero_r();
        calc_int();
        calc_flx();
        calc_new();
        calc_lim();
        if (step % SAVE_STEP == 0) save_vtk(step);
    }
    return 0;
}

void inverse_matr(double a_src[N_BF][N_BF], double am[N_BF][N_BF], int N)
{
	double a[N_BF][N_BF];
    for (int i = 0; i < N_BF; i++) {
        for (int j = 0; j < N_BF; j++) {
            a[i][j] = a_src[i][j];
        }
    }
	double detA = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[0][1] * a[1][2] * a[2][0]
				- a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];

    double m[3][3];
	m[0][0] = a[1][1] * a[2][2] - a[2][1] * a[1][2];
	m[0][1] = a[2][0] * a[1][2] - a[1][0] * a[2][2];
	m[0][2] = a[1][0] * a[2][1] - a[2][0] * a[1][1];
	m[1][0] = a[2][1] * a[0][2] - a[0][1] * a[2][2];
	m[1][1] = a[0][0] * a[2][2] - a[2][0] * a[0][2];
	m[1][2] = a[2][0] * a[0][1] - a[0][0] * a[2][1];
	m[2][0] = a[0][1] * a[1][2] - a[1][1] * a[0][2];
	m[2][1] = a[1][0] * a[0][2] - a[0][0] * a[1][2];
	m[2][2] = a[0][0] * a[1][1] - a[1][0] * a[0][1];

	am[0][0] = m[0][0] / detA;
	am[0][1] = m[1][0] / detA;
	am[0][2] = m[2][0] / detA;
	am[1][0] = m[0][1] / detA;
	am[1][1] = m[1][1] / detA;
	am[1][2] = m[2][1] / detA;
	am[2][0] = m[0][2] / detA;
	am[2][1] = m[1][2] / detA;
	am[2][2] = m[2][2] / detA;
}

void mult_matr_vec(double matr[N_BF][N_BF], double vec[N_BF], double res[N_BF])
{
    for (int i = 0; i < N_BF; i++) {
        res[i] = 0.0;
        for (int j = 0; j < N_BF; j++) {
            res[i] += matr[i][j]*vec[j];
        }
    }
}

void zero_r()
{
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            memset(&r_int[i][j], 0, sizeof(data_t));
        }
    }
}

void init()
{
    HX = (XMAX-XMIN)/NX;
    HY = (YMAX-YMIN)/NY;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            point_t p[N_GP_CELL];
            p[0].x = -1.0/sqrt(3.0);
            p[0].y = -1.0/sqrt(3.0);
            p[1].x =  1.0/sqrt(3.0);
            p[1].y = -1.0/sqrt(3.0);
            p[2].x =  1.0/sqrt(3.0);
            p[2].y =  1.0/sqrt(3.0);
            p[3].x = -1.0/sqrt(3.0);
            p[3].y =  1.0/sqrt(3.0);
            double xmin = XMIN+i*HX;
            double ymin = YMIN+j*HY;
            double xmax = xmin+HX;
            double ymax = ymin+HY;
            for (int i_gp = 0; i_gp < N_GP_CELL; i_gp++) {
                gp_cell[i][j][i_gp].x = 0.5*(xmin+xmax)+p[i_gp].x*(xmax-xmin)*0.5;
                gp_cell[i][j][i_gp].y = 0.5*(ymin+ymax)+p[i_gp].y*(ymax-ymin)*0.5;
                gw_cell[i][j][i_gp] = 1.0;
                gj_cell[i][j] = 0.25*(xmax-xmin)*(ymax-ymin);
            }

            c_cell[i][j].x = 0.5*(xmin+xmax);
            c_cell[i][j].y = 0.5*(ymin+ymax);
        }
    }

    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j < NY; j++) {
            double xmin = XMIN+i*HX;
            double ymin = YMIN+j*HY;
            double ymax = ymin+HY;
            double p[] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
            for (int i_gp = 0; i_gp < N_GP_EDGE; i_gp++) {
                gp_edg_x[i][j][i_gp].x = xmin;
                gp_edg_x[i][j][i_gp].y = 0.5*(ymin+ymax)+p[i_gp]*(ymax-ymin)*0.5;
                gw_edg_x[i][j][i_gp] = 1.0;
                gj_edg_x[i][j] = 0.5*(ymax-ymin);
            }
        }
    }

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j <= NY; j++) {
            double xmin = XMIN+i*HX;
            double ymin = YMIN+j*HY;
            double xmax = xmin+HX;
            double p[] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
            for (int i_gp = 0; i_gp < N_GP_EDGE; i_gp++) {
                gp_edg_y[i][j][i_gp].x = 0.5*(xmin+xmax)+p[i_gp]*(xmax-xmin)*0.5;;
                gp_edg_y[i][j][i_gp].y = ymin;
                gw_edg_y[i][j][i_gp] = 1.0;
                gj_edg_y[i][j] = 0.5*(xmax-xmin);
            }
        }
    }

    double matr[N_BF][N_BF];
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int m = 0; m < N_BF; m++) {
                for (int l = 0; l < N_BF; l++) {
                    matr[m][l] = 0.0;
                    for (int i_gp = 0; i_gp < N_GP_CELL; i_gp++) {
                        matr[m][l] += gw_cell[i][j][i_gp]*bf(m, i, j, gp_cell[i][j][i_gp].x, gp_cell[i][j][i_gp].y)
                            *bf(l, i, j, gp_cell[i][j][i_gp].x, gp_cell[i][j][i_gp].y);
                    }
                    matr[m][l] *= gj_cell[i][j];
                }
            }
            inverse_matr(matr, matr_a[i][j], N_BF);
        }
    }

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            double x = XMIN+(i+0.5)*HX;
            double y = YMIN+(j+0.5)*HY;
            memset(&data[i][j], 0, sizeof(data_t));
            double r, p, u, v, cp, m_mol;
            if (y < -0.005) {
                r = 12.090;
                p = 2.152e+5;
                u = 0.0;
                v = 97.76;
                cp = 1014.16;
                m_mol = 0.02869409;
            }
            else if (y > HY*sin(M_PI*x/0.01)) {
                r = 1.198;
                p = 1.e+5;
                u = 0.0;
                v = 0.0;
                cp = 1014.16;
                m_mol = 0.02869409;
            }
            else {
                r = 6.037;
                p = 1.e+5;
                u = 0.0;
                v = 0.0;
                cp = 1014.16;
                m_mol = 0.02869409;
            }
            double cv = cp-R_GAS/m_mol;
            double gam = cp/cv;
            data[i][j].ro[0] = r;
            data[i][j].ru[0] = r*u;
            data[i][j].rv[0] = r*v;
            data[i][j].re[0] = p/(gam-1.0)+r*(u*u+v*v)*0.5;
        }
    }

}

void calc_int()
{
    double fint[4][N_BF];
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(fint[i_fld], 0, sizeof(double)*N_BF);
            }
            for (int i_gp = 0; i_gp < N_GP_CELL; i_gp++) {
                point_t pt = gp_cell[i][j][i_gp];
                prim_t par;
                cons_to_prim(i, j, pt.x, pt.y, &par);
                double f[4], g[4];
                f[0] = par.r*par.u;
    			f[1] = f[0]*par.u+par.p;
    			f[2] = f[0]*par.v;
    			f[3] = par.u*(par.r*par.e_tot+par.p);

    			g[0] = par.r*par.v;
    			g[1] = g[0]*par.u;
    			g[2] = g[0]*par.v+par.p;
    			g[3] = par.v*(par.r*par.e_tot+par.p);

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < N_BF; k++) {
                        fint[i_fld][k] += gw_cell[i][j][i_gp]*f[i_fld]*bf_dx(k, i, j, pt.x, pt.y);
                        fint[i_fld][k] += gw_cell[i][j][i_gp]*g[i_fld]*bf_dy(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    fint[i_fld][k] *= gj_cell[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    r_int[i][j].flds[i_fld][k] += fint[i_fld][k];
                }
            }


        }
    }
}

void calc_flx()
{
    // X direction
    for (int i = 1; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            double int_m[4][N_BF], int_p[4][N_BF];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_m[i_fld], 0, sizeof(double)*N_BF);
                memset(int_p[i_fld], 0, sizeof(double)*N_BF);
            }
            for (int i_gp = 0; i_gp < N_GP_EDGE; i_gp++) {
                point_t pt = gp_edg_x[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i-1, j, pt.x, pt.y, &par_m);
                cons_to_prim(i,   j, pt.x, pt.y, &par_p);
                double alpha = _MAX_(par_m.cz+sqrt(par_m.u*par_m.u+par_m.v*par_m.v),
                                     par_p.cz+sqrt(par_p.u*par_p.u+par_p.v*par_p.v));
                flx[0] = 0.5*((par_p.r*par_p.u+par_m.r*par_m.u)-alpha*(par_p.r-par_m.r));
                flx[1] = 0.5*((par_p.r*par_p.u*par_p.u+par_p.p+par_m.r*par_m.u*par_m.u+par_m.p)
                       - alpha*(par_p.r*par_p.u-par_m.r*par_m.u));
                flx[2] = 0.5*((par_p.r*par_p.u*par_p.v+par_m.r*par_m.u*par_m.v)
                       - alpha*(par_p.r*par_p.v-par_m.r*par_m.v));
                flx[3] = 0.5*(((par_p.r*par_p.e_tot+par_p.p)*par_p.u+(par_m.r*par_m.e_tot+par_m.p)*par_m.u)
                       - alpha*(par_p.r*par_p.e_tot-par_m.r*par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < N_BF; k++) {
                        int_m[i_fld][k] += gw_edg_x[i][j][i_gp]*flx[i_fld]*bf(k, i-1, j, pt.x, pt.y);
                        int_p[i_fld][k] += gw_edg_x[i][j][i_gp]*flx[i_fld]*bf(k, i,   j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    int_m[i_fld][k] *= gj_edg_x[i][j];
                    int_p[i_fld][k] *= gj_edg_x[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    r_int[i-1][j].flds[i_fld][k] -= int_m[i_fld][k];
                    r_int[i  ][j].flds[i_fld][k] += int_p[i_fld][k];
                }
            }
        }
    }

    // Y direction
    for (int i = 0; i < NX; i++) {
        for (int j = 1; j < NY; j++) {
            double int_m[4][N_BF], int_p[4][N_BF];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_m[i_fld], 0, sizeof(double)*N_BF);
                memset(int_p[i_fld], 0, sizeof(double)*N_BF);
            }
            for (int i_gp = 0; i_gp < N_GP_EDGE; i_gp++) {
                point_t pt = gp_edg_y[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i, j-1, pt.x, pt.y, &par_m);
                cons_to_prim(i, j,   pt.x, pt.y, &par_p);
                double alpha = _MAX_(par_m.cz+sqrt(par_m.u*par_m.u+par_m.v*par_m.v),
                                     par_p.cz+sqrt(par_p.u*par_p.u+par_p.v*par_p.v));
                flx[0] = 0.5*((par_p.r*par_p.v+par_m.r*par_m.v)-alpha*(par_p.r-par_m.r));
                flx[1] = 0.5*((par_p.r*par_p.u*par_p.v+par_m.r*par_m.u*par_m.v)
                       - alpha*(par_p.r*par_p.u-par_m.r*par_m.u));
                flx[2] = 0.5*((par_p.r*par_p.v*par_p.v+par_p.p+par_m.r*par_m.v*par_m.v+par_m.p)
                       - alpha*(par_p.r*par_p.v-par_m.r*par_m.v));
                flx[3] = 0.5*(((par_p.r*par_p.e_tot+par_p.p)*par_p.v+(par_m.r*par_m.e_tot+par_m.p)*par_m.v)
                       - alpha*(par_p.r*par_p.e_tot-par_m.r*par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < N_BF; k++) {
                        int_m[i_fld][k] += gw_edg_y[i][j][i_gp]*flx[i_fld]*bf(k, i, j-1, pt.x, pt.y);
                        int_p[i_fld][k] += gw_edg_y[i][j][i_gp]*flx[i_fld]*bf(k, i, j,   pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    int_m[i_fld][k] *= gj_edg_y[i][j];
                    int_p[i_fld][k] *= gj_edg_y[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    r_int[i][j-1].flds[i_fld][k] -= int_m[i_fld][k];
                    r_int[i][j  ].flds[i_fld][k] += int_p[i_fld][k];
                }
            }
        }
    }

    // Left
    for (int i = 0; i <= 0; i++) {
        for (int j = 0; j < NY; j++) {
            double int_m[4][N_BF], int_p[4][N_BF];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_p[i_fld], 0, sizeof(double)*N_BF);
            }
            for (int i_gp = 0; i_gp < N_GP_EDGE; i_gp++) {
                point_t pt = gp_edg_x[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp-R_GAS/m_mol;
                    double gam = cp/cv;

                    par_m.r =  par_p.r;
                    par_m.u = -par_p.u;
                    par_m.v =  par_p.v;
                    par_m.p =  par_p.p;
                    par_m.e =  par_m.p/(par_m.r*(gam-1.0));
                    par_m.e_tot = par_m.e+(par_m.u*par_m.u+par_m.v*par_m.v)*0.5;
                    par_m.cz = sqrt(gam*par_m.p/par_m.r);
                    par_m.t = par_m.e/cv;

                }
                double alpha = _MAX_(par_m.cz+sqrt(par_m.u*par_m.u+par_m.v*par_m.v),
                                     par_p.cz+sqrt(par_p.u*par_p.u+par_p.v*par_p.v));
                flx[0] = 0.5*((par_p.r*par_p.u+par_m.r*par_m.u)-alpha*(par_p.r-par_m.r));
                flx[1] = 0.5*((par_p.r*par_p.u*par_p.u+par_p.p+par_m.r*par_m.u*par_m.u+par_m.p)
                              - alpha*(par_p.r*par_p.u-par_m.r*par_m.u));
                flx[2] = 0.5*((par_p.r*par_p.u*par_p.v+par_m.r*par_m.u*par_m.v)
                              - alpha*(par_p.r*par_p.v-par_m.r*par_m.v));
                flx[3] = 0.5*(((par_p.r*par_p.e_tot+par_p.p)*par_p.u+(par_m.r*par_m.e_tot+par_m.p)*par_m.u)
                              - alpha*(par_p.r*par_p.e_tot-par_m.r*par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < N_BF; k++) {
                        int_p[i_fld][k] += gw_edg_x[i][j][i_gp]*flx[i_fld]*bf(k, i,   j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    int_p[i_fld][k] *= gj_edg_x[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    r_int[i][j].flds[i_fld][k] += int_p[i_fld][k];
                }
            }
        }
    }

    // Right
    for (int i = NX; i <= NX; i++) {
        for (int j = 0; j < NY; j++) {
            double int_m[4][N_BF], int_p[4][N_BF];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_m[i_fld], 0, sizeof(double)*N_BF);
            }
            for (int i_gp = 0; i_gp < N_GP_EDGE; i_gp++) {
                point_t pt = gp_edg_x[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i-1, j, pt.x, pt.y, &par_m);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp-R_GAS/m_mol;
                    double gam = cp/cv;

                    par_p.r =  par_m.r;
                    par_p.u = -par_m.u;
                    par_p.v =  par_m.v;
                    par_p.p =  par_m.p;
                    par_p.e =  par_p.p/(par_p.r*(gam-1.0));
                    par_p.e_tot = par_p.e+(par_p.u*par_p.u+par_p.v*par_p.v)*0.5;
                    par_p.cz = sqrt(gam*par_p.p/par_p.r);
                    par_p.t = par_p.e/cv;

                }
                double alpha = _MAX_(par_m.cz+sqrt(par_m.u*par_m.u+par_m.v*par_m.v),
                                     par_p.cz+sqrt(par_p.u*par_p.u+par_p.v*par_p.v));
                flx[0] = 0.5*((par_p.r*par_p.u+par_m.r*par_m.u)-alpha*(par_p.r-par_m.r));
                flx[1] = 0.5*((par_p.r*par_p.u*par_p.u+par_p.p+par_m.r*par_m.u*par_m.u+par_m.p)
                              - alpha*(par_p.r*par_p.u-par_m.r*par_m.u));
                flx[2] = 0.5*((par_p.r*par_p.u*par_p.v+par_m.r*par_m.u*par_m.v)
                              - alpha*(par_p.r*par_p.v-par_m.r*par_m.v));
                flx[3] = 0.5*(((par_p.r*par_p.e_tot+par_p.p)*par_p.u+(par_m.r*par_m.e_tot+par_m.p)*par_m.u)
                              - alpha*(par_p.r*par_p.e_tot-par_m.r*par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < N_BF; k++) {
                        int_m[i_fld][k] += gw_edg_x[i][j][i_gp]*flx[i_fld]*bf(k, i-1, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    int_m[i_fld][k] *= gj_edg_x[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    r_int[i-1][j].flds[i_fld][k] -= int_m[i_fld][k];
                }
            }
        }
    }

    // Bottom
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j <= 0; j++) {
            double int_m[4][N_BF], int_p[4][N_BF];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_p[i_fld], 0, sizeof(double)*N_BF);
            }
            for (int i_gp = 0; i_gp < N_GP_EDGE; i_gp++) {
                point_t pt = gp_edg_y[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i, j,   pt.x, pt.y, &par_p);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp-R_GAS/m_mol;
                    double gam = cp/cv;

                    par_m.r = 12.090;
                    par_m.p = 2.152e+5;
                    par_m.u = 0.0;
                    par_m.v = 97.76;
                    par_m.e =  par_m.p/(par_m.r*(gam-1.0));
                    par_m.e_tot = par_m.e+(par_m.u*par_m.u+par_m.v*par_m.v)*0.5;
                    par_m.cz = sqrt(gam*par_m.p/par_m.r);
                    par_m.t = par_m.e/cv;

                }
                double alpha = _MAX_(par_m.cz+sqrt(par_m.u*par_m.u+par_m.v*par_m.v),
                                     par_p.cz+sqrt(par_p.u*par_p.u+par_p.v*par_p.v));
                flx[0] = 0.5*((par_p.r*par_p.v+par_m.r*par_m.v)-alpha*(par_p.r-par_m.r));
                flx[1] = 0.5*((par_p.r*par_p.u*par_p.v+par_m.r*par_m.u*par_m.v)
                              - alpha*(par_p.r*par_p.u-par_m.r*par_m.u));
                flx[2] = 0.5*((par_p.r*par_p.v*par_p.v+par_p.p+par_m.r*par_m.v*par_m.v+par_m.p)
                              - alpha*(par_p.r*par_p.v-par_m.r*par_m.v));
                flx[3] = 0.5*(((par_p.r*par_p.e_tot+par_p.p)*par_p.v+(par_m.r*par_m.e_tot+par_m.p)*par_m.v)
                              - alpha*(par_p.r*par_p.e_tot-par_m.r*par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < N_BF; k++) {
                        int_p[i_fld][k] += gw_edg_y[i][j][i_gp]*flx[i_fld]*bf(k, i, j,   pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    int_p[i_fld][k] *= gj_edg_y[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    r_int[i][j].flds[i_fld][k] += int_p[i_fld][k];
                }
            }
        }
    }


    // Top
    for (int i = 0; i < NX; i++) {
        for (int j = NY; j <= NY; j++) {
            double int_m[4][N_BF], int_p[4][N_BF];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_m[i_fld], 0, sizeof(double)*N_BF);
            }
            for (int i_gp = 0; i_gp < N_GP_EDGE; i_gp++) {
                point_t pt = gp_edg_y[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i, j-1, pt.x, pt.y, &par_m);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp-R_GAS/m_mol;
                    double gam = cp/cv;

                    par_p.r =  par_m.r;
                    par_p.u =  par_m.u;
                    par_p.v = -par_m.v;
                    par_p.p =  par_m.p;
                    par_p.e =  par_p.p/(par_p.r*(gam-1.0));
                    par_p.e_tot = par_p.e+(par_p.u*par_p.u+par_p.v*par_p.v)*0.5;
                    par_p.cz = sqrt(gam*par_p.p/par_p.r);
                    par_p.t = par_p.e/cv;

                }
                double alpha = _MAX_(par_m.cz+sqrt(par_m.u*par_m.u+par_m.v*par_m.v),
                                     par_p.cz+sqrt(par_p.u*par_p.u+par_p.v*par_p.v));
                flx[0] = 0.5*((par_p.r*par_p.v+par_m.r*par_m.v)-alpha*(par_p.r-par_m.r));
                flx[1] = 0.5*((par_p.r*par_p.u*par_p.v+par_m.r*par_m.u*par_m.v)
                              - alpha*(par_p.r*par_p.u-par_m.r*par_m.u));
                flx[2] = 0.5*((par_p.r*par_p.v*par_p.v+par_p.p+par_m.r*par_m.v*par_m.v+par_m.p)
                              - alpha*(par_p.r*par_p.v-par_m.r*par_m.v));
                flx[3] = 0.5*(((par_p.r*par_p.e_tot+par_p.p)*par_p.v+(par_m.r*par_m.e_tot+par_m.p)*par_m.v)
                              - alpha*(par_p.r*par_p.e_tot-par_m.r*par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < N_BF; k++) {
                        int_m[i_fld][k] += gw_edg_y[i][j][i_gp]*flx[i_fld]*bf(k, i, j-1, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    int_m[i_fld][k] *= gj_edg_y[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < N_BF; k++) {
                    r_int[i][j-1].flds[i_fld][k] -= int_m[i_fld][k];
                }
            }
        }
    }
}

void calc_new()
{
    double tmp[N_BF];
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                mult_matr_vec(matr_a[i][j], r_int[i][j].flds[i_fld], tmp);
                for (int k = 0; k < N_BF; k++) {
                    data[i][j].flds[i_fld][k] += TAU*tmp[k];
                }
            }
        }
    }
}

double _minmod(double a, double b, double c) {
    if ((_SIGN_(a) == _SIGN_(b)) && (_SIGN_(b) == _SIGN_(c))) {
        return _SIGN_(a)*_MIN_(_MIN_(fabs(a), fabs(b)), fabs(c));
    }
    else {
        return 0.0;
    }
}

void calc_lim()
{

    for (int i = 1; i < NX-1; i++) {
        for (int j = 0; j < NY; j++) {
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                double fm = data[i-1][j].flds[i_fld][0];
                double fc = data[i  ][j].flds[i_fld][0];
                double fp = data[i+1][j].flds[i_fld][0];
                data[i][j].flds[i_fld][1] = _minmod(data[i][j].flds[i_fld][1], LIM_ALPHA*(fc-fm), LIM_ALPHA*(fp-fc));
            }
        }
    }
    for (int i = 0; i < NX; i++) {
        for (int j = 1; j < NY-1; j++) {
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                double fm = data[i][j-1].flds[i_fld][0];
                double fc = data[i][j  ].flds[i_fld][0];
                double fp = data[i][j+1].flds[i_fld][0];
                data[i][j].flds[i_fld][2] = _minmod(data[i][j].flds[i_fld][2], LIM_ALPHA*(fc-fm), LIM_ALPHA*(fp-fc));
            }
        }
    }
}

double bf    (int i_func, int i, int j, double x, double y)
{
    switch (i_func) {
        case 0:
            return 1.0;
            break;
        case 1:
            return (x-c_cell[i][j].x)/HX;
            break;
        case 2:
            return (y-c_cell[i][j].y)/HY;
            break;
    }
}

double bf_dx  (int i_func, int i, int j, double x, double y)
{
    switch (i_func) {
        case 0:
            return 0.0;
            break;
        case 1:
            return 1.0/HX;
            break;
        case 2:
            return 0.0;
            break;
    }
}

double bf_dy    (int i_func, int i, int j, double x, double y)
{
    switch (i_func) {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 1.0/HY;
            break;
    }
}

double get_fld(int i_fld, int i, int j, double x, double y)
{
    double res = 0.0;
    for (int i_bf = 0; i_bf < N_BF; i_bf++) {
        res += data[i][j].flds[i_fld][i_bf]*bf(i_bf, i, j, x, y);
    }
    return res;
}

void cons_to_prim(int i, int j, double x, double y, prim_t *prim)
{
    double cp = 1014.16;
    double m_mol = 0.02869409;
    double cv = cp-R_GAS/m_mol;
    double gam = cp/cv;

    double ro = get_fld(0, i, j, x, y);
    double ru = get_fld(1, i, j, x, y);
    double rv = get_fld(2, i, j, x, y);
    double re = get_fld(3, i, j, x, y);

    prim->r = ro;
    prim->u = ru/ro;
    prim->v = rv/ro;
    prim->e_tot = re/ro;
    prim->e = prim->e_tot-(prim->u*prim->u+prim->v*prim->v)*0.5;
    prim->p = prim->r*prim->e*(gam-1.0);
    prim->cz = sqrt(gam*prim->p/prim->r);
    prim->t = prim->e/cv;
}

void save_vtk(int num)
{
    char fName[50];
    prim_t pr;
    sprintf(fName, "res_%010d.vtk", num);
    FILE * fp = fopen(fName, "w");
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "GASDIN data file\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", NX+1, NY+1, 1);
    fprintf(fp, "POINTS %d float\n", (NX+1)*(NY+1));
    for (int j = 0; j <= NY; j++) {
        for (int i = 0; i <= NX; i++) {
            double x = XMIN+i*HX;
            double y = YMIN+j*HY;
            fprintf(fp, "%f %f %f\n", x,y,0.0);
        }
    }
    fprintf(fp, "CELL_DATA %d\n", NX*NY);
    fprintf(fp, "SCALARS Density float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            cons_to_prim(i, j, c_cell[i][j].x, c_cell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.r);
        }
    }

    fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            cons_to_prim(i, j, c_cell[i][j].x, c_cell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.p);
        }
    }

    fprintf(fp, "SCALARS TotEnergy float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            cons_to_prim(i, j, c_cell[i][j].x, c_cell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.e_tot);
        }
    }

    fprintf(fp, "SCALARS Energy float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            cons_to_prim(i, j, c_cell[i][j].x, c_cell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.e);
        }
    }

    fprintf(fp, "VECTORS Velosity float\n");
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            cons_to_prim(i, j, c_cell[i][j].x, c_cell[i][j].y, &pr);
            fprintf(fp, "%f %f %f\n", pr.u, pr.v, 0.0);
        }
    }

	fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            cons_to_prim(i, j, c_cell[i][j].x, c_cell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.t);
        }
    }

    fclose(fp);
    printf("File '%s' saved...\n", fName);

}
