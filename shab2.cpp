#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <cassert>
#include <string>
#include <sstream>
#define M_PI 3.14159265358979323846 /* pi */

#include<math.h>
int smoothx(double* q, double* xx, double* xy, int id, int jd, double* a,
    double* b, double* c, int j, double* jac, double* drr, double* drp,
    double* rspec, double* qv, double* dd,
    double epsi, double gamma, double fsmach, double dt)
{

    double* rho, * u_vel, * v_vel, * t_e;

    double eratio, gm1, ggm1, eps, ra, u, v, qq, ss,
        qav, qxx, rr, rp;
    int ib, ie, i, offset, offsetp1, offsetm1, ip, ir, n;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001, gamma);
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    ib = 1;
    ie = id - 1;

    rho = q;
    u_vel = &q[id * jd];
    v_vel = &q[2 * id * jd];
    t_e = &q[3 * id * jd];

    /*     smoothing in xi direction */

    for (i = 0; i < id; i++) {
        offset = id * j + i;
        eps = epsi / jac[offset];
        ra = 1. / rho[offset];
        u = u_vel[offset] * ra;
        v = v_vel[offset] * ra;
        qq = u * u + v * v;
        ss = ggm1 * (t_e[offset] * ra - 0.5 * qq);
        rspec[i] = eps * (fabs(xx[offset] * u + xy[offset] * v) + sqrt((xx[offset] * xx[offset] + xy[offset] * xy[offset]) * ss + 0.01));
        qv[i] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
    }

    for (i = ib; i < ie; i++) {
        ip = i + 1;
        ir = i - 1;
        qxx = qv[ip] - qv[i] * 2. + qv[ir];
        qav = (qv[ip] + qv[i] * 2. + qv[ir]) * .25;
        dd[i] = eratio * fabs(qxx / qav);
    }

    dd[0] = dd[1];
    dd[id - 1] = dd[id - 2];

    for (i = ib; i < ie; i++) {
        ip = i + 1;
        ir = i - 1;
        offset = j * id + i;
        offsetp1 = j * id + i + 1;
        offsetm1 = j * id + i - 1;
        rp = (0.5 * (dd[ip] + dd[i]) + 2.5) * dt * 0.5 * (rspec[ip] + rspec[i]);
        rr = (0.5 * (dd[ir] + dd[i]) + 2.5) * dt * 0.5 * (rspec[ir] + rspec[i]);
        qv[i] = (rr + rp) * jac[offset];
        drr[i] = rr * jac[offsetm1];
        drp[i] = rp * jac[offsetp1];
    }

    for (n = 0; n < 4; n++) {
        for (i = ib; i < ie; i++) {
            offset = (n * 4 + n) * id + i;
            a[offset] -= drr[i];
            b[offset] += qv[i];
            c[offset] -= drp[i];
        }
    }
    return 0;
}

int smoothy(double* q, double* yx, double* yy, int id, int jd, double* a,
    double* b, double* c, int i, double* jac, double* drr, double* drp,
    double* rspec, double* qv, double* dd,
    double epsi, double gamma, double fsmach, double dt)
{

    double* rho, * u_vel, * v_vel, * t_e;

    double eratio, smool, gm1, ggm1, eps, ra, u, v, qq, ss,
        qav, ssfs, qyy, rp, rr;
    int jb, je, j, offset, offsetp1, offsetm1, n,
        jp, jr;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001, gamma);
    smool = 1.0;
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    jb = 1;
    je = jd - 1;

    rho = q;
    u_vel = &q[id * jd];
    v_vel = &q[2 * id * jd];
    t_e = &q[3 * id * jd];

    /*     smoothing in eta direction */

    ssfs = 1. / (0.001 + fsmach * fsmach);
    for (j = 0; j < jd; j++) {
        offset = id * j + i;
        eps = epsi / jac[offset];
        ra = 1. / rho[offset];
        u = u_vel[offset] * ra;
        v = v_vel[offset] * ra;
        qq = u * u + v * v;
        ss = ggm1 * (t_e[offset] * ra - 0.5 * qq) * (1.0 - smool) + smool * qq * ssfs;
        rspec[j] = eps * (fabs(yx[offset] * u + yy[offset] * v) + sqrt((yx[offset] * yx[offset] + yy[offset] * yy[offset]) * ss + 0.01));
        qv[j] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
    }

    for (j = jb; j < je; j++) {
        jp = j + 1;
        jr = j - 1;
        qyy = qv[jp] - qv[j] * 2. + qv[jr];
        qav = (qv[jp] + qv[j] * 2. + qv[jr]) * .25;
        dd[j] = eratio * fabs(qyy / qav);
    }

    dd[0] = dd[1];
    dd[jd - 1] = dd[jd - 2];

    for (j = jb; j < je; j++) {
        jp = j + 1;
        jr = j - 1;
        offset = j * id + i;
        offsetp1 = (j + 1) * id + i;
        offsetm1 = (j - 1) * id + i;
        rp = (0.5 * (dd[jp] + dd[j]) + 2.5) * dt * 0.5 * (rspec[jp] + rspec[j]);
        rr = (0.5 * (dd[jr] + dd[j]) + 2.5) * dt * 0.5 * (rspec[jr] + rspec[j]);
        qv[j] = (rr + rp) * jac[offset];
        drr[j] = rr * jac[offsetm1];
        drp[j] = rp * jac[offsetp1];
    }

    for (n = 0; n < 4; n++) {
        for (j = jb; j < je; j++) {
            offset = (n * 4 + n) * jd + j;
            a[offset] -= drr[j];
            b[offset] += qv[j];
            c[offset] -= drp[j];
        }
    }
    return 0;
}

int smooth(double* q, double* s, double* jac, double* xx, double* xy,
    double* yx, double* yy, int id, int jd, double* s2,
    double* rspec, double* qv, double* dd,
    double epse, double gamma, double fsmach, double dt)
{

    double* rho, * u_vel, * v_vel, * t_e;

    double eratio, smool, gm1, ggm1, cx, cy, eps, ra, u, v, qq, ss, st,
        qav, qxx, ssfs, qyy;
    int ib, ie, jb, je, i, j, offset, offsetp1, offsetm1, ip, ir, n,
        jp, jr;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001, gamma);
    smool = 1.0;
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    ib = 1;
    ie = id - 1;
    jb = 1;
    je = jd - 1;

    cx = 2.;
    cy = 1.;

    rho = q;
    u_vel = &q[id * jd];
    v_vel = &q[2 * id * jd];
    t_e = &q[3 * id * jd];

    /*     smoothing in xi direction */

    for (j = jb; j < je; j++) {
        for (i = 0; i < id; i++) {
            offset = id * j + i;
            eps = epse / jac[offset];
            ra = 1. / rho[offset];
            u = u_vel[offset] * ra;
            v = v_vel[offset] * ra;
            qq = u * u + v * v;
            ss = ggm1 * (t_e[offset] * ra - 0.5 * qq);
            rspec[i] = eps * (fabs(xx[offset] * u + xy[offset] * v) + sqrt((xx[offset] * xx[offset] + xy[offset] * xy[offset]) * ss + 0.01));
            qv[i] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
        }

        for (i = ib; i < ie; i++) {
            ip = i + 1;
            ir = i - 1;
            qxx = qv[ip] - qv[i] * 2. + qv[ir];
            qav = (qv[ip] + qv[i] * 2. + qv[ir]) * .25;
            dd[i] = eratio * fabs(qxx / qav);
        }

        dd[0] = dd[1];
        dd[id - 1] = dd[id - 2];

        for (n = 0; n < 4; n++) {
            for (i = ib; i < ie; i++) {
                offset = (jd * n + j) * id + i;
                offsetp1 = (jd * n + j) * id + i + 1;
                offsetm1 = (jd * n + j) * id + i - 1;
                s2[i] = q[offsetp1] - 2.0 * q[offset] + q[offsetm1];
            }

            s2[0] = s2[1] * -1.;
            s2[id - 1] = s2[id - 2] * -1.;

            for (i = ib; i < ie; i++) {
                ip = i + 1;
                ir = i - 1;
                offset = (jd * n + j) * id + i;
                offsetp1 = (jd * n + j) * id + i + 1;
                offsetm1 = (jd * n + j) * id + i - 1;
                st = ((dd[ip] + dd[i]) * .5 * (q[offsetp1] - q[offset]) - cx / (cx + dd[ip] + dd[i]) * (s2[ip] - s2[i])) * (rspec[ip] + rspec[i]) + ((dd[ir] + dd[i]) * .5 * (q[offsetm1] - q[offset]) - cx / (cx + dd[ir] + dd[i]) * (s2[ir] - s2[i])) * (rspec[ir] + rspec[i]);
                s[offset] += st * .5 * dt;
            }
        }
    }

    /*     smoothing in eta direction */

    ssfs = 1. / (0.001 + fsmach * fsmach);
    for (i = ib; i < ie; i++) {
        for (j = 0; j < jd; j++) {
            offset = id * j + i;
            eps = epse / jac[offset];
            ra = 1. / rho[offset];
            u = u_vel[offset] * ra;
            v = v_vel[offset] * ra;
            qq = u * u + v * v;
            ss = ggm1 * (t_e[offset] * ra - 0.5 * qq) * (1.0 - smool) + smool * qq * ssfs;
            rspec[j] = eps * (fabs(yx[offset] * u + yy[offset] * v) + sqrt((yx[offset] * yx[offset] + yy[offset] * yy[offset]) * ss + 0.01));
            qv[j] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
        }

        for (j = jb; j < je; j++) {
            jp = j + 1;
            jr = j - 1;
            qyy = qv[jp] - qv[j] * 2. + qv[jr];
            qav = (qv[jp] + qv[j] * 2. + qv[jr]) * .25;
            dd[j] = eratio * fabs(qyy / qav);
        }

        dd[0] = dd[1];
        dd[jd - 1] = dd[jd - 2];

        for (n = 0; n < 4; n++) {
            for (j = jb; j < je; j++) {
                offset = (jd * n + j) * id + i;
                offsetp1 = (jd * n + j + 1) * id + i;
                offsetm1 = (jd * n + j - 1) * id + i;
                s2[j] = q[offsetp1] - 2.0 * q[offset] + q[offsetm1];
            }

            s2[0] = s2[1] * -1.;
            s2[jd - 1] = s2[jd - 2] * -1.;

            for (j = jb; j < je; j++) {
                jp = j + 1;
                jr = j - 1;
                offset = (jd * n + j) * id + i;
                offsetp1 = (jd * n + j + 1) * id + i;
                offsetm1 = (jd * n + j - 1) * id + i;
                st = ((dd[jp] + dd[j]) * .5 * (q[offsetp1] - q[offset]) - cy / (cy + dd[jp] + dd[j]) * (s2[jp] - s2[j])) * (rspec[jp] + rspec[j]) + ((dd[jr] + dd[j]) * .5 * (q[offsetm1] - q[offset]) - cy / (cy + dd[jr] + dd[j]) * (s2[jr] - s2[j])) * (rspec[jr] + rspec[j]);
                s[offset] += st * .5 * dt;
            }
        }
    }

    return 0;
}

int btri4s(double* a, double* b, double* c, double* f, int kd, int ks, int ke)
{
    /* Local variables */
    int k, m, n, nd, md;

    float c1, d1, d2, d3, d4, c2, c3, c4, b11, b21, b22, b31, b32, b33,
        b41, b42, b43, b44, u12, u13, u14, u23, u24, u34;


    /*   (A,B,C)F = F, F and B are overloaded, solution in F */

    md = 4;
    nd = 4;

    /*   Part 1. Forward block sweep */

    for (k = ks; k <= ke; k++)
    {

        /*      Step 1. Construct L in B */

        if (k != ks)
        {
            for (m = 0; m < md; m++)
            {
                for (n = 0; n < nd; n++)
                {
                    b[k + kd * (m + md * n)] = b[k + kd * (m + md * n)]
                        - a[k + kd * (m + md * 0)] * b[k - 1 + kd * (0 + md * n)]
                        - a[k + kd * (m + md * 1)] * b[k - 1 + kd * (1 + md * n)]
                        - a[k + kd * (m + md * 2)] * b[k - 1 + kd * (2 + md * n)]
                        - a[k + kd * (m + md * 3)] * b[k - 1 + kd * (3 + md * n)];
                }
            }
        }

        /*      Step 2. Compute L inverse (block matrix) */

        /*          A. Decompose L into L and U */

        b11 = 1. / b[k + kd * (0 + md * 0)];
        u12 = b[k + kd * (0 + md * 1)] * b11;
        u13 = b[k + kd * (0 + md * 2)] * b11;
        u14 = b[k + kd * (0 + md * 3)] * b11;
        b21 = b[k + kd * (1 + md * 0)];
        b22 = 1. / (b[k + kd * (1 + md * 1)] - b21 * u12);
        u23 = (b[k + kd * (1 + md * 2)] - b21 * u13) * b22;
        u24 = (b[k + kd * (1 + md * 3)] - b21 * u14) * b22;
        b31 = b[k + kd * (2 + md * 0)];
        b32 = b[k + kd * (2 + md * 1)] - b31 * u12;
        b33 = 1. / (b[k + kd * (2 + md * 2)] - b31 * u13 - b32 * u23);
        u34 = (b[k + kd * (2 + md * 3)] - b31 * u14 - b32 * u24) * b33;
        b41 = b[k + kd * (3 + md * 0)];
        b42 = b[k + kd * (3 + md * 1)] - b41 * u12;
        b43 = b[k + kd * (3 + md * 2)] - b41 * u13 - b42 * u23;
        b44 = 1. / (b[k + kd * (3 + md * 3)] - b41 * u14 - b42 * u24
            - b43 * u34);

        /*      Step 3. Solve for intermediate vector */

        /*          A. Construct RHS */
        if (k != ks)
        {
            for (m = 0; m < md; m++)
            {
                f[k + kd * m] = f[k + kd * m]
                    - a[k + kd * (m + md * 0)] * f[k - 1 + kd * 0]
                    - a[k + kd * (m + md * 1)] * f[k - 1 + kd * 1]
                    - a[k + kd * (m + md * 2)] * f[k - 1 + kd * 2]
                    - a[k + kd * (m + md * 3)] * f[k - 1 + kd * 3];
            }
        }

        /*          B. Intermediate vector */

        /*          Forward substitution */

        d1 = f[k + kd * 0] * b11;
        d2 = (f[k + kd * 1] - b21 * d1) * b22;
        d3 = (f[k + kd * 2] - b31 * d1 - b32 * d2) * b33;
        d4 = (f[k + kd * 3] - b41 * d1 - b42 * d2 - b43 * d3) * b44;

        /*          Backward substitution */

        f[k + kd * 3] = d4;
        f[k + kd * 2] = d3 - u34 * d4;
        f[k + kd * 1] = d2 - u23 * f[k + kd * 2] - u24 * d4;
        f[k + kd * 0] = d1 - u12 * f[k + kd * 1] - u13 * f[k + kd * 2] - u14 * d4;

        /*      Step 4. Construct U = L ** (-1) * C */
        /*              by columns and store in B */

        if (k != ke)
        {
            for (n = 0; n < nd; n++)
            {

                /*          Forward substitution */

                c1 = c[k + kd * (0 + md * n)] * b11;
                c2 = (c[k + kd * (1 + md * n)] - b21 * c1) * b22;
                c3 = (c[k + kd * (2 + md * n)] - b31 * c1 - b32 * c2) *
                    b33;
                c4 = (c[k + kd * (3 + md * n)] - b41 * c1 - b42 * c2 -
                    b43 * c3) * b44;

                /*          Backward substitution */

                b[k + kd * (3 + md * n)] = c4;
                b[k + kd * (2 + md * n)] = c3 - u34 * c4;
                b[k + kd * (1 + md * n)] = c2 - u23 * b[k + kd * (2 + md * n)] - u24 * c4;
                b[k + kd * (0 + md * n)] = c1 - u12 * b[k + kd * (1 + md * n)]
                    - u13 * b[k + kd * (2 + md * n)] - u14 * c4;
            }
        }
    }

    /*   Part 2. Backward block sweep */

    if (ke == ks)
    {
        return 0;
    }

    for (k = ke - 1; k >= ks; --k)
    {
        for (m = 0; m < md; m++)
        {
            f[k + kd * m] = f[k + kd * m]
                - b[k + kd * (m + md * 0)] * f[k + 1 + kd * 0]
                - b[k + kd * (m + md * 1)] * f[k + 1 + kd * 1]
                - b[k + kd * (m + md * 2)] * f[k + 1 + kd * 2]
                - b[k + kd * (m + md * 3)] * f[k + 1 + kd * 3];
        }
    }

    return 0;

}

// function that reads input parameters
struct Mesh {
    //FSMACH (mach number in infinity), gamma = 1.4, EPSE (smoothing coefficient- around 0.06), deltaT
    double FSMACH;              // = 0.9;
    double p_inf;
    double gamma;               // = 1.4;
    double EPSE;                // = 0.06;
    double deltaT;              // = 1;
    double alpha_inf;               // = 0;
    double rho_inf;

    const int imax;             // = 121;
    const int jmax;             // = 41;

    std::string xFilename;
    std::string yFilename;
    double* xArr;
    double* yArr;
    double* xXiArr;
    double* xEtaArr;
    double* yXiArr;
    double* yEtaArr;
    double* JacobianArr;
    double* invJacobianArr;
    double* XiXArr;
    double* XiYArr;
    double* EtaXArr;
    double* EtaYArr;
    double* uArr;
    double* vArr;
    double* rhoArr;
    double* eArr;
    double* QArr;
    double* SArr;
    double* WArr;
    double* BArr;
    double* AArr;
    double* CArr;
    double* s2;
    double* DArr;
    double* rspec;
    double* qv;
    double* dd;
    double* drr;
    double* drp;

    Mesh(const char* paramFilename = nullptr) : FSMACH(0.9), gamma(1.4), EPSE(0.06), deltaT(1), alpha_inf(0), imax(121), jmax(41), 
        xFilename("x.csv"), yFilename("y.csv") {
        if (paramFilename) {
            parseParams(paramFilename);
        }
        xArr = parseCSV(xFilename.c_str());
        yArr = parseCSV(yFilename.c_str());
        xXiArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        xEtaArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        yXiArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        yEtaArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        JacobianArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        invJacobianArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        XiXArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        XiYArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        EtaXArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        EtaYArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        uArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        vArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        rhoArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        eArr = (double*)malloc((unsigned)(imax * jmax) * sizeof(double));
        QArr = (double*)malloc((unsigned)(imax * jmax * 4) * sizeof(double));
        SArr = (double*)malloc((unsigned)(imax * jmax * 4) * sizeof(double));
        WArr = (double*)malloc((unsigned)(imax * 4) * sizeof(double));
        BArr = (double*)malloc((unsigned)(std::max(imax, jmax) * 16) * sizeof(double));
        AArr = (double*)malloc((unsigned)(std::max(imax, jmax) * 16) * sizeof(double));
        CArr = (double*)malloc((unsigned)(std::max(imax, jmax) * 16) * sizeof(double));
        s2 = (double*)malloc((unsigned)(std::max(imax, jmax)) * sizeof(double));
        DArr = (double*)malloc((unsigned)(std::max(imax, jmax) * 4) * sizeof(double));
        rspec = (double*)malloc((unsigned)(std::max(imax, jmax)) * sizeof(double));
        qv = (double*)malloc((unsigned)(std::max(imax, jmax)) * sizeof(double));
        dd = (double*)malloc((unsigned)(std::max(imax, jmax)) * sizeof(double));
        drr = (double*)malloc((unsigned)(std::max(imax, jmax)) * sizeof(double));
        drp = (double*)malloc((unsigned)(std::max(imax, jmax)) * sizeof(double));
    }

    ~Mesh() {
        delete xArr;
        delete yArr;
        delete xXiArr;
        delete xEtaArr;
        delete yXiArr;
        delete yEtaArr;
        delete JacobianArr;
        delete invJacobianArr;
        delete XiXArr;
        delete XiYArr;
        delete EtaXArr;
        delete EtaYArr;
        delete uArr;
        delete vArr;
        delete rhoArr;
        delete eArr;
        delete QArr;
        delete SArr;
        delete WArr;
        delete BArr;
        delete AArr;
        delete CArr;
        delete s2;
        delete DArr;
        delete rspec;
        delete qv;
        delete dd;
        delete drr;
        delete drp;
    }

    void parseParams(const char* filename) {
        std::ifstream  data(filename);
        if (!data) {
            std::cerr << "failed reading params file" << filename << " you stupid" << std::endl;
        }
        std::cout << "great! you are reading params " << filename << std::endl;
        while (!data.eof()) {
            std::string keyword;
            data >> keyword;
            if (!data.good()) {
                break;
            }
            if (keyword == "FSMACH") { data >> FSMACH; continue; }
            if (keyword == "gamma") { data >> gamma; continue; }
            if (keyword == "EPSE") { data >> EPSE; continue; }
            if (keyword == "deltaT") { data >> deltaT; continue; }
            std::cerr << "unknowm params name " << keyword << std::endl;
            exit(1);
        }
    }

    // function that reads csv files (for the grid)
    double* parseCSV(const char* filename) {
        double* arr = (double*)malloc(jmax * imax * sizeof(double));

        std::ifstream  data(filename);
        std::string line;
        int row = 0;
        while (!data.bad() && std::getline(data, line))
        {
            std::stringstream lineStream(line);
            int col = 0;
            while (col < imax)
            {
                lineStream >> arr[row * imax + col];
                char c = ' ';
                do { lineStream >> c; } while (c != ',' && !lineStream.eof());
                if (col == 5 && row == 9) {
                    std::cout << "row = " << row << ", " << "col = " << col << " => " << arr[row * imax + col] << std::endl;
                }
                col++;
            }
            assert(col == imax);
            row++;
        }
        assert(row == jmax);
        return arr;
    };

    double& x(int i, int j) { return xArr[offset2D(i, j)]; }
    double& y(int i, int j) { return yArr[offset2D(i, j)]; }
    double& xXi(int i, int j) { return xXiArr[offset2D(i, j)]; }
    double& xEta(int i, int j) { return xEtaArr[offset2D(i, j)]; }
    double& yXi(int i, int j) { return yXiArr[offset2D(i, j)]; }
    double& yEta(int i, int j) { return yEtaArr[offset2D(i, j)]; }
    double& Jacobian(int i, int j) { return JacobianArr[offset2D(i, j)]; }
    double& invJacobian(int i, int j) { return invJacobianArr[offset2D(i, j)]; }
    double& XiX(int i, int j) { return XiXArr[offset2D(i, j)]; }
    double& XiY(int i, int j) { return XiYArr[offset2D(i, j)]; }
    double& EtaX(int i, int j) { return EtaXArr[offset2D(i, j)]; }
    double& EtaY(int i, int j) { return EtaYArr[offset2D(i, j)]; }
    double& u(int i, int j) { return uArr[offset2D(i, j)]; }
    double& v(int i, int j) { return vArr[offset2D(i, j)]; }
    double& rho(int i, int j) { return rhoArr[offset2D(i, j)]; }
    double& e(int i, int j) { return eArr[offset2D(i, j)]; }
    double& Q(int i, int j, int k) { return QArr[offset3D(i, j, k)]; }
    double& S(int i, int j, int k) { return SArr[offset3D(i, j, k)]; }
    double& W(int i, int j) { return WArr[offset2D(i, j)]; }
    double& B(int i, int j, int k) { return BArr[offset3D(i, j, k)]; }
    double& A(int i, int j, int k) { return AArr[offset3D(i, j, k)]; }
    double& C(int i, int j, int k) { return CArr[offset3D(i, j, k)]; }
    double& s2(int i, int j) { return s2[offset2D(i, j)]; }
    double& D(int i, int j) { return DArr[offset2D(i, j)]; }



    int offset2D(int i, int j) {
        if (!(i > 0 && i <= jmax)) {
            std::cout << "got bad i = " << i << std::endl;
            assert(true);
        }
        if (!(j > 0 && j <= imax)) {
            std::cout << "got bad j = " << j << std::endl;
            assert(true);
        }
        return (j - 1) * jmax + (i - 1);
    }

    int offset3D(int i, int j, int k) {
        return (k * imax + (j - 1) * jmax + (i - 1));
    }

    void FreeStream() {
        double a;
        double e;
        a = sqrt(gamma * p_inf) / rho_inf;
        e = p_inf / (gamma - 1) + 0.5 * rho_inf * pow(FSMACH * a, 2);
        for (int i = 0; i < imax; i++) {
            for (int j = 0; j < jmax; j++) {
                u(i, j) = FSMACH * a * cos(alpha_inf);
                v(i, j) = FSMACH * a * sin(alpha_inf);
                Q(i, j, 0) = rho_inf;
                Q(i, j, 1) = rho_inf * u(i, j);
                Q(i, j, 2) = rho_inf * v(i, j);
                Q(i, j, 3) = e;
            }
       }
    }

    void MetricCoeffXi(int j) {
        // 2nd order central difference operator
        for (int i = 1; i <= jmax; i++) {
            if (i > 1 && i < jmax) {
                xXi(i, j) = 0.5 * (x(i + 1, j) - x(i - 1, j));
                yXi(i, j) = 0.5 * (y(i + 1, j) - y(i - 1, j));
            }
            else {
                // second order derive
                if (i == 1) {
                    xXi(i, j) = 0.5 * (-x(i + 2, j) + 4 * x(i + 1, j) - 3 * x(i, j));
                    yXi(i, j) = 0.5 * (-y(i + 2, j) + 4 * y(i + 1, j) - 3 * y(i, j));
                }
                else { // i == imax
                    xXi(i, j) = 0.5 * (x(i - 2, j) - 4 * x(i - 1, j) + 3 * x(i, j));
                    yXi(i, j) = 0.5 * (y(i - 2, j) - 4 * y(i - 1, j) + 3 * y(i, j));
                }
            }
        }
    }

    void MetricCoeffEta(int i) {
        // 2nd order central difference operator
        for (int j = 1; j <= imax; j++) {
            if (i > 1 && i < imax) {
                xEta(i, j) = 0.5 * (x(i, j + 1) - x(i, j - 1));
                yEta(i, j) = 0.5 * (y(i, j + 1) - y(i, j - 1));
            }
            else {
                // second order derive
                if (i == 1) {
                    xEta(i, j) = 0.5 * (-x(i, j + 2) + 4 * x(i, j + 1) - 3 * x(i, j));
                    yEta(i, j) = 0.5 * (-y(i, j + 2) + 4 * y(i, j + 1) - 3 * y(i, j));
                }
                else { // i == imax
                    xEta(i, j) = 0.5 * (x(i, j - 2) - 4 * x(i, j - 1) + 3 * x(i, j));
                    yEta(i, j) = 0.5 * (y(i, j - 2) - 4 * y(i, j - 1) + 3 * y(i, j));
                }
            }
        }
    }

    void MetricCoeffFin() {
        for (int i = 1; i <= jmax; i++) {
            MetricCoeffEta(i);
            for (int j = 1; j <= imax; j++) {
                MetricCoeffXi(j);
                Jacobian(i, j) = 1 / (xXi(i, j) * yEta(i, j) - yXi(i, j) * xEta(i, j));
                invJacobian(i, j) = xXi(i, j) * yEta(i, j) - yXi(i, j) * xEta(i, j);
                XiX(i, j) = Jacobian(i, j) * yEta(i, j);
                XiY(i, j) = -Jacobian(i, j) * xEta(i, j);
                EtaX(i, j) = -Jacobian(i, j) * yXi(i, j);
                EtaY(i, j) = Jacobian(i, j) * xXi(i, j);
            }
        }
    }

    void BoundaryConds(int iTEL, int iTEU, int imax, int jmax) {
        double rho01, rhoTEL, rhoTEU, uTEU, uTEL, vTEU, vTEL, p0, pTEL, pTEU, Jacobian0, U;
            // Adiabetic wall: j = jmin, i = itel...iteu
            for (int i = iTEL; i <= iTEU; i++) {
                Q(i, 0, 0) = Q(i, 1, 0);
                rho01 = EtaY(i, 0) * XiX(i, 0) - XiY(i, 0) * EtaX(i, 0);
                U = XiX(i, 1) * Q(i, 1, 1) / rho01 + XiY(i, 1) * Q(i, 1, 2) / rho01;
                Q(i, 0, 1) = EtaY(i, 0) * U * rho01 / Jacobian0;
                Q(i, 0, 2) = -EtaX(i, 0) * U * rho01 / Jacobian0;
                p0 = (gamma - 1) * (Q(i, 1, 3) - 0.5 * rho01 * (pow(Q(i, 1, 1) / rho01, 2) + pow(Q(i, 1, 2) / rho01, 2)));
                Q(i, 0, 3) = p0 / (gamma - 1) + 0.5 * rho01 * (pow(Q(i, 0, 1) / rho01, 2) + pow(Q(i, 0, 2) / rho01, 2));
        }
        // Trailing Edeg
        rhoTEL = Q(iTEL, 0, 0);
        rhoTEU = Q(iTEU, 0, 0);
        uTEL = Q(iTEL, 0, 1) / rhoTEL;
        uTEU = Q(iTEU, 0, 1) / rhoTEU;
        vTEL = Q(iTEL, 0, 2) / rhoTEL;
        vTEU = Q(iTEU, 0, 2) / rhoTEU;
        pTEL = (gamma - 1) * (Q(iTEL, 0, 3) - 0.5 * rhoTEL * (pow(uTEL, 2) + pow(vTEL, 2)));
        pTEU = (gamma - 1) * (Q(iTEU, 0, 3) - 0.5 * rhoTEU * (pow(uTEU, 2) + pow(vTEU, 2)));

        Q(iTEL, 0, 0) = 0.5 * (rhoTEL + rhoTEU);
        Q(iTEL, 0, 1) = 0.5 * (uTEL + uTEU);
        Q(iTEL, 0, 2) = 0.5 * (vTEL + vTEU);
        Q(iTEL, 0, 3) = 0.5 * (pTEL + pTEU);

        for (int k = 0; k < 4; k++) {
            Q(iTEU, 0, k) = Q(iTEL, 0, k);
        }
        // Wake
        for (int i = 0; i < iTEL; i++) {
            for (int k = 0; k < 4; k++) {
                Q(i, 0, k) = 0.5 * (Q(i, 1, k) + Q(imax - 1, 1, k));
                Q(imax - -1 - i, 0, k) = Q(i, 0, k);
            }
        }
        // outflow
        for (int j = 0; j < jmax; j++) {
            for (int k = 0; k < 4; k++) {
                Q(0, j, k) = Q(1, j, k);
                Q(imax - 1, j, k) = Q(imax - 2, j, k);
            }
        }

    }

    void ZeroRHS() {
        for (int i = 0; i < imax; i++) {
            for (int j = 0; j < jmax; j++) {
                for (int k = 0; k < 4; k++) {
                    S(i, j, k) = 0;
                }
            }
        }
    }

    void RHS() {
        //Xi direction
        double UU, VV, p;
        for (int j = 1; j < jmax - 1; j++) {
            for (int i = 0; i < imax; i++) {
                UU = u(i, j) * XiX(i, j) + v(i, j) * XiY(i, j);
                p = (gamma - 1) * (Q(i, j, 3) - 0.5 * Q(i, j, 0) * (pow(Q(i, j, 1) / Q(i, j, 0), 2) + pow(Q(i, j, 2) / Q(i, j, 0), 2)));
                W(i, 0) = Q(i, j, 0) * UU / Jacobian(i, j);
                W(i, 1) = (Q(i, j, 1) * UU + XiX(i, j) * p) / Jacobian(i, j);
                W(i, 2) = (Q(i, j, 2) * UU + XiY(i, j) * p) / Jacobian(i, j);
                W(i, 3) = (Q(i, j, 3) + p) * UU / Jacobian(i, j);
            }
            for (int i = 1; i < imax - 1; i++) {
                for (int k = 0; k < 4; k++) {
                    S(i, j, k) += -0.5*(W(i + 1, k) - W(i - 1, k));
                }
            }
        }
        // Eta direction
        for (int i = 1; i < imax - 1; i++) {
            for (int j = 0; j < jmax; j++) {
                VV = u(i, j) * EtaX(i, j) + v(i, j) * EtaY(i, j);
                p = (gamma -1) * (Q(i, j, 3) - 0.5 * Q(i, j, 0) * (pow(Q(i, j, 1) / Q(i, j, 0), 2) + pow(Q(i, j, 2) / Q(i, j, 0), 2)));
                for (int j = 0; j < jmax; j++) {
                    W(j, 0) = Q(i, j, 0) * VV / Jacobian(i, j);
                    W(j, 1) = (Q(i, j, 1) * VV + EtaX(i, j) * p) / Jacobian(i, j);
                    W(j, 2) = (Q(i, j, 2) * VV + EtaY(i, j) * p) / Jacobian(i, j);
                    W(j, 3) = (Q(i, j, 3) + p) * VV / Jacobian(i, j);
                }
                for (int j = 1; j < jmax - 1; j++) {
                    for (int k = 0; k < 4; k++) {
                        S(i, j, k) += -0.5 * (W(j + 1, k) - W(j - 1, k));
                    }
                }
            }
        }
    }

    void LHSX(int j) {
        double phi, theta, gamma1, gamma2, beta, kx, ky;
        int max_ij;
        int n, m;
        for (int i = 0; i < imax; i++) {
            phi = 0.5 * (gamma - 1) * (pow(u(i, j), 2) + pow(v(i, j), 2));
            kx = XiX(i, j);
            ky = XiY(i, j);
            theta = kx * u(i, j) + ky * v(i, j);
            gamma1 = gamma - 1;
            gamma2 = gamma - 2;
            beta = gamma * Q(i, j, 3) / Q(i, j, 0) - pow(phi, 2);
            max_ij = std::max(imax, jmax);

            B(i, 0, 0) = 0;
            B(i, 0, 1) = kx;
            B(i, 0, 2) = ky;
            B(i, 0, 3) = 0;
            B(i, 1, 0) = kx * pow(phi, 2) - u(i, j) * theta;
            B(i, 1, 1) = theta - kx * gamma2 * u(i, j);
            B(i, 1, 2) = ky * u(i, j) - gamma1 * kx * v(i, j);
            B(i, 1, 3) = kx * gamma1;
            B(i, 2, 0) = ky * pow(phi, 2) - v(i, j) * theta;
            B(i, 2, 1) = kx * v(i, j) - ky * gamma1 * u(i, j);
            B(i, 2, 2) = theta - ky * gamma2 * v(i, j);
            B(i, 2, 3) = ky * gamma1;
            B(i, 3, 0) = theta * (pow(phi, 2) - beta);
            B(i, 3, 1) = kx * beta - gamma1 * u(i, j) * theta;
            B(i, 3, 2) = ky * beta - gamma1 * v(i, j) * theta;
            B(i, 3, 3) = gamma * theta;

            for (n = 0; n < 4; n++) {
                for (m = 0; m < 4; m++) {
                    for (int i = 1; i < imax - 1; i++) {
                        A(i, m, n) = -0.5 * B(i - 1, m, n);
                        C(i, m, n) = 0.5 * B(i + 1, m, n);
                    }
                }
            }
            for (i = 1; i < imax - 1; i++) {
                B(i, 0, 0) = 1.0;
                B(i, 0, 1) = 0;
                B(i, 0, 2) = 0;
                B(i, 0, 3) = 0;
                B(i, 1, 0) = 0;
                B(i, 1, 1) = 1.0;
                B(i, 1, 2) = 0;
                B(i, 1, 3) = 0;
                B(i, 2, 0) = 0;
                B(i, 2, 1) = 0;
                B(i, 2, 2) = 1.0;
                B(i, 2, 3) = 0;
                B(i, 3, 0) = 0;
                B(i, 3, 1) = 0;
                B(i, 3, 2) = 0;
                B(i, 3, 3) = 1.0;
            }
        }
    }

    void LHSY(int i) {
        double phi, theta, gamma1, gamma2, beta, kx, ky;
        int max_ij;
        int n, m;
        for (int j = 0; i < jmax; j++) {
            phi = 0.5 * (gamma - 1) * (pow(u(i, j), 2) + pow(v(i, j), 2));
            kx = EtaX(i, j);
            ky = EtaY(i, j);
            theta = kx * u(i, j) + ky * v(i, j);
            gamma1 = gamma - 1;
            gamma2 = gamma - 2;
            beta = gamma * Q(i, j, 3) / Q(i, j, 0) - pow(phi, 2);
            max_ij = std::max(imax, jmax);

            B(i, 0, 0) = 0;
            B(i, 0, 1) = kx;
            B(i, 0, 2) = ky;
            B(i, 0, 3) = 0;
            B(i, 1, 0) = kx * pow(phi, 2) - u(i, j) * theta;
            B(i, 1, 1) = theta - kx * gamma2 * u(i, j);
            B(i, 1, 2) = ky * u(i, j) - gamma1 * kx * v(i, j);
            B(i, 1, 3) = kx * gamma1;
            B(i, 2, 0) = ky * pow(phi, 2) - v(i, j) * theta;
            B(i, 2, 1) = kx * v(i, j) - ky * gamma1 * u(i, j);
            B(i, 2, 2) = theta - ky * gamma2 * v(i, j);
            B(i, 2, 3) = ky * gamma1;
            B(i, 3, 0) = theta * (pow(phi, 2) - beta);
            B(i, 3, 1) = kx * beta - gamma1 * u(i, j) * theta;
            B(i, 3, 2) = ky * beta - gamma1 * v(i, j) * theta;
            B(i, 3, 3) = gamma * theta;

            for (n = 0; n < 4; n++) {
                for (m = 0; m < 4; m++) {
                    for (int j = 1; j < jmax - 1; j++) {
                        A(j, m, n) = -0.5 * B(j - 1, m, n);
                        C(j, m, n) = 0.5 * B(j + 1, m, n);
                    }
                }
            }
            for (j = 1; j < jmax - 1; j++) {
                B(j, 0, 0) = 1.0;
                B(j, 0, 1) = 0;
                B(j, 0, 2) = 0;
                B(j, 0, 3) = 0;
                B(j, 1, 0) = 0;
                B(j, 1, 1) = 1.0;
                B(j, 1, 2) = 0;
                B(j, 1, 3) = 0;
                B(j, 2, 0) = 0;
                B(j, 2, 1) = 0;
                B(j, 2, 2) = 1.0;
                B(j, 2, 3) = 0;
                B(j, 3, 0) = 0;
                B(j, 3, 1) = 0;
                B(j, 3, 2) = 0;
                B(j, 3, 3) = 1.0;
            }
        }
    }

    void advance() {
        // updates the solution
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                for (int k = 0; k < 4; k++) {
                    Q(i, j, k) += S(i, j, k) * Jacobian(i, j);
                }
            }
        }
    }

    void step() {
        ZeroRHS();
        RHS();
        smooth(Q, S, Jacobian, XiX, XiY,
            EtaX, EtaY, imax, jmax, *s2,
            *rspec, *qv, *dd, 
            espec, gamma, FSMACH, deltaT);
        
        // Xi 
        for (int j = 0; j < jmax; j++) {
            LHSX(j);
            smoothx();
            for (int k = 0; k < 4; k++) {
                for (int i = 0; i < imax; i++) {
                    D(i, k) = S(i, j, k);
                }
            }
            btri4s(A, B, C, D, imax + 1, 0, imax - 1);
            for (int k = 0; k < 4; k++) {
                for (int i = 0; i < imax; i++) {
                    S(i, j, k) = D(i, k);
                }
            }
        }
        // Eta 
        for (int i = 0; i < imax; i++) {
            LHSY(i);
            smoothy();
            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < jmax; j++) {
                    D(j, k) = S(i, j, k);
                }
            }
            btri4s(A, B, C, D, jmax + 1, 0, jmax - 1);
            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < jmax; j++) {
                    S(i, j, k) = D(j, k);
                }
            }
        }
        advance();
    }

    void solve() {

    }
 };



// main
    int main(int argc, const char* argv[]) {
        if (argc > 1) {
            for (int i = 1; i < argc; i++) {
                Mesh mesh(argv[i]);
                mesh.solve();
            }
        }
        else {
            Mesh mesh;
            mesh.solve();
        }
        return 0;
    };