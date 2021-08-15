#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <cassert>
#include <string>
#include <sstream>
#define M_PI 3.14159265358979323846 /* pi */

#include<math.h>
int smoothx(float* q, float* xx, float* xy, int id, int jd, float* a,
    float* b, float* c, int j, float* jac, float* drr, float* drp,
    float* rspec, float* qv, float* dd,
    float epsi, float gamma, float fsmach, float dt)
{

    float* rho, * u_vel, * v_vel, * t_e;

    float eratio, gm1, ggm1, eps, ra, u, v, qq, ss,
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

int smoothy(float* q, float* yx, float* yy, int id, int jd, float* a,
    float* b, float* c, int i, float* jac, float* drr, float* drp,
    float* rspec, float* qv, float* dd,
    float epsi, float gamma, float fsmach, float dt)
{

    float* rho, * u_vel, * v_vel, * t_e;

    float eratio, smool, gm1, ggm1, eps, ra, u, v, qq, ss,
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

int smooth(float* q, float* s, float* jac, float* xx, float* xy,
    float* yx, float* yy, int id, int jd, float* s2,
    float* rspec, float* qv, float* dd,
    float epse, float gamma, float fsmach, float dt)
{

    float* rho, * u_vel, * v_vel, * t_e;

    float eratio, smool, gm1, ggm1, cx, cy, eps, ra, u, v, qq, ss, st,
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

    Mesh(const char* paramFilename = nullptr) : FSMACH(0.9), gamma(1.4), EPSE(0.06), deltaT(1), alpha_inf(0), imax(121), jmax(41), 
        xFilename("x.csv"), yFilename("y.csv") {
        if (paramFilename) {
            parseParams(paramFilename);
        }
        xArr = parseCSV(xFilename.c_str());
        yArr = parseCSV(yFilename.c_str());
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

    void FreeStream(int i, int j) {
        double a;
        double e;
        a = sqrt(gamma * p_inf) / rho_inf;
        e = p_inf / (gamma - 1) + 0.5 * rho_inf * pow(FSMACH * a, 2);
        for (i = 0; i < imax; i++) {
            for (j = 0; j < jmax; j++) {
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


    void BoundaryConds(int i, int j, int k, int iTEL, int iTEU, int imax, int jmax) {
        double rho01, rhoTEL, rhoTEU, uTEU, uTEL, vTEU, vTEL, p0, pTEL, pTEU, Jacobian0, U;
            // Adiabetic wall: j = jmin, i = itel...iteu
            for (i = iTEL; i <= iTEU; i++) {
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

        for (k = 0; k < 4; k++) {
            Q(iTEU, 0, k) = Q(iTEL, 0, k);
        }
        // Wake
        for (i = 0; i < iTEL; i++) {
            for (k = 0; k < 4; k++) {
                Q(i, 0, k) = 0.5 * (Q(i, 1, k) + Q(imax - 1, 1, k));
                Q(imax - -1 - i, 0, k) = Q(i, 0, k);
            }
        }
        // outflow
        for (j = 0; j < jmax; j++) {
            for (k = 0; k < 4; k++) {
                Q(0, j, k) = Q(1, j, k);
                Q(imax - 1, j, k) = Q(imax - 2, j, k);
            }
        }

    }

    void ZeroRHS(int i, int j, int k) {
        for (i = 0; i < imax; i++) {
            for (j = 0; j < jmax; j++) {
                for (k = 0; k < 4; k++) {
                    S(i, j, k) = 0;
                }
            }
        }
    }
    void RHS(int i, int j, int k) {
        //Xi direction
        double UU, VV, p;
        for (j = 1; j < jmax - 1; j++) {
            for (i = 0; i < imax; i++) {
                UU = u(i, j) * XiX(i, j) + v(i, j) * XiY(i, j);
                p = (gamma - 1) * (Q(i, j, 3) - 0.5 * Q(i, j, 0) * (pow(Q(i, j, 1) / Q(i, j, 0), 2) + pow(Q(i, j, 2) / Q(i, j, 0), 2)));
                W(i, 0) = Q(i, j, 0) * UU / Jacobian(i, j);
                W(i, 1) = (Q(i, j, 1) * UU + XiX(i, j) * p) / Jacobian(i, j);
                W(i, 2) = (Q(i, j, 2) * UU + XiY(i, j) * p) / Jacobian(i, j);
                W(i, 3) = (Q(i, j, 3) + p) * UU / Jacobian(i, j);
            }
            for (i = 1; i < imax - 1; i++) {
                for (k = 0; k < 4; k++) {
                    S(i, j, k) += -0.5*(W(i + 1, k) - W(i - 1, k));
                }
            }
        }
        // Eta direction
        for (i = 1; i < imax - 1; i++) {
            for (j = 0; j < jmax; j++) {
                VV = 
                W(j,0) = 
            }
        }
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