#include "gr15.h"

real get_initial_timestep(const real &t, const std::vector<real> &xInteg0, const ForceParameters &forceParams, IntegrationParameters &integParams, const Constants &consts){
    
    real dt;
    if (integParams.dt0!=0.0){
        dt = abs(integParams.dt0);
        if (integParams.tf < integParams.t0){
            dt *= -1.0;
            integParams.dtMax = -abs(integParams.dtMax);
            integParams.dtMin = -abs(integParams.dtMin);
        }
        return dt;
    }
    int order = 15;
    real absMaxPos0, absMaxAcc0, absMaxAcc1Minus0;
    real dtTemp0, dtTemp1;
    std::vector<real> posInteg0(3*integParams.nInteg, 0.0);
    std::vector<real> accInteg1Minus0(3*integParams.nInteg, 0.0);
    std::vector<real> xIntegNext(6*integParams.nInteg, 0.0);

    std::vector<real> accInteg0 = get_state_der(t, xInteg0, forceParams, integParams, consts);
    for (size_t i=0; i<integParams.nInteg; i++){
        for (size_t j=0; j<3; j++){
            posInteg0[3*i+j] = xInteg0[6*i+j];
        }
    }
    vabs_max(posInteg0, absMaxPos0);
    vabs_max(accInteg0, absMaxAcc0);
    if (absMaxPos0 < 1.0e-5 || absMaxAcc0 < 1e-5){
        dtTemp0 = 1.0e-6;
    } else {
        dtTemp0 = 0.01 * (absMaxPos0 / absMaxAcc0);
    }
    // propagate xInteg0 to xIntegNext using an Euler step and a timestep of dtTemp0
    for (size_t i=0; i<integParams.nInteg; i++){
        for (size_t j=0; j<3; j++){
            xIntegNext[6*i+j] = xInteg0[6*i+j] + dtTemp0 * xInteg0[6*i+j+3];
            xIntegNext[6*i+j+3] = xInteg0[6*i+j+3] + dtTemp0 * accInteg0[3*i+j];
        }
    }
    std::vector<real> accInteg1 = get_state_der(t+dtTemp0, xIntegNext, forceParams, integParams, consts);
    vsub(accInteg1, accInteg0, accInteg1Minus0);
    vabs_max(accInteg1Minus0, absMaxAcc1Minus0);
    if (fmax(absMaxAcc0, absMaxAcc1Minus0) <= 1e-15) {
        dtTemp1 = fmax(1.0e-6, dtTemp0 * 1e-3);
    } else {
        dtTemp1 = pow(0.01 / fmax(absMaxAcc0, absMaxAcc1Minus0), 1.0/(order+1));
    }
    dt = fmin(100*dtTemp0, dtTemp1);
    if (abs(integParams.tf-integParams.t0) < dt){
        dt = abs(integParams.tf - integParams.t0);
    }
    if (integParams.tf < integParams.t0){
        dt *= -1.0;
        integParams.dtMax = -abs(integParams.dtMax);
        integParams.dtMin = -abs(integParams.dtMin);
    }
    return dt;
}

void approx_xInteg(const std::vector<real> &xInteg0, const std::vector<real> &accInteg0, std::vector<real> &xIntegNext, const real &dt, const real &h, const std::vector< std::vector<real> > &b, const size_t &nInteg){
    for (size_t i=0; i<nInteg; i++){
        xIntegNext[6*i]   = xInteg0[6*i]   + dt * h * (xInteg0[6*i+3] + dt * h * (accInteg0[3*i]   + h * (b[0][3*i]   / 0.3e1 + h * (b[1][3*i]   / 0.6e1 + h * (b[2][3*i]   / 0.10e2 + h * (b[3][3*i]   / 0.15e2 + h * (b[4][3*i]   / 0.21e2 + h * (b[5][3*i]   / 0.28e2 + h * b[6][3*i]   / 0.36e2))))))) / 0.2e1);
        xIntegNext[6*i+1] = xInteg0[6*i+1] + dt * h * (xInteg0[6*i+4] + dt * h * (accInteg0[3*i+1] + h * (b[0][3*i+1] / 0.3e1 + h * (b[1][3*i+1] / 0.6e1 + h * (b[2][3*i+1] / 0.10e2 + h * (b[3][3*i+1] / 0.15e2 + h * (b[4][3*i+1] / 0.21e2 + h * (b[5][3*i+1] / 0.28e2 + h * b[6][3*i+1] / 0.36e2))))))) / 0.2e1);
        xIntegNext[6*i+2] = xInteg0[6*i+2] + dt * h * (xInteg0[6*i+5] + dt * h * (accInteg0[3*i+2] + h * (b[0][3*i+2] / 0.3e1 + h * (b[1][3*i+2] / 0.6e1 + h * (b[2][3*i+2] / 0.10e2 + h * (b[3][3*i+2] / 0.15e2 + h * (b[4][3*i+2] / 0.21e2 + h * (b[5][3*i+2] / 0.28e2 + h * b[6][3*i+2] / 0.36e2))))))) / 0.2e1);

        xIntegNext[6*i+3] = xInteg0[6*i+3] + dt * h * (accInteg0[3*i]   + h * (b[0][3*i]   / 0.2e1 + h * (b[1][3*i]   / 0.3e1 + h * (b[2][3*i]   / 0.4e1 + h * (b[3][3*i]   / 0.5e1 + h * (b[4][3*i]   / 0.6e1 + h * (b[5][3*i]   / 0.7e1 + h * b[6][3*i]   / 0.8e1)))))));
        xIntegNext[6*i+4] = xInteg0[6*i+4] + dt * h * (accInteg0[3*i+1] + h * (b[0][3*i+1] / 0.2e1 + h * (b[1][3*i+1] / 0.3e1 + h * (b[2][3*i+1] / 0.4e1 + h * (b[3][3*i+1] / 0.5e1 + h * (b[4][3*i+1] / 0.6e1 + h * (b[5][3*i+1] / 0.7e1 + h * b[6][3*i+1] / 0.8e1)))))));
        xIntegNext[6*i+5] = xInteg0[6*i+5] + dt * h * (accInteg0[3*i+2] + h * (b[0][3*i+2] / 0.2e1 + h * (b[1][3*i+2] / 0.3e1 + h * (b[2][3*i+2] / 0.4e1 + h * (b[3][3*i+2] / 0.5e1 + h * (b[4][3*i+2] / 0.6e1 + h * (b[5][3*i+2] / 0.7e1 + h * b[6][3*i+2] / 0.8e1)))))));
    }
}

void compute_g_and_b(const std::vector< std::vector<real> > &AccIntegArr, const size_t &hIdx, std::vector< std::vector<real> > &g, std::vector< std::vector<real> > &b, const size_t &dim){
    const std::vector<real> Acc1 = AccIntegArr[0];
    const std::vector<real> Acc2 = AccIntegArr[1];
    const std::vector<real> Acc3 = AccIntegArr[2];
    const std::vector<real> Acc4 = AccIntegArr[3];
    const std::vector<real> Acc5 = AccIntegArr[4];
    const std::vector<real> Acc6 = AccIntegArr[5];
    const std::vector<real> Acc7 = AccIntegArr[6];
    const std::vector<real> Acc8 = AccIntegArr[7];

    for (size_t i=0; i<dim; i++){
        if (hIdx == 1) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
        } else if (hIdx == 2) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
        } else if (hIdx == 3) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
            g[2][i] = (((Acc4[i] - Acc1[i]) * rMat[3][0] - g[0][i]) * rMat[3][1] - g[1][i]) * rMat[3][2];
        } else if (hIdx == 4) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
            g[2][i] = (((Acc4[i] - Acc1[i]) * rMat[3][0] - g[0][i]) * rMat[3][1] - g[1][i]) * rMat[3][2];
            g[3][i] = ((((Acc5[i] - Acc1[i]) * rMat[4][0] - g[0][i]) * rMat[4][1] - g[1][i]) * rMat[4][2] - g[2][i]) * rMat[4][3];
        } else if (hIdx == 5) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
            g[2][i] = (((Acc4[i] - Acc1[i]) * rMat[3][0] - g[0][i]) * rMat[3][1] - g[1][i]) * rMat[3][2];
            g[3][i] = ((((Acc5[i] - Acc1[i]) * rMat[4][0] - g[0][i]) * rMat[4][1] - g[1][i]) * rMat[4][2] - g[2][i]) * rMat[4][3];
            g[4][i] = (((((Acc6[i] - Acc1[i]) * rMat[5][0] - g[0][i]) * rMat[5][1] - g[1][i]) * rMat[5][2] - g[2][i]) * rMat[5][3] - g[3][i]) * rMat[5][4];
        } else if (hIdx == 6) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
            g[2][i] = (((Acc4[i] - Acc1[i]) * rMat[3][0] - g[0][i]) * rMat[3][1] - g[1][i]) * rMat[3][2];
            g[3][i] = ((((Acc5[i] - Acc1[i]) * rMat[4][0] - g[0][i]) * rMat[4][1] - g[1][i]) * rMat[4][2] - g[2][i]) * rMat[4][3];
            g[4][i] = (((((Acc6[i] - Acc1[i]) * rMat[5][0] - g[0][i]) * rMat[5][1] - g[1][i]) * rMat[5][2] - g[2][i]) * rMat[5][3] - g[3][i]) * rMat[5][4];
            g[5][i] = ((((((Acc7[i] - Acc1[i]) * rMat[6][0] - g[0][i]) * rMat[6][1] - g[1][i]) * rMat[6][2] - g[2][i]) * rMat[6][3] - g[3][i]) * rMat[6][4] - g[4][i]) * rMat[6][5];
        } else if (hIdx == 7) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
            g[2][i] = (((Acc4[i] - Acc1[i]) * rMat[3][0] - g[0][i]) * rMat[3][1] - g[1][i]) * rMat[3][2];
            g[3][i] = ((((Acc5[i] - Acc1[i]) * rMat[4][0] - g[0][i]) * rMat[4][1] - g[1][i]) * rMat[4][2] - g[2][i]) * rMat[4][3];
            g[4][i] = (((((Acc6[i] - Acc1[i]) * rMat[5][0] - g[0][i]) * rMat[5][1] - g[1][i]) * rMat[5][2] - g[2][i]) * rMat[5][3] - g[3][i]) * rMat[5][4];
            g[5][i] = ((((((Acc7[i] - Acc1[i]) * rMat[6][0] - g[0][i]) * rMat[6][1] - g[1][i]) * rMat[6][2] - g[2][i]) * rMat[6][3] - g[3][i]) * rMat[6][4] - g[4][i]) * rMat[6][5];
            g[6][i] = (((((((Acc8[i] - Acc1[i]) * rMat[7][0] - g[0][i]) * rMat[7][1] - g[1][i]) * rMat[7][2] - g[2][i]) * rMat[7][3] - g[3][i]) * rMat[7][4] - g[4][i]) * rMat[7][5] - g[5][i]) * rMat[7][6];
        }
    }
    for (size_t i=0; i<dim; i++){
        if (hIdx == 1) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
        } else if (hIdx == 2) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
        } else if (hIdx == 3) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
            b[2][i] =                                         + cMat[2][2]*g[2][i] + cMat[3][2]*g[3][i] + cMat[4][2]*g[4][i] + cMat[5][2]*g[5][i] + cMat[6][2]*g[6][i];
        } else if (hIdx == 4) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
            b[2][i] =                                         + cMat[2][2]*g[2][i] + cMat[3][2]*g[3][i] + cMat[4][2]*g[4][i] + cMat[5][2]*g[5][i] + cMat[6][2]*g[6][i];
            b[3][i] =                                                              + cMat[3][3]*g[3][i] + cMat[4][3]*g[4][i] + cMat[5][3]*g[5][i] + cMat[6][3]*g[6][i];
        } else if (hIdx == 5) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
            b[2][i] =                                         + cMat[2][2]*g[2][i] + cMat[3][2]*g[3][i] + cMat[4][2]*g[4][i] + cMat[5][2]*g[5][i] + cMat[6][2]*g[6][i];
            b[3][i] =                                                              + cMat[3][3]*g[3][i] + cMat[4][3]*g[4][i] + cMat[5][3]*g[5][i] + cMat[6][3]*g[6][i];
            b[4][i] =                                                                                   + cMat[4][4]*g[4][i] + cMat[5][4]*g[5][i] + cMat[6][4]*g[6][i];
        } else if (hIdx == 6) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
            b[2][i] =                                         + cMat[2][2]*g[2][i] + cMat[3][2]*g[3][i] + cMat[4][2]*g[4][i] + cMat[5][2]*g[5][i] + cMat[6][2]*g[6][i];
            b[3][i] =                                                              + cMat[3][3]*g[3][i] + cMat[4][3]*g[4][i] + cMat[5][3]*g[5][i] + cMat[6][3]*g[6][i];
            b[4][i] =                                                                                   + cMat[4][4]*g[4][i] + cMat[5][4]*g[5][i] + cMat[6][4]*g[6][i];
            b[5][i] =                                                                                                        + cMat[5][5]*g[5][i] + cMat[6][5]*g[6][i];
        } else if (hIdx == 7) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
            b[2][i] =                                         + cMat[2][2]*g[2][i] + cMat[3][2]*g[3][i] + cMat[4][2]*g[4][i] + cMat[5][2]*g[5][i] + cMat[6][2]*g[6][i];
            b[3][i] =                                                              + cMat[3][3]*g[3][i] + cMat[4][3]*g[4][i] + cMat[5][3]*g[5][i] + cMat[6][3]*g[6][i];
            b[4][i] =                                                                                   + cMat[4][4]*g[4][i] + cMat[5][4]*g[5][i] + cMat[6][4]*g[6][i];
            b[5][i] =                                                                                                        + cMat[5][5]*g[5][i] + cMat[6][5]*g[6][i];
            b[6][i] =                                                                                                                             + cMat[6][6]*g[6][i];
        }
    }
}

void refine_b(std::vector< std::vector<real> > &b, std::vector< std::vector<real> > &e, const real &dtRatio, const size_t &dim, const size_t &timestepCounter){
    std::vector< std::vector<real> > bDiff(7, std::vector<real>(dim, 0.0L));
    if (timestepCounter > 1){
        for (size_t i = 0; i < dim; i++){
            bDiff[0][i] = b[0][i] - e[0][i];
            bDiff[1][i] = b[1][i] - e[1][i];
            bDiff[2][i] = b[2][i] - e[2][i];
            bDiff[3][i] = b[3][i] - e[3][i];
            bDiff[4][i] = b[4][i] - e[4][i];
            bDiff[5][i] = b[5][i] - e[5][i];
            bDiff[6][i] = b[6][i] - e[6][i];
        }
    }

    real q = dtRatio;
    real q2 = q * q;
    real q3 = q2 * q;
    real q4 = q2 * q2;
    real q5 = q2 * q3;
    real q6 = q3 * q3;
    real q7 = q2 * q5;

    for (size_t i = 0; i < dim; i++) {
        e[0][i] = q  * (b[6][i] * 7.0  + b[5][i] * 6.0  + b[4][i] * 5.0  + b[3][i] * 4.0 + b[2][i] * 3.0 + b[1][i] * 2.0 + b[0][i]);
        e[1][i] = q2 * (b[6][i] * 21.0 + b[5][i] * 15.0 + b[4][i] * 10.0 + b[3][i] * 6.0 + b[2][i] * 3.0 + b[1][i]);
        e[2][i] = q3 * (b[6][i] * 35.0 + b[5][i] * 20.0 + b[4][i] * 10.0 + b[3][i] * 4.0 + b[2][i]);
        e[3][i] = q4 * (b[6][i] * 35.0 + b[5][i] * 15.0 + b[4][i] * 5.0  + b[3][i]);
        e[4][i] = q5 * (b[6][i] * 21.0 + b[5][i] * 6.0  + b[4][i]);
        e[5][i] = q6 * (b[6][i] * 7.0  + b[5][i]);
        e[6][i] = q7 * (b[6][i]);
    }

    for (size_t i = 0; i < dim; i++) {
        b[0][i] = e[0][i] + bDiff[0][i];
        b[1][i] = e[1][i] + bDiff[1][i];
        b[2][i] = e[2][i] + bDiff[2][i];
        b[3][i] = e[3][i] + bDiff[3][i];
        b[4][i] = e[4][i] + bDiff[4][i];
        b[5][i] = e[5][i] + bDiff[5][i];
        b[6][i] = e[6][i] + bDiff[6][i];
    }
}

void gr15(real t, std::vector<real> xInteg0, Simulation &sim){
    ForceParameters &forceParams = sim.forceParams;
    IntegrationParameters &integParams = sim.integParams;
    size_t dim = 3*integParams.nInteg;
    Constants &consts = sim.consts;

    real dt = get_initial_timestep(t, xInteg0, forceParams, integParams, consts);
    std::vector<real> accInteg0 = get_state_der(t, xInteg0, forceParams, integParams, consts);
    if (dt > 0){
        std::sort(sim.tEval.begin(), sim.tEval.end()); // sort sim.tEval into ascending order
    }
    else if (dt < 0){
        std::sort(sim.tEval.begin(), sim.tEval.end(), std::greater<real>()); // sort sim.tEval into descending order
    }
    std::vector<real> xInteg(2*dim, 0.0);
    std::vector< std::vector<real> > b_old(7, std::vector<real>(dim, 0.0));
    std::vector< std::vector<real> > b(7, std::vector<real>(dim, 0.0));
    std::vector< std::vector<real> > g(7, std::vector<real>(dim, 0.0));
    std::vector< std::vector<real> > e(7, std::vector<real>(dim, 0.0));
    std::vector< std::vector<real> > accIntegArr(hVec.size(), std::vector<real>(dim, 0.0));
    std::vector<real> b6Tilde(dim, 0.0);
    real b6TildeMax, accIntegArr7Max;
    real b6TildeEstim, b6Max, accIntegNextMax;
    real relError, dtReq;
    
    size_t PCmaxIter = 12;
    int maxLoops = 100;
    int loopCounter = 0;
    int keepStepping = 1;
    int oneStepDone = 0;
    if (integParams.t0 == integParams.tf){
        keepStepping = 0;
    }
    while (keepStepping){
        while (!oneStepDone){
            for (size_t PCidx = 1; PCidx < PCmaxIter; PCidx++){
                xInteg = xInteg0;
                accIntegArr[0] = accInteg0;
                compute_g_and_b(accIntegArr, hVec[0], g, b, dim);
                for (size_t hIdx = 1; hIdx < hVec.size(); hIdx++) {
                    approx_xInteg(xInteg0, accInteg0, xInteg, dt, hVec[hIdx], b, integParams.nInteg);
                    accIntegArr[hIdx] = get_state_der(t + hVec[hIdx]*dt, xInteg, forceParams, integParams, consts);
                    compute_g_and_b(accIntegArr, hIdx, g, b, dim);
                }
                for (size_t i = 0; i < dim; i++){
                    b6Tilde[i] = b[6][i] - b_old[6][i];
                }
                vabs_max(b6Tilde, b6TildeMax);
                vabs_max(accIntegArr[7], accIntegArr7Max);
                if (b6TildeMax / accIntegArr7Max < integParams.tolPC){
                    break;
                }
                b_old = b;
            }
            approx_xInteg(xInteg0, accInteg0, xInteg, dt, 1.0, b, integParams.nInteg);
            std::vector<real> accIntegNext = get_state_der(t+dt, xInteg, forceParams, integParams, consts);

            vabs_max(b[6], b6Max);
            vabs_max(accIntegNext, accIntegNextMax);
            b6TildeEstim = b6Max / accIntegNextMax;
            if (integParams.adaptiveTimestep){
                relError = pow(b6TildeEstim/integParams.tolInteg, 1.0L/7.0L);
            }
            else{
                relError = pow(b6TildeEstim/integParams.tolInteg, 0.0L);
            }
            dtReq = dt/relError;
            
            if (relError <= 1 || loopCounter > maxLoops){
                oneStepDone = 1;
                integParams.timestepCounter += 1;
                // // start interpolation call
                std::vector< std::vector<real> > xIntegForInterp(hVec.size(), std::vector<real>(2*dim, 0.0));
                std::vector<real> tVecForInterp (hVec.size(), 0.0);
                xIntegForInterp[0] = xInteg0;
                tVecForInterp[0] = t;
                for (size_t hIdx = 1; hIdx < hVec.size(); hIdx++) {
                    approx_xInteg(xInteg0, accInteg0, xIntegForInterp[hIdx], dt, hVec[hIdx], b, integParams.nInteg);
                    tVecForInterp[hIdx] = t + hVec[hIdx]*dt;
                }
                // xIntegForInterp[hVec.size()] = xInteg;
                // tVecForInterp[hVec.size()] = t + dt;
                interpolate(t+dt, tVecForInterp, xIntegForInterp, sim);
                // // end interpolation call
                t += dt;
                if ((integParams.tf > integParams.t0 && t >= integParams.tf) || (integParams.tf < integParams.t0 && t <= integParams.tf)){
                    keepStepping = 0;
                }
                xInteg0 = xInteg;
                accInteg0 = accIntegNext;
                b_old = b;
                refine_b(b, e, dtReq/dt, dim, integParams.timestepCounter);
                loopCounter = 0;
            }
            else{
                loopCounter += 1;
            }
            if (dtReq/dt > 1.0/integParams.dtChangeFactor){
                dt /= integParams.dtChangeFactor;
            } else if (abs(dtReq) < 1.0e-12L){
                dt *= integParams.dtChangeFactor;
            } else {
                dt = dtReq;
            }
            if ( (integParams.tf > integParams.t0 && dt > integParams.dtMax) || (integParams.tf < integParams.t0 && dt < integParams.dtMax)){
                dt = integParams.dtMax;
            }
            if ( (integParams.tf > integParams.t0 && dt < integParams.dtMin) || (integParams.tf < integParams.t0 && dt > integParams.dtMin)){
                dt = integParams.dtMin;
            }
            if ( (integParams.tf > integParams.t0 && t+dt > integParams.tf) || (integParams.tf < integParams.t0 && t+dt < integParams.tf)){
                dt = integParams.tf-t;
            }
        }
    sim.t = t;
    sim.xInteg = xInteg;
    oneStepDone = 0;
    }
}
