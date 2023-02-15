#include "interpolate.h"

void interpolate(const real &t, const real &dt, const std::vector<real> &xInteg0, const std::vector<real> &accInteg0, const std::vector< std::vector<real> > &b, const std::vector<real> &hVec, Simulation &sim){
    std::vector<real> tVecForInterp (hVec.size(), 0.0);
    std::vector< std::vector<real> > xIntegForInterp(hVec.size(), std::vector<real>(xInteg0.size(), 0.0));
    static std::vector<real> tVecForInterpPrev(tVecForInterp.size(), 0.0);
    static std::vector< std::vector<real> > xIntegForInterpPrev(xIntegForInterp.size(), std::vector<real>(xIntegForInterp[0].size(), 0.0));
    
    tVecForInterp[0] = t;
    xIntegForInterp[0] = xInteg0;
    for (size_t hIdx = 1; hIdx < hVec.size(); hIdx++) {
        tVecForInterp[hIdx] = t + hVec[hIdx]*dt;
        approx_xInteg(xInteg0, accInteg0, xIntegForInterp[hIdx], dt, hVec[hIdx], b, sim.integParams.nInteg);
    }
    // tVecForInterp.push_back(t+dt);
    // xIntegForInterp.push_back(xInteg);
    
    one_timestep_interpolation(t+dt, tVecForInterp, xIntegForInterp, tVecForInterpPrev, xIntegForInterpPrev, sim);
    tVecForInterpPrev = tVecForInterp;
    xIntegForInterpPrev = xIntegForInterp;
}

void get_coeffs(const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &xIntegForInterp, std::vector< std::vector<real> > &coeffs){
    size_t tLen = tVecForInterp.size();
    size_t numStates = xIntegForInterp[0].size();
    std::vector< std::vector< std::vector<real> > > c(numStates, std::vector< std::vector<real> >(tLen, std::vector<real>(tLen, 0.0)));
    for (size_t i = 0; i < numStates; i++){
        for (size_t j = 0; j < tLen; j++){
            c[i][j][0] = xIntegForInterp[j][i];
        }
        for (size_t j = 1; j < tLen; j++){
            for (size_t k = 0; k < tLen-j; k++){
                c[i][k][j] = (c[i][k+1][j-1] - c[i][k][j-1])/(tVecForInterp[k+j]-tVecForInterp[k]);
            }
        }
    }

    for (size_t i = 0; i < numStates; i++){
        for (size_t j = 0; j < tLen; j++){
            coeffs[i][j] = c[i][0][j];
        }
    }
}

void evaluate_one_interpolation(const real &tInterp, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &coeffs, std::vector<real> &xInterp){
    size_t tLen = tVecForInterp.size();
    size_t numStates = xInterp.size();
    size_t n = tLen-1;
    for (size_t i = 0; i < numStates; i++){
        xInterp[i] = coeffs[i][n];
    }
    for (size_t i = 1; i < n+1; i++){
        for (size_t j = 0; j < numStates; j++){
            xInterp[j] = coeffs[j][n-i] + (tInterp-tVecForInterp[n-i])*xInterp[j];
        }
    }
}

void one_timestep_interpolation(const real &tNext, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &xIntegForInterp, const std::vector<real> &tVecForInterpPrev, const std::vector< std::vector<real> > &xIntegForInterpPrev, Simulation &sim){
    size_t tLen = tVecForInterp.size();
    size_t numStates = xIntegForInterp[0].size();
    std::vector< std::vector<real> > coeffs(numStates, std::vector<real>(tLen, 0.0));
    std::vector< std::vector<real> > coeffsPrev(numStates, std::vector<real>(tLen, 0.0));

    get_coeffs(tVecForInterp, xIntegForInterp, coeffs);
    get_coeffs(tVecForInterpPrev, xIntegForInterpPrev, coeffsPrev);

    static size_t interpIdx = 0;
    // std::cout << "interpIdx = " << interpIdx << std::endl;
    bool forwardIntegrate = tVecForInterp[0] < tVecForInterp[tLen-1];
    bool backwardIntegrate = tVecForInterp[0] > tVecForInterp[tLen-1];
    while ( interpIdx < sim.tEval.size()
            &&( (forwardIntegrate && (sim.tEval[interpIdx] == tVecForInterp[0] || (sim.tEval[interpIdx] > tVecForInterp[0] && sim.tEval[interpIdx] <= tNext)))
            ||  (forwardIntegrate && sim.tEval[interpIdx] <= sim.integParams.t0 && sim.tEval[interpIdx]+sim.tEvalMargin >= sim.integParams.t0) || (forwardIntegrate && sim.tEval[interpIdx] >= sim.integParams.tf && sim.tEval[interpIdx]-sim.tEvalMargin <= sim.integParams.tf)
            ||  (backwardIntegrate && (sim.tEval[interpIdx] == tVecForInterp[0] || (sim.tEval[interpIdx] < tVecForInterp[0] && sim.tEval[interpIdx] >= tNext)))
            ||  (backwardIntegrate && sim.tEval[interpIdx] >= sim.integParams.t0 && sim.tEval[interpIdx]-sim.tEvalMargin <= sim.integParams.t0) || (backwardIntegrate && sim.tEval[interpIdx] <= sim.integParams.tf && sim.tEval[interpIdx]+sim.tEvalMargin >= sim.integParams.tf) )
        ){
        real tInterpGeom;
        if (!sim.evalApparentState){
            for (size_t i = 0; i < tLen; i++){
                tInterpGeom = sim.tEval[interpIdx];
                if (tInterpGeom == tVecForInterp[i]){
                    sim.xIntegEval.push_back(xIntegForInterp[i]);
                    interpIdx++;
                    // std::cout << "exactly interpolated tInterpGeom = " << tInterpGeom << std::endl;
                    continue;
                }
            }
        }
        tInterpGeom = sim.tEval[interpIdx];
        // std::cout << "tInterpGeom = " << tInterpGeom << std::endl;
        std::vector<real> xInterpGeom(numStates, 0.0);
        evaluate_one_interpolation(tInterpGeom, tVecForInterp, coeffs, xInterpGeom);
        if (sim.evalApparentState){
            real lightTimeTemp;
            std::vector<real> xRelativeTemp(6, 0.0);
            real distRelativeTemp;
            std::vector<real> xInterpApparentTemp(numStates, 0.0);
            std::vector<real> lightTime(sim.integParams.nInteg, 0.0);
            std::vector<real> xInterpApparent(numStates, 0.0);
            for (size_t i = 0; i < sim.integParams.nInteg; i++){
                for (size_t j = 0; j < 6; j++){
                    xRelativeTemp[j] = xInterpGeom[6*i+j] - sim.xObserver[interpIdx][j];
                }
                vnorm({xRelativeTemp[0], xRelativeTemp[1], xRelativeTemp[2]}, distRelativeTemp);
                lightTimeTemp = distRelativeTemp/sim.consts.clight;
                if (sim.convergedLightTime){
                    real lightTimeTol = 1e-16/86400.0L;
                    real lightTimeTempPrev = 0.0L;
                    size_t maxIter = 20;
                    size_t iter = 0;
                    // keep iterating until max iterations or light time tolerance is met
                    while (iter < maxIter && fabs(lightTimeTemp - lightTimeTempPrev) > lightTimeTol){
                        evaluate_one_interpolation(tInterpGeom-lightTimeTemp, tVecForInterp, coeffs, xInterpApparentTemp);
                        for (size_t j = 0; j < 6; j++){
                            xRelativeTemp[j] = xInterpApparentTemp[6*i+j] - sim.xObserver[interpIdx][j];
                        }
                        vnorm({xRelativeTemp[0], xRelativeTemp[1], xRelativeTemp[2]}, distRelativeTemp);
                        lightTimeTempPrev = lightTimeTemp;
                        lightTimeTemp = distRelativeTemp/sim.consts.clight;
                        iter++;
                    }
                }
                evaluate_one_interpolation(tInterpGeom-lightTimeTemp, tVecForInterp, coeffs, xInterpApparentTemp);
                lightTime[i] = lightTimeTemp;
                for (size_t j = 0; j < 6; j++){
                    xInterpApparent[6*i+j] = xInterpApparentTemp[6*i+j] - sim.xObserver[interpIdx][j];
                }
            }
            sim.xIntegEval.push_back(xInterpApparent);
        } else {
            sim.xIntegEval.push_back(xInterpGeom);
        }
        interpIdx++;
    }
}
