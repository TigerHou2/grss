#include "simulation.h"

/**
 * @param[in] spiceId SPICE ID of the body.
 * @param[in] tMjdTDB Time of observation in Modified Julian Date (TDB).
 * @param[out] baseBodyFrame Name of the body-fixed frame.
 */
void get_baseBodyFrame(const int &spiceId, const real &tMjdTDB,
                       std::string &baseBodyFrame) {
    switch (spiceId) {
        case 10:
            baseBodyFrame = "IAU_SUN";
            break;
        case 1:
        case 199:
            baseBodyFrame = "IAU_MERCURY";
            break;
        case 2:
        case 299:
            baseBodyFrame = "IAU_VENUS";
            break;
        case 301:
            baseBodyFrame = "IAU_MOON";
            break;
        case 399:
            baseBodyFrame = "ITRF93";
            // High precision frame is not defined
            // before 1962 JAN 20 00:00:41.184 TDB or
            // after 2099 AUG 27 00:01:09.182 TDB
            if (tMjdTDB < 37684.0004767L || tMjdTDB > 87942.0008007L) {
                baseBodyFrame = "IAU_EARTH";
            }
            break;
        case 499:
            baseBodyFrame = "IAU_MARS";
            break;
        case 599:
            baseBodyFrame = "IAU_JUPITER";
            break;
        case 699:
            baseBodyFrame = "IAU_SATURN";
            break;
        case 799:
            baseBodyFrame = "IAU_URANUS";
            break;
        case 899:
            baseBodyFrame = "IAU_NEPTUNE";
            break;
        case 999:
            baseBodyFrame = "IAU_PLUTO";
            break;
        default:
            std::cout << "Given base body: " << spiceId << std::endl;
            throw std::invalid_argument("Given base body not supported");
            break;
    }
}

/** 
 * @param[in] tObsMjd Observation time in Modified Julian Date.
 * @param[in] observerInfo Observer information (base body SPICE ID, lon [rad], lat [rad], alt [m])
 * @param[in] propSim PropSimulation object for querying Ephemeris of base body.
 * @param[in] tObsInUTC Flag to indicate if the observation time is in UTC (true) or TDB (false).
 * @param[out] observerState Output observer state (position [DU], velocity [DU/TU]).
 */
void get_observer_state(const real &tObsMjd,
                        const std::vector<real> &observerInfo,
                        PropSimulation *propSim, const bool tObsInUTC,
                        std::vector<real> &observerState) {
    int baseBody = observerInfo[0];
    if ((int)observerInfo[0] == 500) baseBody = 399;
    if (baseBody == 0) {
        observerState[0] = 0.0L;
        observerState[1] = 0.0L;
        observerState[2] = 0.0L;
        observerState[3] = 0.0L;
        observerState[4] = 0.0L;
        observerState[5] = 0.0L;
        return;
    }
    real tObsMjdTDB = tObsMjd;
    if (tObsInUTC) {
        tObsMjdTDB += delta_et_utc(tObsMjd)/86400.0L;
    }
    double baseBodyState[9];
    get_spk_state(baseBody, tObsMjdTDB, propSim->spkEphem, baseBodyState);
    if ((int)observerInfo[0] == 500) {
        observerState[0] = (real) baseBodyState[0] + observerInfo[1];
        observerState[1] = (real) baseBodyState[1] + observerInfo[2];
        observerState[2] = (real) baseBodyState[2] + observerInfo[3];
        observerState[3] = (real) baseBodyState[3] + observerInfo[4];
        observerState[4] = (real) baseBodyState[4] + observerInfo[5];
        observerState[5] = (real) baseBodyState[5] + observerInfo[6];
        return;
    }
    std::string baseBodyFrame;
    get_baseBodyFrame((int)observerInfo[0], tObsMjdTDB, baseBodyFrame);
    std::vector<std::vector<real>> rotMat(6, std::vector<real>(6));
    get_pck_rotMat(baseBodyFrame, "J2000", tObsMjdTDB, propSim->pckEphem, rotMat);
    real lon = observerInfo[1];
    real lat = observerInfo[2];
    real rho = observerInfo[3];
    real bodyFixedX = rho * cos(lat) * cos(lon);
    real bodyFixedY = rho * cos(lat) * sin(lon);
    real bodyFixedZ = rho * sin(lat);
    std::vector<real> bodyFixedState = {bodyFixedX, bodyFixedY, bodyFixedZ,
                                            0.0,        0.0,        0.0};
    mat_vec_mul(rotMat, bodyFixedState, observerState);
    observerState[0] += baseBodyState[0];
    observerState[1] += baseBodyState[1];
    observerState[2] += baseBodyState[2];
    observerState[3] += baseBodyState[3];
    observerState[4] += baseBodyState[4];
    observerState[5] += baseBodyState[5];
}

/** 
 * @param[in] J2 Oblateness coefficient.
 * @param[in] poleRA Right ascension of the pole.
 * @param[in] poleDec Declination of the pole.
*/
void Body::set_J2(real J2, real poleRA, real poleDec) {
    this->J2 = J2;
    if (this->J2 != 0.0L) {
        this->isJ2 = true;
    } else {
        this->isJ2 = false;
    }
    this->poleRA = poleRA * DEG2RAD;
    this->poleDec = poleDec * DEG2RAD;
}

/** 
 * @param[in] poleRA Right ascension of the pole.
 * @param[in] poleDec Declination of the pole.
 * @param[in] nZon Degree of the zonal coefficients.
 * @param[in] nTes Order of the tesseral coefficients.
 * @param[in] J Zonal coefficient vector.
 * @param[in] C Sectoral coefficient array.
 * @param[in] S Tesseral coefficient array.
*/
void Body::set_harmonics(real poleRA, real poleDec, int nZon, int nTes,
                         std::vector<real> J, std::vector<std::vector<real>> C,
                         std::vector<std::vector<real>> S) {
    this->isHarmonic = true;
    if (this->isJ2) {
        throw std::invalid_argument(
            "Body::set_harmonics: Cannot set both J2 and spherical harmonics.");
    }
    this->poleRA = poleRA * DEG2RAD;
    this->poleDec = poleDec * DEG2RAD;
    if (nZon <= 0 || nTes <= 0) {
        throw std::invalid_argument(
            "Body::set_harmonics: The degree and order of the spherical harmonics "
            "must be positive.");
    }
    this->nZon = (size_t)nZon;
    this->nTes = (size_t)nTes;
    this->J = J;
    this->C = C;
    this->S = S;
}

/** 
 * @param[in] name Name of the body.
 * @param[in] spiceId SPICE ID of the body.
 * @param[in] t0 Initial time [MJD TDB].
 * @param[in] mass Mass of the body [kg].
 * @param[in] radius Radius of the body [m].
 */
SpiceBody::SpiceBody(std::string name, int spiceId, real t0, real mass,
                     real radius) {
    this->name = name;
    this->spiceId = spiceId;
    if (this->spiceId > 1000000) {
        this->caTol = 0.05;
    }
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius;
    this->isNongrav = false;
    this->isPPN = false;
    this->isMajor = false;
    // this->pos = {0.0L, 0.0L, 0.0L};
    this->pos[0] = 0.0L;
    this->pos[1] = 0.0L;
    this->pos[2] = 0.0L;
    // this->vel = {0.0L, 0.0L, 0.0L};
    this->vel[0] = 0.0L;
    this->vel[1] = 0.0L;
    this->vel[2] = 0.0L;
    // this->acc = {0.0L, 0.0L, 0.0L};
    this->acc[0] = 0.0L;
    this->acc[1] = 0.0L;
    this->acc[2] = 0.0L;
}

/**
 * @param[in] name Name of the body.
 * @param[in] t0 Initial time.
 * @param[in] mass Mass of the body [kg].
 * @param[in] radius Radius of the body [m].
 * @param[in] cometaryState Initial heliocentric ecliptic cometary state of the body.
 * @param[in] ngParams Nongravitational parameters for the body.
 */
IntegBody::IntegBody(std::string name, real t0, real mass, real radius,
                     std::vector<real> cometaryState,
                     NongravParameters ngParams) {
    this->name = name;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius;
    this->caTol = 0.0;
    std::vector<real> cartesianStateEclip(6);
    std::vector<real> cartesianPos(3);
    std::vector<real> cartesianVel(3);
    this->isCometary = true;
    this->initState = cometaryState;
    this->initCart = std::vector<real>(6, std::numeric_limits<real>::quiet_NaN());

    cometary_to_cartesian(t0, cometaryState, cartesianStateEclip);
    // rotate to eme2000
    std::vector<std::vector<real>> eclipToEquatorial(3, std::vector<real>(3));
    rot_mat_x(EARTH_OBLIQUITY, eclipToEquatorial);
    mat_vec_mul(eclipToEquatorial,
                {cartesianStateEclip[0], cartesianStateEclip[1],
                 cartesianStateEclip[2]},
                cartesianPos);
    mat_vec_mul(eclipToEquatorial,
                {cartesianStateEclip[3], cartesianStateEclip[4],
                 cartesianStateEclip[5]},
                cartesianVel);
    this->pos[0] = cartesianPos[0];
    this->pos[1] = cartesianPos[1];
    this->pos[2] = cartesianPos[2];
    this->vel[0] = cartesianVel[0];
    this->vel[1] = cartesianVel[1];
    this->vel[2] = cartesianVel[2];
    this->acc[0] = 0.0L;
    this->acc[1] = 0.0L;
    this->acc[2] = 0.0L;
    this->isNongrav = false;
    if (ngParams.a1 != 0.0L || ngParams.a2 != 0.0L || ngParams.a3 != 0.0L ||
            ngParams.a1Est || ngParams.a2Est || ngParams.a3Est) {
        this->ngParams.a1 = ngParams.a1;
        this->ngParams.a2 = ngParams.a2;
        this->ngParams.a3 = ngParams.a3;
        this->ngParams.a1Est = ngParams.a1Est;
        this->ngParams.a2Est = ngParams.a2Est;
        this->ngParams.a3Est = ngParams.a3Est;
        this->ngParams.alpha = ngParams.alpha;
        this->ngParams.k = ngParams.k;
        this->ngParams.m = ngParams.m;
        this->ngParams.n = ngParams.n;
        this->ngParams.r0_au = ngParams.r0_au;
        this->isNongrav = true;
    }
    this->isPPN = false;
    this->isMajor = false;
}

/**
 * @param[in] name Name of the body.
 * @param[in] t0 Initial time.
 * @param[in] mass Mass of the body [kg].
 * @param[in] radius Radius of the body [m].
 * @param[in] pos Initial barycentric position of the body [AU].
 * @param[in] vel Initial barycentric velocity of the body [AU/day].
 * @param[in] ngParams Nongravitational parameters for the body.
 */
IntegBody::IntegBody(std::string name, real t0, real mass, real radius,
                     std::vector<real> pos, std::vector<real> vel,
                     NongravParameters ngParams) {
    this->name = name;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius;
    this->caTol = 0.0;
    this->isCometary = false;
    this->initState = {pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]};
    this->initCart = initState;
    this->pos[0] = pos[0];
    this->pos[1] = pos[1];
    this->pos[2] = pos[2];
    this->vel[0] = vel[0];
    this->vel[1] = vel[1];
    this->vel[2] = vel[2];
    this->acc[0] = 0.0L;
    this->acc[1] = 0.0L;
    this->acc[2] = 0.0L;
    this->isNongrav = false;
    if (ngParams.a1 != 0.0L || ngParams.a2 != 0.0L || ngParams.a3 != 0.0L ||
            ngParams.a1Est || ngParams.a2Est || ngParams.a3Est) {
        this->ngParams.a1 = ngParams.a1;
        this->ngParams.a2 = ngParams.a2;
        this->ngParams.a3 = ngParams.a3;
        this->ngParams.a1Est = ngParams.a1Est;
        this->ngParams.a2Est = ngParams.a2Est;
        this->ngParams.a3Est = ngParams.a3Est;
        this->ngParams.alpha = ngParams.alpha;
        this->ngParams.k = ngParams.k;
        this->ngParams.m = ngParams.m;
        this->ngParams.n = ngParams.n;
        this->ngParams.r0_au = ngParams.r0_au;
        this->isNongrav = true;
    }
    this->isPPN = false;
    this->isMajor = false;
}

void IntegBody::prepare_stm(){
    int stmSize = 36;
    size_t numParams = 0;
    if (this->isNongrav) {
        if (ngParams.a1Est) {
            stmSize += 6;
            numParams++;
        }
        if (ngParams.a2Est) {
            stmSize += 6;
            numParams++;
        }
        if (ngParams.a3Est) {
            stmSize += 6;
            numParams++;
        }
    }
    this->stm = std::vector<real>(stmSize, 0.0L);
    for (size_t i = 0; i < 6; i++) {
        this->stm[6 * i + i] = 1.0L;
    }
    this->n2Derivs += (size_t) stmSize/2;
    this->propStm = true;
    // build 6x(6+numParams) state transition matrix first
    if (this->isCometary){
        std::vector<std::vector<real>> partialsEclip(6, std::vector<real>(6, 0.0L));
        get_cartesian_partials(this->t0, this->initState, "com2cart", partialsEclip);
        std::vector<std::vector<real>> bigRotMat(6, std::vector<real>(6, 0.0L));
        bigRotMat[0][0] = bigRotMat[3][3] = 1.0L;
        bigRotMat[1][1] = bigRotMat[4][4] = bigRotMat[2][2] = bigRotMat[5][5] = cos(EARTH_OBLIQUITY);
        bigRotMat[1][2] = bigRotMat[4][5] = -sin(EARTH_OBLIQUITY);
        bigRotMat[2][1] = bigRotMat[5][4] = sin(EARTH_OBLIQUITY);
        std::vector<std::vector<real>> partials(6, std::vector<real>(6, 0.0L));
        mat_mat_mul(bigRotMat, partialsEclip, partials);
        this->dCartdState = partials;
        for (size_t i = 0; i < 36; i++) {
            this->stm[i] = partials[i/6][i%6];
        }
        // add extra entries to each row for nongravitational parameters
        for (size_t i = 0; i < 6; i++) {
            for (size_t j = 0; j < numParams; j++) {
                this->dCartdState[i].push_back(0.0L);
            }
        }
    } else {
        // state transition matrix is the identity matrix in this case
        this->dCartdState = std::vector<std::vector<real>>(6, std::vector<real>(6+numParams, 0.0L));
        for (size_t i = 0; i < 6; i++) {
            this->dCartdState[i][i] = 1.0L;
        }
    }
    // add extra nongravitational parameter blocks at bottom for (6+numParams)x(6+numParams) state transition matrix
    for (size_t i = 0; i < numParams; i++) {
        this->dCartdState.push_back(std::vector<real>(6+numParams, 0.0L));
    }
    // set diagonals of these extra blocks to 1
    for (size_t i = 6; i < 6+numParams; i++) {
        this->dCartdState[i][i] = 1.0L;
    }
}

/**
 * @param[in] propSim PropSimulation object.
 * @param[in] t Time of the event.
 * @param[inout] xInteg State of the body.
 */
void Event::apply_impulsive(PropSimulation *propSim, const real &t, std::vector<real>& xInteg) {
    if (t != this->t) {
        throw std::runtime_error(
            "Event::apply_impulsive: Integration time does "
            "not match event time. Cannot apply impulse.");
    }
    size_t velStartIdx = this->xIntegIndex + 3;
    bool forwardProp = propSim->integParams.tf > propSim->integParams.t0;
    real propDir = forwardProp ? 1.0L : -1.0L;
    for (size_t i = 0; i < 3; i++) {
        xInteg[velStartIdx + i] += propDir * this->multiplier * this->deltaV[i];
    }
    if (propSim->integBodies[this->bodyIndex].propStm && this->eventEst) {
        IntegBody *body = &propSim->integBodies[this->bodyIndex];
        const size_t numNongravs = body->ngParams.a1Est + body->ngParams.a2Est + body->ngParams.a3Est;
        size_t DeltaVStmIdx = this->xIntegIndex + 6 + 36 + 6*numNongravs;
        for (size_t i = 0; i < propSim->eventMngr.impulsiveEvents.size(); i++) {
            const bool sameBody = propSim->eventMngr.impulsiveEvents[i].bodyIndex == this->bodyIndex;
            const bool hasStarted = propSim->eventMngr.impulsiveEvents[i].hasStarted;
            const bool deltaVEst = propSim->eventMngr.impulsiveEvents[i].deltaVEst;
            const bool multiplierEst = propSim->eventMngr.impulsiveEvents[i].multiplierEst;
            if (sameBody && hasStarted && multiplierEst) {
                DeltaVStmIdx += 6*1; // only multiplier
            } else if (sameBody && hasStarted && deltaVEst) {
                DeltaVStmIdx += 6*3; // full deltaV vector
            }
        }
        for (size_t i = 0; i < propSim->eventMngr.continuousEvents.size(); i++) {
            const bool sameBody = propSim->eventMngr.continuousEvents[i].bodyIndex == this->bodyIndex;
            const bool hasStarted = propSim->eventMngr.continuousEvents[i].hasStarted;
            const bool expAccel0Est = propSim->eventMngr.continuousEvents[i].expAccel0Est;
            const bool tauEst = propSim->eventMngr.continuousEvents[i].tauEst;
            if (sameBody && hasStarted && expAccel0Est) {
                DeltaVStmIdx += 6*3; // expAccel0 vector
            }
            if (sameBody && hasStarted && tauEst) {
                DeltaVStmIdx += 6*1; // time constant
            }
        }
        if (this->multiplierEst){
            xInteg[DeltaVStmIdx+3] = propDir * this->deltaV[0];
            xInteg[DeltaVStmIdx+4] = propDir * this->deltaV[1];
            xInteg[DeltaVStmIdx+5] = propDir * this->deltaV[2];
        } else {
            xInteg[DeltaVStmIdx+3] = propDir;
            xInteg[DeltaVStmIdx+10] = propDir;
            xInteg[DeltaVStmIdx+17] = propDir;
        }
    }
    this->hasStarted = true;
}

/**
 * @param[in] name Name of the simulation.
 * @param[in] t0 Initial time.
 * @param[in] defaultSpiceBodies Default SPICE bodies to load.
 *              (0=empty sim with DE440 ephemeris,
 *               430/421=DE430 Sun+planets+Moon+Pluto+big16 asteroids,
 *               440/441=DE440 Sun+planets+Moon+Pluto+big16 asteroids)
 * @param[in] DEkernelPath Path to SPICE metakernel
 */
PropSimulation::PropSimulation(std::string name, real t0,
                               const int defaultSpiceBodies,
                               std::string DEkernelPath) {
    this->name = name;
    this->DEkernelPath = DEkernelPath;
    this->integParams.t0 = t0;

    this->integParams.nInteg = 0;
    this->integParams.nSpice = 0;
    this->integParams.nTotal = 0;
    this->integParams.n2Derivs = 0;
    this->integParams.timestepCounter = 0;

    std::string kernel_sb, kernel_mb;
    switch (defaultSpiceBodies) {
        case 0: {
            kernel_sb = DEkernelPath + "sb441-n16s.bsp";
            kernel_mb = DEkernelPath + "de440.bsp";
            break;
        }
        // DE430 or DE431
        case 430:
        case 431: {
            kernel_sb = DEkernelPath + "sb431-n16s.bsp";
            kernel_mb = DEkernelPath + "de430.bsp";
            real G = 6.6743e-11L /
                (149597870700.0L * 149597870700.0L * 149597870700.0L) *
                86400.0L * 86400.0L;  // default kg au^3 / day^2
            // add planets and planetary bodies from DE431 header
            // (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de431_tech-comments.txt)
            SpiceBody Sun("Sun", 10, this->integParams.t0,
                          2.959122082855911e-4L / G, 6.96e8L);
            SpiceBody MercuryBarycenter("Mercury Barycenter", 1,
                                        this->integParams.t0,
                                        4.91248045036476e-11L / G, 2440.53e3L);
            SpiceBody VenusBarycenter("Venus Barycenter", 2,
                                      this->integParams.t0,
                                      7.24345233264412e-10L / G, 6051.8e3L);
            SpiceBody Earth("Earth", 399, this->integParams.t0,
                            8.887692445125634e-10L / G, 6378136.3L);
            SpiceBody Moon("Moon", 301, this->integParams.t0,
                           1.093189450742374e-11L / G, 1738.1e3L);
            SpiceBody MarsBarycenter("Mars Barycenter", 4, this->integParams.t0,
                                     9.54954869555077e-11L / G, 3396.19e3L);
            SpiceBody JupiterBarycenter("Jupiter Barycenter", 5,
                                        this->integParams.t0,
                                        2.82534584083387e-07L / G, 0.0L);
            SpiceBody SaturnBarycenter("Saturn Barycenter", 6,
                                       this->integParams.t0,
                                       8.45970607324503e-08L / G, 0.0L);
            SpiceBody UranusBarycenter("Uranus Barycenter", 7,
                                       this->integParams.t0,
                                       1.29202482578296e-08L / G, 0.0L);
            SpiceBody NeptuneBarycenter("Neptune Barycenter", 8,
                                        this->integParams.t0,
                                        1.52435734788511e-08L / G, 0.0L);
            SpiceBody PlutoBarycenter("Pluto Barycenter", 9,
                                      this->integParams.t0,
                                      2.17844105197418e-12L / G, 0.0L);
            Sun.isPPN = true;
            Sun.isMajor = true;
            MercuryBarycenter.isPPN = true;
            MercuryBarycenter.isMajor = true;
            VenusBarycenter.isPPN = true;
            VenusBarycenter.isMajor = true;
            Earth.isPPN = true;
            Earth.isMajor = true;
            Moon.isPPN = true;
            Moon.isMajor = true;
            MarsBarycenter.isPPN = true;
            MarsBarycenter.isMajor = true;
            JupiterBarycenter.isPPN = true;
            JupiterBarycenter.isMajor = true;
            SaturnBarycenter.isPPN = true;
            SaturnBarycenter.isMajor = true;
            UranusBarycenter.isPPN = true;
            UranusBarycenter.isMajor = true;
            NeptuneBarycenter.isPPN = true;
            NeptuneBarycenter.isMajor = true;
            PlutoBarycenter.isPPN = true;
            PlutoBarycenter.isMajor = true;
            Sun.set_J2(
                2.1106088532726840e-7L, 286.13L,
                63.87L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            Earth.set_J2(
                0.00108262545L, 0.0L,
                90.0L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            Moon.set_J2(
                2.0321568464952570e-4L, 269.9949L,
                66.5392L);  // (fixed representation ONLY accurate at J2000 - see Table 2 in link for details.) https://doi.org/10.1007/s10569-010-9320-4
            Sun.caTol = 0.25;
            MercuryBarycenter.caTol = 0.1;
            VenusBarycenter.caTol = 0.1;
            Earth.caTol = 0.1;
            Moon.caTol = 0.05;
            MarsBarycenter.caTol = 0.1;
            JupiterBarycenter.caTol = 0.25;
            SaturnBarycenter.caTol = 0.25;
            UranusBarycenter.caTol = 0.25;
            NeptuneBarycenter.caTol = 0.25;
            PlutoBarycenter.caTol = 0.1;
            add_spice_body(Sun);
            add_spice_body(MercuryBarycenter);
            add_spice_body(VenusBarycenter);
            add_spice_body(Earth);
            add_spice_body(Moon);
            add_spice_body(MarsBarycenter);
            add_spice_body(JupiterBarycenter);
            add_spice_body(SaturnBarycenter);
            add_spice_body(UranusBarycenter);
            add_spice_body(NeptuneBarycenter);
            add_spice_body(PlutoBarycenter);

            // add DE431 big16 asteroids from JPL sb431-big16s.bsp, mass
            // parameters from DE431 header
            // (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de431_tech-comments.txt)
            SpiceBody Ceres("Ceres", 2000001, this->integParams.t0,
                            1.400476556172344e-13L / G, 0.0L);
            SpiceBody Vesta("Vesta", 2000004, this->integParams.t0,
                            3.85475018780881e-14L / G, 0.0L);
            SpiceBody Pallas("Pallas", 2000002, this->integParams.t0,
                             3.104448198938713e-14L / G, 0.0L);
            SpiceBody Hygiea("Hygiea", 2000010, this->integParams.t0,
                             1.235800787294125e-14L / G, 0.0L);
            SpiceBody Euphrosyne("Euphrosyne", 2000031, this->integParams.t0,
                                 6.343280473648602e-15L / G, 0.0L);
            SpiceBody Interamnia("Interamnia", 2000704, this->integParams.t0,
                                 5.256168678493662e-15L / G, 0.0L);
            SpiceBody Davida("Davida", 2000511, this->integParams.t0,
                             5.198126979457498e-15L / G, 0.0L);
            SpiceBody Eunomia("Eunomia", 2000015, this->integParams.t0,
                              4.678307418350905e-15L / G, 0.0L);
            SpiceBody Juno("Juno", 2000003, this->integParams.t0,
                           3.617538317147937e-15L / G, 0.0L);
            SpiceBody Psyche("Psyche", 2000016, this->integParams.t0,
                             3.411586826193812e-15L / G, 0.0L);
            SpiceBody Cybele("Cybele", 2000065, this->integParams.t0,
                             3.180659282652541e-15L / G, 0.0L);
            SpiceBody Thisbe("Thisbe", 2000088, this->integParams.t0,
                             2.577114127311047e-15L / G, 0.0L);
            SpiceBody Doris("Doris", 2000048, this->integParams.t0,
                            2.531091726015068e-15L / G, 0.0L);
            SpiceBody Europa("Europa", 2000052, this->integParams.t0,
                             2.476788101255867e-15L / G, 0.0L);
            SpiceBody Patientia("Patientia", 2000451, this->integParams.t0,
                                2.295559390637462e-15L / G, 0.0L);
            SpiceBody Sylvia("Sylvia", 2000087, this->integParams.t0,
                             2.199295173574073e-15L / G, 0.0L);
            Ceres.caTol = 0.05;
            Vesta.caTol = 0.05;
            Pallas.caTol = 0.05;
            add_spice_body(Ceres);
            add_spice_body(Vesta);
            add_spice_body(Pallas);
            add_spice_body(Hygiea);
            add_spice_body(Euphrosyne);
            add_spice_body(Interamnia);
            add_spice_body(Davida);
            add_spice_body(Eunomia);
            add_spice_body(Juno);
            add_spice_body(Psyche);
            add_spice_body(Cybele);
            add_spice_body(Thisbe);
            add_spice_body(Doris);
            add_spice_body(Europa);
            add_spice_body(Patientia);
            add_spice_body(Sylvia);
            break;
        }
        // DE440 or DE441
        case 440:
        case 441: {
            kernel_sb = DEkernelPath + "sb441-n16s.bsp";
            kernel_mb = DEkernelPath + "de440.bsp";
            if (defaultSpiceBodies == 441) {
                std::cout
                    << "WARNING: Choosing DE441 will load long-term ephemeris "
                       "with a significantly higher memory footprint. "
                       "For almost all applications, DE440 is sufficient. "
                       "Be sure that you need the coverage provided by DE441."
                    << std::endl;
                kernel_sb = DEkernelPath + "sb441-n16.bsp";
                kernel_mb = DEkernelPath + "de441.bsp";
            }
            real G = 6.6743e-11L /
                (149597870700.0L * 149597870700.0L * 149597870700.0L) *
                86400.0L * 86400.0L;  // default kg au^3 / day^2
            // add planets and planetary bodies from DE441 header
            // (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_tech-comments.txt)
            SpiceBody Sun("Sun", 10, this->integParams.t0,
                          2.9591220828411956e-04L / G, 6.96e8L);
            SpiceBody MercuryBarycenter(
                "Mercury Barycenter", 1, this->integParams.t0,
                4.9125001948893182e-11L / G, 2440.53e3L);
            SpiceBody VenusBarycenter("Venus Barycenter", 2,
                                      this->integParams.t0,
                                      7.2434523326441187e-10L / G, 6051.8e3L);
            SpiceBody Earth("Earth", 399, this->integParams.t0,
                            8.8876924467071033e-10L / G, 6378136.6L);
            SpiceBody Moon("Moon", 301, this->integParams.t0,
                           1.0931894624024351e-11L / G, 1738.1e3L);
            SpiceBody MarsBarycenter("Mars Barycenter", 4, this->integParams.t0,
                                     9.5495488297258119e-11L / G, 3396.19e3L);
            SpiceBody JupiterBarycenter("Jupiter Barycenter", 5,
                                        this->integParams.t0,
                                        2.8253458252257917e-07L / G, 0.0L);
            SpiceBody SaturnBarycenter("Saturn Barycenter", 6,
                                       this->integParams.t0,
                                       8.4597059933762903e-08L / G, 0.0L);
            SpiceBody UranusBarycenter("Uranus Barycenter", 7,
                                       this->integParams.t0,
                                       1.2920265649682399e-08L / G, 0.0L);
            SpiceBody NeptuneBarycenter("Neptune Barycenter", 8,
                                        this->integParams.t0,
                                        1.5243573478851939e-08L / G, 0.0L);
            SpiceBody PlutoBarycenter("Pluto Barycenter", 9,
                                      this->integParams.t0,
                                      2.1750964648933581e-12L / G, 0.0L);
            Sun.isPPN = true;
            Sun.isMajor = true;
            MercuryBarycenter.isPPN = true;
            MercuryBarycenter.isMajor = true;
            VenusBarycenter.isPPN = true;
            VenusBarycenter.isMajor = true;
            Earth.isPPN = true;
            Earth.isMajor = true;
            Moon.isPPN = true;
            Moon.isMajor = true;
            MarsBarycenter.isPPN = true;
            MarsBarycenter.isMajor = true;
            JupiterBarycenter.isPPN = true;
            JupiterBarycenter.isMajor = true;
            SaturnBarycenter.isPPN = true;
            SaturnBarycenter.isMajor = true;
            UranusBarycenter.isPPN = true;
            UranusBarycenter.isMajor = true;
            NeptuneBarycenter.isPPN = true;
            NeptuneBarycenter.isMajor = true;
            PlutoBarycenter.isPPN = true;
            PlutoBarycenter.isMajor = true;
            Sun.set_J2(
                2.1961391516529825e-07L, 286.13L,
                63.87L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            Earth.set_J2(
                1.0826253900000000e-03L, 0.0L,
                90.0L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            // Moon.set_J2(
            //     2.0321568464952570e-4L, 269.9949L,
            //     66.5392L);  // Table 2 in link for details: https://doi.org/10.1007/s10569-010-9320-4
            // because of our implementation of Legendre polynomials, we need to pad harmonic coefficients with zeros so that indexing is correct
            // std::vector<real> moonJ = {
            //     0.0L,
            //     0.0L,
            //     2.0321568464952570e-4L, 
            //     8.4597026974594570e-6L,
            //    -9.7044138365700000e-6L,
            //     7.4221608384052890e-7L,
            //    -1.3767531350969900e-5L};
            // std::vector<std::vector<real>> moonC = {
            //     {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L},
            //     {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L},
            //     {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L},
            //     {0.0L,  2.8480741195592860e-5L,  4.8449420619770600e-6L,  1.6756178134114570e-6L,  0.0L, 0.0L, 0.0L},
            //     {0.0L, -5.7048697319733210e-6L, -1.5912271792977430e-6L, -8.0678881596778210e-8L, -1.2692158612216040e-7L,  0.0L, 0.0L},
            //     {0.0L, -8.6629769308983560e-7L,  7.1199537967353330e-7L,  1.5399750424904520e-8L,  2.1444704319218450e-8L,  7.6596153884006140e-8L,  0.0L},
            //     {0.0L,  1.2024363601545920e-6L, -5.4703897324156850e-7L, -6.8785612757292010e-8L,  1.2915580402925160e-9L,  1.1737698784460500e-9L, -1.0913395178881540e-9L}};
            // std::vector<std::vector<real>> moonS = {
            //     {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L},
            //     {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L},
            //     {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L},
            //     {0.0L,  5.8915551555318640e-6L,  1.6844743962783900e-6L, -2.4742714379805760e-7L,  0.0L, 0.0L, 0.0L},
            //     {0.0L,  1.5789202789245720e-6L, -1.5153915796731720e-6L, -8.0349266627431070e-7L,  8.2964257754075220e-8L,   0.0L, 0.0L},
            //     {0.0L, -3.5272289393243820e-6L,  1.7107886673430380e-7L,  2.8736257616334340e-7L,  5.2652110720146800e-10L, -6.7824035473995330e-9L,  0.0L},
            //     {0.0L, -2.0453507141252220e-6L, -2.6966834353574270e-7L, -7.1063745295915780e-8L, -1.5361616966632300e-8L,  -8.3465073195142520e-9L,  1.6844213702632920e-9L}};
            std::vector<real> moonJ = {
                0.0L,
                -0.0L,
                9.088124807048e-05L,
                3.19744528042e-06L,
                -3.234795995974e-06L,
                2.237862962417e-07L,
                -3.818428425858e-06L,
                -5.593411497042e-06L,
                -2.346826891883e-06L,
                3.530906115765e-06L,
                1.069295168104e-06L,
                8.852069967396e-07L,
                2.012261355675e-06L,
                -2.564870398577e-07L,
                -4.922527542979e-07L,
                2.280906452514e-07L,
                -5.568503510668e-07L,
                1.033029690169e-06L,
                5.203178351862e-07L,
                -2.156514659605e-07L,
                -3.355607188205e-07L,
                5.601367460853e-08L,
                8.414190265392e-08L,
                3.854355010907e-07L,
                -3.989774022946e-09L,
                2.157607659606e-07L,
                -2.838753756729e-08L,
                6.35484364141e-07L,
                -8.388624419412e-07L,
                6.02391868835e-07L,
                -1.013251448423e-07L,
                -5.429394298733e-07L,
                3.222792355569e-07L,
                -1.080678334159e-07L,
                -1.662297092652e-08L,
                -1.294211920192e-07L,
                -1.673671057062e-07L,
                -2.240662871701e-07L,
                -4.577578937312e-08L,
                3.34287149932e-07L,
                -3.06095510128e-08L,
            };
            std::vector<std::vector<real>> moonC = {
                {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 2.351746492573e-10L, 3.467408607293e-05L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 2.636795509644e-05L, 1.417149107237e-05L, 1.227505227804e-05L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -6.013467637411e-06L, -7.116169094765e-06L, -1.349974758249e-06L, -6.007023544009e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -1.011602594693e-06L, 4.399518982025e-06L, 4.661575916633e-07L, 2.754160545924e-06L, 3.110818615395e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.528258612943e-06L, -4.397301231474e-06L, -3.317525403334e-06L, 3.412141599126e-07L, 1.454413555888e-06L, -4.684211790942e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 7.471664075656e-06L, -6.501523120284e-07L, 5.99421656161e-07L, -8.437047109434e-07L, -2.068084929814e-07L, -1.065340249576e-06L, -1.820256920974e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 4.169368576832e-09L, 3.009314074362e-06L, -1.889048770847e-06L, 3.408687787676e-06L, -1.248066862873e-06L, -1.660417087399e-06L, -1.509668339879e-06L, -2.485685629145e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.866977479281e-06L, 1.927802812038e-06L, -1.992424249344e-06L, -1.884438321323e-06L, -1.562498503898e-06L, -2.127247181721e-06L, -3.914775916875e-06L, -1.311969441463e-06L, -9.383999626693e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 8.41602097663e-07L, 3.572446773007e-07L, 4.841890206877e-07L, -3.572972046925e-06L, 6.99695568131e-07L, -1.272846830237e-07L, -3.99864583838e-06L, -3.55945407182e-06L, -4.753141318719e-06L, 9.479086746857e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -1.94605183845e-08L, 6.995183265342e-07L, 4.084141672369e-07L, -1.086778059564e-06L, 5.033765025768e-08L, 5.243131841919e-07L, -8.450655201886e-08L, -2.095262600571e-06L, -2.334956345061e-06L, -4.680541848292e-06L, -2.888827942133e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -7.068742240042e-07L, -1.052083643652e-07L, 8.55703198874e-07L, 8.130117985304e-07L, 2.052037807514e-10L, 1.089078743533e-06L, 2.139998468072e-06L, 3.985559854863e-07L, -1.345190104344e-06L, -3.07548291113e-06L, -9.818185320143e-07L, 3.027757498068e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.354460596788e-06L, -3.425239492101e-06L, -4.333433404851e-07L, 7.772611476568e-07L, -1.186344628059e-06L, -2.76918839766e-08L, -1.621162213111e-07L, -2.92268559879e-07L, -1.699728664781e-07L, -6.742710165203e-07L, -1.305704137695e-06L, -5.743621652207e-07L, 2.499379485047e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 5.65426394844e-07L, 3.644201662568e-07L, 5.578918008017e-07L, -3.905263791569e-07L, -7.823944826099e-07L, -6.096144345758e-07L, -5.154294731647e-07L, 3.061434396201e-07L, 8.050002610362e-08L, -7.782150514274e-07L, -1.928506836209e-06L, -1.448531525908e-06L, -3.214009044084e-07L, -6.135005506804e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -7.975286301595e-07L, -5.401817264903e-08L, -1.37824413033e-06L, -9.560533835178e-07L, 7.743099756597e-09L, -6.443036337316e-08L, 1.125282382651e-06L, 1.315092334932e-06L, -8.52185835898e-08L, -6.013405127931e-08L, -1.426274722839e-06L, -1.767194004831e-06L, -4.049321212081e-07L, 8.234518677626e-07L, 5.483004282487e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.155919019965e-07L, 1.455232627417e-06L, -1.942731771602e-07L, 6.859523352135e-07L, 7.63423122156e-07L, 7.298127159692e-07L, -2.807030679633e-07L, 2.371112099329e-08L, -7.806793280878e-07L, 2.272112129901e-07L, -3.314720765542e-07L, -1.294179856103e-06L, -1.798844599359e-07L, -6.798745684977e-07L, -3.647704432693e-07L, -9.053671721027e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 5.020203035772e-07L, -2.90725536744e-07L, -2.677029785696e-08L, 1.285116411201e-06L, 1.906694406275e-07L, 7.212398186701e-07L, -1.298775653101e-06L, -2.114532843114e-07L, 3.28203095655e-07L, 4.162934093892e-07L, 1.143965328344e-06L, 9.594776762101e-07L, -3.971860167035e-08L, 4.031162045124e-07L, -6.434950649787e-08L, -3.169632477025e-07L, -3.807578982455e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 3.604168210186e-07L, -3.166161312163e-07L, 9.952234784658e-07L, -8.030666066195e-07L, -5.193042716873e-07L, -1.701307917456e-06L, -2.183009712261e-07L, 4.44591078384e-07L, 5.507266084873e-07L, 3.581956577172e-07L, 4.212978863805e-07L, 8.308545262477e-07L, -1.215831435917e-06L, 6.112330835446e-08L, 2.43220611875e-07L, 8.988119283635e-07L, 1.230015704844e-06L, -2.097047215516e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -5.912414924324e-07L, 1.730032734268e-07L, -6.319853354347e-07L, -9.93486748579e-07L, -2.584465868218e-07L, 1.587099926537e-07L, 1.223449136041e-06L, 6.51193595312e-07L, -6.53875478659e-08L, -7.824230071427e-07L, 9.199853228705e-09L, -1.104366886538e-07L, 7.449839699536e-07L, -1.260595067688e-07L, 8.387886846523e-08L, -5.517925047684e-08L, -1.068986532109e-06L, 1.123829856611e-06L, 6.020542775342e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -1.970524875472e-07L, 8.12799423525e-07L, 6.352115726265e-07L, 4.635540820733e-07L, 2.189821221788e-07L, -1.913179816892e-07L, -6.759398002022e-07L, 2.345666656895e-07L, -3.846565414087e-07L, 3.259859516279e-07L, 2.452891638676e-07L, 1.311401534009e-07L, 8.30471361439e-07L, -1.462339943485e-07L, -3.247069066381e-07L, 9.283544937345e-07L, 7.932690737593e-07L, 8.449720341132e-07L, -2.172424955841e-07L, -1.352338241923e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 3.774308358991e-08L, -1.363839776102e-06L, 3.452768501303e-07L, -7.638070575592e-07L, 3.333231702919e-07L, -6.229457497484e-08L, -3.753367279108e-07L, -8.112958446564e-07L, 6.075614239204e-08L, 4.792892268783e-07L, -3.1497268846e-07L, -2.937670155695e-07L, -2.713801718399e-07L, -1.590905655069e-07L, -3.768494751399e-07L, 1.135223674956e-06L, 8.691153673739e-07L, -1.115487402416e-07L, 1.872114258695e-06L, -4.321492184803e-07L, -3.341158957945e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -3.028755435545e-07L, -3.629235379091e-07L, 5.907282247701e-07L, 3.559481082779e-07L, 4.76698555017e-07L, -2.359536187304e-07L, 3.693191420541e-07L, 1.75020144185e-07L, 7.384647227329e-07L, -3.856865873692e-07L, -3.651356410422e-07L, -2.245225200427e-09L, -3.394126109466e-07L, -9.273232373439e-07L, -2.127150658291e-07L, -6.187348189655e-07L, -7.09649726023e-07L, 5.18019066114e-07L, -4.981089499017e-07L, 6.666281680379e-07L, -1.686492753703e-07L, 3.561470159161e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 2.835896408964e-07L, 4.572665687412e-07L, -6.957931697086e-07L, 4.113716965752e-07L, -2.263619442551e-07L, 4.02166746077e-07L, 4.155030183778e-07L, 3.4876828744e-07L, 1.014257500512e-06L, -6.990140627732e-07L, 2.963144368382e-08L, 3.164223870302e-07L, -7.794685415185e-07L, 7.523151505349e-07L, 2.736241696122e-07L, -5.725485774862e-07L, -3.306379268087e-07L, -3.368904588589e-08L, -4.343098201125e-09L, -1.993014501753e-08L, -5.912061681155e-07L, -9.859864477429e-07L, -1.061370338095e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 3.390385136769e-07L, 4.213812484417e-07L, -1.04690397223e-07L, 5.201617164354e-07L, 4.707346046077e-08L, -3.801458303109e-07L, 3.511100788204e-07L, 5.363182834755e-07L, -6.108208182022e-07L, -4.291169936679e-07L, 2.068941538289e-07L, 1.204887765201e-06L, 3.326046820092e-07L, 8.169144600841e-07L, -1.837499477974e-07L, 2.619307633942e-07L, 5.849810150378e-07L, -9.9560585766e-07L, 4.64187926164e-07L, -3.492952299506e-07L, 2.985474035568e-07L, 1.206920055379e-07L, -1.769404529862e-08L, 9.708785680645e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -4.856612326985e-08L, -6.613634818171e-07L, 6.743636549992e-07L, -5.790562860868e-07L, 2.006091722542e-07L, -2.538754741364e-07L, 6.908558703488e-07L, -9.987414937857e-07L, -3.555481537757e-07L, -2.236217583802e-08L, -2.590970273395e-08L, 1.279323280421e-07L, 3.179845421562e-07L, -6.398529031328e-07L, 1.474221525484e-07L, 2.204943590248e-07L, 3.119674181768e-08L, -8.457427399928e-08L, -8.187399871193e-07L, -4.563800361606e-08L, 6.825278183982e-07L, -3.386910242799e-07L, 9.260030914262e-07L, 9.225008020585e-07L, -1.316505167745e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.639405423885e-07L, 6.585239730648e-07L, 1.932719970775e-07L, 3.465877061304e-07L, 6.030059905384e-08L, -5.199124378934e-07L, 1.15423317525e-08L, -1.141569499357e-09L, 5.844553853316e-08L, -3.097486607555e-07L, 1.635456466728e-07L, -4.246322291764e-07L, 1.320975033515e-07L, -1.204957014134e-07L, 4.413205227482e-08L, -1.028958209046e-07L, -3.584979702101e-07L, 1.051281979747e-08L, -1.941276816526e-07L, -2.601229914951e-07L, -5.301560066899e-07L, -9.320667699728e-07L, 1.034710158199e-07L, -7.25128255457e-07L, 7.446804621392e-08L, 1.0973321842e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 6.163246033951e-07L, -2.049868493341e-07L, 5.601074438448e-07L, -4.00719770898e-07L, -6.904474312355e-07L, 6.475012349976e-09L, 2.573384530733e-07L, 4.12571469733e-08L, -6.17581297644e-07L, 6.558090152813e-08L, 3.453504082105e-07L, -4.400926492876e-07L, 2.356873847179e-07L, 4.592990662147e-07L, 6.181690038861e-08L, 1.513968700232e-07L, 1.881175179844e-07L, 2.055281180856e-07L, 6.839344751758e-08L, 3.941816742625e-07L, 4.419579571303e-07L, 7.905759927657e-07L, -6.487163671831e-08L, 7.645874668034e-07L, -3.737481954021e-07L, -4.488858781848e-07L, -4.229476090756e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 3.286764182251e-08L, -3.379945631701e-08L, -3.093708793202e-07L, -9.16863200467e-07L, 1.48055262858e-07L, 4.680854865667e-07L, -3.764489087056e-07L, -1.411192792247e-07L, -4.297187540611e-07L, -3.481604074913e-07L, 1.869361099809e-07L, 1.144267729633e-07L, 3.91768668552e-07L, 3.618289021574e-07L, -4.414140464284e-07L, -2.387235418067e-07L, -7.974229728158e-08L, -8.597938702178e-08L, -1.629407241524e-07L, 2.799707706537e-07L, -2.576011106534e-07L, -6.204542471532e-07L, -4.750675407146e-09L, 4.065281331427e-07L, 4.813585956337e-07L, 1.334041959922e-07L, -1.596403080344e-07L, 1.338710805338e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -1.595902297753e-07L, -2.032704704877e-07L, -1.966906244805e-07L, -3.026314589709e-07L, 2.196226357699e-07L, -1.277211420108e-07L, -1.076111980842e-07L, 4.441189394194e-07L, 6.016659821326e-08L, 6.210015560552e-07L, -7.84582114152e-08L, -6.437377679767e-07L, 3.622952091543e-07L, 1.881809027494e-08L, -1.994365435076e-07L, 2.676200695387e-07L, -4.044400430658e-07L, -1.314893037868e-07L, 1.090779063386e-08L, -6.74990271558e-09L, 6.295383367731e-08L, -2.440660676362e-07L, 1.576505487583e-07L, 1.556737392332e-08L, -1.3931765224e-07L, 2.643207120599e-07L, 2.061000686885e-07L, -1.151697040428e-07L, 1.842897528242e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -2.944806342383e-07L, -2.164408281118e-07L, -3.616709171735e-07L, -6.226153976562e-08L, 1.857938037775e-07L, 2.4001762416e-08L, -3.448247530555e-07L, -3.853670733082e-07L, -2.676142129194e-07L, -1.147785064062e-07L, 7.085322278412e-09L, -2.73800079521e-07L, 1.114358420441e-07L, 4.997250919797e-07L, -3.408136160889e-07L, 2.355955042331e-07L, 3.276888141725e-07L, 2.961864439237e-07L, -1.617130169838e-07L, -1.544594983541e-07L, -2.536365073624e-07L, -1.397852721084e-08L, -4.38234795544e-07L, 1.346252365618e-07L, 1.578037099937e-07L, 5.512390731824e-08L, -1.221456324711e-07L, -1.378661064142e-07L, 1.692730547791e-07L, -3.241990681921e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 8.304183422004e-09L, 2.466662005158e-08L, -6.043715717861e-08L, 4.317752365093e-07L, -3.994563144302e-07L, 2.072085801081e-07L, 3.603249206912e-07L, 4.914145167184e-07L, -4.118540899566e-07L, 7.913209577978e-07L, 1.324337359124e-07L, -1.460120361694e-07L, 3.453594453011e-07L, 5.021436950256e-08L, -4.042141245692e-07L, -7.262407233629e-09L, -1.497439509134e-07L, 2.824671732706e-07L, -6.023122203834e-09L, -3.631991536621e-08L, -2.80224312823e-07L, 2.987661559268e-07L, 1.086568285506e-06L, 1.515747923514e-07L, -6.110666730027e-07L, -3.259384983049e-07L, 5.47041326435e-08L, 5.204424777774e-08L, 4.403742514956e-07L, 1.152856858309e-08L, 1.965100644866e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -2.841313891525e-07L, -6.314600994058e-07L, 9.956305039002e-08L, 1.716138890666e-08L, -1.748538122623e-07L, 2.080848285772e-07L, -1.795784172654e-07L, 3.644992867629e-08L, -1.149970226755e-07L, 1.452671518142e-07L, -1.224363438307e-07L, 3.243095708245e-07L, -9.946424621774e-08L, 1.849980665985e-07L, -3.71093061671e-07L, -9.766882893012e-08L, -4.109513816319e-07L, 3.958918500978e-07L, 1.558715576152e-08L, 3.643494175144e-07L, 3.490224836617e-07L, -4.523972033178e-07L, 1.747306287398e-07L, -1.154723965902e-07L, -1.127750150001e-07L, 5.650968517176e-07L, 1.180633857846e-07L, 4.269933888287e-07L, 1.761023943193e-07L, -3.154556680814e-07L, 3.430622426277e-07L, 1.625251864049e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.840000385064e-07L, 6.716259248662e-07L, 1.018104651955e-07L, 1.45226565844e-07L, 2.677950359236e-07L, -2.195707246077e-07L, 6.033335129823e-07L, 2.433137228845e-07L, 3.609398018201e-07L, -2.587252703994e-07L, 1.042435954498e-07L, 2.836595081553e-07L, -3.837477881591e-07L, -2.706577782075e-07L, -1.816237945769e-07L, 3.840969598269e-07L, 2.760230368454e-07L, -3.639044584903e-08L, 2.814024115261e-07L, 5.204562236718e-08L, -6.163690233555e-07L, -3.460239575159e-08L, 1.106203583824e-07L, -2.622263473138e-07L, 3.254129142283e-07L, 3.133079770763e-07L, 6.265497710812e-08L, -5.367543564807e-07L, -1.462665081894e-07L, -1.509151375444e-07L, -2.606062376547e-07L, -5.049067932428e-07L, 1.778650345322e-08L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 2.821075064704e-07L, 2.204609001406e-07L, -3.543457193257e-07L, 1.938518248592e-07L, 2.081103667849e-07L, 2.591523672961e-07L, 2.014870768163e-07L, -3.949689929297e-07L, 2.295124013566e-07L, 4.513107718146e-08L, 1.507421756571e-07L, 7.441648974769e-08L, -4.215658091637e-08L, -8.439966288462e-08L, -8.22913958492e-08L, 2.908753605209e-07L, 8.247527158938e-08L, -2.303959681161e-07L, -3.683069431367e-07L, -2.440688996187e-07L, 5.249279801868e-07L, 6.275905012396e-07L, -4.041287587356e-07L, -8.046112879803e-08L, 8.335791607322e-08L, -3.167032572514e-07L, -1.323183811286e-08L, -1.803958350509e-07L, 1.948883248412e-07L, -1.18864594199e-07L, 3.071884151071e-07L, -4.041355610794e-08L, -2.033619011149e-07L, 9.468602446773e-09L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 3.08651945893e-07L, -3.380244138245e-07L, 5.052522197209e-08L, 2.953237511278e-07L, 1.01981182845e-07L, -1.351143459128e-07L, -6.497676247145e-09L, 6.944972011106e-08L, 4.419121434856e-07L, -4.056994828217e-07L, -3.819284083219e-07L, 1.765795229026e-07L, -3.945354940515e-07L, 2.866383969217e-07L, 9.19980646963e-08L, 2.257671272459e-07L, 2.850602419426e-07L, -2.72031938076e-07L, -4.399411129254e-08L, -8.596806226926e-08L, 5.087103660808e-07L, -1.871243561882e-07L, -8.658671210836e-08L, -7.473634626882e-09L, 1.219766952361e-08L, -2.178455893829e-07L, 1.127443238304e-07L, 1.777165457753e-07L, 5.295313802421e-07L, 1.130014575145e-07L, 3.097786946665e-07L, -4.087367399657e-07L, 3.073354436315e-07L, 3.403369978317e-07L, -6.203652984923e-08L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -2.87499658204e-07L, 2.367010434623e-07L, 6.775016707304e-07L, -7.530002290689e-08L, -9.194554700051e-09L, -3.023330217233e-07L, -5.141008980657e-08L, -1.099079179674e-09L, -1.782981398325e-07L, 2.938509462004e-08L, 1.520336172235e-07L, -2.11774127483e-07L, -4.275923255779e-07L, 7.66576872491e-08L, 8.665137172898e-08L, -1.613174132406e-07L, -4.966567484407e-07L, -3.71520548425e-08L, 2.861205192288e-07L, -3.135101121129e-07L, 4.529518111558e-07L, -8.909228915519e-08L, -3.182676758167e-07L, 4.077332078719e-08L, -2.647035550822e-07L, 8.482586524913e-08L, 2.93107833983e-07L, 6.807487078475e-08L, -1.402900469201e-07L, -1.59171967976e-07L, -2.119155328096e-07L, -4.628848921532e-08L, 4.747981345985e-07L, -7.740560148733e-08L, 3.162467405698e-07L, 2.178276488934e-08L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 2.417538713476e-08L, -2.491248706017e-07L, -2.507631488025e-07L, -3.480796446282e-08L, 2.735367201607e-07L, 4.355893660955e-08L, 4.41770216548e-07L, -8.507978882833e-08L, -2.030218724412e-07L, 8.216279574697e-07L, 4.091009427739e-08L, -1.24426880443e-07L, 3.682872723138e-08L, 1.921334507518e-07L, 9.016909967937e-08L, -1.091044826627e-07L, 9.347912364753e-08L, 1.6899827287e-07L, -1.97499543334e-07L, -2.318086686862e-07L, -2.041035474606e-07L, 1.617808291856e-08L, 2.982790833133e-07L, -1.714617967475e-07L, -1.236453360901e-07L, 7.375722297842e-08L, -2.266847253858e-07L, -1.788810154232e-07L, -1.115380601188e-07L, -1.223490589359e-07L, -4.050829603485e-08L, 2.210540632798e-07L, 1.211853302505e-07L, -7.495097594831e-08L, 3.89106465138e-07L, -4.645242933507e-07L, 1.995588035727e-07L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 3.687936048555e-08L, 1.327264642699e-08L, -9.281077466455e-08L, -5.436838041953e-08L, 1.917527887537e-07L, 4.860268978087e-08L, -4.463964761623e-07L, -5.542568279733e-07L, 1.65286891255e-07L, 1.816001385249e-07L, -1.322750463105e-07L, -9.612102178565e-08L, -8.253192977785e-08L, 8.835064815704e-09L, -2.939277365377e-08L, 3.725899207366e-08L, -4.984457218139e-08L, 2.743873866582e-07L, 4.664476362031e-07L, -2.773034157495e-07L, 8.213253392241e-08L, -1.170497482922e-07L, 2.061221132139e-08L, 7.25990488422e-08L, -5.326711593721e-09L, -5.653593487278e-08L, 1.062339224208e-07L, 2.986024496247e-07L, 4.530749265587e-08L, -4.723404681276e-07L, -6.654523443189e-08L, 1.060055359561e-07L, -1.800153020971e-07L, -1.291807945674e-07L, -5.062547832979e-08L, -3.826626225491e-07L, -1.664575995012e-08L, -3.203343390791e-07L, 0.0L, 0.0L, },
                {0.0L, -1.686558705431e-07L, -1.567707944034e-07L, 2.371818417063e-07L, 2.284942752381e-08L, -5.89015234993e-09L, 1.711356677333e-07L, -6.937357081146e-08L, 1.844435077992e-07L, -2.365934054387e-08L, -2.471380210442e-07L, -2.707358752071e-08L, -1.167776409852e-07L, 1.338442232596e-07L, 2.994036917555e-07L, -7.799686103213e-08L, 2.304236762608e-07L, -2.439663543089e-07L, 8.99189995809e-08L, 5.203905584616e-08L, -4.820642258629e-07L, 2.370120653456e-07L, 1.53216793987e-08L, 3.974105294128e-07L, -3.114601073054e-07L, -1.642643267979e-07L, 8.291260338137e-08L, 1.441876120204e-07L, -2.408243540694e-08L, -4.989899797356e-08L, 1.047286402491e-07L, 2.867200618931e-07L, 3.930665188016e-07L, -2.743489161464e-07L, 9.444378836274e-08L, 1.683073984106e-07L, -3.762042719394e-07L, 2.461104641394e-07L, -2.777099373105e-07L, 2.165258793652e-07L, 0.0L, },
                {0.0L, 5.298536176321e-07L, -8.664146374574e-08L, 2.473922278404e-07L, -2.754135132048e-07L, -1.04517037428e-07L, -1.208724192291e-07L, -2.991605003572e-07L, 1.158243346084e-07L, -6.103929018276e-08L, 1.66052130446e-08L, 1.937911735148e-08L, 1.169914335933e-07L, 4.872774302249e-08L, -1.190877551044e-08L, -1.374181798308e-07L, -1.432257892093e-07L, 4.485793715105e-08L, -1.669399018713e-07L, 6.314033255248e-08L, -3.480779173925e-08L, 3.909780803098e-07L, -1.234984419028e-07L, 1.394926425286e-08L, 1.66388578795e-07L, 2.787455281081e-07L, -2.858069371344e-07L, -6.909759076387e-08L, -1.975642621008e-07L, -1.933450919083e-07L, -6.680511878855e-08L, -2.680955816735e-07L, -3.693410553192e-07L, -5.709262542542e-07L, 2.145555827746e-07L, -2.139824580743e-07L, -9.73178102322e-08L, -1.617642921496e-07L, 8.322883168651e-08L, 1.060281133896e-07L, -7.717580348722e-09L, },
            };
            std::vector<std::vector<real>> moonS = {
                {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.044474062376e-09L, -1.339776662195e-10L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 5.454539983454e-06L, 4.877974232126e-06L, -1.774317683429e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.664316234498e-06L, -6.777045803795e-06L, -1.344499123017e-05L, 3.926536974166e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -4.118918173072e-06L, 1.057141493807e-06L, 8.698899369844e-06L, 6.762259094131e-08L, -2.754537967015e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -2.59959663689e-06L, -2.16769838031e-06L, -3.427429908688e-06L, -4.058046543825e-06L, -1.034177169942e-05L, 7.229882841312e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -1.197358462935e-07L, 2.411111240194e-06L, 2.357332156782e-06L, 7.565103968368e-07L, 1.069310900234e-06L, 1.100462730689e-06L, -1.600061209739e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.098037428797e-06L, 1.930565327482e-06L, 9.54467857119e-07L, -5.282456437967e-07L, 2.918558502819e-06L, -2.114684853858e-06L, 3.268875515845e-06L, 2.116401754605e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 8.104011483527e-08L, -1.387556589354e-06L, 2.201755361894e-06L, -1.425802105006e-06L, -3.524681297892e-06L, -3.002638828765e-06L, -1.068314188477e-07L, -2.203396882848e-06L, 2.488182623272e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -9.540706794646e-07L, -2.651140380086e-07L, 6.688288283126e-07L, 1.578917095792e-06L, -3.145818444791e-07L, -2.095301146882e-06L, -9.107479035398e-07L, 2.848852925879e-06L, -5.15709694709e-08L, -1.719363865228e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 6.142459427231e-07L, -1.968039101834e-06L, 7.3626740905e-07L, 2.190865808919e-06L, 2.588022664298e-06L, -1.367017856084e-06L, -2.841080837348e-06L, -3.829208951302e-07L, 2.526419021932e-07L, -1.735740282035e-07L, -1.746342814782e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.647787714916e-06L, 1.507815517001e-06L, -4.05895618219e-06L, -3.446065746044e-07L, 8.936239514436e-07L, 5.452460917158e-07L, -9.67550678751e-07L, -2.20654449047e-06L, 1.094999000556e-06L, -1.273937739806e-06L, -3.284553467758e-07L, 1.246948397101e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -7.666473040453e-08L, 3.839634842546e-07L, -1.456309283997e-06L, -5.819236782287e-07L, 1.021336387082e-06L, 1.384223888561e-06L, -6.726606407854e-07L, -1.210513854683e-06L, 5.274962761747e-07L, 4.528921005315e-07L, 4.390632645253e-07L, -1.522336427881e-06L, -2.668728525327e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -6.821045531959e-07L, 3.552093494384e-07L, 1.679963027112e-07L, -2.582060879249e-06L, -2.160029285717e-07L, 1.862282308997e-06L, 1.161513049457e-06L, 4.414611335583e-07L, -1.6380502802e-06L, 5.930174421348e-07L, 1.937210949587e-06L, -5.394799988645e-07L, -1.017786062391e-06L, 1.071082146428e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 3.536182255415e-07L, 7.44655359357e-07L, -6.307094540459e-07L, -1.323373444078e-06L, -8.933914895019e-07L, -3.306734482261e-07L, 1.489417555027e-06L, 1.88656871524e-06L, -1.011633317443e-06L, 4.227357479522e-07L, -1.502438889732e-06L, -4.717880217535e-07L, 3.095282915492e-07L, 5.614579324705e-07L, 7.778728243735e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 8.252521364808e-07L, 7.654915689896e-09L, 5.233336049571e-07L, 3.12441562348e-07L, -6.484841123361e-07L, -7.28070250666e-07L, -5.803978614532e-07L, 3.200903605372e-07L, -8.233740097447e-07L, 3.355487191747e-08L, 7.721814417669e-07L, 9.158381607381e-07L, -1.629951910605e-07L, 1.776100258325e-07L, -1.015801567309e-06L, -9.187880805939e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 7.625806000147e-08L, -8.725823953726e-07L, -5.022504201313e-07L, 2.037548108394e-06L, -6.042884581559e-07L, -1.245859787129e-06L, -1.556268740031e-06L, -4.704433232921e-07L, 1.371639609494e-06L, 6.137159977419e-07L, -2.019999915838e-07L, 1.795595357475e-07L, -2.222894053828e-07L, 3.045656127311e-07L, -5.045794793982e-07L, -2.515414449896e-07L, 1.504669144869e-06L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 2.683984901673e-08L, -7.261171114228e-07L, 4.326731275642e-07L, 1.091313292446e-06L, 7.318700860829e-08L, 7.338987413744e-07L, -1.133723378003e-06L, 2.868069850685e-07L, -6.374033796797e-07L, -5.501063826309e-07L, -1.297980427538e-06L, 2.3476604878e-07L, -1.001865063468e-07L, -1.132048001672e-06L, -1.084448923955e-06L, -2.13580890541e-07L, -4.71239003349e-07L, -2.851264159973e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.313156308342e-07L, -3.717786212158e-07L, 9.766063752082e-07L, -1.053904338461e-06L, 1.483264507519e-06L, 3.800645165356e-07L, 4.9115804117e-07L, 8.666830547178e-07L, -4.222245276742e-07L, -2.259927200653e-07L, -1.752645937746e-07L, -2.871721695615e-07L, 2.231847216298e-06L, 5.39892526904e-07L, 8.251565487773e-08L, 6.894550448986e-07L, -1.868706390978e-06L, -2.166372909266e-08L, 3.653256650455e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -3.825714830644e-07L, 3.41229266847e-09L, -2.65491690579e-07L, -4.888520402021e-07L, -3.031436186534e-07L, 3.299406194466e-07L, -7.31491013526e-07L, 6.92823732073e-08L, 1.098718563014e-07L, 4.890299780408e-08L, 1.129429187955e-06L, -9.832406165536e-07L, 6.176687475653e-07L, 8.283482102639e-08L, 9.198513136095e-08L, -6.115845329403e-07L, -1.185035249911e-06L, -3.418036284974e-08L, 3.051553162035e-07L, -3.587712729262e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 5.530573059903e-07L, -3.909911550754e-07L, -9.283116489533e-07L, 7.418801815701e-07L, 9.014896504071e-08L, -4.021903722654e-07L, 6.192316402257e-07L, -3.254440883349e-08L, -5.41629120388e-07L, -5.734850203776e-08L, -1.730793967632e-07L, 4.039477516997e-07L, -1.76884462041e-07L, 3.486754910762e-07L, 1.926492958209e-07L, -2.919836294116e-07L, 1.46678343981e-06L, -6.064004749891e-07L, 1.489300742765e-07L, 2.544076698055e-06L, 3.240729316779e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 2.514713825708e-07L, -9.016490558287e-08L, 6.582280522973e-07L, 7.156229651902e-08L, -7.716438589418e-07L, 1.839210594086e-07L, 6.64028465675e-07L, -4.071334091255e-07L, 1.188605286958e-07L, 6.046400565132e-07L, -7.102555792955e-07L, 4.021022695535e-08L, 2.015747707424e-07L, 1.18901155243e-06L, -1.768907341255e-07L, 4.793566788058e-07L, -7.659208159485e-08L, -2.666970716475e-07L, 3.254091989824e-07L, -5.140453051984e-07L, 7.622999223576e-08L, 6.564732766735e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -8.344223573882e-07L, -8.720162772709e-08L, 3.190707733568e-07L, -6.001262482625e-07L, 7.308599518535e-07L, 9.624794445702e-08L, -7.767477862986e-07L, -3.505148568725e-07L, 7.236582539038e-08L, 1.724835013205e-07L, 1.583208588105e-07L, -4.733169296864e-07L, 2.260816077464e-07L, -4.321216480685e-07L, 6.543988167744e-08L, 2.636038652699e-07L, -3.834567094311e-07L, 3.045899015026e-07L, -6.097388529209e-08L, 3.810293475185e-07L, 3.041881945015e-07L, -4.889028658922e-07L, 2.494189772105e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -2.729034532285e-07L, 3.630446013164e-07L, -9.526698803286e-07L, 1.966601169573e-07L, -1.01380129791e-07L, 2.973255869331e-08L, 1.465890204928e-07L, 2.126609697074e-07L, -7.016696185965e-08L, 1.348424701022e-07L, -3.809187771401e-07L, -4.163232534338e-08L, -3.834689542754e-07L, -6.480080800445e-07L, 3.807235058171e-07L, -2.137899149198e-07L, 2.515893767991e-07L, -6.533770332598e-07L, -6.31344872182e-07L, -1.68760940854e-07L, 1.57735219217e-07L, 7.239025922275e-07L, -3.568917116157e-07L, -9.445272356012e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 5.291851584188e-07L, 5.391468770915e-08L, -1.457329715102e-07L, -7.099666970322e-08L, 8.101712812197e-07L, 3.345087345497e-07L, 9.665545512157e-07L, 9.957575835315e-07L, 9.835031183462e-08L, 7.126165668393e-07L, -3.992671286612e-07L, 4.431774644267e-09L, 2.737872510877e-07L, 6.848378037851e-07L, 3.708045396001e-08L, -5.942461613879e-07L, -4.867898874374e-07L, 1.324585936755e-07L, -1.334047366906e-07L, -4.007682012942e-07L, -2.907680982473e-07L, -7.880917504188e-07L, -3.588007478087e-07L, -3.066367132013e-08L, 5.601500097409e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 2.222336870325e-07L, -1.962622913373e-07L, 6.64473067708e-07L, 2.447331386851e-07L, -7.765566444202e-07L, -1.350487538877e-07L, -4.393223211152e-07L, 1.521586415309e-07L, 7.009025394695e-07L, -1.092377130461e-08L, 3.232900189515e-07L, -7.900953431015e-07L, 4.266593359953e-07L, 8.898778789525e-08L, -2.542898061917e-07L, 4.0056264357e-07L, 8.402912110609e-08L, 2.997083091876e-07L, 2.114133399463e-07L, 6.826550930673e-07L, -5.377462676796e-07L, 8.594766412428e-07L, -6.273303618222e-08L, -1.171161887714e-08L, -3.400624243204e-08L, 4.780189217452e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 2.29285392446e-08L, -4.292190020953e-07L, 8.740218337209e-08L, 3.213599371303e-07L, 1.391027811042e-08L, 4.631886574389e-07L, -3.585575379525e-07L, -7.587956431698e-08L, 1.002209325107e-07L, -6.36132221614e-07L, -4.384866238147e-07L, -5.205172066273e-07L, -1.492465718418e-07L, -7.281935606007e-07L, 4.102812185938e-07L, 3.910584018969e-07L, 6.809294239008e-07L, -6.787968462406e-07L, -1.923158005742e-07L, -5.21870581523e-07L, -4.587511698186e-07L, 8.75628015862e-07L, -6.418419491224e-07L, -1.64207822735e-07L, 5.77789722028e-07L, -7.206039392975e-07L, -2.770623564877e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 7.93918175382e-07L, 7.104729910688e-07L, -8.397260216268e-08L, 6.436936827206e-08L, 4.895017819616e-08L, 3.445603155883e-08L, 2.630327137446e-07L, 3.033023452887e-08L, -9.158822011114e-08L, -9.021925816672e-07L, -4.182500433986e-08L, -5.158638405244e-08L, -3.461137935007e-07L, -2.946964330896e-07L, 3.42954848577e-07L, -8.449626742855e-07L, -1.887915078794e-07L, 6.46138901191e-07L, -3.7486577711e-08L, -5.683766632748e-07L, 7.202595910908e-08L, -9.252631544963e-08L, 2.450128816412e-07L, 7.536665127763e-09L, -5.169168196915e-07L, 9.087008001661e-07L, 6.720722284915e-07L, -4.42115114823e-09L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -4.383562931228e-07L, -7.004605889154e-07L, 2.659397958359e-07L, 1.473540897279e-07L, 1.690101167548e-08L, -1.88112219107e-07L, -3.852915599005e-07L, -2.266408376259e-07L, -6.360342177816e-07L, 3.447057105217e-07L, 4.354528259873e-07L, -4.918695901789e-07L, 3.842989807918e-07L, 8.077061811645e-07L, 6.885060477046e-08L, -8.701199765469e-08L, -3.072597394555e-07L, -5.157505606126e-09L, -7.61005571234e-08L, 3.850018487241e-07L, 2.933433726825e-07L, -4.748970116642e-07L, 3.299120415832e-07L, 3.719812319233e-08L, 2.249319547257e-07L, -3.946803752199e-07L, 1.912208336303e-07L, 2.493256361544e-07L, -1.46040329807e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -4.29125617833e-09L, -2.322808564937e-07L, 3.514294674176e-07L, -2.222791613673e-07L, -3.126401031567e-07L, 3.477860718182e-07L, -3.079964896133e-07L, -3.359042411346e-07L, 5.051666141115e-08L, -1.008749529556e-07L, -2.348209554301e-07L, 2.386419134008e-07L, 3.914484680254e-07L, -2.057105425827e-07L, 1.659496122874e-07L, -2.693493169004e-07L, 3.503390241802e-07L, 3.205336249258e-08L, 5.729579188325e-08L, -2.556248870676e-08L, -1.951987992034e-08L, 3.929165698255e-07L, 1.438976874542e-07L, 4.530258630329e-08L, 2.46114343289e-07L, -3.70543662118e-07L, 3.390712764326e-07L, -6.169911123294e-07L, -1.404467228199e-07L, 9.582805242295e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -4.194911784593e-07L, 1.260974818464e-07L, -5.413168681722e-07L, 2.035392967677e-07L, -1.781932995919e-07L, 1.042441630579e-07L, -2.188262974551e-07L, -2.652269841736e-07L, -1.281660268794e-07L, -2.381376768425e-07L, 5.000901451383e-07L, 2.977322004162e-07L, 3.608535599251e-07L, 1.904806808337e-08L, -1.364191230043e-08L, -3.000068520955e-07L, -1.614323141714e-08L, -3.653972347049e-07L, -5.144311466918e-08L, -1.572705703698e-07L, -9.611297734068e-08L, -1.620955449039e-08L, 6.208412270959e-07L, -1.032532540641e-07L, -1.392406410171e-07L, -4.234271377724e-07L, 2.454576319149e-07L, -2.887007074567e-07L, 2.540631082712e-07L, 1.560869445876e-07L, -1.928240513299e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 2.126467537629e-07L, 7.126675863991e-08L, -2.514259553883e-07L, -1.092462396596e-07L, -5.158320212966e-07L, -2.757490695386e-07L, -5.347380353121e-07L, 6.719830544143e-08L, 4.224079980322e-07L, 1.988150044197e-07L, -1.334580442632e-07L, 4.232516912204e-08L, -1.910945348649e-07L, -2.994802064769e-07L, 2.813844077705e-07L, 3.603516922205e-08L, -1.505556827427e-07L, 8.568391803426e-08L, 2.623102696602e-07L, -4.29509021138e-07L, 1.610733235787e-07L, 7.175601227161e-08L, -1.520761757626e-07L, 1.108420840711e-07L, 2.467523432591e-07L, 6.017903929537e-08L, 4.042425714452e-07L, -1.732713454176e-07L, 4.947047703168e-07L, -4.12296811588e-07L, -2.656228490679e-07L, -5.422820784761e-08L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -5.460656018599e-08L, -3.327451150041e-07L, 5.823803295664e-07L, 1.533638526264e-07L, 1.134934361193e-07L, -1.252766200641e-07L, 9.719689958218e-08L, 3.515645969136e-07L, 2.68869409436e-07L, 1.349179872582e-07L, 1.851030173599e-07L, -1.020764006673e-08L, -2.881505221128e-07L, 1.484266235999e-07L, 6.78189216613e-07L, 3.059281062431e-07L, -3.91558561736e-07L, -3.284010144349e-07L, 2.569473460313e-07L, -1.355673113805e-07L, -2.238144500253e-09L, -2.631179391852e-07L, 1.070067160548e-07L, 2.252685432701e-07L, 3.521040817802e-07L, 3.51470670728e-07L, -9.351103960262e-08L, -3.955517005484e-07L, -1.435758995715e-07L, 3.395887219015e-07L, 2.87872927185e-07L, 3.822172663903e-07L, 6.856304792339e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, -2.39664277919e-07L, -2.464657645197e-08L, -1.093059012276e-07L, -1.90180833553e-07L, 1.196074057703e-07L, 1.612081086426e-07L, 3.159933538399e-07L, 2.001879090814e-07L, -2.046941397447e-07L, -5.764423322558e-07L, -3.591247909357e-07L, 4.154565887978e-07L, 1.879588351651e-07L, -1.210992630495e-07L, 2.660211923449e-08L, -3.758636989957e-07L, -1.915441994942e-07L, 1.471213448277e-07L, 2.107882098529e-07L, -1.475229060039e-07L, -1.311081455416e-07L, 4.377468907859e-08L, -1.373105901331e-07L, -1.632045355309e-07L, 3.311864406837e-07L, 3.501291767722e-08L, -1.541191321421e-07L, 3.477268046638e-07L, -6.094202077302e-08L, -2.420271032126e-09L, 5.667232113088e-08L, 6.824508815334e-08L, 4.176239383879e-07L, -4.377353586955e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 4.722453814787e-07L, 2.135743194085e-07L, 4.76210702342e-07L, -5.835612979894e-07L, 2.023010729576e-07L, -4.598119886563e-10L, -4.50693934552e-08L, 4.525489478466e-07L, 4.311356850144e-08L, -1.775781085813e-07L, 1.56553154736e-07L, 5.631467884047e-08L, -1.501991716605e-07L, -2.197263120772e-07L, -2.931122743692e-07L, 4.10489044412e-08L, -2.566690488652e-07L, 2.674250364416e-07L, 1.316079356725e-07L, 1.614644136051e-07L, 1.034102727846e-07L, -5.428304919648e-07L, -3.800361851438e-07L, -1.65495544583e-07L, 7.261309432491e-10L, -5.673026784372e-07L, -6.820138128797e-08L, 2.009786491523e-07L, -8.618195130827e-08L, 4.74981941912e-07L, 3.156639610424e-07L, 2.151486449859e-07L, -1.825718534019e-07L, -1.478488844463e-07L, 2.34678714783e-07L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 4.891874964385e-07L, -1.250287883638e-07L, -3.72118694153e-07L, -1.308091915335e-07L, 7.94828655103e-08L, 1.748523072116e-07L, 1.44215745105e-07L, 2.914110822072e-07L, -2.609798555396e-07L, 1.11790328765e-07L, 5.338423585597e-08L, 1.617589231834e-07L, -7.613129901404e-08L, -2.262898579315e-07L, -1.293800613211e-07L, 5.457588873549e-07L, -8.877908917685e-08L, 1.494824339083e-07L, 3.612959333098e-07L, -1.742492031228e-07L, -8.835968772083e-08L, 1.771974705077e-07L, 3.556563228484e-07L, 9.768877213234e-08L, 3.413947540595e-07L, 5.010338209208e-08L, 5.411405068125e-07L, 3.116157742247e-07L, 1.494368725133e-07L, 1.584728624616e-07L, -3.523767646797e-07L, -2.54697717473e-07L, -3.681367699627e-07L, 1.870041521666e-07L, 4.571506836222e-07L, 2.22396888922e-07L, 0.0L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.053990915421e-07L, -3.779754000684e-07L, 3.089027970266e-07L, 5.109830084639e-08L, 4.585346751007e-07L, -1.609089365007e-07L, -9.969212378631e-08L, -2.51736882417e-07L, -3.347052129209e-07L, -8.994535246947e-08L, -1.326591181253e-07L, -1.88105821136e-07L, 1.984141287611e-07L, -4.076621873926e-07L, 1.633090271434e-07L, 3.395977387943e-07L, -5.714667960032e-07L, 4.16677683888e-07L, -2.745587452498e-07L, -1.144751449903e-08L, -4.259026022179e-08L, -5.24311907523e-07L, 4.711661725878e-07L, 1.125926348188e-08L, 6.04156828375e-08L, -1.609516524757e-07L, -3.625125620208e-07L, -4.795688558048e-07L, 6.814931971859e-08L, -1.523247980856e-07L, 1.954250438277e-07L, 4.124829572539e-07L, 2.276760397848e-07L, 2.918889218245e-07L, 1.944303661509e-07L, -1.883294263995e-07L, -1.14011041218e-07L, 0.0L, 0.0L, 0.0L, },
                {0.0L, 1.351250452538e-07L, 1.956626082291e-08L, 3.099851094477e-07L, 2.323458610498e-07L, -1.000683435823e-07L, 2.472791425524e-07L, -4.857309065529e-08L, -4.577762320372e-08L, -9.083087186688e-08L, 4.010320993992e-08L, -5.345160500579e-08L, 2.141885769205e-07L, -3.283039318928e-07L, 1.448829318667e-08L, 6.511456840293e-07L, -6.670909042756e-08L, -3.028310762569e-07L, -6.850176962823e-08L, 1.514396174737e-07L, 1.545144174681e-07L, -2.71159265009e-07L, 2.620487586347e-07L, 3.162515312535e-07L, -3.699523841614e-07L, 2.295421502024e-08L, -3.534468576474e-08L, -5.96971349401e-08L, 1.918165057874e-07L, 6.911936186668e-08L, -8.617876623366e-09L, 2.906034817143e-07L, 7.759793234428e-08L, -8.494146464567e-08L, -3.026649722361e-07L, -1.890901393434e-07L, 2.238488062264e-07L, 3.153481315397e-08L, -3.035140583812e-08L, 0.0L, 0.0L, },
                {0.0L, -1.532975761088e-07L, -4.83758085411e-08L, -3.035549753884e-08L, -8.435719162789e-08L, -2.26400624678e-07L, -2.248804716919e-07L, 1.193613709939e-07L, -3.147493082332e-08L, -1.390612188784e-07L, 4.506531638294e-08L, 2.043486875048e-07L, 1.62788124182e-07L, 1.609718750729e-07L, 2.446930909348e-07L, -4.25670054758e-07L, -4.600761182478e-07L, -1.932549197265e-07L, 2.494094979886e-07L, -2.369909867542e-07L, -7.936060114598e-08L, 3.343589721657e-08L, 2.824811329321e-07L, 5.491443173647e-10L, -3.002526287364e-07L, 2.157850858548e-07L, -2.189556816239e-07L, -8.474195391328e-08L, 3.635736209889e-07L, -6.586735635505e-09L, -2.934978816695e-07L, 1.920283324574e-07L, -3.135361793815e-09L, -4.259444918322e-09L, -1.813275258487e-07L, 5.339876823957e-07L, 4.775466776249e-07L, -2.76980444284e-07L, 1.385227497451e-07L, -1.02829784935e-07L, 0.0L, },
                {0.0L, 1.077057702761e-07L, 2.315010980541e-07L, 2.165906896888e-07L, -7.998331784095e-09L, -1.910378670181e-07L, 2.06616620689e-07L, -5.475178330136e-08L, 4.599127945747e-09L, -1.655513302748e-07L, 2.367528240854e-08L, -1.919174169068e-07L, -6.138335412944e-08L, 1.101287205646e-07L, 1.442220062286e-07L, -4.236549831124e-07L, -2.13021116611e-07L, 1.887453849318e-07L, 2.291603112464e-07L, 2.498209560982e-07L, -1.888852504442e-07L, -4.050865824876e-08L, -8.905967363098e-08L, -6.611225698635e-08L, 5.313821405955e-08L, 6.972630382721e-08L, -3.670433583564e-08L, -7.597272083162e-08L, 2.548810577039e-07L, -4.596247685622e-08L, -1.688063163977e-07L, -1.123332382724e-07L, -2.167757509774e-07L, -2.685218825747e-07L, 9.899936915079e-08L, 1.663428860602e-07L, -1.411829390477e-07L, 8.761604542582e-08L, 1.60827946814e-07L, 2.571010304203e-07L, 1.033111214101e-08L, },
            };
            Moon.set_harmonics(269.9949L, 66.5392L, 
                40, 40, moonJ, moonC, moonS);
            Sun.caTol = 0.25;
            MercuryBarycenter.caTol = 0.1;
            VenusBarycenter.caTol = 0.1;
            Earth.caTol = 0.1;
            Moon.caTol = 0.05;
            MarsBarycenter.caTol = 0.1;
            JupiterBarycenter.caTol = 0.25;
            SaturnBarycenter.caTol = 0.25;
            UranusBarycenter.caTol = 0.25;
            NeptuneBarycenter.caTol = 0.25;
            PlutoBarycenter.caTol = 0.1;
            add_spice_body(Sun);
            add_spice_body(MercuryBarycenter);
            add_spice_body(VenusBarycenter);
            add_spice_body(Earth);
            add_spice_body(Moon);
            add_spice_body(MarsBarycenter);
            add_spice_body(JupiterBarycenter);
            add_spice_body(SaturnBarycenter);
            add_spice_body(UranusBarycenter);
            add_spice_body(NeptuneBarycenter);
            add_spice_body(PlutoBarycenter);

            // add DE441 big16 asteroids from JPL SSD IOM 392R-21-005
            // (ftp://ssd.jpl.nasa.gov/pub/eph/small_bodies/asteroids_de441/SB441_IOM392R-21-005_perturbers.pdf)
            SpiceBody Ceres("Ceres", 2000001, this->integParams.t0,
                            1.3964518123081070e-13L / G, 0.0L);
            SpiceBody Vesta("Vesta", 2000004, this->integParams.t0,
                            3.8548000225257904e-14L / G, 0.0L);
            SpiceBody Pallas("Pallas", 2000002, this->integParams.t0,
                             3.0471146330043200e-14L / G, 0.0L);
            SpiceBody Hygiea("Hygiea", 2000010, this->integParams.t0,
                             1.2542530761640810e-14L / G, 0.0L);
            SpiceBody Davida("Davida", 2000511, this->integParams.t0,
                             8.6836253492286545e-15L / G, 0.0L);
            SpiceBody Interamnia("Interamnia", 2000704, this->integParams.t0,
                                 6.3110343420878887e-15L / G, 0.0L);
            SpiceBody Europa("Europa", 2000052, this->integParams.t0,
                             5.9824315264869841e-15L / G, 0.0L);
            SpiceBody Sylvia("Sylvia", 2000087, this->integParams.t0,
                             4.8345606546105521e-15L / G, 0.0L);
            SpiceBody Eunomia("Eunomia", 2000015, this->integParams.t0,
                              4.5107799051436795e-15L / G, 0.0L);
            SpiceBody Juno("Juno", 2000003, this->integParams.t0,
                           4.2823439677995011e-15L / G, 0.0L);
            SpiceBody Psyche("Psyche", 2000016, this->integParams.t0,
                             3.5445002842488978e-15L / G, 0.0L);
            SpiceBody Camilla("Camilla", 2000107, this->integParams.t0,
                              3.2191392075878588e-15L / G, 0.0L);
            SpiceBody Thisbe("Thisbe", 2000088, this->integParams.t0,
                             2.6529436610356353e-15L / G, 0.0L);
            SpiceBody Iris("Iris", 2000007, this->integParams.t0,
                           2.5416014973471498e-15L / G, 0.0L);
            SpiceBody Euphrosyne("Euphrosyne", 2000031, this->integParams.t0,
                                 2.4067012218937576e-15L / G, 0.0L);
            SpiceBody Cybele("Cybele", 2000065, this->integParams.t0,
                             2.0917175955133682e-15L / G, 0.0L);
            Ceres.caTol = 0.05;
            Vesta.caTol = 0.05;
            Pallas.caTol = 0.05;
            add_spice_body(Ceres);
            add_spice_body(Vesta);
            add_spice_body(Pallas);
            add_spice_body(Hygiea);
            add_spice_body(Davida);
            add_spice_body(Interamnia);
            add_spice_body(Europa);
            add_spice_body(Sylvia);
            add_spice_body(Eunomia);
            add_spice_body(Juno);
            add_spice_body(Psyche);
            add_spice_body(Camilla);
            add_spice_body(Thisbe);
            add_spice_body(Iris);
            add_spice_body(Euphrosyne);
            add_spice_body(Cybele);
            break;
        }
        default:
            throw std::invalid_argument(
                "The defaultSpiceBodies argument is only defined for no "
                "preloaded SPICE bodies (case 0) or DE431 (case 431) or DE441 "
                "(case 441).");
            break;
    }
    this->spkEphem.mbPath = kernel_mb;
    this->spkEphem.sbPath = kernel_sb;
    this->pckEphem.histPckPath = DEkernelPath + "earth_historic.bpc";
    this->pckEphem.latestPckPath = DEkernelPath + "earth_latest.bpc";
    this->pckEphem.predictPckPath = DEkernelPath + "earth_predict.bpc";
}

/**
 * @param[in] name Name of the simulation.
 * @param[in] simRef Reference simulation to copy.
 */
PropSimulation::PropSimulation(std::string name, const PropSimulation& simRef) {
    this->name = name;
    this->DEkernelPath = simRef.DEkernelPath;
    this->spkEphem.mbPath = simRef.spkEphem.mbPath;
    this->spkEphem.sbPath = simRef.spkEphem.sbPath;
    this->spkEphem.mb = nullptr;
    this->spkEphem.sb = nullptr;
    this->pckEphem.histPckPath = simRef.pckEphem.histPckPath;
    this->pckEphem.latestPckPath = simRef.pckEphem.latestPckPath;
    this->pckEphem.predictPckPath = simRef.pckEphem.predictPckPath;
    this->pckEphem.histPck = nullptr;
    this->pckEphem.latestPck = nullptr;
    this->pckEphem.predictPck = nullptr;
    this->consts = simRef.consts;
    this->integParams = simRef.integParams;
    this->integParams.nInteg = 0;
    this->integParams.n2Derivs = 0;
    this->integParams.nTotal = simRef.integParams.nSpice;
    this->integParams.timestepCounter = 0;
    this->spiceBodies = simRef.spiceBodies;
    this->tEvalUTC = simRef.tEvalUTC;
    this->evalApparentState = simRef.evalApparentState;
    this->evalMeasurements = simRef.evalMeasurements;
    this->convergedLightTime = simRef.convergedLightTime;
    this->observerInfo = simRef.observerInfo;
    this->xObserver = simRef.xObserver;
    this->tEvalMargin = simRef.tEvalMargin;
    this->tEval = simRef.tEval;
    this->obsType = simRef.obsType;
    this->isPreprocessed = false;
}

/**
 * @param[in] tEval Vector of times at which to evaluate the integrated state.
 * @param[in] observerInfo Observer information array.
 */
void PropSimulation::prepare_for_evaluation(
    std::vector<real>& tEval, std::vector<std::vector<real>>& observerInfo) {
    const bool forwardProp = this->integParams.t0 <= this->integParams.tf;
    const bool backwardProp = this->integParams.t0 >= this->integParams.tf;
    if (forwardProp && backwardProp) {
        throw std::invalid_argument(
            "The initial and final times must be different.");
    }
    // sort tEval into ascending order or descending order based on the
    // integration direction (not during orbit fits)
    if (observerInfo.size() == 0) {
        std::vector<size_t> tEvalSortedIdx(tEval.size());
        sort_vector(tEval, forwardProp, tEvalSortedIdx);
    }
    if (forwardProp) {
        int removeCounter = 0;
        while (tEval[0] < this->integParams.t0 - this->tEvalMargin) {
            // remove any tEval values that are before the integration start
            // time
            tEval.erase(tEval.begin());
            if (observerInfo.size() != 0)
                observerInfo.erase(observerInfo.begin());
            removeCounter++;
        }
        while (tEval.back() > this->integParams.tf + this->tEvalMargin) {
            // remove any tEval values that are after the integration end time
            tEval.pop_back();
            if (observerInfo.size() != 0) observerInfo.pop_back();
            removeCounter++;
        }
        if (removeCounter > 0) {
            std::cout << "WARNING: " << removeCounter
                      << " tEval and observerInfo value(s) were removed "
                         "because they were outside the interpolation range, "
                         "i.e., integration range with a margin of "
                      << this->tEvalMargin << " day(s)." << std::endl;
        }
    } else if (backwardProp) {
        int removeCounter = 0;
        while (tEval[0] > this->integParams.t0 + this->tEvalMargin) {
            // remove any tEval values that are after the integration start time
            tEval.erase(tEval.begin());
            if (observerInfo.size() != 0)
                observerInfo.erase(observerInfo.begin());
            removeCounter++;
        }
        while (tEval.back() < this->integParams.tf - this->tEvalMargin) {
            // remove any tEval values that are before the integration end time
            tEval.pop_back();
            if (observerInfo.size() != 0) observerInfo.pop_back();
            removeCounter++;
        }
        if (removeCounter > 0) {
            std::cout << "WARNING: " << removeCounter
                      << " tEval and observerInfo value(s) were removed "
                         "because they were outside the interpolation range, "
                         "i.e., integration range with a margin of "
                      << this->tEvalMargin << " day(s)." << std::endl;
        }
    }

    if (observerInfo.size() == 0) {
        observerInfo = std::vector<std::vector<real>>(
            tEval.size(), std::vector<real>(4, 0.0L));
    }

    if (this->observerInfo.size() == 0) {
        this->observerInfo = observerInfo;
    } else if (this->observerInfo.size() != 0) {
        for (size_t i = 0; i < observerInfo.size(); i++) {
            this->observerInfo.push_back(observerInfo[i]);
        }
    }

    std::vector<std::vector<real>> xObserver = std::vector<std::vector<real>>(
        tEval.size(), std::vector<real>(6, 0.0L));
    std::vector<int> obsType = std::vector<int>(tEval.size(), 0);
    this->map_ephemeris();
    if (tEval.size() != 0) {
        for (size_t i = 0; i < tEval.size(); i++) {
            if (observerInfo[i].size() == 4 || observerInfo[i].size() == 7) {
                obsType[i] = 0;
                // obsType = 3 is optical gaia measurement, set externally from Python when needed
            } else if (observerInfo[i].size() == 9) {
                obsType[i] = 1;
            } else if (observerInfo[i].size() == 10) {
                obsType[i] = 2;
            } else {
                throw std::invalid_argument(
                    "The observerInfo vector must have 4/7 (optical), 9 (radar "
                    "delay), or 10 elements (radar doppler).");
            }
            get_observer_state(tEval[i], observerInfo[i], this,
                            this->tEvalUTC, xObserver[i]);
        }
    }
    this->unmap_ephemeris();

    if (this->tEval.size() == 0) {
        this->tEval = tEval;
        this->xObserver = xObserver;
        this->obsType = obsType;
    } else if (this->tEval.size() != 0) {
        for (size_t i = 0; i < tEval.size(); i++) {
            this->tEval.push_back(tEval[i]);
            this->xObserver.push_back(xObserver[i]);
            this->obsType.push_back(obsType[i]);
        }
    }
}

void PropSimulation::map_ephemeris(){
    if (this->spkEphem.mb == nullptr){
        this->spkEphem.mb = spk_init(this->spkEphem.mbPath);
    }
    if (this->spkEphem.sb == nullptr){
        this->spkEphem.sb = spk_init(this->spkEphem.sbPath);
    }
    if (this->pckEphem.histPck == nullptr){
        this->pckEphem.histPck = pck_init(this->pckEphem.histPckPath);
    }
    if (this->pckEphem.latestPck == nullptr){
        this->pckEphem.latestPck = pck_init(this->pckEphem.latestPckPath);
    }
    if (this->pckEphem.predictPck == nullptr){
        this->pckEphem.predictPck = pck_init(this->pckEphem.predictPckPath);
    }
}

void PropSimulation::unmap_ephemeris(){
    if (this->spkEphem.mb != nullptr){
        spk_free(this->spkEphem.mb);
        this->spkEphem.mb = nullptr;
    }
    if (this->spkEphem.sb != nullptr){
        spk_free(this->spkEphem.sb);
        this->spkEphem.sb = nullptr;
    }
    if (this->pckEphem.histPck != nullptr){
        pck_free(this->pckEphem.histPck);
        this->pckEphem.histPck = nullptr;
    }
    if (this->pckEphem.latestPck != nullptr){
        pck_free(this->pckEphem.latestPck);
        this->pckEphem.latestPck = nullptr;
    }
    if (this->pckEphem.predictPck != nullptr){
        pck_free(this->pckEphem.predictPck);
        this->pckEphem.predictPck = nullptr;
    }
}

/**
 * @param[in] t Time at which to get the state.
 * @param[in] bodyName Name of the body.
 * @return std::vector<real> State of the body at time t.
 */
std::vector<real> PropSimulation::get_spiceBody_state(const real t, const std::string &bodyName) {
    int spiceId = -1;
    for (size_t i = 0; i < this->spiceBodies.size(); i++){
        if (this->spiceBodies[i].name == bodyName){
            spiceId = this->spiceBodies[i].spiceId;
            break;
        }
    }
    if (spiceId == -1){
        throw std::invalid_argument("SPICE Body with name " + bodyName +
                                        " does not exist in simulation " +
                                        this->name);
    }
    if (this->spkEphem.mb == nullptr || this->spkEphem.sb == nullptr){
        throw std::invalid_argument(
            "get_spiceBody_state: Ephemeris kernels are not loaded. Memory map "
            "the ephemeris using PropSimulation.map_ephemeris() method first.");
    }
    double spiceState[9];
    get_spk_state(spiceId, t, this->spkEphem, spiceState);
    std::vector<real> state = {spiceState[0], spiceState[1], spiceState[2],
                               spiceState[3], spiceState[4], spiceState[5]};
    return state;
}

/**
 * @param[in] body SpiceBody object to add to the simulation.
 */
void PropSimulation::add_spice_body(SpiceBody body) {
    // check if body already exists. if so, throw error
    for (size_t i = 0; i < this->spiceBodies.size(); i++) {
        if (this->spiceBodies[i].name == body.name) {
            throw std::invalid_argument("SPICE Body with name " + body.name +
                                        " already exists in simulation " +
                                        this->name);
        }
    }
    body.radius /= this->consts.du2m;
    this->spiceBodies.push_back(body);
    this->integParams.nSpice++;
    this->integParams.nTotal++;
}

/**
 * @param[in] body IntegBody object to add to the simulation.
 */
void PropSimulation::add_integ_body(IntegBody body) {
    // check if body already exists. if so, throw error
    for (size_t i = 0; i < this->integBodies.size(); i++) {
        if (this->integBodies[i].name == body.name) {
            throw std::invalid_argument(
                "Integration body with name " + body.name +
                " already exists in simulation " + this->name);
        }
    }
    if (body.t0 != this->integParams.t0) {
        throw std::invalid_argument(
            "Integration body " + body.name + " has initial time MJD " +
            std::to_string(body.t0) +
            " TDB which is different from the simulation initial time: MJD " +
            std::to_string(this->integParams.t0) + " TDB.");
    }
    body.radius /= this->consts.du2m;
    this->integBodies.push_back(body);
    this->integParams.nInteg++;
    this->integParams.nTotal++;
    this->integParams.n2Derivs += body.n2Derivs;
}

/**
 * @param[in] name Name of the body to remove.
 */
void PropSimulation::remove_body(std::string name) {
    for (size_t i = 0; i < this->spiceBodies.size(); i++) {
        if (this->spiceBodies[i].name == name) {
            this->spiceBodies.erase(this->spiceBodies.begin() + i);
            this->integParams.nSpice--;
            this->integParams.nTotal--;
            return;
        }
    }
    for (size_t i = 0; i < this->integBodies.size(); i++) {
        if (this->integBodies[i].name == name) {
            this->integBodies.erase(this->integBodies.begin() + i);
            this->integParams.nInteg--;
            this->integParams.nTotal--;
            return;
        }
    }
    std::cout << "Error: Body " << name << " not found." << std::endl;
}

static size_t event_preprocess(PropSimulation *propSim, const std::string &eventBodyName,
                        const real &tEvent) {
    // check if tEvent is valid
    const bool forwardProp = propSim->integParams.tf > propSim->integParams.t0;
    const bool backwardProp = propSim->integParams.tf < propSim->integParams.t0;
    if ((forwardProp &&
         (tEvent < propSim->integParams.t0 || tEvent >= propSim->integParams.tf)) ||
        (backwardProp &&
         (tEvent > propSim->integParams.t0 || tEvent <= propSim->integParams.tf))) {
        throw std::invalid_argument("Event time " + std::to_string(tEvent) +
                                    " is not within simulation time bounds.");
    }
    // check if body exists
    bool bodyExists = false;
    size_t bodyIndex;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        if (propSim->integBodies[i].name == eventBodyName) {
            bodyExists = true;
            bodyIndex = i;
            break;
        }
    }
    if (!bodyExists) {
        throw std::invalid_argument("Integration body with name " + eventBodyName +
                                    " does not exist in simulation " +
                                    propSim->name);
    }
    return bodyIndex;
}

static void event_stm_handling(PropSimulation *propSim, const Event &event){
    int numEventParams = -1;
    if (!event.isContinuous && !event.deltaVEst && event.multiplierEst) {
        // only multiplier is estimated for an impulsive delta-V
        numEventParams = 1;
    } else if (!event.isContinuous && event.deltaVEst && !event.multiplierEst) {
        // only delta-V vector is estimated for an impulsive delta-V
        numEventParams = 3;
    } else if (event.isContinuous && !event.expAccel0Est && event.tauEst) {
        // only time constant is estimated for a continuous event
        numEventParams = 1;
    } else if (event.isContinuous && event.expAccel0Est && !event.tauEst) {
        // only exponential initial acceleration is estimated for a continuous event
        numEventParams = 3;
    } else if (event.isContinuous && event.expAccel0Est && event.tauEst) {
        // exponential initial acceleration and time constant are estimated for a continuous event
        numEventParams = 4;
    } else {
        throw std::invalid_argument(
            "event_stm_handling: Invalid event estimation configuration for analytic STM "
            "propagation.");
    }
    propSim->integBodies[event.bodyIndex].n2Derivs += 3*numEventParams;
    propSim->integParams.n2Derivs += 3*numEventParams;

    std::vector<real> extraVec = std::vector<real>(6*numEventParams, 0.0);
    // add extra vector at the end of stm vector
    for (size_t i = 0; i < extraVec.size(); i++) {
        propSim->integBodies[event.bodyIndex].stm.push_back(extraVec[i]);
    }
}

static void event_postprocess(PropSimulation *propSim, Event &event) {
    // modify stm if estimating Delta-V
    if (propSim->integBodies[event.bodyIndex].propStm && event.eventEst) {
        event_stm_handling(propSim, event);
    }
    // assign event xIntegIndex
    event.xIntegIndex = 0;
    for (size_t i = 0; i < event.bodyIndex; i++) {
        event.xIntegIndex += 2*propSim->integBodies[i].n2Derivs;
    }
    // add event to either continuous or impulse events based on isContinuous
    std::vector<Event> *eventsList;
    if (event.isContinuous) {
        eventsList = &propSim->eventMngr.continuousEvents;
    } else {
        eventsList = &propSim->eventMngr.impulsiveEvents;
    }
    if (eventsList->size() == 0) {
        eventsList->push_back(event);
    } else {
        for (size_t i = 0; i < eventsList->size(); i++) {
            if (event.t < eventsList->at(i).t) {
                eventsList->insert(eventsList->begin() + i, event);
                break;
            } else if (i == eventsList->size() - 1) {
                eventsList->push_back(event);
                break;
            }
        }
    }
}

/**
 * @param[in] event Event object to add to the simulation.
 */
void PropSimulation::add_event(Event event) {
    size_t bodyIndex = event_preprocess(this, event.bodyName, event.t);
    event.bodyIndex = bodyIndex;
    event.hasStarted = false;
    if (event.isContinuous) {
        // make sure exponential accelerations are not nan and time constant is postive
        if (!std::isfinite(event.expAccel0[0]) || !std::isfinite(event.expAccel0[1]) ||
            !std::isfinite(event.expAccel0[2]) || !std::isfinite(event.tau) || event.tau <= 0.0) {
            throw std::invalid_argument(
                "add_event: Exponential acceleration must be finite and time "
                "constant must be finite and positive.");
        }
        // make sure deltaV is nan and multiplier is nan
        if (std::isfinite(event.deltaV[0]) || std::isfinite(event.deltaV[1]) ||
            std::isfinite(event.deltaV[2]) || std::isfinite(event.multiplier)) {
            throw std::invalid_argument("add_event: Delta-V and multiplier must be nan.");
        }
        this->eventMngr.nConEvents++;
    } else {
        // make sure deltaV is not nan and multiplier is not nan
        if (!std::isfinite(event.deltaV[0]) || !std::isfinite(event.deltaV[1]) ||
            !std::isfinite(event.deltaV[2]) || !std::isfinite(event.multiplier)) {
            throw std::invalid_argument("add_event: Delta-V and multiplier must be finite.");
        }
        // make sure exponential acceleration is nan
        if (std::isfinite(event.expAccel0[0]) || std::isfinite(event.expAccel0[1]) ||
            std::isfinite(event.expAccel0[2]) || std::isfinite(event.tau)) {
            throw std::invalid_argument(
                "add_event: Exponential acceleration must be zero.");
        }
        this->eventMngr.nImpEvents++;
    }
    event_postprocess(this, event);
}

/**
 * @param[in] du2m Distance unit conversion factor (default: AU to meters).
 * @param[in] tu2s Time unit conversion factor (default: days to seconds).
 * @param[in] G Gravitational constant (default: kg AU^3/day^2).
 * @param[in] clight Speed of light (default: AU/day).
  */
void PropSimulation::set_sim_constants(real du2m, real tu2s, real G,
                                       real clight) {
    this->consts.du2m = du2m;
    this->consts.tu2s = tu2s;
    this->consts.duptu2mps = du2m / tu2s;
    this->consts.G = G;
    this->consts.clight = clight;
    this->consts.j2000Jd = 2451545.0;
    this->consts.JdMinusMjd = 2400000.5;
}

/**
 * @param[in] tf Final time.
 * @param[in] tEval Vector of times at which to evaluate the integrated state (default: empty).
 * @param[in] tEvalUTC Flag to indicate if the evaluation times are in UTC (default: false).
 * @param[in] evalApparentState Flag to indicate if the apparent state is evaluated (default: false).
 * @param[in] convergedLightTime Flag to indicate if converged light time needs to be computed (default: false).
 * @param[in] observerInfo Observer information array (default: empty).
 * @param[in] adaptiveTimestep Flag to use adaptive timestep (default: true).
 * @param[in] dt0 Initial timestep (default: 0.0).
 * @param[in] dtMin Minimum timestep (default: 0.005).
 * @param[in] dtChangeFactor Maximum factor by which to change timestep (default: 0.25).
 * @param[in] tolInteg Tolerance for integrator (default: 1e-11).
 * @param[in] tolPC Tolerance for Gauss-Radau predictor-corrector loop (default: 1e-16).
 */
void PropSimulation::set_integration_parameters(
    real tf, std::vector<real> tEval, bool tEvalUTC, bool evalApparentState,
    bool convergedLightTime, std::vector<std::vector<real>> observerInfo,
    bool adaptiveTimestep, real dt0, real dtMin,
    real dtChangeFactor, real tolInteg, real tolPC) {
    this->integParams.tf = tf;
    this->tEvalUTC = tEvalUTC;
    this->evalApparentState = evalApparentState;
    this->convergedLightTime = convergedLightTime;
    if (tEval.size() != 0) {
        prepare_for_evaluation(tEval, observerInfo);
    }
    this->integParams.dt0 = dt0;
    this->integParams.dtMin = dtMin;
    this->integParams.dtChangeFactor = dtChangeFactor;
    this->integParams.adaptiveTimestep = adaptiveTimestep;
    this->integParams.tolPC = tolPC;
    this->integParams.tolInteg = tolInteg;
}

/**
 * @return std::vector<real> Vector of simulation constants.
 */
std::vector<real> PropSimulation::get_sim_constants() {
    std::vector<real> constants = {
        this->consts.du2m,      this->consts.tu2s,   this->consts.duptu2mps,
        this->consts.G,         this->consts.clight, this->consts.j2000Jd,
        this->consts.JdMinusMjd};
    return constants;
}

/**
 * @return std::vector<real> Vector of simulation integration parameters.
 */
std::vector<real> PropSimulation::get_integration_parameters() {
    std::vector<real> integration_parameters = {
        (real)this->integParams.nInteg,
        (real)this->integParams.nSpice,
        (real)this->integParams.nTotal,
        this->integParams.t0,
        this->integParams.tf,
        (real)this->integParams.adaptiveTimestep,
        this->integParams.dt0,
        this->integParams.dtMin,
        this->integParams.dtChangeFactor,
        this->integParams.tolInteg,
        this->integParams.tolPC};
    return integration_parameters;
}

void PropSimulation::preprocess() {
    if (!this->isPreprocessed) {
        this->t = this->integParams.t0;
        for (size_t i = 0; i < this->integParams.nInteg; i++) {
            if (this->integBodies[i].isCometary) {
                double sunState[9];
                get_spk_state(10, this->integBodies[i].t0, this->spkEphem, sunState);
                this->integBodies[i].pos[0] += sunState[0];
                this->integBodies[i].pos[1] += sunState[1];
                this->integBodies[i].pos[2] += sunState[2];
                this->integBodies[i].vel[0] += sunState[3];
                this->integBodies[i].vel[1] += sunState[4];
                this->integBodies[i].vel[2] += sunState[5];
                this->integBodies[i].initCart[0] = this->integBodies[i].pos[0];
                this->integBodies[i].initCart[1] = this->integBodies[i].pos[1];
                this->integBodies[i].initCart[2] = this->integBodies[i].pos[2];
                this->integBodies[i].initCart[3] = this->integBodies[i].vel[0];
                this->integBodies[i].initCart[4] = this->integBodies[i].vel[1];
                this->integBodies[i].initCart[5] = this->integBodies[i].vel[2];
            }
            for (size_t j = 0; j < 3; j++) {
                this->xInteg.push_back(this->integBodies[i].pos[j]);
            }
            for (size_t j = 0; j < 3; j++) {
                this->xInteg.push_back(this->integBodies[i].vel[j]);
            }
            if (this->integBodies[i].propStm) {
                for (size_t j = 0; j < this->integBodies[i].stm.size(); j++) {
                    this->xInteg.push_back(this->integBodies[i].stm[j]);
                }
            }
        }
        bool backwardProp = this->integParams.t0 > this->integParams.tf;
        if (backwardProp) {
            std::reverse(this->eventMngr.impulsiveEvents.begin(), this->eventMngr.impulsiveEvents.end());
            std::reverse(this->eventMngr.continuousEvents.begin(), this->eventMngr.continuousEvents.end());
        }
        this->eventMngr.nextImpEventIdx = 0;
        this->eventMngr.nextConEventIdx = 0;
        if (this->eventMngr.nImpEvents > 0) {
            this->eventMngr.tNextImpEvent = this->eventMngr.impulsiveEvents[0].t;
        }
        if (this->eventMngr.nConEvents > 0) {
            this->eventMngr.tNextConEvent = this->eventMngr.continuousEvents[0].t;
            this->eventMngr.allConEventDone = false;
        } else {
            this->eventMngr.allConEventDone = true;
        }
        this->isPreprocessed = true;
    }
}
/**
 * @brief Extend the simulation to a new final time.
 * 
 * @param[in] tf New final time.
 * @param[in] tEvalNew New vector of times at which to evaluate the integrated state.
 * @param[in] observerInfoNew New vector of observer information.
 */
void PropSimulation::extend(real tf, std::vector<real> tEvalNew,
                            std::vector<std::vector<real>> observerInfoNew) {
    std::cout << "WARNING: The extend() method is under development and may "
                 "not work properly with the interpolation/observable "
                 "computation routines."
              << std::endl;

    // empty existing vectors from previous integration
    this->caParams.clear();
    this->interpParams.t0 = 0.0;
    this->interpParams.dt0 = 0.0;
    this->interpParams.xInteg0.clear();
    this->interpParams.b0.clear();
    this->interpParams.accInteg0.clear();
    this->interpParams.tStack.clear();
    this->interpParams.xIntegStack.clear();
    this->interpParams.bStack.clear();
    this->interpParams.accIntegStack.clear();
    this->interpParams.interpIdx = 0;
    this->xObserver.clear();
    this->observerInfo.clear();
    this->tEval.clear();
    this->obsType.clear();
    this->lightTimeEval.clear();
    this->xIntegEval.clear();
    this->opticalObs.clear();
    this->opticalPartials.clear();
    this->opticalObsCorr.clear();
    this->radarObs.clear();
    this->radarPartials.clear();

    // first prepare for integration and then integrate
    this->integParams.t0 = this->t;
    this->set_integration_parameters(tf, tEvalNew, this->tEvalUTC,
                                     this->evalApparentState,
                                     this->convergedLightTime, observerInfoNew);
    this->integrate();
}

/**
 * @param[in] filename Name of the file to save the simulation to.
 */
void PropSimulation::save(std::string filename) {
    struct stat fileExistsBuffer;
    if (stat(filename.c_str(), &fileExistsBuffer) == 0) {
        throw std::invalid_argument("Cannot save PropSimulation '" + this->name +
                                    "' to file " + filename +
                                    " (File already exists).");
    }
    auto timeWidth = std::setw(8);
    auto SpiceIdWidth = std::setw(7);
    auto floatWidth = std::setw(12);
    auto timeFloatPrec = std::setprecision(8);
    auto doubleWidth = std::setw(20);
    auto doublePrec = std::setprecision(12);
    int maxChars = 121;
    std::string tab = "    ";
    std::string halfTab = "  ";
    std::string subsectionFull = std::string(maxChars, '-');
    std::string sectionFull = std::string(maxChars, '=');
    std::string headerSectionHalf = std::string((int)(maxChars-13)/2, '=');
    std::ofstream file(filename, std::ios::out);
    // print header
    file << sectionFull << std::endl;
    #if defined(GRSS_VERSION)
        file << headerSectionHalf << " GRSS v" << GRSS_VERSION <<" " << headerSectionHalf << std::endl;
    #else
        file << headerSectionHalf << " GRSS vINFTY " << headerSectionHalf << std::endl;
    #endif
    file << sectionFull << std::endl;

    time_t now = time(nullptr);
    tm *utc = gmtime(&now);
    const char *format = "%a %b %d %Y %H:%M:%S UTC";
    char buffer[80];
    strftime(buffer, 80, format, utc);
    std::string timeStr(buffer);
    file << std::string((int)(maxChars-timeStr.size())/2, ' ') << timeStr << std::string((int)(maxChars-timeStr.size())/2, ' ') << std::endl;

    file << std::endl;
    std::string nextSubsection = this->name;
    file << std::string((int)(maxChars-nextSubsection.size())/2, '-') << nextSubsection << std::string((int)(maxChars-nextSubsection.size())/2, '-') << std::endl;
    file << subsectionFull << std::endl;
    file << "Integration from MJD " << timeWidth << std::fixed << timeFloatPrec << this->integParams.t0
         << " to MJD " << timeWidth << std::fixed << timeFloatPrec << this->integParams.tf << " [TDB]"
         << std::endl;
    file << "Main-body kernel path: " << this->spkEphem.mbPath << std::endl;
    file << "Small-body kernel path: " << this->spkEphem.sbPath << std::endl;

    file << std::endl;
    nextSubsection = std::to_string(this->integParams.nInteg) + " Integration bodies";
    file << std::string((int)(maxChars-nextSubsection.size())/2, '-') << nextSubsection << std::string((int)(maxChars-nextSubsection.size())/2, '-') << std::endl;
    file << subsectionFull << std::endl;
    size_t starti = 0;
    for (size_t i = 0; i < this->integBodies.size(); i++) {
        file << this->integBodies[i].name << std::endl;
        file << "Initial time : MJD " << timeWidth << std::fixed << timeFloatPrec << this->integBodies[i].t0
             << " [TDB]" << std::endl;
        file << "Radius       : " << floatWidth << std::fixed << timeFloatPrec << this->integBodies[i].radius*this->consts.du2m << " [m]" << std::endl;
        file << "Mass         : " << floatWidth << std::fixed << timeFloatPrec << this->integBodies[i].mass << " [kg]" << std::endl;
        file << "Initial state:" << std::endl;
        if (this->integBodies[i].isCometary) {
            file << " Cometary, heliocentric IAU76/J2000 ecliptic:" << std::endl;
            file << doubleWidth << "ecc" << doubleWidth << "peri. dist. [AU]" << doubleWidth << "peri. t. [MJD TDB]"
                    << doubleWidth << "l. asc. node [rad]" << doubleWidth << "arg. peri. [rad]" << doubleWidth << "inc. [rad]" << std::endl;
            for (size_t j = 0; j < 6; j++) {
                file << doubleWidth << std::scientific << doublePrec << this->integBodies[i].initState[j];
            }
            file << std::endl;
        }
        file << " Cartesian, J2000 barycentric:" << std::endl;
        file << doubleWidth << "x [AU]" << doubleWidth << "y [AU]" << doubleWidth << "z [AU]"
                << doubleWidth << "vx [AU/day]" << doubleWidth << "vy [AU/day]" << doubleWidth << "vz [AU/day]" << std::endl;
        for (size_t j = 0; j < 6; j++) {
            file << doubleWidth << std::scientific << doublePrec << this->integBodies[i].initCart[j];
        }
        file << std::endl;
        file << "Final state:" << std::endl;
        file << " Cartesian, J2000 barycentric:" << std::endl;
        file << doubleWidth << "x [AU]" << doubleWidth << "y [AU]" << doubleWidth << "z [AU]"
                << doubleWidth << "vx [AU/day]" << doubleWidth << "vy [AU/day]" << doubleWidth << "vz [AU/day]" << std::endl;
        for (size_t j = 0; j < 6; j++) {
            file << doubleWidth << std::scientific << doublePrec << this->xInteg[starti + j];
        }
        file << std::endl;
        starti += 6;
        if (this->integBodies[i].propStm) {
            size_t numParams = (this->integBodies[i].n2Derivs - 21)/3;
            if (this->integBodies[i].isCometary) {
                file << "STM (final Cartesian state w.r.t. initial Cometary state + any params):" << std::endl;
                file << doubleWidth << "ecc" << doubleWidth << "peri. dist. [AU]" << doubleWidth << "peri. t. [MJD TDB]"
                        << doubleWidth << "l. asc. node [rad]" << doubleWidth << "arg. peri. [rad]" << doubleWidth << "inc. [rad]";
            } else {
                file << "STM (final Cartesian state w.r.t. initial Cartesian state + any params):" << std::endl;
                file << doubleWidth << "x [AU]" << doubleWidth << "y [AU]" << doubleWidth << "z [AU]"
                        << doubleWidth << "vx [AU/day]" << doubleWidth << "vy [AU/day]" << doubleWidth << "vz [AU/day]";
            }
            if (this->integBodies[i].ngParams.a1Est) file << doubleWidth << "A1 [AU/day^2]";
            if (this->integBodies[i].ngParams.a2Est) file << doubleWidth << "A2 [AU/day^2]";
            if (this->integBodies[i].ngParams.a3Est) file << doubleWidth << "A3 [AU/day^2]";
            file << std::endl;
            std::vector<real> stmFinalFlat = std::vector<real>(36+6*numParams, 0.0);
            for (size_t j = 0; j < 36+6*numParams; j++) {
                stmFinalFlat[j] = this->xInteg[starti + j];
            }
            std::vector<std::vector<real>> stmFinal = reconstruct_stm(stmFinalFlat);
            for (size_t j = 0; j < stmFinal.size(); j++) {
                for (size_t k = 0; k < stmFinal[0].size(); k++) {
                    file << doubleWidth << std::scientific << doublePrec << stmFinal[j][k];
                }
                file << std::endl;
            }
            if (this->integBodies[i].isCometary) {
                file << "STM (initial Cartesian state + any params w.r.t. initial Cometary state + any params):" << std::endl;
                file << doubleWidth << "ecc" << doubleWidth << "peri. dist. [AU]" << doubleWidth << "peri. t. [MJD TDB]"
                        << doubleWidth << "l. asc. node [rad]" << doubleWidth << "arg. peri. [rad]" << doubleWidth << "inc. [rad]";
                if (this->integBodies[i].ngParams.a1Est) file << doubleWidth << "A1 [AU/day^2]";
                if (this->integBodies[i].ngParams.a2Est) file << doubleWidth << "A2 [AU/day^2]";
                if (this->integBodies[i].ngParams.a3Est) file << doubleWidth << "A3 [AU/day^2]";
                file << std::endl;
                for (size_t j = 0; j < this->integBodies[i].dCartdState.size(); j++) {
                    for (size_t k = 0; k < this->integBodies[i].dCartdState[0].size(); k++) {
                        file << doubleWidth << std::scientific << doublePrec << this->integBodies[i].dCartdState[j][k];
                    }
                    file << std::endl;
                }
            }
        }
        if (i!=this->integBodies.size()-1) file << std::endl << subsectionFull << std::endl;
    }

    file << std::endl;
    nextSubsection = std::to_string(this->integParams.nSpice) + " SPICE bodies";
    file << std::string((int)(maxChars-nextSubsection.size())/2, '-') << nextSubsection << std::string((int)(maxChars-nextSubsection.size())/2, '-') << std::endl;
    file << subsectionFull << std::endl;
    for (size_t i = 0; i < this->spiceBodies.size(); i++) {
        file << SpiceIdWidth << this->spiceBodies[i].spiceId << halfTab
             << this->spiceBodies[i].name << std::endl;
    }

    if (this->impactParams.size() != 0) {
        file << std::endl;
        nextSubsection = std::to_string(this->impactParams.size()) + " Impacts detected";
        file << std::string((int)(maxChars-nextSubsection.size())/2, '-') << nextSubsection << std::string((int)(maxChars-nextSubsection.size())/2, '-') << std::endl;
        file << subsectionFull << std::endl;
        for (size_t i = 0; i < this->impactParams.size(); i++) {
            ImpactParameters imp = this->impactParams[i];
            file << "MJD " << timeWidth << std::fixed << imp.t << " TDB:" << std::endl;
            file << "  " << imp.flybyBody << " impacted " << imp.centralBody << " with a relative velocity of " << imp.vel << " AU/d." << std::endl;
            file << "  Impact location: " << std::endl;
            file << "    Longitude: " << imp.lon*RAD2DEG << " deg" << std::endl;
            file << "    Latitude: " << imp.lat*RAD2DEG << " deg" << std::endl;
            file << "    Altitude: " << imp.alt*this->consts.du2m/1.0e3 << " km" << std::endl;
        }
    }

    if (this->caParams.size() != 0) {
        file << std::endl;
        nextSubsection = std::to_string(this->caParams.size()) + " Close approaches detected";
        file << std::string((int)(maxChars-nextSubsection.size())/2, '-') << nextSubsection << std::string((int)(maxChars-nextSubsection.size())/2, '-') << std::endl;
        file << subsectionFull << std::endl;
        for (size_t i = 0; i < this->caParams.size(); i++) {
            CloseApproachParameters ca = this->caParams[i];
            file << "MJD " << timeWidth << std::fixed << ca.t << " TDB:" << std::endl;
            file << "  " << ca.flybyBody << " approached " << ca.centralBody << " at " << ca.dist << " AU." << std::endl;
            file << "  Relative Velocity: " << ca.vel << " AU/d. V-infinity: " << ca.vInf << " AU/d." << std::endl;
            file << "  Gravitational focusing factor: " << ca.gravFocusFactor << ". Impact: " << std::boolalpha << ca.impact << std::endl;
        }
    }

    // insert machine readable data for close approach and impact parameters
    file.precision(18);
    if (this->impactParams.size() > 0) {
        file << std::endl;
        nextSubsection = "Impact data";
        file << std::string((int)(maxChars-nextSubsection.size())/2, '-') << nextSubsection << std::string((int)(maxChars-nextSubsection.size())/2, '-') << std::endl;
        file << subsectionFull << std::endl;
        file << "$$IMPACT_START" << std::endl;
        file << "t | xRel | tMap | xRelMap | dist | vel | vInf | flybyBody | flybyBodyIdx | centralBody | centralBodyIdx | centralBodySpiceId | impact | tPeri | tLin | bVec | bMag | gravFocusFactor | kizner_x | kizner_y | kizner_z | opik_x | opik_y | opik_z | scaled_x | scaled_y | scaled_z | mtp_x | mtp_y | mtp_z | xRelBodyFixed | lon | lat | alt" << std::endl;
        for (size_t i = 0; i < this->impactParams.size(); i++) {
            ImpactParameters imp = this->impactParams[i];
            file << imp.t << " | ";
            file << "[";
            for (size_t j = 0; j < imp.xRel.size(); j++) {
                file << imp.xRel[j] << ",";
            }
            file << "] | ";
            file << imp.tMap << " | ";
            file << "[";
            for (size_t j = 0; j < imp.xRelMap.size(); j++) {
                file << imp.xRelMap[j] << ",";
            }
            file << "] | ";
            file << imp.dist << " | ";
            file << imp.vel << " | ";
            file << imp.vInf << " | ";
            file << imp.flybyBody << " | ";
            file << imp.flybyBodyIdx << " | ";
            file << imp.centralBody << " | ";
            file << imp.centralBodyIdx << " | ";
            file << imp.centralBodySpiceId << " | ";
            file << imp.impact << " | ";
            file << imp.tPeri << " | ";
            file << imp.tLin << " | ";
            file << "[";
            for (size_t j = 0; j < imp.bVec.size(); j++) {
                file << imp.bVec[j] << ",";
            }
            file << "] | ";
            file << imp.bMag << " | ";
            file << imp.gravFocusFactor << " | ";
            file << imp.kizner.x << " | ";
            file << imp.kizner.y << " | ";
            file << imp.kizner.z << " | ";
            file << imp.opik.x << " | ";
            file << imp.opik.y << " | ";
            file << imp.opik.z << " | ";
            file << imp.scaled.x << " | ";
            file << imp.scaled.y << " | ";
            file << imp.scaled.z << " | ";
            file << imp.mtp.x << " | ";
            file << imp.mtp.y << " | ";
            file << imp.mtp.z << " | ";
            file << "[";
            for (size_t j = 0; j < imp.xRelBodyFixed.size(); j++) {
                file << imp.xRelBodyFixed[j] << ",";
            }
            file << "] | ";
            file << imp.lon << " | ";
            file << imp.lat << " | ";
            file << imp.alt;
            file << std::endl;
        }
        file << "$$IMPACT_END" << std::endl;
    }
    if (this->caParams.size() > 0) {
        file << std::endl;
        nextSubsection = "Close approach data";
        file << std::string((int)(maxChars-nextSubsection.size())/2, '-') << nextSubsection << std::string((int)(maxChars-nextSubsection.size())/2, '-') << std::endl;
        file << subsectionFull << std::endl;
        file << "$$CA_START" << std::endl;
        file << "t | xRel | tMap | xRelMap | dist | vel | vInf | flybyBody | flybyBodyIdx | centralBody | centralBodyIdx | centralBodySpiceId | impact | tPeri | tLin | bVec | bMag | gravFocusFactor | kizner_x | kizner_y | kizner_z | kizner_dx | kizner_dy | opik_x | opik_y | opik_z | opik_dx | opik_dy | scaled_x | scaled_y | scaled_z | scaled_dx | scaled_dy | mtp_x | mtp_y | mtp_z | mtp_dx | mtp_dy" << std::endl;
        for (size_t i = 0; i < this->caParams.size(); i++) {
            CloseApproachParameters ca = this->caParams[i];
            file << ca.t << " | ";
            file << "[";
            for (size_t j = 0; j < ca.xRel.size(); j++) {
                file << ca.xRel[j] << ",";
            }
            file << "] | ";
            file << ca.tMap << " | ";
            file << "[";
            for (size_t j = 0; j < ca.xRelMap.size(); j++) {
                file << ca.xRelMap[j] << ",";
            }
            file << "] | ";
            file << ca.dist << " | ";
            file << ca.vel << " | ";
            file << ca.vInf << " | ";
            file << ca.flybyBody << " | ";
            file << ca.flybyBodyIdx << " | ";
            file << ca.centralBody << " | ";
            file << ca.centralBodyIdx << " | ";
            file << ca.centralBodySpiceId << " | ";
            file << ca.impact << " | ";
            file << ca.tPeri << " | ";
            file << ca.tLin << " | ";
            file << "[";
            for (size_t j = 0; j < ca.bVec.size(); j++) {
                file << ca.bVec[j] << ",";
            }
            file << "] | ";
            file << ca.bMag << " | ";
            file << ca.gravFocusFactor << " | ";
            file << ca.kizner.x << " | ";
            file << ca.kizner.y << " | ";
            file << ca.kizner.z << " | ";
            file << "[";
            for (size_t j = 0; j < ca.kizner.dx.size(); j++) {
                file << ca.kizner.dx[j] << ",";
            }
            file << "] | ";
            file << "[";
            for (size_t j = 0; j < ca.kizner.dy.size(); j++) {
                file << ca.kizner.dy[j] << ",";
            }
            file << "] | ";
            file << ca.opik.x << " | ";
            file << ca.opik.y << " | ";
            file << ca.opik.z << " | ";
            file << "[";
            for (size_t j = 0; j < ca.opik.dx.size(); j++) {
                file << ca.opik.dx[j] << ",";
            }
            file << "] | ";
            file << "[";
            for (size_t j = 0; j < ca.opik.dy.size(); j++) {
                file << ca.opik.dy[j] << ",";
            }
            file << "] | ";
            file << ca.scaled.x << " | ";
            file << ca.scaled.y << " | ";
            file << ca.scaled.z << " | ";
            file << "[";
            for (size_t j = 0; j < ca.scaled.dx.size(); j++) {
                file << ca.scaled.dx[j] << ",";
            }
            file << "] | ";
            file << "[";
            for (size_t j = 0; j < ca.scaled.dy.size(); j++) {
                file << ca.scaled.dy[j] << ",";
            }
            file << "] | ";
            file << ca.mtp.x << " | ";
            file << ca.mtp.y << " | ";
            file << ca.mtp.z << " | ";
            file << "[";
            for (size_t j = 0; j < ca.mtp.dx.size(); j++) {
                file << ca.mtp.dx[j] << ",";
            }
            file << "] | ";
            file << "[";
            for (size_t j = 0; j < ca.mtp.dy.size(); j++) {
                file << ca.mtp.dy[j] << ",";
            }
            file << "]";
            file << std::endl;
        }
        file << "$$CA_END" << std::endl;
    }

    file << std::endl;
    file << sectionFull << std::endl;
    file << headerSectionHalf << " END OF FILE " << headerSectionHalf << std::endl;
    file << sectionFull << std::endl;
    file.close();
}
