#ifndef GR15_H
#define GR15_H
#include "utilities.h"
#include "simulation.h"

const std::vector<real> hVec = {
    0.0,
    0.0562625605369221464656522,
    0.1802406917368923649875799,
    0.3526247171131696373739078,
    0.5471536263305553830014486,
    0.7342101772154105315232106,
    0.8853209468390957680903598,
    0.9775206135612875018911745
    };
const real hSize = 8;
const std::vector< std::vector<real> > rMat = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {17.773808914078000840752659565672904106978971632681, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {5.5481367185372165056928216140765061758579336941398, 8.0659386483818866885371256689687154412267416180207, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {2.8358760786444386782520104428042437400879003147949, 3.3742499769626352599420358188267460448330087696743, 5.8010015592640614823286778893918880155743979164251, 0.0, 0.0, 0.0, 0.0, 0.0},
    {1.8276402675175978297946077587371204385651628457154, 2.0371118353585847827949159161566554921841792590404, 2.7254422118082262837742722003491334729711450288807, 5.1406241058109342286363199091504437929335189668304, 0.0, 0.0, 0.0, 0.0},
    {1.3620078160624694969370006292445650994197371928318, 1.4750402175604115479218482480167404024740127431358, 1.8051535801402512604391147435448679586574414080693, 2.6206449263870350811541816031933074696730227729812, 5.34597689987110751412149096322778980457703366603548, 0.0, 0.0, 0.0},
    {1.1295338753367899027322861542728593509768148769105, 1.2061876660584456166252036299646227791474203527801, 1.4182782637347391537713783674858328433713640692518, 1.8772424961868100972169920283109658335427446084411, 2.9571160172904557478071040204245556508352776929762, 6.6176620137024244874471284891193925737033291491748, 0.0, 0.0},
    {1.0229963298234867458386119071939636779024159134103, 1.0854721939386423840467243172568913862030118679827, 1.2542646222818777659905422465868249586862369725826, 1.6002665494908162609916716949161150366323259154408, 2.3235983002196942228325345451091668073608955835034, 4.1099757783445590862385761824068782144723082633980, 10.846026190236844684706431007823415424143683137181, 0.0}
    };
const std::vector< std::vector<real> > cMat = {
    {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {-0.562625605369221464656522e-1, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.1014080283006362998648180399549641417413495311078e-1, -0.2365032522738145114532321e0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {-0.35758977292516175949344589284567187362040464593728e-2, 0.9353769525946206589574845561035371499343547051116e-1, -0.5891279693869841488271399e0, 1, 0, 0, 0, 0},
    {0.19565654099472210769005672379668610648179838140913e-2, -0.54755386889068686440808430671055022602028382584495e-1, 0.41588120008230686168862193041156933067050816537030e0, -0.11362815957175395318285885e1, 1, 0, 0, 0},
    {-0.14365302363708915424459554194153247134438571962198e-2, 0.42158527721268707707297347813203202980228135395858e-1, -0.36009959650205681228976647408968845289781580280782e0, 0.12501507118406910258505441186857527694077565516084e1, -0.18704917729329500633517991e1, 1, 0, 0},
    {0.12717903090268677492943117622964220889484666147501e-2, -0.38760357915906770369904626849901899108502158354383e-1, 0.36096224345284598322533983078129066420907893718190e0, -0.14668842084004269643701553461378480148761655599754e1, 0.29061362593084293014237914371173946705384212479246e1, -0.27558127197720458314421589e1, 1, 0},
    {-0.12432012432012432012432013849038719237133940238163e-2, 0.39160839160839160839160841227582657239289159887563e-1, -0.39160839160839160839160841545895262429018228668896e0, 0.17948717948717948717948719027866738711862551337629e1, -0.43076923076923076923076925231853900723503338586335e1, 0.56000000000000000000000001961129300233768803845526e1, -0.37333333333333333333333334e1, 1}
    };

void gr15(real t, std::vector<real> xInteg, Simulation &sim);

void Simulation::integrate(){
    gr15(this->t, this->xInteg, *this);
}
#endif
