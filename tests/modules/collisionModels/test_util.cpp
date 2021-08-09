#include "CollisionModel_util.hpp"
#include "catch.hpp"

TEST_CASE( "Test ion parameter estimation functions", "[CollisionModels][Util]") {

    SECTION( "Estimation function for collision diameter from mass should be correct") {
        //CollisionModel::StatisticalDiffusionModel sds(100000,273.15,28,);

        //n=6, m=91, d=0.541621(est), Ko=1.89928(est), MFPo=3.2037e-005, Vo=0.252097 [*]
        double diamFromMass_nm = CollisionModel::util::estimateCollisionDiameterFromMass(91);
        CHECK(Approx(diamFromMass_nm)== 0.541621);

        //n=4, m=55, d=0.457934(est), Ko=2.34714(est), MFPo=4.39026e-005, Vo=0.32427 [*]
        diamFromMass_nm = CollisionModel::util::estimateCollisionDiameterFromMass(55);
        CHECK(Approx(diamFromMass_nm)== 0.457934);
    }

    SECTION( "Estimation function for mobility should be correct") {

        //n=6, m=91, d=0.541621(est), Ko=1.89928(est), MFPo=3.2037e-005, Vo=0.252097 [*]
        double massIon_amu = 91;
        double diameterIon_nm = 0.541621;
        double massGas_amu = 28.94515; //SDS default value for air
        double diameterGas_nm = 0.366;

        //CollisionModel::StatisticalDiffusionModel sds(100000,273.15,massGas_amu,diameterGas);

        double mobility = CollisionModel::util::estimateMobility(
                massIon_amu, diameterIon_nm,massGas_amu,diameterGas_nm);

        CHECK(Approx(mobility)== 1.89928);
    }
}