#include <gtest/gtest.h>
#include "../engine/linalg/specialfunctions.h"

#define NEAR_THRESH 1e-7

TEST(FStatToPVal, SanityCheck) {

    // Test against 
    double pval = alglib::fdistribution(3000, 1, 5.0);
    ASSERT_NEAR(pval,0.65475315127310253,NEAR_THRESH);
    
    pval = alglib::fdistribution(1,1,5.0);
    ASSERT_NEAR(pval,0.73227952719876888,NEAR_THRESH);
    
    pval = alglib::fdistribution(1, 3000, 5.0);
    ASSERT_NEAR(pval,0.97457942728192826,NEAR_THRESH);
    
}

