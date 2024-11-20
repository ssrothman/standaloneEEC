#include <stdio.h>
#include "test_proj.h"
#include "test_res3.h"
#include "test_res4.h"

constexpr unsigned NPART = 50;
constexpr unsigned NTEST = 10;

constexpr const char* PASSSTR = "\x1B[32mpass\033[0m";
constexpr const char* FAILSTR = "\x1B[31mfail\033[0m";

int main() noexcept {

    printf("testing proj...\n");
    printf("\tsecond order: %s\n", 
            test_proj(NTEST, NPART, 2) ? PASSSTR : FAILSTR);
    printf("\tthird order : %s\n",  
            test_proj(NTEST, NPART, 3) ? PASSSTR : FAILSTR);
    printf("\tfourth order: %s\n", 
            test_proj(NTEST, NPART, 4) ? PASSSTR : FAILSTR);
    //printf("\tfifth order : %s\n",  
    //        test_proj(NTEST, NPART, 5) ? PASSSTR : FAILSTR);
    //printf("\tsixth order : %s\n",  
    //        test_proj(NTEST, NPART, 6) ? PASSSTR : FAILSTR);
    //printf("\n");

    printf("testing resolved...\n");
    printf("\tres3: %s\n",
            test_res3(NTEST, NPART) ? PASSSTR : FAILSTR);
    printf("\tres4: %s\n",
            test_res4(NTEST, NPART) ? PASSSTR : FAILSTR);
    printf("\n");

    printf("testing PU...\n");
    printf("\tproj2: %s\n",
            test_proj_PU(NTEST, NPART, 2) ? PASSSTR : FAILSTR);
    printf("\tproj3: %s\n",
            test_proj_PU(NTEST, NPART, 3) ? PASSSTR : FAILSTR);
    printf("\tproj4: %s\n",
            test_proj_PU(NTEST, NPART, 4) ? PASSSTR : FAILSTR);
    //printf("\tproj5: %s\n",
    //        test_proj_PU(NTEST, NPART, 5) ? PASSSTR : FAILSTR);
    //printf("\tproj6: %s\n",
    //        test_proj_PU(NTEST, NPART, 6) ? PASSSTR : FAILSTR);
    printf("\tres3: %s\n",
            test_res3_PU(NTEST, NPART) ? PASSSTR : FAILSTR);
    printf("\tres4: %s\n",
            test_res4_PU(NTEST, NPART) ? PASSSTR : FAILSTR);
    printf("\n");

    printf("testing transfer...\n");
    printf("\tproj2: %s\n",
            test_proj_transfer(NTEST, NPART, 2) ? PASSSTR : FAILSTR);
    printf("\tproj3: %s\n",
            test_proj_transfer(NTEST, NPART, 3) ? PASSSTR : FAILSTR);
    printf("\tproj4: %s\n",
            test_proj_transfer(NTEST, NPART, 4) ? PASSSTR : FAILSTR);
    //printf("\tproj5: %s\n",
    //        test_proj_transfer(NTEST, NPART, 5) ? PASSSTR : FAILSTR);
    //printf("\tproj6: %s\n",
    //        test_proj_transfer(NTEST, NPART, 6) ? PASSSTR : FAILSTR);
    printf("\tres3: %s\n",
            test_res3_transfer(NTEST, NPART) ? PASSSTR : FAILSTR);
    return 0;
}
