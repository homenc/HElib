set(SRC
BenesNetwork.cpp
blockMatmul1D.cpp
blockMatmul.cpp
bluestein.cpp
CModulus.cpp
Ctxt.cpp
debugging.cpp
DoubleCRT.cpp
EncryptedArray.cpp
EvalMap.cpp
extractDigits.cpp
FHEContext.cpp
FHE.cpp
hypercube.cpp
IndexSet.cpp
KeySwitching.cpp
matching.cpp
matmul1D.cpp
matmul.cpp
NumbTh.cpp
OptimizePermutations.cpp
PAlgebra.cpp
PermNetwork.cpp
permutations.cpp
polyEval.cpp
powerful.cpp
recryption.cpp
replicate.cpp
rotations.cpp
timing.cpp
)
if(FHE_THREADS)
    add_definitions(-DFHE_THREADS)
    set(RUN_LIB fhe ntl gmp pthread)
else(FHE_THREADS)
    set(RUN_LIB fhe ntl gmp)
endif(FHE_THREADS)

add_library(fhe STATIC SHARED ${SRC})

if (ENABLE_TEST)
    file(GLOB TEST_FILES RELATIVE "${CMAKE_SOURCE_DIR}/src" "Test_*.cpp")
    foreach(TEST_FILE ${TEST_FILES})
        add_executable(${TEST_FILE}_exe ${TEST_FILE})
        target_link_libraries(${TEST_FILE}_exe ${RUN_LIB})
    endforeach(TEST_FILE)
    set(OTHER_PRG 
        MulTime.cpp
        cgauss.cpp
        params.cpp
        eqtesting.cpp)
    foreach(PRG ${OTHER_PRG})
        add_executable(${PRG}_exe ${PRG})
        target_link_libraries(${PRG}_exe ${RUN_LIB})
    endforeach(PRG)
endif(ENABLE_TEST)
