add_executable(TEST_AES simpleAES.cpp homAES.cpp Test_AES.cpp)
target_link_libraries(TEST_AES PUBLIC helib)
