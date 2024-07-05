# Executables
add_executable(VerifyQintSym VerifyQintSym.cpp)
add_executable(VerifyQtriSym VerifyQtriSym.cpp)

# Linking
target_link_libraries(VerifyQintSym cxsc cxsc_toolbox)
target_link_libraries(VerifyQtriSym cxsc cxsc_toolbox)

#---------------------------------------------------------------------------------------#
