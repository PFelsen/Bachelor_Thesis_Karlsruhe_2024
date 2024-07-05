# modified CMake file to make cxsc compile together with mpp

include_directories(
    # where cxscconf.h is
    ${PROJECT_MPP_DIR}/cxsc

    # the C-XSC sources
    ${PROJECT_MPP_DIR}/cxsc/src
    ${PROJECT_MPP_DIR}/cxsc/src/fi_lib
    ${PROJECT_MPP_DIR}/cxsc/src/rts
    ${PROJECT_MPP_DIR}/cxsc/src/asm)

include(src/CMakeLists.txt)
include(src/fi_lib/CMakeLists.txt)
include(src/rts/CMakeLists.txt)

add_library(cxsc STATIC 
    ${cxsc_src}
    ${rts_src}
    ${fi_lib_src})