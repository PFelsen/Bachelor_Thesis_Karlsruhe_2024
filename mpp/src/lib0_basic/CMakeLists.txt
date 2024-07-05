# Clean up duplicates and evaluate values
list(REMOVE_DUPLICATES SET_OPTIONS)
foreach(option IN LISTS SET_OPTIONS)
    list(APPEND PRINT_OPTION_ENTRIES_ARGS ${option} ${${option}})
endforeach()

# Generate PrintOptionEntries.cpp
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/generated/logging/PrintOptionEntries.cpp
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/logging/gen_build_info.py
        -o ${CMAKE_CURRENT_BINARY_DIR}/generated/logging/PrintOptionEntries.cpp
        ${PRINT_OPTION_ENTRIES_ARGS}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/logging/gen_build_info.py
    COMMENT "Generate PrintOptionEntries.cpp"
    VERBATIM
    COMMAND_EXPAND_LISTS)

set(basic_src
        Point.cpp
        utility/TimeDate.cpp
        utility/Config.cpp
        utility/ctools.cpp
        utility/SaveLoad.cpp
        utility/Assertion.cpp
        utility/M_IOFiles.cpp
        utility/MemoryLogger.cpp
        logging/Logging.cpp
        parallel/Buffer.cpp
        parallel/ExchangeBuffer.cpp
        parallel/Parallel.cpp
        parallel/Communicator.cpp
        utility/StringUtil.cpp
        utility/Functools.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/generated/logging/PrintOptionEntries.cpp
        )

if (USE_CXSC)
    set(basic_src ${basic_src}
            intervalarithmetic/IAIO.cpp
            intervalarithmetic/IAInterval.cpp
            intervalarithmetic/IACInterval.cpp
            parallel/IAMpiOperators.cpp
            )
endif ()

add_library(LIB_BASIC SHARED ${basic_src})
target_link_libraries(LIB_BASIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} stdc++fs)

if (TIRPC_FOUND)
    target_link_libraries(LIB_BASIC ${TIRPC_LIBRARIES})
endif()


if (USE_CXSC)
    target_link_libraries(LIB_BASIC cxsc)
endif ()
