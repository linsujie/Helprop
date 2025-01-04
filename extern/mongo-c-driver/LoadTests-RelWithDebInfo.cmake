      set (TEST_LIBMONGOC_EXE [[/home/linsj/works/solar/HelProp/mongo-c-driver/src/libmongoc/test-libmongoc]])
      set (SRC_ROOT [[/home/linsj/works/solar/HelProp/mongo-c-driver]])
      set (IS_MULTICONF 0)
      if (NOT IS_MULTICONF OR CTEST_CONFIGURATION_TYPE STREQUAL "RelWithDebInfo")
         # We are not in multi-conf, or the current config matches our config.
         include ("${SRC_ROOT}/build/cmake/LoadTests.cmake")
      elseif (NOT CTEST_CONFIGURATION_TYPE)
         # We are in multi-conf, but no '-C' config was specified
         message (WARNING "Specify a --build-config when using CTest with a multi-config build")
      else ()
         # Do nothing. Not our config.
      endif ()
   