aux_source_directoy(. FUNC_DIR_SRC)
add_library(func ${FUNC_DIR_SRC})
ADD_EXECUTABLE(shock main.cc particle.cc import.cc)