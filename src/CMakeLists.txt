
include_directories ("${PROJECT_SOURCE_DIR}/src" ${BROOM_INCLUDES})

# Add other common functions here

add_library( broom STATIC weights.cpp physics.cpp Quadrature.cpp Solve.cpp DSASolve.cpp diffSolve.cpp ProblemSpecifications.cpp GridFuncBroom.cpp SwifferProblems.cpp MeshLoader.cc Distortions.cc )

install( TARGETS broom DESTINATION lib )
install( FILES physics.hpp Quadrature.hpp Solve.hpp DSASolve.hpp diffSolve.hpp ProblemSpecifications.hpp GridFuncBroom.hpp SwifferProblems.hpp DESTINATION include )

add_executable(Broom Broom.cpp)
add_executable(Swiffer Swiffer.cpp)

target_link_libraries( Broom broom ${BROOM_LIBS} )
target_link_libraries( Swiffer broom ${BROOM_LIBS} )

install( TARGETS Broom Swiffer DESTINATION bin )
