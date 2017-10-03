
SET(TEST_LIBS
  ${Boost_LIBRARIES}
  libminipart
)

SET(TESTS
  hypergraph
  io
)

FOREACH(TEST_BASE ${TESTS})
  SET(TEST ${TEST_BASE}_test)
  ADD_EXECUTABLE(${TEST} ${TEST}.cc)
  target_link_libraries(${TEST} ${TEST_LIBS})
  ADD_TEST(${TEST} ${TEST})
ENDFOREACH()

