
find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
    message("Doxygen not found, I will not generate/install the documentation")
else()
   configure_file(doxyfilerc Doxyfile)

   set(DOXYGEN_INPUT ${CMAKE_BINARY_DIR}/docs/Doxyfile)
   set(DOXYGEN_OUTPUT ${CMAKE_BINARY_DIR}/docs/html/index.html)

   add_custom_target(doc ALL
     COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
     COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT} > /dev/null
     COMMAND ${CMAKE_COMMAND} -E echo "Done."
     )

  install(DIRECTORY ${CMAKE_BINARY_DIR}/docs/html
          DESTINATION ${CMAKE_BINARY_DIR}/docs)

endif()
