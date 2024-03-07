INCLUDE(CMakeParseArguments)

# Borrowed from Kokkos Kernels (https://github.com/kokkos/kokkos-kernels)
FUNCTION(DIVERG_ADD_OPTION SUFFIX DEFAULT TYPE DOCSTRING)
  CMAKE_PARSE_ARGUMENTS(OPT
    ""
    ""
    "VALID_ENTRIES" #if this is a list variable, the valid values in the list
    ${ARGN}
  )

  SET(CAMEL_NAME DiverG_${SUFFIX})
  STRING(TOUPPER ${CAMEL_NAME} UC_NAME)

  # Make sure this appears in the cache with the appropriate DOCSTRING
  SET(${CAMEL_NAME} ${DEFAULT} CACHE ${TYPE} ${DOCSTRING})

  #I don't love doing it this way because it's N^2 in number options, but cest la vie
  FOREACH(opt ${DIVERG_GIVEN_VARIABLES})
    STRING(TOUPPER ${opt} OPT_UC)
    IF ("${OPT_UC}" STREQUAL "${UC_NAME}")
      IF (NOT "${opt}" STREQUAL "${CAMEL_NAME}")
        MESSAGE(FATAL_ERROR "Matching option found for ${CAMEL_NAME} with the wrong case ${opt}. Please delete your CMakeCache.txt and change option to -D${CAMEL_NAME}=${${opt}}. This is now enforced to avoid hard-to-debug CMake cache inconsistencies.")
      ENDIF()
    ENDIF()
  ENDFOREACH()


  #okay, great, we passed the validation test - use the default
  IF (DEFINED ${CAMEL_NAME})
    IF (OPT_VALID_ENTRIES)
      STRING(TOUPPER   "${OPT_VALID_ENTRIES}" OPT_VALID_ENTRIES_UC)
      FOREACH(entry ${${CAMEL_NAME}})
        STRING(TOUPPER ${entry} ENTRY_UC)
        IF (NOT ${ENTRY_UC} IN_LIST OPT_VALID_ENTRIES_UC)
          MESSAGE(FATAL_ERROR "Given entry ${entry} in list for option ${SUFFIX}. "
                  "Valid case-insensitive values are any of ${OPT_VALID_ENTRIES}")
        ENDIF()
      ENDFOREACH()
      STRING(TOUPPER "${${CAMEL_NAME}}" GIVEN_ENTRIES_UC)
      SET(${UC_NAME} ${GIVEN_ENTRIES_UC} PARENT_SCOPE)
    ELSE()
      SET(${UC_NAME} ${${CAMEL_NAME}} PARENT_SCOPE)
    ENDIF()
  ELSE()
    SET(${UC_NAME} ${DEFAULT} PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

FUNCTION(DIVERG_IS_ENABLED)
  CMAKE_PARSE_ARGUMENTS(PARSE
    ""
    "OUTPUT_VARIABLE"
    "COMPONENTS"
    ${ARGN})

  IF (DIVERG_ENABLED_COMPONENTS STREQUAL "ALL")
    SET(${PARSE_OUTPUT_VARIABLE} TRUE PARENT_SCOPE)
  ELSEIF(PARSE_COMPONENTS)
    SET(ENABLED TRUE)
    FOREACH(comp ${PARSE_COMPONENTS})
      STRING(TOUPPER ${comp} COMP_UC)
      # make sure this is in the list of enabled components
      IF(NOT "${COMP_UC}" IN_LIST DIVERG_ENABLED_COMPONENTS)
        # if not in the list, one or more components is missing
        SET(ENABLED FALSE)
      ENDIF()
    ENDFOREACH()
    SET(${PARSE_OUTPUT_VARIABLE} ${ENABLED} PARENT_SCOPE)
  ELSE()
    # we did not enable all components and no components
    # were given as part of this - we consider this enabled
    SET(${PARSE_OUTPUT_VARIABLE} TRUE PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()
