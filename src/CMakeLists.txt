configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/online/sts_config.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/online/sts_config.h"
  )

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/online")

set(STS_CPP_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/online/bpp_shim.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/gsl.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/util.cc

  # ONLINE
  ${CMAKE_CURRENT_SOURCE_DIR}/online/beagle_tree_likelihood.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/branch_length_prior.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/composite_tree_likelihood.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/guided_online_add_sequence_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/likelihood_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/multiplier_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/multiplier_proposal.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/multiplier_smc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/node_slider_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/node_slider_smc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/online_add_sequence_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/online_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/online_smc_init.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/online_util.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/tree_particle.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/uniform_online_add_sequence_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/attachment_likelihood.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/lcfit_online_add_sequence_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/lcfit_rejection_sampler.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/tripod_optimizer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online/uniform_length_online_add_sequence_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../lib/lcfit/lcfit_cpp_src/lcfit_cpp.cc
  )
add_library(sts-static STATIC ${STS_CPP_FILES})

set(STS_PHYLO_LIBS
  ${BPP_LIBRARIES}
  smctc
  ${GSL_LIBRARIES}
  smctc
  jsoncpp
  lcfit-static
  ${HMS_BEAGLE_LIBRARIES})

# Export
set(STS_PHYLO_LIBS ${STS_PHYLO_LIBS} PARENT_SCOPE)

add_executable(sts-online ${CMAKE_CURRENT_SOURCE_DIR}/online/sts_online.cc)
target_link_libraries(sts-online sts-static ${STS_PHYLO_LIBS})
install(TARGETS sts-online
        RUNTIME DESTINATION bin)
