configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/sts_config.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/sts_config.h"
  )
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/online")

set(STS_CPP_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/base_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bpp_shim.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/child_swap_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/delta_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/edge.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/exponential_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/forest_likelihood.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/gamma_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/gsl.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/log/json_logger.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/log/sampler.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/metropolis_hastings_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/node.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/online_calculator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rooted_merge.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/smc_init.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/smc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/state.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/uniform_bl_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/uniform_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/util.cc

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

add_executable(sts ${CMAKE_CURRENT_SOURCE_DIR}/sts.cc)
target_link_libraries(sts sts-static ${STS_PHYLO_LIBS})
install(TARGETS sts
        RUNTIME DESTINATION bin)

add_executable(sts-online ${CMAKE_CURRENT_SOURCE_DIR}/online/sts_online.cc)
target_link_libraries(sts-online sts-static ${STS_PHYLO_LIBS})
install(TARGETS sts-online
        RUNTIME DESTINATION bin)
