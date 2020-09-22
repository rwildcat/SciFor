# -----------------------------------------------------------------------------
# CMake for RANDLIB as internal lib for SciFor
# Ramon Solano <ramon.solano@gmail.com>
# 
# RANDLIB90 Retrieved from:
#	https://biostatistics.mdanderson.org/SoftwareDownload/SingleSoftware/Index/27
#
# Sep/16/2020
#    initial version
# -----------------------------------------------------------------------------

# -- main target sources
set (SRC ecuyer_cote_mod.f90 random_beta_mod.f90 random_binomial_mod.f90 random_chisq_mod.f90 random_exponential_mod.f90 random_f_mod.f90 random_gamma_mod.f90 random_multinomial_mod.f90 random_multivariate_normal_mod.f90 random_nc_chisq_mod.f90 random_nc_f_mod.f90 random_negative_binomial_mod.f90 random_normal_mod.f90 random_permutation_mod.f90 random_poisson_mod.f90 random_standard_exponential_mod.f90 random_standard_gamma_mod.f90 random_standard_normal_mod.f90 random_standard_uniform_mod.f90 random_uniform_integer_mod.f90 random_uniform_mod.f90 user_set_generator.f90)

# -- main target, to be linked, not referred
add_library(randlib OBJECT ${SRC})