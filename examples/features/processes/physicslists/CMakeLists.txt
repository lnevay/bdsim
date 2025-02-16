# high precision tests take a while to load neutron data, so labelled long
simple_testing(physics-all-particles          "--file=all_particles.gmad"          ${OVERLAP_CHECK})
simple_testing(physics-annihi-to-mumu         "--file=annihi_to_mumu.gmad"         ${OVERLAP_CHECK})
simple_testing(physics-charge-exchange        "--file=charge_exchange.gmad"        ${OVERLAP_CHECK})
simple_testing(physics-decay                  "--file=decay.gmad"                  ${OVERLAP_CHECK})
simple_testing(physics-decay-radioactive      "--file=decay_radioactive.gmad"      ${OVERLAP_CHECK})
simple_testing(physics-em                     "--file=em.gmad"                     ${OVERLAP_CHECK})
simple_testing(physics-em-extra               "--file=em_extra.gmad"               ${OVERLAP_CHECK})
simple_testing(physics-em-extra-options       "--file=em_extra_options.gmad"       ${OVERLAP_CHECK})
simple_testing(physics-em-penelope            "--file=em_penelope.gmad"            ${OVERLAP_CHECK})
simple_testing(physics-em-livermore           "--file=em_livermore.gmad"           ${OVERLAP_CHECK})
simple_testing(physics-em-livermore-polarised "--file=em_livermore_polarised.gmad" ${OVERLAP_CHECK})
simple_testing(physics-em-low-ep              "--file=em_low_ep.gmad"              ${OVERLAP_CHECK})
simple_testing(physics-em-ss                  "--file=em_ss.gmad"                  ${OVERLAP_CHECK})
simple_testing(physics-em-wvi                 "--file=em_wvi.gmad"                 ${OVERLAP_CHECK})
simple_testing(physics-em-1                   "--file=em_1.gmad"                   ${OVERLAP_CHECK})
simple_testing(physics-em-2                   "--file=em_2.gmad"                   ${OVERLAP_CHECK})
simple_testing(physics-em-3                   "--file=em_3.gmad"                   ${OVERLAP_CHECK})
simple_testing(physics-em-4                   "--file=em_4.gmad"                   ${OVERLAP_CHECK})
simple_testing(physics-ftfp-bert              "--file=ftfp_bert.gmad"              ${OVERLAP_CHECK})
simple_testing_long(physics-ftfpberthp        "--file=ftfp_bert_hp.gmad"           ${OVERLAP_CHECK})
simple_testing(physics-g4-ftfp-bert           "--file=g4FTFP_BERT.gmad"            ${OVERLAP_CHECK})
simple_testing(physics-g4-qbbc                "--file=g4QBBC.gmad"                 ${OVERLAP_CHECK})
simple_testing(physics-gamma-to-mumu          "--file=gamma_to_mumu.gmad"          ${OVERLAP_CHECK})
simple_testing(physics-hadronic-elastic       "--file=hadronic_elastic.gmad"       ${OVERLAP_CHECK})
simple_testing(physics-hadronic-elastic-d     "--file=hadronic_elastic_d.gmad"     ${OVERLAP_CHECK})
simple_testing(physics-hadronic-elastic-h     "--file=hadronic_elastic_h.gmad"     ${OVERLAP_CHECK})
simple_testing(physics-hadronic-elastic-hp    "--file=hadronic_elastic_hp.gmad"    ${OVERLAP_CHECK})
simple_testing(physics-hadronic-elastic-xs    "--file=hadronic_elastic_xs.gmad"    ${OVERLAP_CHECK})
simple_testing(physics-ion                    "--file=ion.gmad"                    ${OVERLAP_CHECK})
simple_testing(physics-ion-elastic            "--file=ion_elastic.gmad"            ${OVERLAP_CHECK})
simple_testing(physics-ion-elastic-qmd        "--file=ion_elastic_qmd.gmad"        ${OVERLAP_CHECK})
simple_testing(physics-ion-em-dissociation    "--file=ion_em_dissociation.gmad"    ${OVERLAP_CHECK})
simple_testing(physics-ionisation             "--file=ionisation.gmad"             ${OVERLAP_CHECK})
simple_testing(physics-muon                   "--file=muon.gmad"                   ${OVERLAP_CHECK})
simple_testing(physics-muon-inelastic         "--file=muon_inelastic.gmad"         ${OVERLAP_CHECK})
simple_testing(physics-neutron-tracking-cut   "--file=neutron_tracking_cut.gmad"   ${OVERLAP_CHECK})
simple_testing(physics-optical                "--file=optical.gmad"                ${OVERLAP_CHECK})
simple_testing(physics-qgsp-bert              "--file=qgsp_bert.gmad"              ${OVERLAP_CHECK})
simple_testing_long(physics-qgsp-bert-hp      "--file=qgsp_bert_hp.gmad"           ${OVERLAP_CHECK})
simple_testing(physics-qgsp-bic               "--file=qgsp_bic.gmad"               ${OVERLAP_CHECK})
simple_testing_long(physics-qgsp-bic-hp       "--file=qgsp_bic_hp.gmad"            ${OVERLAP_CHECK})
simple_testing_long(physics-shielding         "--file=shielding.gmad"              ${OVERLAP_CHECK})
simple_testing(physics-stopping               "--file=stopping.gmad"               ${OVERLAP_CHECK})
simple_testing(physics-synch-rad              "--file=synch_rad.gmad"              ${OVERLAP_CHECK})
simple_testing(physics-verbose                "--file=em_verbose.gmad"             ${OVERLAP_CHECK})

# slight modifications to geant4 reference physics lists
simple_testing(physics-g4-bdsim-range         "--file=g4FTFP_BERT_bdsim_range.gmad"      ${OVERLAP_CHECK})
simple_testing(physics-g4-bdsim-minimum-ek    "--file=g4FTFP_BERT_bdsim_minimum_ek.gmad" ${OVERLAP_CHECK})

if (G4_MINOR_VERSION GREATER 1)
  simple_testing(physics-alias-spindecay   "--file=alias_cerenkov.gmad"     ${OVERLAP_CHECK})
  simple_testing(physics-cherenkov         "--file=cherenkov.gmad"          ${OVERLAP_CHECK})
  simple_testing(physics-em-gs             "--file=em_gs.gmad"              ${OVERLAP_CHECK})
endif()

if (G4_MINOR_VERSION GREATER 2)
  simple_testing(physics-decay-spin        "--file=decay_spin.gmad"         ${OVERLAP_CHECK})
endif()

if (G4_MINOR_VERSION GREATER 3)
  simple_testing(physics-decay-muonic-atom "--file=decay_muonic_atom.gmad"  ${OVERLAP_CHECK})
  simple_testing(physics-dna               "--file=dna.gmad"                ${OVERLAP_CHECK})
  simple_testing(physics-dna-1             "--file=dna_1.gmad"              ${OVERLAP_CHECK})
  simple_testing(physics-dna-2             "--file=dna_2.gmad"              ${OVERLAP_CHECK})
  simple_testing(physics-dna-3             "--file=dna_3.gmad"              ${OVERLAP_CHECK})
  simple_testing(physics-dna-4             "--file=dna_4.gmad"              ${OVERLAP_CHECK})
  simple_testing(physics-dna-5             "--file=dna_5.gmad"              ${OVERLAP_CHECK})
  simple_testing(physics-dna-6             "--file=dna_6.gmad"              ${OVERLAP_CHECK})
  simple_testing(physics-dna-7             "--file=dna_7.gmad"              ${OVERLAP_CHECK})
  simple_testing(physics-radioactivation        "--file=radioactivation.gmad"        ${OVERLAP_CHECK})
endif()

# test invalid physics list
simple_fail(physics-invalid          "--file=gibberish.gmad"        ${OVERLAP_CHECK})

# these require the optional G4PARTICLEHPDATA set from the Geant4 website
if(DEFINED ENV{G4PARTICLEHPDATA})
  message(STATUS "Tests will use G4PARTICLEHPDATA")
  if (G4_MINOR_VERSION GREATER 2 OR (G4_MINOR_VERSION EQUAL 2 AND G4_PATCH_LEVEL GREATER 2))
    simple_testing_long(physics-ion-php         "--file=ion_php.gmad"               ${OVERLAP_CHECK})
    simple_testing_long(physics-alias-ionphp         "--file=alias_ionphp.gmad"          ${OVERLAP_CHECK})
  endif()
  simple_testing(physics-ion-inclxx             "--file=ion_inclxx.gmad"            ${OVERLAP_CHECK})
  simple_testing(physics-alias-ion-inclxx       "--file=alias_ioninclxx.gmad"       ${OVERLAP_CHECK})
  simple_testing(physics-ion-binary             "--file=ion_binary.gmad"            ${OVERLAP_CHECK})
  simple_testing(physics-alias-ion-binary       "--file=alias_ionbinary.gmad"       ${OVERLAP_CHECK})
else()
  if (G4_MINOR_VERSION GREATER 2 OR (G4_MINOR_VERSION EQUAL 2 AND G4_PATCH_LEVEL GREATER 2))
    simple_fail(physics-ion-php-no-hp           "--file=ion_php.gmad"               ${OVERLAP_CHECK})
    simple_fail(physics-alias-ionphp-no-hp      "--file=alias_ionphp.gmad"          ${OVERLAP_CHECK})
  endif()
  simple_fail(physics-ion-inclxx-no-hp          "--file=ion_inclxx.gmad"            ${OVERLAP_CHECK})
  simple_fail(physics-alias-ion-inclxx-no-hp    "--file=alias_ioninclxx.gmad"       ${OVERLAP_CHECK})
  simple_fail(physics-ion-binary-no-hp          "--file=ion_binary.gmad"            ${OVERLAP_CHECK})
  simple_fail(physics-alias-ion-binary-no-hp    "--file=alias_ionbinary.gmad"       ${OVERLAP_CHECK})
endif()

# these require the optional G4LENDDATA set from the Geant4 website
if(DEFINED ENV{G4LENDDATA})
  message(STATUS "Tests will use G4LENDDATA")
  simple_testing(physics-hadronic-elastic-lend  "--file=hadronic_elastic_lend.gmad" ${OVERLAP_CHECK})
  if (G4_MINOR_VERSION GREATER 3)
    simple_testing(physics-em-extra-lend        "--file=em_extra_lend.gmad"         ${OVERLAP_CHECK})
    simple_testing_long(physics-shielding-lend  "--file=shielding_lend.gmad"        ${OVERLAP_CHECK})
  endif()
else()
  simple_fail(physics-hadronic-elastic-lend-no-lend "--file=hadronic_elastic_lend.gmad" ${OVERLAP_CHECK})
  if (G4_MINOR_VERSION GREATER 3)
    simple_fail(physics-em-extra-lend-no-lend   "--file=em_extra_lend.gmad"         ${OVERLAP_CHECK})
    simple_fail(physics-shielding-lend-no-lend  "--file=shielding_lend.gmad"        ${OVERLAP_CHECK})
  endif()
endif()

if (G4_MINOR_VERSION GREATER 1 AND G4_MAJOR_VERSION GREATER 10)
# (REMOVED)  simple_testing(physics-xray-reflection  "--file=xray_reflection.gmad"   ${OVERLAP_CHECK})
endif()

#aliases
simple_testing(physics-alias-cerenkov    "--file=alias_cerenkov.gmad"     ${OVERLAP_CHECK})
simple_testing(physics-alias-emlow       "--file=alias_emlow.gmad"        ${OVERLAP_CHECK})
simple_testing_long(physics-alias-hadronic-hp "--file=alias_hadronic_hp.gmad"  ${OVERLAP_CHECK})
simple_testing(physics-alias-synchrad    "--file=alias_synchrad.gmad"     ${OVERLAP_CHECK})


#incompatible physics lists
simple_fail(physics-g4-incorrect      "--file=g4_incorrect.gmad"       ${OVERLAP_CHECK})
simple_fail(physics-incompatible-1    "--file=incompatible1.gmad"      ${OVERLAP_CHECK})
simple_fail(physics-incompatible-2    "--file=incompatible2.gmad"      ${OVERLAP_CHECK})
simple_fail(physics-incompatible-3    "--file=incompatible3.gmad"      ${OVERLAP_CHECK})

