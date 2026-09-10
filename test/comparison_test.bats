#!/usr/bin/env bats

# Run before every test.
setup() {
    # All Lua scripts expect to be executed relative to test/cases.
    cd "$BATS_TEST_DIRNAME/cases"

    # Create directories for generated test results.
    mkdir -p "$BATS_TEST_DIRNAME/output/automated"
    mkdir -p "$BATS_TEST_DIRNAME/output/comparison"
}


# Run an original Lua script and check that it succeeds.
run_original() {
    local script="$1"

    run ugshell -ex "$BATS_TEST_DIRNAME/original/$script"
    [ "$status" -eq 0 ]
}


# Run a plugin_Util test script and check that it succeeds.
run_plugin_test() {
    local script="$1"

    run ugshell -ex "$BATS_TEST_DIRNAME/cases/$script"
    [ "$status" -eq 0 ]
}


@test "MGStats custom: original and plugin_Util produce identical solution" {

    # Run both implementations.
    run_original "laplace_mgstats_custom.lua"

    # Save the freshly generated original solution.
    mv "laplace_mgstats_custom_3d.vec" \
       "$BATS_TEST_DIRNAME/output/comparison/laplace_mgstats_custom_original.vec"

    run_plugin_test "laplace_mgstats_custom_test.lua"

    # Compare the two freshly generated solutions exactly.
    diff \
        "$BATS_TEST_DIRNAME/output/comparison/laplace_mgstats_custom_original.vec" \
        "$BATS_TEST_DIRNAME/output/automated/laplace_mgstats_custom_test_3d.vec"
}


@test "MGStats standard: original and plugin_Util produce identical solution" {

    run_original "laplace_mgstats_standard.lua"

    mv "laplace_mgstats_standard_3d.vec" \
       "$BATS_TEST_DIRNAME/output/comparison/laplace_mgstats_standard_original.vec"

    run_plugin_test "laplace_mgstats_standard_test.lua"

    diff \
        "$BATS_TEST_DIRNAME/output/comparison/laplace_mgstats_standard_original.vec" \
        "$BATS_TEST_DIRNAME/output/automated/laplace_mgstats_standard_test_3d.vec"
}


@test "Poisson: original and plugin_Util produce identical solution" {

    run_original "poisson.lua"

    mv "poisson_2d.vec" \
       "$BATS_TEST_DIRNAME/output/comparison/poisson_original.vec"

    run_plugin_test "poisson_test.lua"

    diff \
        "$BATS_TEST_DIRNAME/output/comparison/poisson_original.vec" \
        "$BATS_TEST_DIRNAME/output/automated/poisson_test_2d.vec"
}


@test "Cooler: original and plugin_Util produce identical solution" {

    run_original "cooler.lua"

    run_plugin_test "cooler_test.lua"

    diff \
        "sol_cooler_t0020.vtu" \
        "$BATS_TEST_DIRNAME/output/automated/cooler_test_t0020.vtu"
}

@test "Henry: original and plugin_Util produce identical solution" {

    run_original "henry_stat.lua"

    run_plugin_test "henry_stat_test.lua"

    diff \
        "HenryStat-Boussinesq-GL7.vtu" \
        "$BATS_TEST_DIRNAME/output/automated/HenryStat_test-Boussinesq-GL7.vtu"
}
