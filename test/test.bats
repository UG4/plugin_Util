# This function runs before every test.
setup() {
    # Switch to the directory containing the Lua test cases.
    # This ensures that all relative paths inside the Lua files are resolved correctly.
    cd "$BATS_TEST_DIRNAME/cases"

    # Ensure that the directory for automatically generated files exists.
    mkdir -p "$BATS_TEST_DIRNAME/output/automated"
}

# This helper executes one Lua test in serial or parallel mode.
run_ugshell_test() {
    # The first argument passed to this function is the Lua filename.
    local lua_script="$1"

    # Use MPI when UG4 was compiled with parallel support.
    if [[ "$UG_CMAKE_PARALLEL" == "ON" ]]; then
        echo "run parallel test: $lua_script"
        mpirun --allow-run-as-root -n 1 ugshell -ex "$lua_script"
    else
        echo "run serial test: $lua_script"
        ugshell -ex "$lua_script"
    fi
}

@test "basic functionality preconditioner" {
    run_ugshell_test "precond_test.lua"
}

@test "MGStats with standard configuration" {
    run_ugshell_test "laplace_mgstats_standard_test.lua"
}

@test "MGStats with custom configuration" {
    run_ugshell_test "laplace_mgstats_custom_test.lua"
}

@test "basic functionality poisson" {
    run_ugshell_test "poisson_test.lua"
}

@test "basic functionality henry_stat" {
    run_ugshell_test "henry_stat_test.lua"
}

@test "basic functionality cooler" {
    run_ugshell_test "cooler_test.lua"
}
