@test "basic functionality precond" {
    if [[ "$UG_CMAKE_PARALLEL" == "ON" ]]; then
        echo "run parallel test"
        mpirun --allow-run-as-root -n 1 ugshell -ex ../test/precond/test_case/precond.lua
    else
        echo "run serial test"
        ugshell -ex ../test/precond/test_case/precond.lua
    fi
}
@test "basic functionality Newton_Solver" {
    if [[ "$UG_CMAKE_PARALLEL" == "ON" ]]; then
        echo "run parallel test"
        mpirun --allow-run-as-root -n 1 ugshell -ex ../test/precond/test_case/precond_debug_test.lua
    else
        echo "run serial test"
        ugshell -ex ../test/precond/test_case/precond_debug_test.lua
    fi
}
@test "basic functionality laplace" {
    if [[ "$UG_CMAKE_PARALLEL" == "ON" ]]; then
        echo "run parallel test"
        mpirun --allow-run-as-root -n 1 ugshell -ex ../test/laplace/test_case/laplace_test.lua
    else
        echo "run serial test"
        ugshell -ex ../test/laplace/test_case/laplace_test.lua
    fi
}
@test "basic functionality poisson" {
    if [[ "$UG_CMAKE_PARALLEL" == "ON" ]]; then
        echo "run parallel test"
        mpirun --allow-run-as-root -n 1 ugshell -ex ../test/poisson/test_case/poisson_test.lua
    else
        echo "run serial test"
        ugshell -ex ../test/poisson/test_case/poisson_test.lua
    fi
}
@test "basic functionality henry_stat" {
    if [[ "$UG_CMAKE_PARALLEL" == "ON" ]]; then
        echo "run parallel test"
        mpirun --allow-run-as-root -n 1 ugshell -ex ../test/henry_stat/test_case/henry_stat_test.lua
    else
        echo "run serial test"
        ugshell -ex ../test/henry_stat/test_case/henry_stat_test.lua
    fi
}
@test "basic functionality cooler" {
    if [[ "$UG_CMAKE_PARALLEL" == "ON" ]]; then
        echo "run parallel test"
        mpirun --allow-run-as-root -n 1 ugshell -ex ../test/cooler/test_case/cooler_test.lua
    else
        echo "run serial test"
        ugshell -ex ../test/cooler/test_case/cooler_test.lua
    fi
}