include(FetchContent)
FetchContent_Declare(
    ipc-toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG "913f3886a7c1c41e94ef6dcc2bda6215da92edbf"
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ipc-toolkit)