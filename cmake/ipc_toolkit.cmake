include(FetchContent)
FetchContent_Declare(
    ipc-toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG "e01609eab1d06c30164d8d6b486f294535a088de"
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ipc-toolkit)