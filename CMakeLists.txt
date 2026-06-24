cmake_minimum_required(VERSION 3.12.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(FASTLCKIN VERSION 0.7.0 LANGUAGES CXX)

set(PROJECT_AUTHOR "Shujia Huang")
set(PROJECT_DESCRIPTION "fastlckin: Fast Maximum Likelihood Kinship Estimation from Low-Coverage Sequencing Data")

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h
)

add_compile_options(-O3 -fPIC)

# ── Static build option ─────────────────────────────────────────────
# Linux  : Bundle libstdc++/libgcc + htslib + compression libs statically;
#          glibc stays dynamic (portable to any glibc >= build-host version).
# macOS  : Best-effort static; only Apple system frameworks remain dynamic.
if(STATIC_BUILD)
    if(APPLE)
        message(STATUS "STATIC_BUILD=ON (macOS) — linking third-party libs statically via explicit .a paths")
    else()
        message(STATUS "STATIC_BUILD=ON (Linux/glibc) — bundling libstdc++/libgcc/htslib/compression libs, glibc dynamic")
    endif()
endif()

# For static builds:
# --disable-libcurl  : we only read local files, no remote URLs.
# --without-libdeflate: htslib's ./configure auto-detects libdeflate when
#                       present on the host (e.g. ubuntu-22.04 GH runners),
#                       and then emits libhts.a with libdeflate_* references.
#                       Since we do not link -ldeflate, we explicitly disable it.
#                       htslib falls back to zlib transparently with no
#                       functional loss.
if(STATIC_BUILD)
    set(HTSLIB_CONFIGURE_FLAGS --disable-libcurl --without-libdeflate)
else()
    set(HTSLIB_CONFIGURE_FLAGS "")
endif()

# ── htslib ──────────────────────────────────────────────────────────
include(ProcessorCount)
ProcessorCount(NPROC)
if(NPROC EQUAL 0)
    set(NPROC 2)
endif()

set(HTSLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/htslib)
set(HTSLIB_A   ${HTSLIB_DIR}/libhts.a)

# Stamp file encodes configure flags so toggling STATIC_BUILD triggers rebuild
set(HTSLIB_STAMP ${CMAKE_CURRENT_BINARY_DIR}/htslib_flags.stamp)
string(REPLACE ";" " " _HTSLIB_FLAGS_STR "${HTSLIB_CONFIGURE_FLAGS}")
file(WRITE ${HTSLIB_STAMP}.in "${_HTSLIB_FLAGS_STR}\n")
configure_file(${HTSLIB_STAMP}.in ${HTSLIB_STAMP} COPYONLY)

add_custom_command(
    OUTPUT  ${HTSLIB_A}
    COMMAND ${CMAKE_COMMAND} -E echo "Configuring htslib with: ${HTSLIB_CONFIGURE_FLAGS}"
    COMMAND make distclean || true
    COMMAND autoreconf -i
    COMMAND ./configure ${HTSLIB_CONFIGURE_FLAGS}
    COMMAND make -j${NPROC} libhts.a
    WORKING_DIRECTORY ${HTSLIB_DIR}
    DEPENDS ${HTSLIB_STAMP}
    COMMENT "Building htslib (parallel, j=${NPROC})"
    VERBATIM
)
add_custom_target(Buildhts DEPENDS ${HTSLIB_A})

include_directories(${HTSLIB_DIR})
include_directories(${HTSLIB_DIR}/htslib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# ── System libs ─────────────────────────────────────────────────────
set(COMMON_SYS_LIB pthread z bz2 m lzma curl)
if(NOT APPLE)
    list(APPEND COMMON_SYS_LIB ssl crypto)
endif()

set(HTSLIB ${HTSLIB_A})

# ── fastlckin executable ───────────────────────────────────────────
file(GLOB_RECURSE SRC_FL CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)
# Exclude test files
list(FILTER SRC_FL EXCLUDE REGEX ".*test.*")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
add_executable(fastlckin ${SRC_FL})
if(APPLE)
    target_link_options(fastlckin PRIVATE -Wl,-no_compact_unwind)
endif()
add_dependencies(fastlckin Buildhts)

# ── Link libraries (static vs dynamic) ─────────────────────────────
if(STATIC_BUILD AND APPLE)
    # macOS best-effort static: link zlib/bz2/lzma as .a from Homebrew;
    # Apple system frameworks stay dynamic.
    find_library(ZLIB_STATIC  libz.a
        HINTS /opt/homebrew/opt/zlib/lib /usr/local/opt/zlib/lib
        REQUIRED NO_DEFAULT_PATH)
    find_library(BZ2_STATIC   libbz2.a
        HINTS /opt/homebrew/opt/bzip2/lib /usr/local/opt/bzip2/lib
        REQUIRED NO_DEFAULT_PATH)
    find_library(LZMA_STATIC  liblzma.a
        HINTS /opt/homebrew/opt/xz/lib /usr/local/opt/xz/lib
        REQUIRED NO_DEFAULT_PATH)
    message(STATUS "Static libs: zlib=${ZLIB_STATIC} bz2=${BZ2_STATIC} lzma=${LZMA_STATIC}")
    message(STATUS "Dynamic (system): curl")
    target_link_libraries(fastlckin
        ${HTSLIB_A}
        curl
        ${LZMA_STATIC}
        ${BZ2_STATIC}
        ${ZLIB_STATIC}
        pthread m
    )
elseif(STATIC_BUILD AND NOT APPLE)
    # Linux partial-static: bundle libstdc++/libgcc + compression libs;
    # glibc stays dynamic for portability.
    find_library(ZLIB_STATIC  NAMES libz.a    REQUIRED)
    find_library(BZ2_STATIC   NAMES libbz2.a  REQUIRED)
    find_library(LZMA_STATIC  NAMES liblzma.a REQUIRED)
    message(STATUS "Linux static: zlib=${ZLIB_STATIC} bz2=${BZ2_STATIC} lzma=${LZMA_STATIC}")
    target_link_options(fastlckin PRIVATE -static-libstdc++ -static-libgcc)
    target_link_libraries(fastlckin
        ${HTSLIB_A}
        ${LZMA_STATIC}
        ${BZ2_STATIC}
        ${ZLIB_STATIC}
        pthread m
    )
else()
    target_link_libraries(fastlckin ${HTSLIB} ${COMMON_SYS_LIB})
endif()

# ── Unit tests ──────────────────────────────────────────────────────
enable_testing()
add_subdirectory(tests)
