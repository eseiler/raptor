# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

project (update_cpm_packages LANGUAGES NONE)

# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------

function (replace_if_unique_match old_string new_string output_var input_string)
    string (FIND ${input_string} ${old_string} FIND1)
    string (FIND ${input_string} ${old_string} FIND2 REVERSE)
    if (FIND1 STREQUAL FIND2)
        string (REPLACE "${old_string}" "${new_string}" ${output_var} ${input_string})
    else ()
        message (WARNING "UPDATE FAILED: Found string \"${old_string}\" twice.")
    endif ()
    return (PROPAGATE ${output_var})
endfunction ()

function (check_if_commit repo tag result)
    set (${result} FALSE)
    if (tag)
        string (LENGTH ${tag} TAG_LENGTH)
        if (TAG_LENGTH EQUAL 40)
            set (${result} TRUE)
        else ()
            execute_process (COMMAND "gh" "api" "/repos/${repo}/branches/${tag}"
                             RESULT_VARIABLE IS_BRANCH
                             OUTPUT_QUIET ERROR_QUIET
            )
            if (NOT IS_BRANCH EQUAL 0)
                set (${result} TRUE)
            endif ()
        endif ()
    endif ()
    return (PROPAGATE ${result})
endfunction ()

# ----------------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------------

if (NOT PACKAGE_LOCK)
    message (FATAL_ERROR "Please pass the path to the CPM package lock file: -DPACKAGE_LOCK=<path>")
endif ()
if (NOT EXISTS "${PACKAGE_LOCK}")
    message (FATAL_ERROR "File does not exist: \"${PACKAGE_LOCK}\"")
endif ()
if (NOT IS_READABLE "${PACKAGE_LOCK}")
    message (FATAL_ERROR "File exists, but cannot be read: \"${PACKAGE_LOCK}\"")
endif ()
if (NOT IS_WRITABLE "${PACKAGE_LOCK}")
    message (FATAL_ERROR "File exists, but cannot be written to: \"${PACKAGE_LOCK}\"")
endif ()

file (DOWNLOAD https://github.com/cpm-cmake/CPM.cmake/releases/latest/download/CPM.cmake
      ${CMAKE_CURRENT_BINARY_DIR}/CPM.cmake
)
include (${CMAKE_CURRENT_BINARY_DIR}/CPM.cmake)

set (PACKAGE_LOCK "${CMAKE_CURRENT_LIST_DIR}/../cmake/package-lock.cmake")
file (READ "${PACKAGE_LOCK}" CONTENT)
CPMUsePackageLock ("${PACKAGE_LOCK}")

set (PACKAGE_LIST "")
# cmake-format: off
block (SCOPE_FOR VARIABLES PROPAGATE PACKAGE_LIST)
    set (TMP_PACKAGE_LOCK "${CMAKE_CURRENT_BINARY_DIR}/package-lock.cmake")
    configure_file ("${PACKAGE_LOCK}" "${TMP_PACKAGE_LOCK}")
    file (READ "${TMP_PACKAGE_LOCK}" STR)
    string (REGEX MATCHALL "\\(([^()]+)\\)" PACKAGE_LIST ${STR})
endblock ()
# cmake-format: on

# ----------------------------------------------------------------------------
# Do the thing
# ----------------------------------------------------------------------------

foreach (PACKAGE ${PACKAGE_LIST})
    set (REPO "")
    set (TAG "")
    set (VERSION "")
    set (BRANCH "")
    set (TAG_IS_COMMIT FALSE)

    if (PACKAGE MATCHES "GITHUB_REPOSITORY +([^ \t\r\n]+)")
        set (REPO "${CMAKE_MATCH_1}")

        if (PACKAGE MATCHES "GIT_TAG +([^ \t\r\n]+)[^#\r\n]*([ #]*)([^ \t\r\n]*)")
            set (TAG "${CMAKE_MATCH_1}")
            set (BRANCH "${CMAKE_MATCH_3}")
        elseif (PACKAGE MATCHES "VERSION +([^ \t\r\n]+)")
            set (VERSION "${CMAKE_MATCH_1}")
        endif ()

        check_if_commit (${REPO} "${TAG}" TAG_IS_COMMIT)
    elseif (PACKAGE MATCHES "URL +([^ \t\r\n]+)([ #]*)([^ \t\r\n]*)")
        set (URL "${CMAKE_MATCH_1}")
        set (BRANCH "${CMAKE_MATCH_3}")

        if (URL MATCHES "https://github.com/([^/]+)/([^/]+)/archive/(.*)")
            set (REPO "${CMAKE_MATCH_1}/${CMAKE_MATCH_2}")
            set (SUBURL "${CMAKE_MATCH_3}")
            if (SUBURL MATCHES "refs/heads/(.*)\.(tar\.gz|zip)")
                set (BRANCH "${CMAKE_MATCH_1}")
            elseif (SUBURL MATCHES "refs/tags/(.*)\.(tar\.gz|zip)")
                set (TAG "${CMAKE_MATCH_1}")
            elseif (SUBURL MATCHES "([a-zA-Z0-9]*)\.(tar\.gz|zip)")
                set (TAG "${CMAKE_MATCH_1}")
                set (TAG_IS_COMMIT TRUE)
            endif ()
        endif ()
    endif ()

    if (VERSION)
        execute_process (COMMAND "gh" "api" "/repos/${REPO}/releases/latest" "-q" ".tag_name"
                         OUTPUT_VARIABLE LATEST_RELEASE
                         OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY
        )
        if (LATEST_RELEASE MATCHES "[vV]*(.*)")
            set (LATEST_RELEASE "${CMAKE_MATCH_1}")
        else ()
            message (NOTICE "-- Skipping ${REPO}\n" #
                            "    Unknown version string for latest release: \"${LATEST_RELEASE}\""
            )
        endif ()
        if (VERSION STREQUAL LATEST_RELEASE)
            message (STATUS "${REPO} is up-to-date\n" #
                            "    Current release: \"${VERSION}\"\n" #
                            "    Latest release : \"${LATEST_RELEASE}\""
            )
        else ()
            message (STATUS "Updating ${REPO}\n" #
                            "    Current release: \"${VERSION}\"\n" #
                            "    Latest release : \"${LATEST_RELEASE}\""
            )
            replace_if_unique_match ("${VERSION}" "${LATEST_RELEASE}" CONTENT ${CONTENT})
        endif ()
    elseif (BRANCH)
        execute_process (COMMAND "gh" "api" "/repos/${REPO}/branches/${BRANCH}" "-q" ".commit .sha"
                         OUTPUT_VARIABLE LATEST_COMMIT
                         OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY
        )
        if (NOT TAG)
            message (NOTICE "-- Skipping ${REPO}\n" #
                            "    Branch: \"${BRANCH}\"\n" #
                            "    URL points to a branch archive. Nothing to update."
            )
        elseif (TAG STREQUAL LATEST_COMMIT)
            message (STATUS "${REPO} is up-to-date\n" #
                            "    Branch: \"${BRANCH}\"\n" #
                            "    Current commit: \"${TAG}\"\n" #
                            "    Latest commit : \"${LATEST_COMMIT}\""
            )
        else ()
            message (STATUS "Updating ${REPO}\n" #
                            "    Branch: \"${BRANCH}\"\n" #
                            "    Current commit: \"${TAG}\"\n" #
                            "    Latest commit : \"${LATEST_COMMIT}\""
            )
            replace_if_unique_match ("${TAG}" "${LATEST_COMMIT}" CONTENT ${CONTENT})
        endif ()
    elseif (TAG AND NOT TAG_IS_COMMIT)
        message (NOTICE "-- Skipping ${REPO}\n" #
                        "    Tag: \"${TAG}\"\n" #
                        "    Current tag/URL points to a tag.\n" #
                        "    It is not generally possible to tell what would be considered a newer tag."
        )
    elseif (TAG AND TAG_IS_COMMIT)
        message (NOTICE "-- Skipping ${REPO}\n" #
                        "    Tag: \"${TAG}\"\n" #
                        "    Current tag/URL points to a commit.\n" #
                        "    Append \"# <branch>\" to the line to configure which branch should be queried for the newest commit.\n" #
                        "    For example, \"# main\"."
        )
    endif ()
endforeach ()

file (WRITE "${PACKAGE_LOCK}" ${CONTENT})
