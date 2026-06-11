#!/bin/bash

function banner() {
    echo
    echo "================================================================================"
    echo "$*"
    echo "================================================================================"
    echo
}

#####################################################################
# Takes two parameters, a "name" and a "command".
# Runs the command and prints out whether it succeeded or failed, and
# also tracks a list of failed steps in $failures.
#####################################################################
function run() {
    local name=$1
    local cmd=$2

    banner "Running $name [$cmd]"
    set +e
    $cmd
    exit_code=$?
    set -e

    if [[ $exit_code == 0 ]]; then
        echo "Passed $name: [$cmd]"
    else
        echo "Failed $name: [$cmd]"
        if [ -z "$failures" ]; then
            failures="$failures $name"
        else
            failures="$failures, $name"
        fi
    fi
}

#####################################################################
# Regenerates the `fqtk demux` usage embedded in the README and fails
# if it leaves the working tree dirty (i.e. the README was out of
# date). Mirrors the README sync check enforced in CI.
#####################################################################
function check_docs() {
    cargo build --locked || return 1
    bash .github/scripts/update-docs.sh || return 1
    git diff --exit-code -- README.md
}

run "Formatting"  "cargo fmt --check --all"
run "Clippy"      "cargo clippy --locked --all-features --all-targets -- -D warnings"
run "Unit Tests"  "cargo test --locked"
run "Docs"        "check_docs"

if [ -z "$failures" ]; then
    banner "Checks Passed"
else
    banner "Checks Failed with failures in: $failures"
    exit 1
fi
