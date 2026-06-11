#!/bin/bash

set -euo pipefail

# Regenerates the per-tool usage blocks in README.md from each subcommand's `--help` output.
# Each tool's block lives between a matching pair of marker comments:
#   <!-- start usage:TOOL --> ... <!-- end usage:TOOL -->
# GitHub Actions verifies the README is up to date, so run this whenever a tool's arguments or
# docstring change and commit the result.

tools=(demux shard subsample)

for tool in "${tools[@]}"; do
    usage="usage-${tool}.txt"
    {
        echo "<!-- start usage:${tool} -->"
        echo '````console'
        echo
        ./target/debug/fqtk "${tool}" --help | sed -e 's_^[ ]*$__g'
        echo '````'
        echo "<!-- end usage:${tool} -->"
    } > "${usage}"

    # Replace everything between (and including) this tool's markers with the freshly generated
    # block.  Outside the marker range lines pass through unchanged.
    sed -e "/<!-- start usage:${tool} -->/,/<!-- end usage:${tool} -->/!b" \
        -e "/<!-- end usage:${tool} -->/!d;r ${usage}" \
        -e d \
        README.md > README.md.new
    mv README.md.new README.md
    rm "${usage}"
done
