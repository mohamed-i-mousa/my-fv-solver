#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
# SPDX-License-Identifier: Apache-2.0
#
# Regression check: run the solver on a fixed-iteration steady case and 
# compare the pinned total drag coefficient against a stored reference.
# A fixed iteration count (never convergence-gated) plus a globalSum-reduced
# force integral makes the value deterministic and rank-count-independent.
#
# Usage:
#   checkRegression.sh <solver> <caseFile> <forcesFile> <expectedCd> <relTol>
# -----------------------------------------------------------------------------

set -euo pipefail

solver="$1"
caseFile="$2"
forcesFile="$3"
expected="$4"
relTol="$5"

# The writers never create directories, so the output folder must exist first
outputDir="$(dirname "$forcesFile")"
solverLog="$outputDir/solver.log"

mkdir -p "$outputDir"
rm -f "$forcesFile" "$solverLog"

status=0
"$solver" "$caseFile" > "$solverLog" 2>&1 || status=$?

if [ "$status" -ne 0 ]; then
    echo "regression: solver exited $status (log: $solverLog)" >&2
    tail -n 40 "$solverLog" >&2
    exit 1
fi

if [ ! -f "$forcesFile" ]; then
    echo "regression: forces file not produced: $forcesFile" >&2
    tail -n 40 "$solverLog" >&2
    exit 1
fi

# Total drag coefficient is the fourth column of the line beginning with "Cd"
actual="$(awk '/^Cd /{print $4; exit}' "$forcesFile")"

if [ -z "$actual" ]; then
    echo "regression: could not parse Cd from $forcesFile" >&2
    exit 1
fi

# Relative-error comparison (awk handles the scientific-notation arithmetic)
awk -v a="$actual" -v e="$expected" -v t="$relTol" 'BEGIN {
    d = a - e; if (d < 0) d = -d;
    m = e;     if (m < 0) m = -m;
    rel = (m > 0) ? d / m : d;
    printf "regression Cd: actual=%s expected=%s relErr=%g tol=%s\n", a, e, rel, t;
    exit (rel <= t) ? 0 : 1;
}'
