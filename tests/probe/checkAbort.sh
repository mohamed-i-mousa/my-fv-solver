#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
# SPDX-License-Identifier: Apache-2.0
#
# Death check: the probe must print "FATAL ERROR" AND terminate non-zero.
#
# Usage:
#   checkAbort.sh <probe> <trigger>
# -----------------------------------------------------------------------------

set -uo pipefail

probe="$1"
trigger="$2"

status=0
output="$("$probe" "$trigger" 2>&1)" || status=$?

if [ "$status" -eq 0 ]; then
    echo "death: $trigger did not abort (exit 0)" >&2
    echo "$output" >&2
    exit 1
fi

if ! printf '%s\n' "$output" | grep -q 'FATAL ERROR'; then
    echo "death: $trigger aborted ($status) without a FATAL ERROR message" >&2
    echo "$output" >&2
    exit 1
fi

echo "death: $trigger aborted with status $status and a FATAL ERROR message"