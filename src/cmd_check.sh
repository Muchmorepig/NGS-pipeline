#!/bin/bash
command_exists() {
    command -v $1 >/dev/null 2>&1 || {
        echo >&2 "Require '$1' but it's not detected.  Please load '$1'"
        exit 1
    }
}
