#!/bin/bash
ABS_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$ABS_SCRIPT_DIR"

#python3 -m pip install -U --no-deps .


python3 -m pip install -e  .
