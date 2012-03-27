#!/usr/bin/env bash

source env.sh

puts "Bash options:"
set -o errexit
set -o errtrace
set -o

puts "args: $@"


whereami
check-initial
run-md
setup-msmbuilder
assign
check-result
package
cleanup
