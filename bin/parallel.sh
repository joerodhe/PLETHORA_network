#!/bin/bash

###############################################################################
# Description:
#   This script defines a function named `parallel()` that allows you to run 
#   multiple processes in parallel using a custom argument parsing mechanism. 
#   It optionally activates a conda environment, queues up a specified command 
#   (`$func`) with an array of arguments, and manages concurrency based on 
#   a user-defined number of threads.
###############################################################################

set -eu  # Exit on error (-e) and treat unset variables as errors (-u)

###############################################################################
# Function: parallel()
# --------------------
# Description:
#   Parses command-line options to determine a command ($func), number of threads
#   ($threads), whether conda activation is requested ($conda), the environment
#   to activate ($env), and the set of arguments ($args). It then executes
#   `$func` in parallel for each item in $args, respecting the max concurrency.
#
# Options recognized:
#   -f, --func    : The command or function to run.
#   -t, --threads : Number of parallel threads/processes allowed.
#   -c, --conda   : Boolean (true/false) indicating if a conda environment 
#                   should be activated.
#   -e, --env     : The name of the conda environment to activate (default: "base").
#   -a, --args    : The remainder of arguments, forming an array to be processed.
#
# Usage Example:
#   parallel -f echo -t 3 -c true -e myenv -a "Hello" "World" "Test"
#
# Behavior:
#   1) Parses the options and saves them to local variables.
#   2) If conda is true, activates the specified environment.
#   3) Checks that $func is in the PATH, and that $threads is valid.
#   4) Queues $func calls in parallel, up to $threads concurrently.
###############################################################################
parallel() {

    # Declare local variables
    local func
    local -a args
    local threads=1
    local conda=false
    local env=base

    # We skip using getopt for parsing to keep the example straightforward:
    # Instead, we manually parse arguments in a while-loop.
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -f|--func)
                func="$2"
                shift 2
                ;;
            -t|--threads)
                threads="$2"
                shift 2
                ;;
            -c|--conda)
                conda="$2"
                shift 2
                ;;
            -e|--env)
                env="$2"
                shift 2
                ;;
            -a|--args)
                shift
                args=("$@")  # Collect the remainder of arguments
                break        # Exit the loop
                ;;
            *)
                echo "$1 unrecognized option." >&2
                exit 1
                ;;
        esac
    done

    echo
    # Display the parsed configuration
    echo "Command: $func"
    echo "Number of processes to be queued: ${#args[@]}"
    echo "threads: $threads"
    echo "conda: $conda"
    echo "env: $env"

    ###########################################################################
    # Optional: Conda environment activation
    ###########################################################################
    if $conda; then
        echo "Activating $env environment"
        # Evaluate shell hook and activate environment
        eval "$(conda shell.bash hook)"
        source activate "$env"
    fi

    ###########################################################################
    # Validation checks
    ###########################################################################
    # Ensure the provided command ($func) is in the PATH
    if ! command -v "$func" &> /dev/null; then
        echo "Command not found: $func"
        exit 1
    fi

    # Ensure number of threads is valid
    if [[ $threads -lt 1 ]]; then
        echo "-t, --threads option not valid: $threads"
        exit 1
    fi

    echo
    echo "Starting..."
    echo

    # Array to store process IDs (PIDs) of running tasks
    local -a pids

    ###########################################################################
    # Main loop: Spawn processes in parallel, respecting the $threads limit
    ###########################################################################
    for arg in "${args[@]}"; do
        queued=true
        while $queued; do
            # Check if we still have a free thread
            if [[ $threads -gt 0 ]]; then
                parallel=true
            else
                parallel=false
            fi

            # If a free thread is available, run the command in background
            if $parallel; then
                "$func" "$arg" &
                pid=$!
                pids+=($pid)

                echo "Queued the following process:"
                echo "  $func $arg"
                echo "  PID: $pid"

                # Decrement the thread counter to reflect one busy slot
                (( threads -= 1 ))
                queued=false
            fi  

            # Check running processes, reclaim threads for any processes that ended
            for i in "${!pids[@]}"; do
                pid=${pids[$i]}
                if ! ps -p "$pid" > /dev/null; then
                    unset pids[$i]  # Remove from array
                    (( threads += 1 ))
                    echo
                    echo "*******************"
                    echo "PID $pid finished."
                    echo "*******************"
                    echo
                fi
            done
        done
    done
    
    ###########################################################################
    # Wait for all background processes to finish
    ###########################################################################
    wait
}

###############################################################################
# If needed, invoke the function with all passed command-line arguments
# parallel "$@"
#
# Example usage:
#   ./script.sh -f echo -t 2 -c false -a "Hello" "World" "123"
###############################################################################
