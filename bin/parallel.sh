#!/bin/bash

set -eu

# Define the parse_arguments function
parallel() {

    local func
    local -a args
    local threads=1
    local conda=false
    local env=base

    # Parse arguments using getopt
    #ARGS=$(getopt -o "f:t:c:e:" --long "func:,threads:,conda:,env::" -n "parallel" -- "$@")

    # Check for getopt parsing errors (getopt)
    #if [ $? -ne 0 ]; then
    #    echo "Error: Failed to parse arguments." >&2
    #    exit 1
    #fi

    # Reset positional parameters to the parsed arguments (getopt)
    #eval set -- "$ARGS"

    # Process parsed options
    #while true; do
    while [[ $# -gt 0 ]]; do
        case "$1" in
             -f | --func)
                func="$2"
                shift 2
                ;;
            -t | --threads)
                threads="$2"
                shift 2
                ;;
            -c | --conda)
                conda="$2"
                shift 2
                ;;
            -e | --env)
                env="$2"
                shift 2
                ;;
            -a | --args)
                shift
                args=("$@")
                break
                ;;
            *)
                echo "$1 unrecognized option." >&2
                exit 1
                ;;
        esac
    done

    echo
    # Output parsed arguments
    echo "Command: $func"
    echo "Number of processes to be queued: ${#args[@]}"
    echo "threads: $threads"
    echo "conda: $conda"
    echo "env: $env"

    if $conda; then
        echo "Activating $env environment"
        eval "$(conda shell.bash hook)"
        source activate "$env"
    fi

    if ! command -v $func &> /dev/null; then
        echo "Command not found."
        exit 1
    fi

    if [ $threads -lt 1 ]; then
        echo "-t, --threads option not valid."
        exit 1
    fi

    echo; echo "Starting..."; echo 

    local -a pids
    
    for arg in "${args[@]}"; do
        queued=true
        while $queued; do
            if [ $threads -gt 0 ]; then
                    parallel=true
                else
                    parallel=false
            fi
            if $parallel; then
                $func $arg &
                pid=$!
                pids+=($pid)
                echo "Queued the following process:"
                echo "  $func $arg"
                echo "  PID: $pid"
                (( threads -= 1 ))
                queued=false
            fi  
            for i in ${!pids[@]}; do
                pid=${pids[$i]}
                if ! ps -p $pid > /dev/null; then
                    unset pids[$i]
                    (( threads += 1 ))
                    echo; echo "*******************"
                    echo "PID $pid finished."
                    echo "*******************"; echo
                fi
            done
        done
    done
    
    wait
}

# Call the function with command-line arguments
# parallel "$@"