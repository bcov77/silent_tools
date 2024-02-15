#!/bin/bash

unset SILENT_PARAMS
unset SILENT_AT
unset SILENT_FLAGS
unset SILENT_FLAGS_NO_J


POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -p)
      SILENT_PARAMS="$SILENT_PARAMS -extra_res_fa $2"
      SILENT_FLAGS="$SILENT_FLAGS -p $2"
      SILENT_FLAGS_NO_J="$SILENT_FLAGS_NO_J -p $2"
      shift # past argument
      shift # past value
      ;;
    -j)
      SILENT_J="$2"
      SILENT_FLAGS="$SILENT_FLAGS -j$2"
      shift # past argument
      shift # past value
      ;;
    -j*)
      SILENT_J="${1:2:${#1}-2}"
      SILENT_FLAGS="$SILENT_FLAGS $1"
      shift # past value
      ;;
    @*)
      SILENT_AT="$SILENT_AT $1"
      SILENT_FLAGS="$SILENT_FLAGS $1"
      SILENT_FLAGS_NO_J="$SILENT_FLAGS_NO_J $1"
      shift # past value
      ;;
    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done

export SILENT_PARAMS
export SILENT_AT
export SILENT_FLAGS
export SILENT_FLAGS_NO_J

set -- "${POSITIONAL[@]}" # restore positional parameters

