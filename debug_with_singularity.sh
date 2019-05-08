#!/usr/bin/env bash
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

if ! [ "$SCRIPTPATH" -ef "$PWD" ]; then
    echo "This script just works if ran from $SCRIPTPATH"
    exit 1
fi

PANDORA_DEBUG_IMG="pandora_debug.img"
if [ -f "$PANDORA_DEBUG_IMG" ]; then
    echo "Using already created $PANDORA_DEBUG_IMG"
else
    echo "Creating $PANDORA_DEBUG_IMG..."
    sudo singularity build --writable pandora_debug.img Singularity_debug.pandora
fi

echo "You can proceed from here (e.g. mkdir build_debug && cd build_debug && cmake -DCMAKE_BUILD_TYPE=DEBUG .. && make -j4 && ctest -VV)"
echo "Use gdb for command-line debug or gdbserver for IDE debug (gdbserver :<port> /path/to/pandora <args>)"
sudo singularity exec -B $SCRIPTPATH:/pandora pandora_debug.img bash /pandora_start.sh