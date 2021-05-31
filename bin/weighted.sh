#!/bin/bash
options=$(getopt -o brg --long tmpdir:,maxrss:,maxtmp:,timeout:,task: -- "$@")
INPUT="$( cd "$( dirname "$1" )" && pwd )/$( basename "$1" )"
cd "$(dirname "$0")"

TMPDIR=""
MAXRSS=""
MAXTMP=""
TIMEOUT=""
TASK=""
eval set -- "$options"
while true; do
    case "$1" in
    --tmpdir)
        TMPDIR=$2
        shift 2
        ;;
    --maxrss)
        MAXRSS=$2
        shift 2
        ;;
    --maxtmp)
        MAXTMP=$2
        shift 2
        ;;
    --timeout)
        TIMEOUT=$2
        shift 2
        ;;
    --task)
		TASK=$2
        shift 2
		;;
    --)
        shift
        break
        ;;
    esac
done

CS=$((500*MAXRSS-500))

echo "c o Input file: $INPUT"
echo "c o Tmp dir: $TMPDIR"
echo "c o Max RAM: $MAXRSS"
echo "c o Max TMP: $MAXTMP"
echo "c o Timeout: $TIMEOUT"
echo "c o Task: $TASK"
echo "c o Cs: $CS"

./sharpSATW -decot 120 -decow 100 -tmpdir "$TMPDIR" -cs "$CS" "$INPUT"
