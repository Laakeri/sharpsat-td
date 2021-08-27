#!/bin/bash
rm -rf build/
mkdir build
cd build
if [ "$1" = "static" ]; then
	cmake -DSTATIC=ON ..
else
	cmake ..
fi
make
cd ..
cp build/sharpSAT bin/sharpSAT
cp build/flow_cutter_pace17 bin/flow_cutter_pace17