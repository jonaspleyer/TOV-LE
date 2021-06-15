OUTPUT_FILE=SolverBIN
INPUT_FILE=Solver.cpp

OPTIONS=(-Wall -lstdc++ -lm -ldl)

INCLUDED_LIBRARIES=(libmongocxx)

echo "Compiling $INPUT_FILE"
echo "Include flags $(pkg-config --cflags --libs ${INCLUDED_LIBRARIES[@]})"
time gcc $INPUT_FILE -o $OUTPUT_FILE ${OPTIONS[@]} $(pkg-config --cflags --libs ${INCLUDED_LIBRARIES[@]})
echo "Written to $OUTPUT_FILE"
