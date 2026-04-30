#!/usr/bin/env bash

set -e  # keep errors, but remove -u for now

LOG_DIR="logs"
mkdir -p "$LOG_DIR"

NP=16 # Number of computer cores

echo "========================================"
echo "Stage 0: Setup + Decompose"
echo "========================================"

./cleanCase.sh || true

echo "Running blockMesh..."
blockMesh | tee "$LOG_DIR/blockMesh.log"

echo "Running setFields..."
setFields | tee "$LOG_DIR/setFields.log"

echo "Running decomposePar..."
decomposePar -force | tee "$LOG_DIR/decompose.log"

echo "========================================"
echo "Stage 1: 0 -> 10 sec (fine logging)"
echo "========================================"

foamDictionary system/controlDict -entry startFrom -set startTime
foamDictionary system/controlDict -entry startTime -set 0
foamDictionary system/controlDict -entry endTime -set 10
foamDictionary system/controlDict -entry writeInterval -set 0.005

echo "Running Stage 1 solver..."

mpirun -np $NP foamRun -solver incompressibleVoF -parallel \
| tee "$LOG_DIR/stage1.log" \
| awk '
BEGIN {
    t_start = 0;
    t_end = 10;
    weight = 0.6;

    wall_start = systime();
}
/^Time =/ {
    t = $3;
    gsub("s","",t);

    wall_now = systime();
    wall_elapsed = wall_now - wall_start;

    sim_elapsed = t - t_start;

    if (sim_elapsed > 0) {
        speed = wall_elapsed / sim_elapsed;
    } else {
        speed = 0;
    }

    stage_progress = sim_elapsed / (t_end - t_start);
    global_progress = weight * stage_progress;

    remaining_sim = (t_end - t);
    eta = remaining_sim * speed;

    printf("\r[Global %5.1f%%] [Stage1 %5.1f%%] Time = %.3f | Wall = %ds | Speed = %.2f s/sim-s | ETA = %ds",
           100*global_progress,
           100*stage_progress,
           t,
           wall_elapsed,
           speed,
           eta);

    fflush();
}
END { print "" }'

echo "========================================"
echo "Stage 2: 10 sec -> end (coarse logging)"
echo "========================================"

foamDictionary system/controlDict -entry startFrom -set latestTime
foamDictionary system/controlDict -entry endTime -set 100
foamDictionary system/controlDict -entry writeInterval -set 0.02
foamDictionary system/controlDict -entry "functions.forces.writeInterval" -set 0.02
foamDictionary system/controlDict -entry "functions.liquidCOM.writeInterval" -set 0.02
foamDictionary system/controlDict -entry "functions.liquidCOM.executeInterval" -set 0.02

echo "Running Stage 2 solver..."

mpirun -np $NP foamRun -solver incompressibleVoF -parallel \
| tee "$LOG_DIR/stage2.log" \
| awk '
BEGIN {
    t_start = 10;
    t_end = 100;

    weight_prev = 0.6;
    weight = 0.4;

    wall_start = systime();
}
/^Time =/ {
    t = $3;
    gsub("s","",t);

    wall_now = systime();
    wall_elapsed = wall_now - wall_start;

    sim_elapsed = t - t_start;

    if (sim_elapsed > 0) {
        speed = wall_elapsed / sim_elapsed;
    } else {
        speed = 0;
    }

    stage_progress = sim_elapsed / (t_end - t_start);
    global_progress = weight_prev + weight * stage_progress;

    remaining_sim = (t_end - t);
    eta = remaining_sim * speed;

    printf("\r[Global %5.1f%%] [Stage2 %5.1f%%] Time = %.3f | Wall = %ds | Speed = %.2f s/sim-s | ETA = %ds",
           100*global_progress,
           100*stage_progress,
           t,
           wall_elapsed,
           speed,
           eta);

    fflush();
}
END { print "" }'

echo "========================================"
echo "Merging data"
echo "========================================"

python3 Merge_data.py | tee "$LOG_DIR/merge.log"

echo "DONE"
