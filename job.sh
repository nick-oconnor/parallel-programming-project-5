#!/bin/bash
> results.txt
for NODES in 64 128 256; do
    for TASKS in 1 4 16 64; do
        for RANKS_PER_FILE in 1 4 8 32; do
            for COMPACT in 0 1; do
                THREADS=$((64/$TASKS))
                PARTITION="small"
                if [ $NODES -gt 64 ]; then
                    PARTITION="medium"
                fi
                echo "Running $TASKS tasks with $THREADS threads on $NODES nodes"
                srun --partition=$PARTITION --time=15 --overcommit --runjob-opts="--mapping TEDCBA" --nodes=$NODES --ntasks-per-node=$TASKS ./a.out $THREADS $RANKS_PER_FILE $COMPACT >> results.txt
                if [ $? -ne 0 ]
                then
                    echo "Error: Non-zero return code"
                    exit 0
                fi
                echo "Success"
            done
        done
    done
done
