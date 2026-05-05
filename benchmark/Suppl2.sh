#!/bin/bash

INPUT_FILES="*.fasta"
OUTPUT_CSV="benchmark_results.csv"
TMP_TIME="tmp_time.txt"

THREADS=(1 2 4 8 16 32)
KMERS=(21 26 31)
PARABOLA_S=(5000 1000 100 10)
MASH_S=(1000 10000 30000 100000)

echo "Tool,Threads,K-mer,SketchSize,Elapsed_Time(s),MaxResident_Memory(KB)" > "$OUTPUT_CSV"

for p in "${THREADS[@]}"; do
    for k in "${KMERS[@]}"; do
        for s in "${PARABOLA_S[@]}"; do
            echo -n "Running Parabola | -p $p | -k $k | -s $s ... "

            /usr/bin/time -f "%e,%M" -o "$TMP_TIME" \
                ../parabola/parabola sketch -p "$p" -k "$k" -s "$s" $INPUT_FILES > /dev/null 2>&1

            RESULT=$(cat "$TMP_TIME")

            echo "parabola,$p,$k,$s,$RESULT" >> "$OUTPUT_CSV"
        done
    done
done

for p in "${THREADS[@]}"; do
        for s in "${PARABOLA_S[@]}"; do
            echo -n "Running Skani | -t $p | -c $s ... "
	    rm -rf tmp/
            
            /usr/bin/time -f "%e,%M" -o "$TMP_TIME" \
                skani sketch -t "$p" -c "$s" $INPUT_FILES  -o tmp > /dev/null 2>&1
            
            RESULT=$(cat "$TMP_TIME")
            
            echo "skani,$p,-,$s,$RESULT" >> "$OUTPUT_CSV"
        done
done


for p in "${THREADS[@]}"; do
    for k in "${KMERS[@]}"; do
        for s in "${MASH_S[@]}"; do
            echo -n "Running Mash | -p $p | -k $k | -s $s ... "

            /usr/bin/time -f "%e,%M" -o "$TMP_TIME" \
                mash sketch -p "$p" -k "$k" -s "$s" $INPUT_FILES > /dev/null 2>&1

            RESULT=$(cat "$TMP_TIME")

            echo "mash,$p,$k,$s,$RESULT" >> "$OUTPUT_CSV"
        done
    done
done

rm -f "$TMP_TIME"