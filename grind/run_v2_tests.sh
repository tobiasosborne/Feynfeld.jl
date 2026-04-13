#!/bin/bash
# Sequential v2 test runner with incremental output
cd /home/tobiasosborne/Projects/Feynfeld.jl
exec > grind/v2_test_results.txt 2>&1
echo "=== v2 test suite (post-Phase-12d fix) ==="
for f in test/v2/test_*.jl; do
    line=$(julia --project=. "$f" 2>&1 | grep -E '^(Test Summary|ERROR|Some tests)' | tr '\n' ' ')
    echo "$(basename $f): $line"
done
echo "=== DONE ==="
