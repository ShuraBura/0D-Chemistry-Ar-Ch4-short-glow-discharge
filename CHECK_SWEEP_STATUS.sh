#!/bin/bash
# Quick status check for CSB parallel sweep

echo "=================================="
echo "CSB SWEEP STATUS CHECK"
echo "=================================="

# Check if sweep is running
PROC_COUNT=$(ps aux | grep -E "python3.*sweep_csb_parallel" | grep -v grep | wc -l)
echo "Running processes: $PROC_COUNT"

if [ $PROC_COUNT -gt 0 ]; then
    echo "Status: ✅ RUNNING"
    echo ""
    echo "Worker processes:"
    ps aux | grep -E "python3.*sweep_csb_parallel" | grep -v grep | head -5
else
    echo "Status: ⚠️  NOT RUNNING (completed or crashed)"
fi

echo ""
echo "=================================="
echo "RESULTS CHECK"
echo "=================================="

# Check for results file
if [ -f csb_sweep_results_*.json ]; then
    RESULTS_FILE=$(ls -t csb_sweep_results_*.json | head -1)
    echo "Results file: $RESULTS_FILE"
    
    # Count completed runs
    COMPLETED=$(python3 -c "import json; data=json.load(open('$RESULTS_FILE')); print(data.get('completed_runs', 0))" 2>/dev/null || echo "0")
    echo "Completed runs: $COMPLETED / 36"
    
    # Show file size
    ls -lh $RESULTS_FILE
else
    echo "No results file yet (still on first simulations)"
fi

echo ""
echo "=================================="
echo "To view full results when complete:"
echo "  python3 -m json.tool csb_sweep_results_*.json | less"
echo "=================================="
