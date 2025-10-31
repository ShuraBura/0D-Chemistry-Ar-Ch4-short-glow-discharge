#!/bin/bash
echo "=================================="
echo "OVERNIGHT TEST STATUS CHECK"
echo "=================================="
echo "Current time: $(date)"
echo ""

# Check if running
if ps aux | grep -v grep | grep "test_single_overnight" > /dev/null; then
    echo "Status: ✅ RUNNING"
    echo ""
    ps aux | grep "test_single_overnight" | grep -v grep
    echo ""
    echo "Log tail:"
    tail -10 single_overnight.log
else
    echo "Status: ⏹️  STOPPED (completed or crashed)"
fi

echo ""
echo "=================================="
echo "RESULTS CHECK"
echo "=================================="

if [ -f single_test_result.json ]; then
    echo "✓ Results file exists!"
    echo ""
    python3 -c "
import json
with open('single_test_result.json') as f:
    data = json.load(f)
    if data.get('success'):
        print(f\"✓ TEST SUCCESSFUL!\")
        print(f\"Runtime: {data['timing']['integration_s']/60:.1f} minutes\")
        print(f\"\nResults:\")
        print(f\"  H:  {data['results']['H']:.2e} (target: {data['targets']['H']:.2e})\")
        print(f\"  CH: {data['results']['CH']:.2e} (target: {data['targets']['CH']:.2e})\")
        print(f\"  C2: {data['results']['C2']:.2e} (target: {data['targets']['C2']:.2e})\")
        print(f\"\nTotal error: {data['errors']['total']:.3f}\")
        print(f\"\nSweep estimates:\")
        runtime = data['timing']['integration_s']
        print(f\"  12 runs: {runtime*12/3600:.1f} hours\")
        print(f\"  36 runs: {runtime*36/3600:.1f} hours\")
    else:
        print(f\"✗ TEST FAILED\")
        print(f\"Message: {data.get('message', data.get('exception', 'Unknown error'))}\")
" 2>/dev/null || echo "Error reading results file"
else
    echo "No results file yet (still running or crashed)"
fi

echo ""
echo "=================================="
