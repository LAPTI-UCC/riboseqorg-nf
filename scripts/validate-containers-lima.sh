#!/usr/bin/env bash
# Wrapper script for running container validation with Lima/Apptainer on Mac

set -e

echo "========================================="
echo "Container Validation with Lima"
echo "========================================="
echo ""

# Check if limactl is installed
if ! command -v limactl &>/dev/null; then
    echo "❌ Error: limactl not found"
    echo ""
    echo "Lima is required to run Singularity on Mac."
    echo "Install it with: brew install lima"
    echo ""
    echo "Then run this script again."
    exit 1
fi

# Check if apptainer VM exists
if ! limactl list | grep -q apptainer; then
    echo "📦 Creating apptainer VM with Rosetta (this may take a few minutes)..."
    echo ""
    echo "Note: Using Rosetta to enable x86_64 container compatibility on Mac"
    limactl start --name=apptainer --vm-type=vz --rosetta template:apptainer
    echo ""
    echo "✅ Apptainer VM created successfully!"
    echo ""
fi

# Start the VM if not running
echo "🚀 Starting Lima apptainer VM..."
limactl start apptainer 2>/dev/null || true

# Wait a moment for VM to be ready
sleep 2

# Get the current directory
CURRENT_DIR=$(pwd)

echo "🧪 Running validation script inside Lima VM..."
echo ""

# Run the validation script inside the VM
# The VM has access to your home directory at the same path
limactl shell apptainer bash -c "cd '$CURRENT_DIR' && ./scripts/validate-containers.sh"

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Validation completed successfully!"
else
    echo "❌ Validation failed with exit code $EXIT_CODE"
fi

echo ""
echo "💡 Tip: To save resources, stop the VM when not in use:"
echo "   limactl stop apptainer"
echo ""
echo "   To start it again later:"
echo "   limactl start apptainer"

exit $EXIT_CODE
