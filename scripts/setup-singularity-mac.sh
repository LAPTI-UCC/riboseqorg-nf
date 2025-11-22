#!/usr/bin/env bash
# Setup script for Singularity/Apptainer on Mac using Lima

echo "========================================="
echo "Singularity/Apptainer Setup for Mac"
echo "========================================="
echo ""

# Check if limactl is installed
if ! command -v limactl &>/dev/null; then
    echo "❌ Lima not found"
    echo ""
    echo "Installing Lima via Homebrew..."
    if command -v brew &>/dev/null; then
        brew install lima
        echo "✅ Lima installed!"
    else
        echo "❌ Homebrew not found. Please install Homebrew first:"
        echo "   /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
        exit 1
    fi
    echo ""
fi

# Check if apptainer VM exists
if limactl list 2>/dev/null | grep -q apptainer; then
    echo "✅ Apptainer VM already exists"

    # Check if it's running
    if limactl list 2>/dev/null | grep apptainer | grep -q Running; then
        echo "✅ Apptainer VM is running"
    else
        echo "🚀 Starting Apptainer VM..."
        limactl start apptainer
        echo "✅ Apptainer VM started"
    fi
else
    echo "📦 Creating Apptainer VM with Rosetta (this may take 5-10 minutes)..."
    echo ""
    echo "Note: Using Rosetta to enable x86_64 container compatibility on Apple Silicon"
    limactl start --name=apptainer --vm-type=vz --rosetta template:apptainer
    echo ""
    echo "✅ Apptainer VM created and started!"
fi

echo ""
echo "🧪 Testing Apptainer..."
if limactl shell apptainer apptainer --version &>/dev/null; then
    VERSION=$(limactl shell apptainer apptainer --version)
    echo "✅ Apptainer works! Version: $VERSION"
else
    echo "❌ Apptainer test failed"
    exit 1
fi

echo ""
echo "========================================="
echo "Setting up shell alias"
echo "========================================="
echo ""

# Detect shell
SHELL_RC=""
if [ -n "$ZSH_VERSION" ]; then
    SHELL_RC="$HOME/.zshrc"
elif [ -n "$BASH_VERSION" ]; then
    SHELL_RC="$HOME/.bashrc"
else
    echo "⚠️  Could not detect shell type"
    echo "Add this line to your shell config manually:"
    echo '  alias singularity="limactl shell apptainer apptainer"'
    echo ""
    exit 0
fi

# Check if alias already exists
if grep -q "alias singularity=" "$SHELL_RC" 2>/dev/null; then
    echo "✅ Singularity alias already exists in $SHELL_RC"
else
    echo "Adding singularity alias to $SHELL_RC..."
    echo "" >> "$SHELL_RC"
    echo "# Singularity/Apptainer via Lima (added by riboseqorg-nf setup)" >> "$SHELL_RC"
    echo 'alias singularity="limactl shell apptainer apptainer"' >> "$SHELL_RC"
    echo "✅ Alias added!"
fi

echo ""
echo "========================================="
echo "Setup Complete!"
echo "========================================="
echo ""
echo "To start using Singularity, either:"
echo ""
echo "1. Reload your shell:"
echo "   source $SHELL_RC"
echo ""
echo "2. Or open a new terminal window"
echo ""
echo "Then you can use singularity commands like:"
echo "   singularity --version"
echo "   singularity pull docker://biocontainers/samtools:1.20--h50ea8bc_0"
echo ""
echo "To run validation:"
echo "   ./scripts/validate-containers.sh  # Will now detect Singularity!"
echo ""
echo "Or use the Lima-specific wrapper:"
echo "   ./scripts/validate-containers-lima.sh"
echo ""
echo "To manage the VM:"
echo "   limactl start apptainer    # Start VM"
echo "   limactl stop apptainer     # Stop VM (saves resources)"
echo "   limactl list               # List all VMs"
echo ""
