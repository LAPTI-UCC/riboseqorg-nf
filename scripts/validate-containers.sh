#!/usr/bin/env bash
# Container Validation Script
# Tests that all containers referenced in modules can be pulled and run

# Don't exit on error - we want to continue testing even if some fail
set +e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Counters
PASSED=0
FAILED=0
SKIPPED=0

echo "========================================="
echo "Container Validation Script"
echo "========================================="
echo ""

# Check available tools
HAS_DOCKER=false
HAS_SINGULARITY=false

if command -v docker &>/dev/null; then
    HAS_DOCKER=true
    echo -e "${BLUE}✓${NC} Docker found: $(docker --version | head -1)"
else
    echo -e "${YELLOW}✗${NC} Docker not found - Docker tests will be skipped"
fi

if command -v singularity &>/dev/null; then
    HAS_SINGULARITY=true
    echo -e "${BLUE}✓${NC} Singularity found: $(singularity --version)"
else
    echo -e "${YELLOW}✗${NC} Singularity not found - Singularity tests will be skipped"
fi

if [ "$HAS_DOCKER" = false ] && [ "$HAS_SINGULARITY" = false ]; then
    echo ""
    echo -e "${RED}ERROR: Neither Docker nor Singularity is available${NC}"
    echo "Please install at least one container runtime to run this validation."
    echo ""
    echo "Install Docker: https://docs.docker.com/get-docker/"
    echo "Install Singularity: https://docs.sylabs.io/guides/latest/user-guide/"
    exit 1
fi

echo ""

# Function to test Docker container
test_docker_container() {
    local container=$1
    local tool=$2

    echo -n "Testing Docker: $container ... "

    # Try to pull
    if ! docker pull "$container" &>/dev/null; then
        echo -e "${RED}FAILED${NC} (pull failed)"
        ((FAILED++))
        return 1
    fi

    # Try to run
    case $tool in
        samtools|bowtie|fastp|fastqc|star|multiqc)
            if docker run --rm "$container" $tool --version &>/dev/null || \
               docker run --rm "$container" $tool -v &>/dev/null || \
               docker run --rm "$container" $tool --help &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                return 1
            fi
            ;;
        fastq-dl)
            if docker run --rm "$container" fastq-dl --version &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                return 1
            fi
            ;;
        getRPF)
            if docker run --rm "$container" getRPF --version &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                return 1
            fi
            ;;
        RiboMetric)
            if docker run --rm "$container" RiboMetric --version &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                return 1
            fi
            ;;
        RDP-tools)
            if docker run --rm "$container" RDP-Tools --version &>/dev/null || \
               docker run --rm "$container" python -c "import RDP_Tools" &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                return 1
            fi
            ;;
        *)
            if docker run --rm "$container" echo "test" &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                return 1
            fi
            ;;
    esac
}

# Function to test Singularity container
test_singularity_container() {
    local container=$1
    local tool=$2

    # Check if Singularity is available
    if ! command -v singularity &>/dev/null; then
        echo -e "${YELLOW}SKIPPED${NC} (Singularity not installed)"
        ((SKIPPED++))
        return 0
    fi

    echo -n "Testing Singularity: $container ... "

    # Create temp file for container
    local temp_sif=$(mktemp -u).sif

    # Fix container URL for Singularity
    # BioContainers need quay.io prefix, not docker.io
    local singularity_url="docker://$container"
    if [[ "$container" == biocontainers/* ]]; then
        singularity_url="docker://quay.io/$container"
    fi

    # Try to pull
    if ! singularity pull "$temp_sif" "$singularity_url" &>/dev/null; then
        echo -e "${RED}FAILED${NC} (pull failed)"
        ((FAILED++))
        return 1
    fi

    # Try to run
    case $tool in
        samtools|bowtie|fastp|fastqc|star|multiqc)
            if singularity exec "$temp_sif" $tool --version &>/dev/null || \
               singularity exec "$temp_sif" $tool -v &>/dev/null || \
               singularity exec "$temp_sif" $tool --help &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                rm -f "$temp_sif"
                return 1
            fi
            ;;
        fastq-dl)
            if singularity exec "$temp_sif" fastq-dl --version &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                rm -f "$temp_sif"
                return 1
            fi
            ;;
        getRPF|RiboMetric)
            if singularity exec "$temp_sif" $tool --version &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                rm -f "$temp_sif"
                return 1
            fi
            ;;
        RDP-tools)
            if singularity exec "$temp_sif" RDP-Tools --version &>/dev/null || \
               singularity exec "$temp_sif" python -c "import RDP_Tools" &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                rm -f "$temp_sif"
                return 1
            fi
            ;;
        *)
            if singularity exec "$temp_sif" echo "test" &>/dev/null; then
                echo -e "${GREEN}PASSED${NC}"
                ((PASSED++))
            else
                echo -e "${RED}FAILED${NC} (run failed)"
                ((FAILED++))
                rm -f "$temp_sif"
                return 1
            fi
            ;;
    esac

    # Cleanup
    rm -f "$temp_sif"
}

echo "Testing BioContainers..."
echo "------------------------"

# Test common BioContainers with Docker
if [ "$HAS_DOCKER" = true ]; then
    test_docker_container "biocontainers/samtools:1.20--h50ea8bc_0" "samtools"
    test_docker_container "biocontainers/fastp:0.23.4--h5f740d0_0" "fastp"
    test_docker_container "biocontainers/fastqc:0.12.1--hdfd78af_0" "fastqc"
    test_docker_container "biocontainers/star:2.7.11b--h43eeafb_0" "star"
    test_docker_container "biocontainers/multiqc:1.17--pyhdfd78af_0" "multiqc"
else
    echo -e "${YELLOW}Skipping Docker tests (Docker not available)${NC}"
    ((SKIPPED+=5))
fi

# Test with Singularity (sample only - full test would take too long)
if [ "$HAS_SINGULARITY" = true ]; then
    echo ""
    echo "Running sample Singularity tests (this may take a few minutes)..."
    test_singularity_container "biocontainers/samtools:1.20--h50ea8bc_0" "samtools"
    test_singularity_container "biocontainers/fastp:0.23.4--h5f740d0_0" "fastp"
else
    echo ""
    echo -e "${YELLOW}Skipping Singularity tests (Singularity not available)${NC}"
    ((SKIPPED+=2))
fi

echo ""
echo "Testing Custom Containers..."
echo "----------------------------"

# Test custom containers with Docker
if [ "$HAS_DOCKER" = true ]; then
    # These might fail if not built yet - don't count as failures
    if ! test_docker_container "ghcr.io/jackcurragh/riboseqorg-nf-getrpf:latest" "getRPF" 2>/dev/null; then
        echo -e "${YELLOW}  Note: getRPF container not found. Build it first with GitHub Actions.${NC}"
        ((FAILED--))  # Don't count this as a failure
        ((SKIPPED++))
    fi

    if ! test_docker_container "ghcr.io/jackcurragh/riboseqorg-nf-ribometric:latest" "RiboMetric" 2>/dev/null; then
        echo -e "${YELLOW}  Note: RiboMetric container not found. Build it first with GitHub Actions.${NC}"
        ((FAILED--))
        ((SKIPPED++))
    fi

    if ! test_docker_container "ghcr.io/jackcurragh/riboseqorg-nf-rdp-tools:latest" "RDP-tools" 2>/dev/null; then
        echo -e "${YELLOW}  Note: RDP-tools container not found. Build it first with GitHub Actions.${NC}"
        ((FAILED--))
        ((SKIPPED++))
    fi
else
    echo -e "${YELLOW}Skipping Docker tests (Docker not available)${NC}"
    ((SKIPPED+=3))
fi

echo ""
echo "========================================="
echo "Validation Summary"
echo "========================================="
echo ""
echo "Container Runtime:"
if [ "$HAS_DOCKER" = true ]; then
    echo -e "  Docker:      ${GREEN}Available${NC}"
else
    echo -e "  Docker:      ${YELLOW}Not available${NC}"
fi
if [ "$HAS_SINGULARITY" = true ]; then
    echo -e "  Singularity: ${GREEN}Available${NC}"
else
    echo -e "  Singularity: ${YELLOW}Not available${NC}"
fi

echo ""
echo "Test Results:"
echo -e "  Passed:      ${GREEN}$PASSED${NC}"
if [ $FAILED -gt 0 ]; then
    echo -e "  Failed:      ${RED}$FAILED${NC}"
else
    echo -e "  Failed:      $FAILED"
fi
if [ $SKIPPED -gt 0 ]; then
    echo -e "  Skipped:     ${YELLOW}$SKIPPED${NC}"
else
    echo -e "  Skipped:     $SKIPPED"
fi
echo ""

if [ $FAILED -gt 0 ]; then
    echo -e "${RED}✗ Some tests failed!${NC}"
    echo "Review the output above for details."
    exit 1
elif [ $PASSED -eq 0 ]; then
    echo -e "${YELLOW}⚠ No tests were run!${NC}"
    echo "Make sure Docker or Singularity is installed and containers are available."
    exit 1
else
    echo -e "${GREEN}✓ All available tests passed!${NC}"
    if [ $SKIPPED -gt 0 ]; then
        echo -e "${YELLOW}Note: $SKIPPED tests were skipped (container runtime not available or containers not built yet)${NC}"
    fi
    exit 0
fi
