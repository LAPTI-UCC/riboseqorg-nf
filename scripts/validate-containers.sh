#!/usr/bin/env bash
# Container Validation Script
# Tests that all containers referenced in modules can be pulled and run

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Counters
PASSED=0
FAILED=0
SKIPPED=0

echo "========================================="
echo "Container Validation Script"
echo "========================================="
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

    # Try to pull
    if ! singularity pull "$temp_sif" "docker://$container" &>/dev/null; then
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

# Check if Docker is available
if ! command -v docker &>/dev/null; then
    echo -e "${YELLOW}Warning: Docker not found. Docker tests will be skipped.${NC}"
    echo ""
fi

# Check if Singularity is available
if ! command -v singularity &>/dev/null; then
    echo -e "${YELLOW}Warning: Singularity not found. Singularity tests will be skipped.${NC}"
    echo ""
fi

echo "Testing BioContainers..."
echo "------------------------"

# Test common BioContainers
if command -v docker &>/dev/null; then
    test_docker_container "biocontainers/samtools:1.20--h50ea8bc_0" "samtools"
    test_docker_container "biocontainers/fastp:0.23.4--h5f740d0_0" "fastp"
    test_docker_container "biocontainers/fastqc:0.12.1--hdfd78af_0" "fastqc"
    test_docker_container "biocontainers/star:2.7.11b--h43eeafb_0" "star"
    test_docker_container "biocontainers/multiqc:1.17--pyhdfd78af_0" "multiqc"
fi

if command -v singularity &>/dev/null; then
    test_singularity_container "biocontainers/samtools:1.20--h50ea8bc_0" "samtools"
    test_singularity_container "biocontainers/fastp:0.23.4--h5f740d0_0" "fastp"
fi

echo ""
echo "Testing Custom Containers..."
echo "----------------------------"

# Test custom containers
if command -v docker &>/dev/null; then
    # These might fail if not built yet
    set +e
    test_docker_container "ghcr.io/jackcurragh/riboseqorg-nf-getrpf:latest" "getRPF" || \
        echo -e "${YELLOW}Note: Custom container not found. Build it first with GitHub Actions.${NC}"
    test_docker_container "ghcr.io/jackcurragh/riboseqorg-nf-ribometric:latest" "RiboMetric" || \
        echo -e "${YELLOW}Note: Custom container not found. Build it first with GitHub Actions.${NC}"
    test_docker_container "ghcr.io/jackcurragh/riboseqorg-nf-rdp-tools:latest" "RDP-tools" || \
        echo -e "${YELLOW}Note: Custom container not found. Build it first with GitHub Actions.${NC}"
    set -e
fi

echo ""
echo "========================================="
echo "Validation Summary"
echo "========================================="
echo -e "Passed:  ${GREEN}$PASSED${NC}"
echo -e "Failed:  ${RED}$FAILED${NC}"
echo -e "Skipped: ${YELLOW}$SKIPPED${NC}"
echo ""

if [ $FAILED -gt 0 ]; then
    echo -e "${RED}Some tests failed!${NC}"
    exit 1
else
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
fi
