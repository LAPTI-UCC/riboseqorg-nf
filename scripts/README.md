# Validation Scripts

This directory contains scripts for validating and testing containers.

## Available Scripts

### 1. validate-containers.sh

Main container validation script that tests both Docker and Singularity containers.

**Usage:**
```bash
./scripts/validate-containers.sh
```

**Features:**
- Automatically detects available container runtimes (Docker/Singularity)
- Tests BioContainers and custom containers
- Gracefully skips unavailable runtimes
- Color-coded output
- Works on Mac with just Docker

**Requirements:**
- At least one of: Docker or Singularity/Apptainer
- Internet connection for pulling containers

---

### 2. validate-containers-lima.sh (Mac only)

Wrapper script that runs validation inside a Lima VM with Apptainer/Singularity support.

**Usage:**
```bash
./scripts/validate-containers-lima.sh
```

**Features:**
- Automatically starts Lima VM if needed
- Runs full validation with both Docker and Singularity
- Handles file path mapping between Mac and VM
- Provides helpful tips for VM management

**Requirements:**
- Lima installed (`brew install lima`)
- Apptainer VM created (script will create if needed)

**When to use:**
- You want to test Singularity containers on Mac
- You want to validate the full workflow before deploying to HPC
- You need to ensure Singularity containers work correctly

---

### 3. setup-singularity-mac.sh (Mac only)

One-time setup script that configures Singularity/Apptainer on your Mac using Lima.

**Usage:**
```bash
./scripts/setup-singularity-mac.sh
```

**What it does:**
1. Checks if Lima is installed (installs if needed via Homebrew)
2. Creates Apptainer VM if it doesn't exist
3. Tests that Apptainer works
4. Adds `singularity` alias to your shell config
5. Provides next steps

**After running:**
```bash
# Reload your shell
source ~/.zshrc  # or ~/.bashrc

# Now singularity commands work!
singularity --version

# And the regular validation script will detect Singularity
./scripts/validate-containers.sh
```

**One-time setup** - only need to run once per Mac.

---

## Quick Start Guide

### For Mac Users

#### Option 1: Docker Only (Fastest)
```bash
# Just run the validation - will use Docker only
./scripts/validate-containers.sh
```

#### Option 2: Docker + Singularity (Most Complete)
```bash
# One-time setup
./scripts/setup-singularity-mac.sh

# Reload shell
source ~/.zshrc

# Now run validation (will test both Docker and Singularity)
./scripts/validate-containers.sh
```

#### Option 3: Use Lima Wrapper Directly
```bash
# No setup needed - wrapper handles everything
./scripts/validate-containers-lima.sh
```

### For Linux/HPC Users

```bash
# Just run validation - Singularity should be natively available
./scripts/validate-containers.sh
```

---

## Understanding the Output

### Successful Run
```
=========================================
Container Validation Script
=========================================

✓ Docker found: Docker version 24.0.5
✓ Singularity found: apptainer version 1.2.4

Testing BioContainers...
------------------------
Testing Docker: biocontainers/samtools:1.20--h50ea8bc_0 ... PASSED
Testing Docker: biocontainers/fastp:0.23.4--h5f740d0_0 ... PASSED

Running sample Singularity tests...
Testing Singularity: biocontainers/samtools:1.20--h50ea8bc_0 ... PASSED

Testing Custom Containers...
----------------------------
Testing Docker: ghcr.io/jackcurragh/riboseqorg-nf-getrpf:latest ... PASSED

=========================================
Validation Summary
=========================================

Container Runtime:
  Docker:      Available
  Singularity: Available

Test Results:
  Passed:      8
  Failed:      0
  Skipped:     0

✓ All available tests passed!
```

### Mac with Docker Only
```
=========================================
Container Validation Script
=========================================

✓ Docker found: Docker version 24.0.5
✗ Singularity not found - Singularity tests will be skipped

Testing BioContainers...
------------------------
Testing Docker: biocontainers/samtools:1.20--h50ea8bc_0 ... PASSED

Skipping Singularity tests (Singularity not available)

Testing Custom Containers...
----------------------------
  Note: getRPF container not found. Build it first with GitHub Actions.

Test Results:
  Passed:      5
  Failed:      0
  Skipped:     5

✓ All available tests passed!
Note: 5 tests were skipped (container runtime not available or containers not built yet)
```

---

## Troubleshooting

### "Neither Docker nor Singularity is available"

**Install Docker:**
- Mac: Download from [docker.com](https://www.docker.com/products/docker-desktop)
- Linux: `sudo apt-get install docker.io` or equivalent

**Or install Singularity via Lima (Mac):**
```bash
./scripts/setup-singularity-mac.sh
```

### "VM not found" (Lima/Mac)

Run the setup script:
```bash
./scripts/setup-singularity-mac.sh
```

### Tests fail for custom containers

Custom containers need to be built first:
1. Push your code to GitHub
2. Wait for GitHub Actions to build containers
3. Run validation again

Or build locally:
```bash
cd docker/getRPF
docker build -t getrpf:test .
```

### Lima VM is slow

```bash
# Stop the VM when not in use
limactl stop apptainer

# Or allocate more resources
limactl edit apptainer
# Edit cpus and memory settings
```

---

## Best Practices

1. **Development on Mac:**
   - Use Docker for daily development (faster)
   - Use Singularity for final validation before HPC deployment

2. **Before committing:**
   ```bash
   ./scripts/validate-containers.sh
   ```

3. **Before deploying to HPC:**
   ```bash
   # On Mac
   ./scripts/setup-singularity-mac.sh
   source ~/.zshrc
   ./scripts/validate-containers.sh

   # Should show both Docker and Singularity tests passing
   ```

4. **CI/CD:**
   - GitHub Actions automatically runs validation
   - Check the Actions tab for results

---

## Files Overview

```
scripts/
├── README.md                        # This file
├── validate-containers.sh           # Main validation script
├── validate-containers-lima.sh      # Lima/Mac wrapper
└── setup-singularity-mac.sh         # One-time Mac setup
```

---

## Advanced Usage

### Test Specific Container

Edit `validate-containers.sh` and comment out containers you don't want to test:

```bash
# Test only samtools
test_docker_container "biocontainers/samtools:1.20--h50ea8bc_0" "samtools"
```

### Run in CI/CD

```yaml
- name: Validate containers
  run: ./scripts/validate-containers.sh
```

### Run with verbose output

```bash
# Show all pull output
bash -x ./scripts/validate-containers.sh
```

### Pre-pull containers

```bash
# Pull all containers first (for offline testing)
docker pull biocontainers/samtools:1.20--h50ea8bc_0
docker pull biocontainers/fastp:0.23.4--h5f740d0_0
# ... etc

# Then run validation (will be faster)
./scripts/validate-containers.sh
```

---

## Resources

- [Docker Documentation](https://docs.docker.com/)
- [Singularity Documentation](https://docs.sylabs.io/)
- [Lima Documentation](https://lima-vm.io/)
- [Apptainer Documentation](https://apptainer.org/)
- [Main Container Docs](../docs/CONTAINERS.md)
- [Singularity on Mac Guide](../docs/SINGULARITY_ON_MAC.md)
