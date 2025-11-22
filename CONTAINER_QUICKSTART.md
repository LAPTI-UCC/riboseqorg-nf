# Container Quick Start Guide

This guide will help you get the riboseqorg-nf pipeline running with Singularity containers.

## Prerequisites

- Nextflow installed
- Singularity/Apptainer installed (for HPC) OR Docker (for local)
- GitHub account (for building custom containers)

## Option 1: Use Existing Containers (Fastest)

### Step 1: Run with Singularity

```bash
nextflow run main.nf \
  -profile singularity \
  --input samples.csv \
  --outdir results
```

Nextflow will automatically:
- Pull BioContainers from depot.galaxyproject.org
- Pull custom containers from GitHub Releases
- Cache everything in `./singularity/`

### Step 2: Monitor Progress

Watch for container pulls in the log:
```
Pulling Singularity image docker://biocontainers/samtools:1.20--h50ea8bc_0 [cache /path/to/singularity/...]
```

## Option 2: Build Custom Containers (If Needed)

### Step 1: Trigger GitHub Actions Build

Since you're on Mac and the custom GitHub tools might have been updated:

1. Go to your repository on GitHub
2. Navigate to **Actions** tab
3. Select **"Build and Push Containers"** workflow
4. Click **"Run workflow"** dropdown
5. Select branch (main) and click **"Run workflow"**

This will:
- Build Docker images for getRPF, RiboMetric, and RDP-tools
- Convert to Singularity (.sif) format
- Upload to GitHub Releases

### Step 2: Download Singularity Images (Optional)

If you want to pre-download for offline use:

```bash
# Create singularity cache directory
mkdir -p singularity

# Download from GitHub Releases
cd singularity
wget https://github.com/JackCurragh/riboseqorg-nf/releases/download/containers-latest/getRPF.sif
wget https://github.com/JackCurragh/riboseqorg-nf/releases/download/containers-latest/RiboMetric.sif
wget https://github.com/JackCurragh/riboseqorg-nf/releases/download/containers-latest/RDP-tools.sif
cd ..
```

### Step 3: Run Pipeline

```bash
nextflow run main.nf -profile singularity --input samples.csv --outdir results
```

## Option 3: Docker (Mac/Local Development)

```bash
nextflow run main.nf \
  -profile docker \
  --input samples.csv \
  --outdir results
```

Docker images will be pulled from:
- `biocontainers/*` for standard tools
- `ghcr.io/jackcurragh/riboseqorg-nf-*` for custom tools

## Troubleshooting

### Problem: "Failed to pull singularity image"

**Solution 1**: Check internet connectivity and retry
```bash
nextflow run main.nf -profile singularity -resume
```

**Solution 2**: Pre-pull manually
```bash
singularity pull docker://biocontainers/samtools:1.20--h50ea8bc_0
```

**Solution 3**: Fall back to conda
```bash
nextflow run main.nf -profile conda --input samples.csv --outdir results
```

### Problem: "Cannot access GitHub Container Registry"

**Solution**: Authenticate with GitHub
```bash
# For Docker
echo $GITHUB_TOKEN | docker login ghcr.io -u YOUR_USERNAME --password-stdin

# For Singularity
export SINGULARITY_DOCKER_USERNAME=YOUR_USERNAME
export SINGULARITY_DOCKER_PASSWORD=$GITHUB_TOKEN
```

### Problem: Custom containers are outdated

**Solution**: Rebuild via GitHub Actions (see Option 2 above)

### Problem: Out of disk space

**Solution**: Clean old containers
```bash
# Docker
docker system prune -a

# Singularity
rm -rf work/.singularity
rm -rf singularity/*
```

## Verifying Container Setup

### Test Individual Containers

```bash
# Test standard tool
singularity exec docker://biocontainers/samtools:1.20--h50ea8bc_0 samtools --version

# Test custom tool (Docker)
docker run --rm ghcr.io/jackcurragh/riboseqorg-nf-getrpf:latest getRPF --version

# Test custom tool (Singularity)
singularity exec getRPF.sif getRPF --version
```

### Test Full Pipeline with Containers

```bash
# Use test profile if available
nextflow run main.nf -profile test,singularity

# Or run with small dataset
nextflow run main.nf \
  -profile singularity \
  --input test_samples.csv \
  --outdir test_results \
  --max_cpus 4 \
  --max_memory 8.GB
```

## Container Locations

After running, containers are cached in:

**Singularity**:
```
./singularity/
├── biocontainers-samtools-1.20--h50ea8bc_0.img
├── biocontainers-fastp-0.23.4--h5f740d0_0.img
└── ... (other containers)
```

**Docker**:
```bash
docker images | grep -E "biocontainers|ghcr.io"
```

## Next Steps

1. **Production Use**: Pin container versions in modules instead of using `:latest`
2. **Offline Use**: Download all .sif files to a shared filesystem
3. **Customization**: Modify Dockerfiles in `docker/` directory and rebuild
4. **Documentation**: See [docs/CONTAINERS.md](docs/CONTAINERS.md) for detailed information

## Getting Help

- **Container Issues**: Check [docs/CONTAINERS.md](docs/CONTAINERS.md)
- **Build Issues**: Check [docker/README.md](docker/README.md)
- **Pipeline Issues**: Check main README.md or open an issue
