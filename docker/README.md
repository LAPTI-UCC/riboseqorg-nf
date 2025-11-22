# Custom Container Dockerfiles

This directory contains Dockerfiles for custom tools that require pip installations from GitHub or specific package configurations.

## Directory Structure

```
docker/
├── getRPF/
│   └── Dockerfile        # Container for getRPF tool
├── RiboMetric/
│   └── Dockerfile        # Container for RiboMetric tool
├── RDP-tools/
│   └── Dockerfile        # Container for RiboSeq-DP-Tools
└── README.md            # This file
```

## Containers

### getRPF
- **Purpose**: RNA structure detection and RPF extraction
- **Base**: condaforge/mambaforge
- **Key packages**: Python 3.10, biopython
- **Pip install**: git+https://github.com/JackCurragh/get-RPF.git

### RiboMetric
- **Purpose**: Ribosome profiling quality control metrics
- **Base**: condaforge/mambaforge
- **Key packages**: Python 3.10, biopython, pysam
- **Pip install**: git+https://github.com/JackCurragh/RiboMetric.git

### RDP-tools
- **Purpose**: Ribosome profiling data processing
- **Base**: condaforge/mambaforge
- **Key packages**: Python 3.10.13
- **Pip install**: RiboSeq-DP-Tools==0.1.10 (from PyPI)

## Building Containers

### Automatic Build (via GitHub Actions)

The easiest way to build these containers is through GitHub Actions:

1. Make changes to any Dockerfile in this directory
2. Commit and push to the main branch:
   ```bash
   git add docker/
   git commit -m "Update getRPF container"
   git push origin main
   ```
3. GitHub Actions will automatically:
   - Build multi-arch Docker images (amd64, arm64)
   - Push to GitHub Container Registry (ghcr.io)
   - Convert to Singularity format
   - Create a release with .sif files

See [.github/workflows/build-containers.yml](../.github/workflows/build-containers.yml) for details.

### Manual Build (Mac/Linux)

#### Docker Build

```bash
# Build a specific container
cd docker/getRPF
docker build -t getrpf:dev .

# Test it
docker run -it --rm getrpf:dev getRPF --version

# Push to GitHub Container Registry (requires authentication)
docker tag getrpf:dev ghcr.io/jackcurragh/riboseqorg-nf-getrpf:dev
docker push ghcr.io/jackcurragh/riboseqorg-nf-getrpf:dev
```

#### Singularity Build (Linux only)

On a Linux system with Singularity installed:

```bash
# Pull from Docker Hub and convert
singularity pull getRPF.sif docker://ghcr.io/jackcurragh/riboseqorg-nf-getrpf:latest

# Or build from local Docker daemon
docker build -t getrpf:latest docker/getRPF/
singularity build getRPF.sif docker-daemon://getrpf:latest
```

**Note**: Mac users cannot build Singularity images natively. Use the GitHub Actions workflow or build on a Linux system.

### Manual Workflow Trigger

You can also trigger the build workflow manually from GitHub:

1. Go to [Actions tab](https://github.com/JackCurragh/riboseqorg-nf/actions)
2. Select "Build and Push Containers"
3. Click "Run workflow"
4. Optionally specify a specific container to build (or leave as "all")

## Testing Containers

### Quick Test

```bash
# Docker
docker run -it --rm ghcr.io/jackcurragh/riboseqorg-nf-getrpf:latest getRPF --help

# Singularity (if you have .sif file)
singularity exec getRPF.sif getRPF --help
```

### Test with Nextflow

```bash
# Create a test profile in nextflow.config
nextflow run main.nf -profile test,singularity
```

## Updating Containers

### When to Update

Update containers when:
- The upstream GitHub repository has important fixes
- You need a different Python or package version
- New dependencies are required

### Update Process

1. **Edit the Dockerfile**:
   ```dockerfile
   # Example: Update Python version
   RUN mamba install -y -c conda-forge -c bioconda \
       python=3.11 \  # Changed from 3.10
       pip \
       biopython && \
       mamba clean -a -y
   ```

2. **Test locally** (optional but recommended):
   ```bash
   cd docker/getRPF
   docker build -t getrpf:test .
   docker run -it --rm getrpf:test getRPF --version
   ```

3. **Commit and push**:
   ```bash
   git add docker/getRPF/Dockerfile
   git commit -m "Update getRPF to Python 3.11"
   git push
   ```

4. **Wait for GitHub Actions** to complete the build

5. **Update module if needed**:
   If you want to pin to a specific version, update the module's container directive:
   ```groovy
   container "ghcr.io/jackcurragh/riboseqorg-nf-getrpf:v1.1.0"
   ```

## Troubleshooting

### Build Failures

Check the GitHub Actions logs:
1. Go to the [Actions tab](https://github.com/JackCurragh/riboseqorg-nf/actions)
2. Click on the failed workflow run
3. Expand the failed step to see error messages

Common issues:
- **"gcc not found" or compilation errors**: Install dependencies via conda instead of pip (see [TROUBLESHOOTING.md](TROUBLESHOOTING.md#issue-arm64-build-fails-with-gcc-not-found))
- **Pip install fails**: The GitHub repository might be unavailable or have breaking changes
- **Conda package conflicts**: Adjust package versions in the Dockerfile
- **Docker build timeout**: The base image or packages are too large
- **ARM64 build fails**: Use conda packages which provide pre-compiled binaries for all architectures

**Pro tip:** See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for detailed solutions to common build issues.

### Container Pull Fails in Pipeline

If Nextflow can't pull containers:

1. **Check authentication**:
   ```bash
   echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin
   ```

2. **Pre-pull containers**:
   ```bash
   singularity pull getRPF.sif docker://ghcr.io/jackcurragh/riboseqorg-nf-getrpf:latest
   mv getRPF.sif ./singularity/
   ```

3. **Use conda fallback**:
   ```bash
   nextflow run main.nf -profile conda
   ```

## Resources

- [Docker Documentation](https://docs.docker.com/)
- [Singularity Documentation](https://docs.sylabs.io/)
- [GitHub Container Registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry)
- [BioContainers Best Practices](https://biocontainers-edu.readthedocs.io/)
