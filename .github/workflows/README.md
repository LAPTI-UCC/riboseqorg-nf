# GitHub Actions Workflows

This directory contains automated workflows for the riboseqorg-nf pipeline.

## Available Workflows

### 1. Build and Push Containers ([build-containers.yml](build-containers.yml))

**Purpose**: Automatically build Docker images and convert to Singularity format.

**Triggers**:
- Push to main branch (changes to `docker/` or workflow file)
- Manual dispatch via GitHub Actions UI
- Pull requests (for testing, but doesn't push)

**What it does**:
1. Builds multi-architecture Docker images (linux/amd64, linux/arm64)
2. Pushes to GitHub Container Registry (ghcr.io)
3. Converts Docker images to Singularity (.sif) format
4. Creates GitHub Release with .sif files as downloadable assets
5. Caches layers for faster rebuilds

**Containers built**:
- getRPF (Python + biopython + getRPF from GitHub)
- RiboMetric (Python + pysam + RiboMetric from GitHub)
- RDP-tools (Python + RiboSeq-DP-Tools from PyPI)

**Manual trigger**:
1. Go to [Actions](https://github.com/JackCurragh/riboseqorg-nf/actions)
2. Select "Build and Push Containers"
3. Click "Run workflow"
4. Optionally specify which container to build (default: all)

**Outputs**:
- Docker images: `ghcr.io/jackcurragh/riboseqorg-nf-{tool}:latest`
- Singularity files: Available in [Releases](https://github.com/JackCurragh/riboseqorg-nf/releases)

---

### 2. Validate Containers ([validate-containers.yml](validate-containers.yml))

**Purpose**: Ensure all containers referenced in modules are valid and functional.

**Triggers**:
- Pull requests (changes to `modules/`, `docker/`, or workflow file)
- Push to main branch
- Manual dispatch
- Weekly schedule (every Monday at 6am UTC)

**What it does**:
1. **Extracts** all container URLs from module files
2. **Validates BioContainers** (both Docker and Singularity):
   - Pulls containers from repositories
   - Tests that tools inside actually run
   - Verifies version commands work
3. **Validates Custom Containers**:
   - Tests getRPF, RiboMetric, RDP-tools
   - Checks both Docker and Singularity formats
4. **Checks Module Definitions**:
   - Ensures all processes have container directives
   - Validates container URL format

**Test matrix**:
The workflow dynamically creates test jobs based on containers found in modules, so adding new modules automatically adds validation.

**Weekly monitoring**:
The weekly run helps catch:
- Broken upstream containers
- Deprecated container versions
- Registry availability issues

**Viewing results**:
1. Go to [Actions](https://github.com/JackCurragh/riboseqorg-nf/actions)
2. Select "Validate Containers" workflow
3. View the summary for pass/fail status
4. Click individual jobs for detailed logs

**Local testing**:
You can run similar tests locally:
```bash
./scripts/validate-containers.sh
```

---

## Workflow Best Practices

### For Contributors

**When modifying Docker images**:
1. Update the Dockerfile in `docker/{tool}/`
2. Test locally if possible: `docker build -t test docker/{tool}/`
3. Create a pull request
4. Validation will run automatically
5. After merge, containers will be built and pushed

**When adding new modules**:
1. Include both `conda` and `container` directives
2. Use the standard format:
   ```groovy
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/tool:version--hash' :
       'biocontainers/tool:version--hash' }"
   ```
3. Validation will automatically test the new container

**When containers fail validation**:
1. Check the Actions logs for specific error
2. Common issues:
   - Container doesn't exist in registry
   - Tool command changed (adjust test in workflow)
   - Upstream breaking changes
3. Fix the issue and push again

### For Maintainers

**Approving PRs**:
- Ensure "Validate Containers" check passes
- Review any new container definitions
- Check that custom containers have corresponding Dockerfiles

**Monitoring**:
- Check weekly validation runs
- Update container versions if deprecated
- Address any failing tests promptly

**Releasing**:
- Container builds are automatic on merge to main
- Tag releases for stable versions
- Update module container URIs to use version tags instead of `latest` for releases

---

## Troubleshooting

### Build failures

**"Failed to push to ghcr.io"**
- Check repository permissions
- Ensure `GITHUB_TOKEN` has package write permissions
- Verify GitHub Actions is enabled

**"Singularity build timeout"**
- Container might be too large
- Check if base image is accessible
- Consider using smaller base image

### Validation failures

**"Failed to pull container"**
- Container might not exist yet (for custom containers)
- Check registry URL is correct
- Verify container was successfully built

**"Tool command failed"**
- Tool might have different CLI in new version
- Update the test command in workflow
- Check if tool installation was successful

**"Container definition missing"**
- Add container directive to module
- Follow the standard format shown above

### Weekly validation failures

If weekly validation starts failing:
1. Check if upstream BioContainers changed
2. Look for deprecation notices
3. Update container versions in modules
4. Re-run validation to confirm fix

---

## Adding New Workflows

When adding new workflows:
1. Create `.yml` file in this directory
2. Add clear comments explaining purpose
3. Use appropriate triggers
4. Add to this README
5. Test with manual dispatch first

---

## Resources

- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [Docker Build and Push Action](https://github.com/docker/build-push-action)
- [Singularity Setup Action](https://github.com/eWaterCycle/setup-singularity)
- [BioContainers](https://biocontainers.pro/)
