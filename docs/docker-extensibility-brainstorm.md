# Sunbeam Docker Extensibility Brainstorm

The current Docker images bundle Sunbeam and preload the project conda
environments, which makes the images convenient out of the box but harder to
extend. Below are options to rework both Dockerfiles so they are easier for
users to layer new extensions (`sunbeam extend <repo>`) and additional Python
packages such as Snakemake executor plugins.

## 1. Split the image into base and runtime layers

* **Create a `sunbeam-base` stage** that installs OS-level dependencies,
  Mamba/conda configuration, and copies the core Sunbeam source. Export this
  stage for downstream users so they can `FROM sunbeam-base` to add their own
  packages without repeating the heavy bootstrap steps.
* **Add a thin `sunbeam-runtime` stage** that keeps the current behaviour of
  pre-creating conda environments. Publish both images so users can choose
  between quick-start (`sunbeam-runtime`) and extendable (`sunbeam-base`).
* **Make `slim.Dockerfile` reuse the shared base stage** to avoid divergence and
  to ensure extensions can target the same minimal foundation.

## 2. Externalise extension installation

* **Accept a build argument or bind-mount for the extensions directory.** During
  build, copy only the base workflow. Users can supply their own
  `$SUNBEAM_EXTENSIONS` directory in a derived Dockerfile or at runtime via a
  volume mount.
* **Document an `ENTRYPOINT` wrapper** that runs `sunbeam extend` for any Git
  URLs passed through an environment variable (e.g. `SUNBEAM_AUTO_EXTEND`). This
  allows derived images to add extensions without re-implementing the bootstrap
  steps.

## 3. Provide pip-friendly hooks

* **Install Sunbeam in editable mode** (`pip install -e .`) or install into a
  virtual environment located in `/opt/sunbeam`. This makes it straightforward
  for derived images to `pip install` extra Python packages, including Snakemake
  plugins, without conflicting with system Python.
* **Expose a requirements overlay file** (e.g. copy `requirements.ext.txt` if it
  exists). Derived Dockerfiles can simply drop in their dependency file and run
  `pip install -r requirements.ext.txt` to add Snakemake executor plugins.

## 4. Cache conda environments without locking the image

* **Move the conda environment creation into an optional script** (e.g.
  `/usr/local/bin/prewarm-conda-envs.sh`). The default image can still invoke it
  during build, but derived images can skip it to keep their layer mutable.
* **Add a `SUNBEAM_PREWARM_CONDA=false` build argument** that disables the
  `sunbeam run --conda-create-envs-only` step for users who want a lighter base
  to extend.

## 5. Support runtime extension installs

* **Include a non-root `sunbeam` user with passwordless sudo** for package
  installs at runtime. This way users can run `pip install` or `mamba install`
  inside a running container without rebuilding.
* **Mount `/home/sunbeam/extensions` as a volume by default** in the published
  `docker run` examples so users can persist and share extensions across
  container rebuilds.

---

Each of the ideas above can be combined. A likely minimal change would be to
introduce a shared base stage (Idea 1) plus an optional build argument for
pre-warming conda (Idea 4). That keeps the existing user experience intact while
making it easy for downstream images to extend Sunbeam and add Snakemake
plugins.
