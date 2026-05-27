# Conda Channel Publishing

This recipe is used by `.github/workflows/publish-conda-channel.yml` to publish `cell-gps` into the static channel at `hutaobo/conda-channel`.

The workflow runs when a GitHub release is published, when a `v*` tag is pushed, or manually through `workflow_dispatch` with a version. It builds from the PyPI sdist for the resolved version, copies the noarch package into the channel repository, re-indexes the channel, and pushes the channel update.

Required repository secret: `CONDA_CHANNEL_PAT`, a GitHub token with write access to `hutaobo/conda-channel`.

