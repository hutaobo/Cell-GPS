# scverse template alignment

Cell-GPS was not originally generated from `scverse/cookiecutter-scverse`.
The repository now adopts the parts of that structure that are low risk for an
existing package with published users.

## Adopted

- `src/` package layout with the `cellgps` import namespace.
- Sphinx documentation under `docs/` and ReadTheDocs configuration.
- GitHub Actions for tests, release publishing, and package build checks.
- Structured issue templates for bug reports and feature requests.
- `pyproject.toml` as the central package metadata and tool configuration file.
- Development metadata for pytest, coverage, ruff, pre-commit, and Codecov.
- scverse-style public namespace aliases: `cellgps.pp`, `cellgps.tl`, and `cellgps.pl`.

## Deliberately preserved

- The build backend remains `setuptools.build_meta`; switching to hatchling is
  unnecessary for users and would change the packaging path.
- Python support remains `>=3.9` instead of the template's newer default, because
  the package already advertises Python 3.9 compatibility.
- The PyPI publishing workflow remains `.github/workflows/python-publish.yml`,
  matching the configured trusted publisher.
- Existing `.rst` documentation pages are preserved instead of converting the
  documentation tree to MyST Markdown.

## Future optional work

- Add Codecov as a GitHub App if coverage comments should appear on pull requests.
- Enable pre-commit.ci after the codebase has been formatted and linted in
  a separate maintenance change.
- Consider a future hatch/uv migration only if it is worth changing the local
  developer workflow and release tooling.
