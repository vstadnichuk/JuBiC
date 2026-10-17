# Building and Publishing the Documentation

The HTML site is built with Documenter.jl from the Markdown sources in
`docs/src/`. The detailed HNDP generator and solver guides are maintained in
`examples/HNDP/NETWORK_GENERATION.md` and
`examples/HNDP/HNDP_SOLVER_OVERVIEW.md`. `docs/make.jl` incorporates them and their
JSON examples automatically, so fixes appear in both the repository guides and
the website without maintaining two copies.

## Local build

From the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. docs/make.jl
```

Open `docs/build/index.html` in a browser. Local builds use `.html` links that
work directly from disk. Generated HTML and copied reference sources are
ignored by Git. Building the documentation imports JuBiC but does not run the
optimization examples, so the build does not need a Gurobi license. Running the
examples or solver tests does require the relevant solver setup.

## GitHub Pages

The configured public URL is <https://vstadnichuk.github.io/JuBiC/>. HTML is
published by `.github/workflows/documentation.yml`, using GitHub Pages' artifact
deployment. The Markdown sources and workflow remain part of this repository;
the generated site is served by GitHub Pages rather than as a raw HTML file in
the GitHub code viewer.

One-time setup after reviewing and committing the workflow:

1. Push the documentation changes to `v0.6` or `main`.
2. In the repository, open **Settings → Pages → Build and deployment** and set
   **Source** to **GitHub Actions**.
3. If necessary, allow `main` and `v0.6` in the **github-pages** environment's
   deployment branch rules under **Settings → Environments**.
4. Open **Actions → Documentation** and rerun the workflow (or use **Run
   workflow** on a branch containing it).
5. Confirm the deployment succeeds and the README's HTML documentation link opens.

Pushes to `main` and `v0.6` publish the latest site. Pull requests targeting
these branches build and check the docs without publishing. Forks cannot deploy
to this repository's site. No personal access token or deployment key is needed;
the deployment job uses GitHub's scoped Pages token permissions.

Until the workflow is committed/pushed and Pages enabled, the public URL may
return 404. A successful local build confirms the generated site, not a live
GitHub deployment.

The build treats Documenter warnings as errors. Fix broken internal links,
missing navigation pages, or invalid documentation blocks before publishing.
Code examples shown as ordinary `julia` fences are instructional snippets, not
executed doctests. Optimizer-dependent validation is run separately.
