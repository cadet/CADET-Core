# CADET-Core Docs: Sphinx → MyST/JupyterBook Migration

## Context

- **Source:** `doc/` -- ~160 RST files, Alabaster theme, no notebooks
- **Goal:** Migrate to MyST-Markdown + JupyterBook v2 for modern tooling, notebook support, and better UX
- **Note:** JupyterBook v2 does not use Sphinx; it builds natively with MyST-MD

## Content Structure

```
doc/
├── getting_started/   # Installation, build, tutorials
├── modelling/         # Physics/math models (binding, reactions, crystallization, unit ops)
├── simulation/        # Solver and time integration
├── interface/         # HDF5 file format, inputs/outputs (mirrors modelling/ structure)
├── developer_guide/   # Architecture, testing, CADET-Python integration
├── examples/          # 5 workflow examples (chromatography, CSTR, RTD, etc.)
├── _static/           # CSS, logos
└── literature.bib     # BibTeX references (used in publications.rst)
```

## Risk Items

- **Versioning:** deferred; not in scope for now
- **RST conversion quality:** ~160 files; custom directives (`.. todo::`, raw HTML, includes) need manual fixes after batch conversion

## Styling Guidelines

These apply to all new and converted documentation.

- One sentence per line (improves diff readability in version control)
- No em-dashes; use a comma, semicolon, or restructure the sentence
- Use `@label` for internal cross-references (modern MyST, replaces `:ref:` / `{numref}`)

## Deferred: Future Restructuring

After conversion is stable, split and reorganize into:
- **CADET-Core** -- solver internals, developer guide, simulation
- **Modelling + Interface** -- unified section (currently split but parallel in structure)
- **CADET-Process** -- integrate separate project docs

Not in scope for initial migration.

## Deferred: Notation and Content Updates

Updates to defer until after conversion:
- Rename solid-phase concentration symbol: $q$ → $c^s$ throughout modelling and interface sections
- (add further notation/content todos here)

## Tasks

### 1. Bootstrap JupyterBook
- [ ] Install `jupyter-book` (v2)
- [ ] Add `myst.yml` (JB v2 config, theme: `book-theme`) and table of contents mirroring current `index.rst` toctree
- [ ] Verify `jb build doc/` produces valid output

### 2. Replace Theme
- [ ] Use `book-theme` default styling (no custom theme config needed)
- [ ] Port any necessary overrides from `_static/custom.css`

### 3. Batch Convert RST → MyST-MD
- [ ] Run `rst2myst` on all `.rst` files (preferred over pandoc; Sphinx-directive-aware)
- [ ] Audit conversion output for broken directives, raw HTML, `.. include::`, `.. literalinclude::`
- [ ] Fix cross-references: `:doc:`, `:ref:` → MyST `@label` syntax
- [ ] Verify math blocks (numbered equations, MathJax)

### 4. Bibliography
- [ ] Drop `sphinxcontrib.bibtex`; use native MyST citations (`@key` syntax, see https://mystmd.org/guide/citations)
- [ ] Prefer DOI-based references where possible (MyST resolves them automatically)
- [ ] Migrate `literature.bib` entries to the MyST references format
- [ ] Convert all `.. bibliography::` / `zbibliography` directives and `:cite:` roles to `@key` inline citations

### 5. CI/CD Update
- [ ] Replace `sphinx-build` invocations with `jb build`
- [ ] Update requirements file for docs build

### 6. Final Cleanup
- [ ] Remove `conf.py` and legacy Sphinx config entirely
- [ ] Update `README.md` and `CONTRIBUTING` build instructions
