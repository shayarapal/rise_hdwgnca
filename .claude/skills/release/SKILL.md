---
name: release
description: Cut a release — bump the version, write the CHANGELOG entry, commit, tag, and push to GitHub. Use when the user says "release", "cut a version", "ship this", "push these changes", or asks to update the changelog and version after a batch of work.
---

# Release

Bundles a batch of work into a versioned, tagged, pushed commit. Run this after a
meaningful change lands — not after every edit.

## Before anything else: guard the data

This repo sits next to a ~11 GB dataset. **Never let it reach GitHub.**

```bash
git add -A
# Nothing here may exceed 50 MB. If it does, STOP and fix .gitignore.
git diff --cached --name-only -z | xargs -0 -I{} sh -c \
  '[ -f "{}" ] && [ $(stat -f%z "{}") -gt 52428800 ] && echo "TOO BIG: {}"'
```

`.gitignore` already excludes `*.rds`, `*.tar`, `*.mtx`, `*.h5Seurat`, and the
`DATA_UNZIPPED/`, `DATA_GSE233866/` and `analysis_results/` directories outright — as of
2026-09-07 all three live in `$SHARED_DIR` outside the repo and nothing in them is tracked.
If a new data format appears, add it there rather than committing it.

## 1. Decide the version

Read the current version from the top non-`[Unreleased]` heading in `CHANGELOG.md`, then bump
per [SemVer](https://semver.org/):

| Bump | When |
|---|---|
| **major** (`1.0.0`) | Breaking change to an endpoint contract, request/response shape, or input format |
| **minor** (`0.3.0`) | New endpoint, new supported input, new script or capability |
| **patch** (`0.2.1`) | Bug fix, doc fix, or internal change with no new capability |

If the bump is ambiguous, ask rather than guess.

## 2. Write the CHANGELOG entry

Edit `CHANGELOG.md`. Insert the new section directly below `## [Unreleased]`, following the
existing [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) structure. Use only the
subsections that apply: `### Added`, `### Fixed`, `### Changed`, `### Removed`, `### Notes`.

Rules that matter for this repo:

- **State the actual defect, not the symptom.** "SeuratDisk calls `GetAssayData(slot=)`, made
  defunct in SeuratObject 5.0.0" beats "fixed a loading bug."
- Reference real files and functions by name, so a lab member can find them.
- Anything known-broken goes under `## [Unreleased]` → `### Known issues`. Do not quietly ship
  a problem you know about.
- Date the entry with today's real date (`date +%Y-%m-%d`), not a guess.

## 3. Commit, tag, push

```bash
git add -A
git commit -F - <<'EOF'
<subject: imperative, <=72 chars, no trailing period>

<body: what was broken and why the change was necessary. Explain the reasoning a
reviewer could not reconstruct from the diff. Wrap at ~80 columns.>

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
EOF

git tag -a vX.Y.Z -m "vX.Y.Z — <one-line summary>"
git push origin main
git push origin vX.Y.Z
```

Confirm both pushes succeeded, then report the tag and the commit SHA.

## Do not

- Force-push, amend, or rewrite history on `main`.
- Commit anything from `DATA_UNZIPPED/`, `DATA_GSE233866/` or `analysis_results/` — those
  live outside the repo now and must stay untracked.
- Tag a version that already exists (`git tag` to check first).
- Invent changelog entries for work you did not verify actually happened — read the diff.
