# Notes

Working notes for Claude sessions (and anyone else) on this project. Read these before starting work;
update them when you learn something non-obvious or finish a chunk of work.

- `codebase.md` - how the code is laid out and the conventions/gotchas that aren't obvious from reading it
- `open_issues.md` - known bugs, modeling gaps and decisions still waiting on Brendan
- `log.md` - dated record of what was changed and why (newest first)

## How Brendan wants to work
A primary goal of this project is Brendan's own learning. Be hesitant about mass edits:
- **Physics / engineering** (correlations, models, equations, design choices): explain the issue and the fix,
  point to the lines, and let him make the change unless he asks you to. Teach rather than do.
- **Software / math structure** (refactors, solvers, plumbing, typos, crashes): fine to edit more freely, but
  still keep changes small and say what changed.
When in doubt, describe instead of edit.

Also remove fixed items from `open_issues.md`; e.g. the throat arc now uses `rc_throat` (fixed by Brendan 2026-10-23).

Keep entries short. Put facts the code or git history already records elsewhere; put *why* and *what's
still wrong* here. Remove items once they're resolved.
