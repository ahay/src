---
name: finding-workflow-examples
description: Use when writing an end-to-end Madagascar processing flow for a named geophysical task (NMO, migration, well-tie, denoising, dataset fetch, etc.) — teaches how to discover the right recipe in book/ and adapt it.
---

# Finding and adapting workflow examples from book/

Madagascar ships with `book/` — a corpus of ~1,747 runnable SConstructs covering published research and teaching material. Every common geophysical workflow you are likely asked to write has a canonical example in there. This skill teaches you to find the right one and adapt it, instead of reinventing the flow.

**The catalog is the starting point.** See `CATALOG.md` next to this file. It is grouped into **Workflows** (what-to-do) and **Datasets** (what-to-load). Each entry names one concrete path into `book/`.

## Discovery workflow (three steps)

1. **Read `CATALOG.md`.** Scan the headings, find the entry that matches your task. If nothing fits, go to step 3.
2. **Open the real SConstruct.** The catalog's one-line description is a pointer, not content. The authoritative source is the file itself. Read it top-to-bottom before adapting anything.
3. **Grep fallback when no catalog entry fits.** Useful patterns:
   - By `sf*` program: `grep -rln "sfnmo\b" /Users/jgoai/m8r/src/book/`
   - By imported dataset helper: `grep -rln "import sigsbee" /Users/jgoai/m8r/src/book/`
   - By directory-name keyword: `find /Users/jgoai/m8r/src/book -type d -iname "*kirchhoff*"`
   - By `Fetch()` target: `grep -rln "Fetch.*marmvel" /Users/jgoai/m8r/src/book/`

## Adaptation pattern

Most `book/` SConstructs encode two kinds of knowledge: **structural** (which `sf*` programs chain together, in what order, with what axis conventions) and **parametric** (specific `n1`, `d1`, `nt`, labels, dataset paths). Adapt parametric content; preserve structural content.

1. **Strip any paper-rendering shell.** If the top of the SConstruct is `from rsf.tex import *` and the bottom is `End(color='...')` or similar, that file is a LaTeX paper wrapper, not a processing flow. The real flow lives in a subdirectory of the same folder — check the child directories. Swap in `from rsf.proj import *` for a regular processing SConstruct.
2. **Copy the relevant `Flow()` / `Plot()` / `Result()` / `Fetch()` calls.** Keep them in the same order. The author tuned this order.
3. **Replace parameters.** `n1`, `d1`, `o1`, `nt`, `dt`, velocities, dataset paths, and labels are case-specific. Replace them with the user's values. Do not guess — ask the user if a parameter is not specified.
4. **Resolve `book/Recipes/` imports.** Some SConstructs do `from rsf.recipes.helderiv import Helderiv` or import local helpers like `import awe, wplot`. If the recipe is short, inline the Flow calls it produces. If long, vendor the helper file into the user's project directory.
5. **Check `Fetch()` servers.** `Fetch('foo','marm')` and `Fetch('foo.tgz','freeusp')` use named servers configured in the Madagascar install — they work if the user has the standard Madagascar setup. Flag external servers (`Fetch(..., server='https://...')`) in case the user needs to adjust.
6. **Preserve plot labels and units.** Wrong `label1=`/`unit1=` on a plot is how small bugs hide in the output.

## When CATALOG.md has nothing for the task

Do the grep fallback (step 3 above). If still nothing useful, tell the user the corpus does not appear to have a direct match, describe what you found, and propose building the flow from primitives using the `using-sf-programs` and `writing-rsf-flows` skills. Do not fabricate.

## What this skill is not

- Not a tutorial for the `sf*` programs — see `using-sf-programs`.
- Not a tutorial for the rsf.proj DSL — see `writing-rsf-flows`.
- Not an exhaustive index of book/ — the catalog is curated. If you need broader search, grep `book/` directly.
