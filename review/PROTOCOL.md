# Cross-agent review loop

Two agents review this work in writing, in the repository, so the exchange is
part of the record: **Codex** (OpenAI) reviews, **Claude** (who wrote the code
and the analysis) answers each point and either fixes it or argues against it.

## Files

```
review/prompts/roundN.md     the prompt given to Codex for round N
review/roundN-codex.md       Codex's findings, verbatim
review/roundN-claude.md      Claude's response: accepted / fixed / disputed
```

Nothing is edited after the fact. A finding that turns out to be wrong stays,
with the argument against it; a finding that was right stays with the commit
that fixed it.

## Running a round

Codex runs read-only, so it can read the repository and run commands but cannot
change anything; its output is captured into the round file.

```bash
cd /home/lgomez/MHD
codex exec -s read-only -C /home/lgomez/MHD "$(cat review/prompts/round1.md)" \
  2>&1 | tee review/round1-codex.md
```

Then Claude reads `review/roundN-codex.md`, writes `review/roundN-claude.md`,
applies the fixes it accepts, and commits both files together with the fixes.

To let Codex write its own file instead of piping (it can then also leave
comments next to code), swap the sandbox for `-s workspace-write` and add
"write your findings to review/roundN-codex.md" to the prompt. Read-only is the
default here because the point of the loop is the argument, not the edit.

## What each round is for

| round | question put to the reviewer |
|---|---|
| 1 | Is the science right, and does the code do what the write-up claims? |
| 2 | Whatever round 1 leaves unresolved, plus the fixes made in response |

## Rules for both sides

- Every claim points at a file, a line, or a number that can be recomputed.
- "I could not verify this" is a valid finding; "this looks fine" is not a review.
- The reviewer ranks findings by whether they change a conclusion, not by taste.
- The author answers every finding: fixed (with the commit), accepted but not
  fixed (with the reason), or disputed (with the evidence).
