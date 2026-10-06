# AI disclosure

AI models produced the content of this repository: the reduction, the evaluator and its error analysis, the checkers,
the local certificates, the runs, the figures and the text. The repository owner chose the problem, directed the work
and decided on scope and publication. The owner did not check the mathematics or the code line by line.

**Models**
- Claude Opus 5.5 (Anthropic, via Claude Code) did almost all of the work, as several coordinated agents.
- The independent arb checker (`c/arbcheck/`), the second local certificate (`local2/`) and the first-principles
  check (`firstprinciples/`) were written by separate agent instances. The first two worked from a specification only
  (`SPEC.md` in each directory), the third from the problem statement and the file formats only. Each was told not to
  read the code it checks.
- Reviews were run as separate instances of:
  - Claude Opus 5.5, Claude Fable 5.1 and Claude Sonnet 5 (Anthropic);
  - gpt-6-astra (OpenAI, via the Codex CLI).
- The commits carry a `Co-Authored-By: Claude Opus 5.5` trailer.

**Reviews.** Every component except the first-principles check (`firstprinciples/`, itself a check) was reviewed by at
least one model instance that did not write it.

| Component | Reviews |
|---|---|
| Reduction (`REDUCTION.md`) | a Claude agent; gpt-6-astra; Claude Fable 5.1 |
| Error analysis of the evaluator (`c/ERROR_ANALYSIS_v2.md`) | Claude Opus 5.5; gpt-6-astra; Claude Fable 5.1 |
| arb checker (`c/arbcheck/`) | Claude Fable 5.1, adversarial, with 74 false-claim test cases |
| Coverage check (`c/check_done.py`) and changes to the local certificate | Claude Sonnet 5 |
| Second local certificate (`local2/`), $`d = 6`$ | Claude Fable 5.1 |
| `local2/` generalised to $`d = 5, 6, 7`$, with the joint final step | gpt-6-astra |
| Public repository as a whole | gpt-6-astra (stopped early by a usage limit); Claude Opus 5.5 |

The AI reviews found real errors. All of them were fixed:
- **First evaluator.** gpt-6-astra judged it not rigorous as claimed. It exhibited:
  - an invalid affine enclosure at scale $`2^{-200}`$;
  - a wrong symmetry discard at scale $`2^{-537}`$;
  - a resume merge that could accept a torn checkpoint line;
  - coverage gaps for grid settings other than the one used;
  - reliance on library `log` and `pow` without error bounds.

  The reviewer noted that the two arithmetic counterexamples are at scales that the documented run cannot reach.

  The evaluator was rewritten as `c/smale_bb_v2.c`, with model validity checks, a correctly rounded `log` (CORE-MATH)
  and checkpoints bound to the configuration. The degree-6 run in this repository is the run of that version.
- **Error analysis.** The reviews found 16 write-up gaps in total, from gpt-6-astra and Claude Fable 5.1. Examples are
  the scope of a modulus lemma near overflow, the majorant induction (M), and the stated trust base.
- **Displayed bounds.** In the README, displayed bounds had been rounded in the unsafe direction (gpt-6-astra). A sweep
  of all displayed constants followed. Exactly known values are now written exactly.
- **Reduction.** Missing hypotheses in the corollary ($`i \ne j`$, $`S_i \ne 0`$) and imprecise wording of the equality
  configurations (Claude Fable 5.1).
- **Second local certificate, $`d = 7`$.** Five displayed decimal enclosures in `local2/README.md` had been rounded
  inward, and the run scripts returned exit status 0 after a failure (gpt-6-astra). The certified rational constants
  were unaffected.
- **arb checker.** Non-finite task coordinates were accepted, and a cache reset was in the wrong place, a latent
  hazard (Claude Fable 5.1).

AI reviews are not peer review, and no human expert has checked this work. The reviewing models are not fully
independent of the authoring model: two of the three review families are Claude models. In the terminology of the
Lean community this is a *warrant*, not a human-readable proof (see the note at the top of `README.md`).

**What is checked by software**
- **Branch-and-bound run.** `c/smale_bb_v2.c` covers the compact chart in binary64 ball arithmetic. `c/check_done.py`
  checks in exact rational arithmetic that the finished tasks cover the grid.
- **Subdivision tree.** `c/arbcheck/` re-verifies the claim of every leaf of the exported tree in FLINT/arb ball
  arithmetic. It is independent of the evaluator's floating-point error analysis.
- **Local certificates.** `local/` and `local2/` certify the inequality near the extremal configuration in arb ball
  arithmetic. `local2/` also uses exact arithmetic in $`\mathbb{Q}(\omega)`$.

What remains to be trusted:
- the written arguments in `REDUCTION.md` (Proposition R and its corollary) and the hand-over between the cover and
  the local certificates, which only AI models have checked;
- that the code computes the quantities the written arguments refer to;
- FLINT/arb, python-flint, CORE-MATH, the compilers and the hardware floating-point unit.
