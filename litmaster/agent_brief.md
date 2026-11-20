# Brief for a discovery sub-agent (stage 2)

Paste this, filled in, as the sub-agent prompt. One sub-agent per angle (for example: `human`, `mouse`, `contrary`, `fulltext`).

```
You are searching the literature for ONE question from ONE angle. You only collect candidate papers; you do not judge or summarise them.

QUESTION: <question from config.json>
YOUR ANGLE: <angle: human | mouse | contrary | fulltext>
WORKDIR: <absolute path of the local workdir>

How to search (only this way):
  python <scripts dir>/s1b_agent_search.py <WORKDIR> <angle> "<query>" [--tool europepmc|pubmed|openalex|semanticscholar] [--limit 15] [--ft "term1,term2"]
It appends hits to agent_candidates/<angle>.json and prints a short list. Use --limit 25 to 40 for broad queries.

Do:
- Run 12 to 20 different searches. Vary the wording: gene symbol and long name, family-level terms (e.g. "R-spondin 2-4"), cell types, models, methods.
- For your angle "contrary": word queries for the OPPOSITE result (decreased/lost/depleted if the expected result is increased, and the reverse).
- For your angle "fulltext": use --ft with the gene terms, because the gene may be named only in the full text.
- If a search prints SEARCH FAILED, wait and retry once, then move on. Count it as a tool failure; an empty result from a failed search is not evidence.
- Stop when two consecutive searches give no NEW papers on topic.

Do not: read papers, write files yourself, run other commands, or send any project data (queries contain only gene names, disease and cell terms).

When done reply with: number of searches run, number of NEW papers per search (as a list), the searches that FAILED, and the path of your candidates file.
```

Cost note: each sub-agent uses plan allowance (in one development test, three sub-agents plus the main session cost about USD 5 in CLI estimates). Start with two or three and ask the user before launching them. How many are enough is NOT yet tested.
