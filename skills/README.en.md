# How to use the metacellcnv-bootstrap skill

`metacellcnv-bootstrap.ja.SKILL.md` / `metacellcnv-bootstrap.en.SKILL.md` in this directory are
"skills" (agent instruction files) that tell an LLM coding agent how to take nothing but a
CellRanger output directory and, end to end: build an isolated Python environment and install
dependencies, fetch [metacellcnv](https://github.com/takaho/metacellcnv), auto-resolve the gene
annotation (GTF) and chromosome map from a genome name, run the full pipeline, and report
progress/ETA along the way. The Japanese and English versions are identical in content.

> **Verification status**: this skill has **only been tested in Claude (Claude Code / Cowork)**.
> The instructions below for Codex, Antigravity, and OpenCode describe each tool's *general*
> mechanism for persistent custom instructions — **this specific skill has not actually been
> run through those tools**. These tools change quickly, so check each one's own current docs
> before relying on the details below.

## Claude (verified)

**Claude Code**:
```
~/.claude/skills/metacellcnv-bootstrap/SKILL.md         # personal, available in every project
<project>/.claude/skills/metacellcnv-bootstrap/SKILL.md  # project-scoped
```
Copy either `metacellcnv-bootstrap.ja.SKILL.md` or `.en.SKILL.md` into one of these locations
under the name `SKILL.md` (to keep both languages available, use separate directory names, e.g.
`metacellcnv-bootstrap-ja/` and `metacellcnv-bootstrap-en/`). Once placed, just ask naturally
("get CNV calls from this CellRanger output") and it triggers automatically.

**Cowork (this chat)**: you can paste the downloaded file straight into a conversation, or, if
your organization's settings allow it, save it from the file card as a skill.

## OpenAI Codex CLI (not verified)

Codex supports `AGENTS.md` at the repository root (persistent project instructions, always
loaded), plus a "custom prompts" feature: markdown files under `~/.codex/prompts/` (user-wide)
or `.codex/prompts/` (project-only) each become an invocable slash command.

- **Always-on**: paste the body of `metacellcnv-bootstrap.en.SKILL.md` (the English version is
  recommended — Codex tends to follow English instructions more reliably; you can drop the
  leading `---`-fenced YAML block) into `AGENTS.md`.
- **Invoke on demand**: save the body as `~/.codex/prompts/metacellcnv-bootstrap.md` (or
  `.codex/prompts/metacellcnv-bootstrap.md`) and call it with `/metacellcnv-bootstrap`.

Exact directory names/behavior can change — check `codex --help` or the official docs at
developers.openai.com/codex for the current convention.

## Google Antigravity (not verified)

Antigravity has a project/workspace/user-level "Rules" mechanism for persistent instructions.
Register the body as a Rules file, or, if unsure, simply paste the body of
`metacellcnv-bootstrap.en.SKILL.md` as the first message of a new agent session and ask it to
follow it. This area of Antigravity evolves quickly, so check antigravity.google/docs/rules for
the current mechanism and exact file locations.

## OpenCode (not verified)

OpenCode also supports `AGENTS.md` (project root or global) for persistent instructions, plus
markdown+frontmatter-based "Commands" (`.opencode/command/*.md`, invoked as a slash command) and
"Agents" (`.opencode/agent/*.md`).

- **Always-on**: append the body to `AGENTS.md`.
- **Invoke on demand**: save it as `.opencode/command/metacellcnv-bootstrap.md` and call it with
  `/metacellcnv-bootstrap`.

## A fallback that works with any agent

If none of the above mechanisms is available or clear, the body of
`metacellcnv-bootstrap.en.SKILL.md` (or `.ja.`) is just plain-language instructions — any agent
with shell and file access should be able to follow the same procedure if you simply **paste the
file's contents as the first message of a new task**. The leading YAML frontmatter (the
`---`-fenced `name`/`description` block) is a Claude-specific convention and can safely be
ignored by any tool that doesn't use it.
