# Project Management Guide

This guide covers branch management, task design, agent orientation, and overall project organization for FUSILLI.

---

## 📋 Table of Contents

1. [Branch Management](#branch-management)
2. [Task Design & Specification](#task-design--specification)
3. [Orienting AI Agents](#orienting-ai-agents)
4. [Subagents & Work Organization](#subagents--work-organization)
5. [Overall Organization](#overall-organization)

---

## 🌳 Branch Management

### Strategy

**Main Branches:**
- `main` - Production-ready, stable code
- `develop` - Integration branch for features (optional, use if multiple concurrent features)

**Feature Branches:**
- `feature/description` - New features (e.g., `feature/add-bowtie-alignment`)
- `fix/description` - Bug fixes (e.g., `fix/breakpoint-window-calculation`)
- `refactor/description` - Code refactoring (e.g., `refactor/string-matcher-optimization`)
- `docs/description` - Documentation updates (e.g., `docs/add-performance-tuning-guide`)

**Workflow:**
```bash
# Start new feature
git checkout -b feature/my-feature develop  # or main if no develop branch

# Work and commit frequently
git add .
git commit -m "feat: add initial implementation"

# Keep branch updated
git fetch origin
git rebase origin/develop  # or main

# When ready, push and create PR
git push origin feature/my-feature
```

### Branch Naming Conventions

- Use lowercase with hyphens
- Prefix with type: `feature/`, `fix/`, `refactor/`, `docs/`, `test/`
- Be descriptive but concise: `feature/add-bowtie-alignment` not `feature/new-stuff`

### When to Create Branches

**Create a branch for:**
- Any feature that takes >1 commit
- Bug fixes that might need review
- Refactoring that affects multiple files
- Experimental work you might abandon

**Work directly on main for:**
- Small documentation fixes
- Typo corrections
- Trivial config updates (if you're the only developer)

---

## 📝 Task Design & Specification

### Task Documentation Structure

Create and maintain these files:

#### `tasks/tasks.md` - Active Task List

```markdown
# Tasks

## In Progress
- [ ] Task ID: TASK-001
  - **Title:** Add Bowtie2 alignment option
  - **Status:** In Progress
  - **Assignee:** @username
  - **Priority:** High
  - **Dependencies:** None
  - **Acceptance Criteria:**
    - [ ] Bowtie2 can be selected via config
    - [ ] Alignment results match string matching counts
    - [ ] Tests pass
    - [ ] Documentation updated
  - **Notes:** Need to benchmark performance vs string matching

## Backlog
- [ ] Task ID: TASK-002
  - **Title:** Optimize string matcher for large libraries
  - **Status:** Backlog
  - **Priority:** Medium
  - **Dependencies:** None

## Completed
- [x] Task ID: TASK-000
  - **Title:** Initial pipeline implementation
  - **Completed:** 2024-01-15
```

#### `docs/status.md` - Project Status

```markdown
# Project Status

**Last Updated:** 2024-01-20

## Current Sprint Focus
- Performance optimization
- Adding alignment-based detection

## Recent Completions
- ✅ Reproducibility metadata capture
- ✅ MultiQC integration

## Blockers
- None currently

## Next Steps
1. Complete Bowtie2 integration
2. Benchmark alignment vs string matching
3. Update documentation
```

#### `docs/PRD.md` - Product Requirements Document (for major features)

```markdown
# Product Requirements: [Feature Name]

## Overview
Brief description of what this feature does and why it's needed.

## Goals
- Primary goal
- Secondary goal

## User Stories
- As a [user type], I want [action] so that [benefit]

## Technical Requirements
- Requirement 1
- Requirement 2

## Acceptance Criteria
- [ ] Criterion 1
- [ ] Criterion 2

## Dependencies
- Dependency 1
- Dependency 2

## Out of Scope
- What we're explicitly not doing

## Success Metrics
- Metric 1
- Metric 2
```

### Task Specification Best Practices

1. **Be Specific:** "Add Bowtie2 alignment" not "Improve detection"
2. **Include Context:** Why is this needed? What problem does it solve?
3. **Define Acceptance Criteria:** Clear, testable conditions for completion
4. **List Dependencies:** What needs to be done first?
5. **Estimate Effort:** Small/Medium/Large or time estimate
6. **Link to Issues:** Reference GitHub issues or related tasks

### Task Lifecycle

```
Backlog → In Progress → Review → Completed
                ↓
            Blocked
```

Update `tasks/tasks.md` when status changes.

---

## 🤖 Orienting AI Agents

> **See also:** [AI-Assisted Coding Guide](AI_ASSISTED_CODING.md) for detailed best practices on working with AI assistants.

### Context Files for AI Assistants

Create these files to help AI understand your project:

#### `.cursorrules` (Already exists)
- Project-specific rules and conventions
- Keep this updated as patterns emerge

#### `docs/ARCHITECTURE.md` - System Architecture

```markdown
# Architecture Overview

## System Components

### Workflow Layer (Snakemake)
- **Entry Point:** `workflow/Snakefile`
- **Rule Modules:** `workflow/rules/*.smk`
- **Scripts:** `workflow/scripts/*.py`

### Configuration
- **Main Config:** `config/config.yaml`
- **Schemas:** `workflow/schemas/*.schema.yaml`
- **Validation:** Automatic via Snakemake

### Data Flow
1. FASTQ input → Preprocessing (bbduk, bbmerge)
2. Reference generation (fusion_sequences.py)
3. Detection (string_matcher.py)
4. Aggregation → Output

## Key Design Decisions
- Why Snakemake: Reproducibility, parallelization
- Why string matching: Fast, sufficient for known breakpoints
- Why modular rules: Maintainability

## Module Boundaries
- Rules should not import from scripts directly
- Scripts are standalone and testable
- Config is the single source of truth
```

#### `docs/TECHNICAL.md` - Technical Patterns

```markdown
# Technical Patterns & Conventions

## Code Style
- Python: Black formatting, isort imports
- Maximum line length: 88 characters
- Type hints required for all functions

## Testing
- Unit tests in `tests/`
- Integration tests for full pipeline
- Use pytest fixtures from `conftest.py`

## Error Handling
- Use specific exceptions
- Log errors with context
- Fail fast for configuration errors

## Performance
- Use pyfastx for FASTQ parsing when available
- Parallelize with Snakemake resources
- Profile before optimizing

## Common Patterns
- Config loading: Use `common.smk` functions
- File paths: Use `pathlib.Path`
- Progress: Use `tqdm` for long operations
```

### Agent Orientation Checklist

When starting work with an AI agent, ensure:

1. **Project Context:**
   - Agent has read `README.md`
   - Agent understands the pipeline purpose
   - Agent knows the tech stack (Snakemake, Python)

2. **Current State:**
   - Agent has reviewed `docs/status.md`
   - Agent knows what's in progress
   - Agent understands any blockers

3. **Task Context:**
   - Agent has the specific task from `tasks/tasks.md`
   - Agent understands acceptance criteria
   - Agent knows dependencies

4. **Codebase Context:**
   - Agent understands relevant modules
   - Agent knows existing patterns
   - Agent has reviewed similar implementations

### Providing Context to Agents

**Good prompts:**
```
"Add a new detection method 'bowtie2' following the pattern in
process_strings.smk. See TASK-001 in tasks/tasks.md for requirements."
```

**Bad prompts:**
```
"Add bowtie2"
```

**Include in prompts:**
- Task ID or reference
- Relevant file paths
- Similar existing code patterns
- Constraints or requirements

---

## 🔀 Subagents & Work Organization

### When to Use Subagents

**Use subagents for:**
- Large, independent features that can be developed in parallel
- Different expertise areas (e.g., one for Snakemake rules, one for Python scripts)
- Separate concerns (e.g., one for implementation, one for tests)

**Don't use subagents for:**
- Small tasks (< 1 day)
- Tightly coupled work
- Tasks requiring frequent coordination

### Organizing Work with Subagents

#### Approach 1: Feature-Based Split

```
Main Agent: Orchestrates, integrates
├── Subagent A: Implements core feature (workflow/scripts/feature.py)
├── Subagent B: Adds Snakemake rules (workflow/rules/feature.smk)
└── Subagent C: Writes tests (tests/test_feature.py)
```

#### Approach 2: Layer-Based Split

```
Main Agent: Overall architecture, integration
├── Subagent A: Data processing layer
├── Subagent B: Detection algorithms
└── Subagent C: Output formatting
```

### Subagent Communication

**Shared Context Files:**
- `docs/status.md` - Current work status
- `tasks/tasks.md` - Task assignments
- `docs/ARCHITECTURE.md` - System design

**Coordination:**
- Update `docs/status.md` when starting/blocking/completing
- Mark dependencies clearly in `tasks/tasks.md`
- Use clear commit messages referencing task IDs

### Subagent Task Assignment

In `tasks/tasks.md`, assign tasks:

```markdown
- [ ] Task ID: TASK-001
  - **Title:** Add Bowtie2 alignment
  - **Status:** In Progress
  - **Assignee:** Subagent-A (Core Implementation)
  - **Sub-tasks:**
    - [ ] Subagent-A: Implement alignment script
    - [ ] Subagent-B: Add Snakemake rules
    - [ ] Subagent-C: Write tests
```

---

## 📁 Overall Organization

### Directory Structure

```
fusilli/
├── config/              # Configuration files
│   ├── config.yaml      # Main config
│   ├── samples.csv      # Sample metadata
│   └── fusion_partners.csv
├── docs/                # Documentation
│   ├── ARCHITECTURE.md
│   ├── TECHNICAL.md
│   ├── PROJECT_MANAGEMENT.md (this file)
│   └── status.md
├── tasks/               # Task management
│   └── tasks.md
├── workflow/            # Snakemake pipeline
│   ├── Snakefile
│   ├── rules/           # Rule modules
│   ├── scripts/         # Python scripts
│   └── schemas/         # Config schemas
├── tests/               # Test suite
├── references/          # Reference sequences
├── resources/           # Adapter sequences, etc.
├── results/             # Pipeline outputs (gitignored)
├── logs/                # Pipeline logs (gitignored)
└── stats/               # QC stats (gitignored)
```

### File Organization Principles

1. **Separation of Concerns:**
   - Config: Configuration only
   - Docs: Documentation only
   - Tasks: Task tracking only
   - Code: Implementation only

2. **Modularity:**
   - Small, focused files (< 300 lines)
   - Clear module boundaries
   - Minimal cross-dependencies

3. **Discoverability:**
   - Clear naming conventions
   - Logical directory structure
   - README files where helpful

### Documentation Hierarchy

```
README.md (entry point)
  ├── docs/ARCHITECTURE.md (system design)
  ├── docs/TECHNICAL.md (patterns, conventions)
  ├── docs/PROJECT_MANAGEMENT.md (this file)
  └── docs/status.md (current state)
```

### Keeping Things Updated

**Update when:**
- Starting new work → Update `tasks/tasks.md`, `docs/status.md`
- Completing work → Mark complete, update status
- Architecture changes → Update `docs/ARCHITECTURE.md`
- New patterns emerge → Update `docs/TECHNICAL.md`
- Blockers arise → Update `docs/status.md`

**Review periodically:**
- Weekly: `docs/status.md`, `tasks/tasks.md`
- Monthly: `docs/ARCHITECTURE.md`, `docs/TECHNICAL.md`
- As needed: `README.md`

---

## 🎯 Quick Reference

### Starting a New Task

1. Create/update task in `tasks/tasks.md`
2. Create feature branch: `git checkout -b feature/task-name`
3. Update `docs/status.md` with current focus
4. Orient AI agent with context
5. Work and commit frequently
6. Update task status as you progress
7. Mark complete when done

### Working with AI Agents

1. Provide task ID and reference `tasks/tasks.md`
2. Point to relevant documentation
3. Show similar code patterns
4. Specify acceptance criteria
5. Review AI output critically

### Coordinating Multiple Agents

1. Assign clear tasks in `tasks/tasks.md`
2. Update `docs/status.md` frequently
3. Use clear commit messages with task IDs
4. Communicate blockers immediately
5. Integrate work incrementally

---

## 📚 Additional Resources

- [Git Flow](https://nvie.com/posts/a-successful-git-branching-model/)
- [Conventional Commits](https://www.conventionalcommits.org/)
- [Snakemake Best Practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html)
