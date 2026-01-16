# AI-Assisted Coding Best Practices

This guide covers best practices for working with AI coding assistants (like Cursor's AI) in the FUSILLI project, and how to leverage the project's documentation structure for effective collaboration.

---

## 📋 Table of Contents

1. [Core Principles](#core-principles)
2. [Effective Prompting](#effective-prompting)
3. [Context Management](#context-management)
4. [Iterative Development](#iterative-development)
5. [Code Review & Validation](#code-review--validation)
6. [Using Project Documentation](#using-project-documentation)
7. [Common Patterns](#common-patterns)
8. [Anti-Patterns to Avoid](#anti-patterns-to-avoid)

---

## 🎯 Core Principles

### 1. AI as a Collaborative Tool, Not a Replacement

**Best Practice:**
- Use AI to accelerate development, not replace understanding
- Always review and understand AI-generated code
- Question assumptions and verify logic
- Use AI for repetitive tasks, boilerplate, and pattern matching

**Example:**
```
❌ Bad: "Write a fusion detection pipeline"
✅ Good: "Add a new detection method 'bowtie2' following the pattern in
         process_strings.smk. See TASK-001 in tasks/tasks.md for requirements."
```

### 2. Provide Context, Not Just Instructions

**Best Practice:**
- Reference existing code patterns
- Point to relevant documentation
- Include task IDs and acceptance criteria
- Show similar implementations

**Example:**
```
❌ Bad: "Add error handling"
✅ Good: "Add error handling to fusion_sequences.py following the pattern
         in string_matcher.py (lines 45-60). Handle FileNotFoundError and
         ValueError with context messages as specified in docs/TECHNICAL.md."
```

### 3. Incremental, Focused Changes

**Best Practice:**
- Break large tasks into smaller, focused requests
- Complete and test one change before moving to the next
- Use the task tracking system to manage complexity

**Example:**
```
❌ Bad: "Refactor the entire pipeline to use classes"
✅ Good: "Refactor fusion_sequences.py to extract BreakpointGenerator class.
         Keep the CLI interface unchanged. See TASK-002 subtask 1."
```

---

## 💬 Effective Prompting

### Structure Your Prompts

**Template:**
```
[Context] I'm working on [task/feature] (TASK-XXX)
[Reference] Following the pattern in [file:lines] and [doc section]
[Goal] I need to [specific action]
[Constraints] Must [requirement], should not [limitation]
[Acceptance] Success means [criterion]
```

**Example:**
```
Context: I'm working on adding Bowtie2 alignment (TASK-001)
Reference: Following the pattern in workflow/rules/process_strings.smk
           and the detection method structure in docs/ARCHITECTURE.md
Goal: Create a new rule module process_bowtie.smk that aligns reads and
      detects fusions
Constraints: Must integrate with existing config system, should not break
             string matching method
Acceptance: New method selectable via config, produces same output format,
            tests pass
```

### Use Project-Specific References

**Reference Files:**
- `tasks/tasks.md` - Task requirements and acceptance criteria
- `docs/ARCHITECTURE.md` - System design and extension points
- `docs/TECHNICAL.md` - Code patterns and conventions
- `docs/status.md` - Current project state

**Example:**
```
"Add a new QC metric following the pattern in workflow/rules/qc.smk.
See the 'Adding New QC Tools' section in docs/ARCHITECTURE.md for the
extension pattern. The metric should measure [specific thing] and output
to results/{experiment}/qc_metrics.csv."
```

### Be Specific About Code Style

**Reference Technical Docs:**
```
"Add type hints to all functions in fusion_sequences.py following
docs/TECHNICAL.md conventions. Use Optional[Type] format and include
Google-style docstrings."
```

---

## 📚 Context Management

### How to Orient an AI Agent

**Step 1: Provide Project Context**
```
"Review the README.md to understand this is a Snakemake pipeline for
fusion detection. The main entry point is workflow/Snakefile."
```

**Step 2: Provide Task Context**
```
"See TASK-001 in tasks/tasks.md. I need to implement the Bowtie2 detection
method. Acceptance criteria are listed there."
```

**Step 3: Provide Code Context**
```
"Look at workflow/rules/process_strings.smk to see how string matching
is implemented. I want a similar structure for Bowtie2 alignment."
```

**Step 4: Provide Technical Context**
```
"Follow the patterns in docs/TECHNICAL.md for error handling and logging.
Use pathlib.Path for file paths, and include type hints."
```

### Using the Documentation Hierarchy

```
README.md (entry point)
  ↓
docs/ARCHITECTURE.md (system design)
  ↓
docs/TECHNICAL.md (code patterns)
  ↓
Specific files (implementation details)
```

**Example Prompt:**
```
"I need to add a new preprocessing step. I've reviewed:
- README.md: Understands pipeline purpose
- docs/ARCHITECTURE.md: Knows data flow and module boundaries
- docs/TECHNICAL.md: Knows code style requirements
- workflow/rules/filter_paired.smk: Sees existing preprocessing pattern

Now add a contamination filtering step using bbduk, following the same
pattern as adapter trimming."
```

### Context Files for Different Tasks

**For New Features:**
- `tasks/tasks.md` - Requirements
- `docs/ARCHITECTURE.md` - Extension points
- Similar existing code - Pattern to follow

**For Bug Fixes:**
- Error message/logs
- `docs/TECHNICAL.md` - Error handling patterns
- Related code sections

**For Refactoring:**
- `docs/TECHNICAL.md` - Code organization principles
- `docs/ARCHITECTURE.md` - Module boundaries
- Current implementation

---

## 🔄 Iterative Development

### Workflow with AI

**1. Plan (Human)**
- Define task in `tasks/tasks.md`
- Break into subtasks
- Identify dependencies

**2. Implement (AI-Assisted)**
- Start with one subtask
- Request AI to implement following patterns
- Review and test incrementally

**3. Iterate (Collaborative)**
- Test and validate
- Request refinements
- Move to next subtask

**Example Workflow:**

```
Human: "I've added TASK-001 to tasks/tasks.md. Start with subtask 1:
        Create the Bowtie2 alignment script."

AI: [Implements script]

Human: "Good, but add error handling for missing reference files.
       See docs/TECHNICAL.md error handling section."

AI: [Adds error handling]

Human: "Perfect. Now move to subtask 2: Add the Snakemake rule."
```

### Incremental Validation

**After Each Change:**
1. Review the code
2. Run relevant tests
3. Check against acceptance criteria
4. Update task status

**Example:**
```
"Add the Bowtie2 script. After implementation, verify:
- Type hints are present (docs/TECHNICAL.md)
- Error handling follows patterns
- Can be run standalone
- Tests pass"
```

---

## ✅ Code Review & Validation

### Reviewing AI-Generated Code

**Checklist:**
- [ ] Follows project code style (`docs/TECHNICAL.md`)
- [ ] Matches existing patterns
- [ ] Includes type hints
- [ ] Has appropriate error handling
- [ ] Includes docstrings
- [ ] Handles edge cases
- [ ] Doesn't break existing functionality

### Validation Prompts

**Ask AI to Self-Review:**
```
"Review the code you just generated. Check it against:
1. docs/TECHNICAL.md code style requirements
2. Similar code in [file] for pattern consistency
3. Error handling patterns in docs/TECHNICAL.md
4. Type hints and docstrings

Identify any issues or improvements."
```

**Ask for Test Cases:**
```
"Generate test cases for the function you just created. Include:
- Normal operation
- Edge cases
- Error conditions
Follow the test patterns in tests/test_fusion_sequences.py"
```

### Testing AI-Generated Code

**Always Test:**
```bash
# Run unit tests
pytest tests/test_new_feature.py

# Run integration tests
pytest tests/test_integration.py

# Dry run Snakemake
snakemake -n

# Check linting
black --check workflow/scripts/
isort --check workflow/scripts/
```

**Request AI to Help Debug:**
```
"The test test_bowtie_alignment is failing with [error].
Review the implementation in [file] and the test in [test_file].
What's wrong?"
```

---

## 📖 Using Project Documentation

### For Different Types of Tasks

#### Adding New Features

**Reference These Docs:**
1. `tasks/tasks.md` - Task requirements
2. `docs/ARCHITECTURE.md` - Extension points section
3. `docs/TECHNICAL.md` - Code patterns
4. Similar existing code

**Example Prompt:**
```
"Add a new detection method 'bowtie2'. I've reviewed:
- TASK-001 in tasks/tasks.md for requirements
- 'Extension Points' in docs/ARCHITECTURE.md
- process_strings.smk for the pattern
- docs/TECHNICAL.md for code style

Implement following the same structure as string matching."
```

#### Fixing Bugs

**Reference These Docs:**
1. Error logs
2. `docs/TECHNICAL.md` - Error handling patterns
3. Related code sections
4. `docs/status.md` - Known issues

**Example Prompt:**
```
"Fix the breakpoint generation error. The error is [error message].
Review docs/TECHNICAL.md error handling section and the similar
error handling in string_matcher.py (lines 120-135). Apply the
same pattern."
```

#### Refactoring

**Reference These Docs:**
1. `docs/TECHNICAL.md` - Code organization
2. `docs/ARCHITECTURE.md` - Module boundaries
3. Current implementation

**Example Prompt:**
```
"Refactor fusion_sequences.py to extract a BreakpointGenerator class.
Follow docs/TECHNICAL.md 'Code Organization' section. Keep the CLI
interface unchanged. Maintain all existing functionality."
```

### Documentation as Context

**When Starting a Session:**
```
"I'm working on [feature]. Please review:
- README.md for project overview
- docs/ARCHITECTURE.md for system design
- docs/status.md for current state
- tasks/tasks.md for [TASK-XXX] requirements

Then help me implement [specific part]."
```

**When Continuing Work:**
```
"Continuing TASK-001. Previous work completed [what].
Current status in docs/status.md. Now implement [next part]."
```

---

## 🔧 Common Patterns

### Pattern 1: Adding a New Script

**Prompt Structure:**
```
"Create a new script workflow/scripts/bowtie_matcher.py following the
pattern in string_matcher.py. Requirements:
- Standalone CLI interface (see docs/TECHNICAL.md 'Configuration Loading')
- Type hints for all functions
- Error handling following docs/TECHNICAL.md patterns
- Progress reporting using utils.ProgressReporter
- Output format matches fusion_counts.csv structure"
```

### Pattern 2: Adding a New Rule

**Prompt Structure:**
```
"Add a new rule module workflow/rules/process_bowtie.smk following
process_strings.smk. Requirements:
- Integrate with config system (use DETECTION_METHOD from common.smk)
- Follow rule structure in docs/TECHNICAL.md
- Use helper functions from common.smk
- Include proper logging
- Update get_all_targets() in common.smk if needed"
```

### Pattern 3: Adding Tests

**Prompt Structure:**
```
"Add tests for bowtie_matcher.py in tests/test_bowtie_matcher.py.
Follow patterns in tests/test_string_matcher.py:
- Use fixtures from conftest.py
- Test normal operation, edge cases, error conditions
- Mock file I/O where appropriate
- Aim for 80%+ coverage"
```

### Pattern 4: Updating Documentation

**Prompt Structure:**
```
"Update docs/ARCHITECTURE.md to document the new Bowtie2 detection method.
Add it to:
- System Components section
- Data Flow diagram
- Extension Points section
Follow the existing documentation style."
```

---

## ❌ Anti-Patterns to Avoid

### 1. Vague Requests

```
❌ "Make it better"
✅ "Optimize the string matching loop in string_matcher.py.
    Use set lookup instead of list for partner_ends (see
    docs/TECHNICAL.md Performance section)."
```

### 2. Too Large Scope

```
❌ "Refactor the entire pipeline"
✅ "Refactor fusion_sequences.py to extract BreakpointGenerator class.
    Keep CLI interface unchanged. See TASK-002 subtask 1."
```

### 3. Ignoring Project Patterns

```
❌ "Add a function to process reads" (without context)
✅ "Add a function process_reads() to string_matcher.py following
    the existing pattern. Include type hints and error handling as
    specified in docs/TECHNICAL.md."
```

### 4. Not Reviewing AI Output

```
❌ Accept AI code without review
✅ Review code, run tests, verify against acceptance criteria,
    then integrate
```

### 5. Not Updating Documentation

```
❌ Make changes without updating docs/status.md or tasks/tasks.md
✅ Update task status, document decisions, keep status current
```

### 6. Missing Context

```
❌ "Fix the bug"
✅ "Fix the FileNotFoundError in fusion_sequences.py line 45.
    The issue is [specific problem]. See similar error handling
    in string_matcher.py lines 120-135."
```

---

## 🎓 Learning from AI Suggestions

### When AI Suggests Alternatives

**Evaluate Against:**
- Project patterns (`docs/TECHNICAL.md`)
- Architecture constraints (`docs/ARCHITECTURE.md`)
- Task requirements (`tasks/tasks.md`)

**Example:**
```
AI: "I could use pandas for this, or a dictionary..."

You: "Check docs/TECHNICAL.md. For this use case, we prefer
     dictionaries for performance. Use dict as specified."
```

### When AI Questions Your Approach

**Consider:**
- Is the AI right? Review architecture docs
- Does it fit project patterns?
- Is there a documented reason for current approach?

**Example:**
```
AI: "Why not use classes here?"

You: "Check docs/ARCHITECTURE.md. Scripts are intentionally
     standalone functions for testability. Keep current structure."
```

---

## 🔄 Workflow Integration

### Daily Workflow

**Morning:**
1. Review `docs/status.md` for current state
2. Check `tasks/tasks.md` for priorities
3. Orient AI with current context

**During Work:**
1. Reference docs for patterns
2. Make incremental changes
3. Test frequently
4. Update task status

**End of Day:**
1. Update `docs/status.md`
2. Mark completed tasks
3. Note blockers or issues

### Task Completion Workflow

**When Starting:**
```
"Starting TASK-001. Review:
- tasks/tasks.md for requirements
- docs/ARCHITECTURE.md extension points
- Similar code in [file]

Implement [first subtask]."
```

**When Completing:**
```
"TASK-001 complete. Update:
- tasks/tasks.md: Mark as completed
- docs/status.md: Add to recent completions
- README.md: Update if needed (check requirements)"
```

---

## 📝 Example: Complete Feature Addition

### Scenario: Adding Bowtie2 Detection

**Step 1: Task Definition (Human)**
```markdown
# tasks/tasks.md
- [ ] TASK-001
  - **Title:** Add Bowtie2 alignment detection method
  - **Status:** In Progress
  - **Acceptance Criteria:**
    - [ ] Bowtie2 selectable via config
    - [ ] Produces same output format
    - [ ] Tests pass
    - [ ] Documentation updated
```

**Step 2: Orient AI**
```
"I'm working on TASK-001. Please review:
- tasks/tasks.md for requirements
- docs/ARCHITECTURE.md 'Extension Points' section
- workflow/rules/process_strings.smk for pattern
- docs/TECHNICAL.md for code style

Then create the Bowtie2 alignment script following the string_matcher.py pattern."
```

**Step 3: Implement Script (AI-Assisted)**
```
AI generates: workflow/scripts/bowtie_matcher.py

You review: Check against docs/TECHNICAL.md
- Type hints? ✓
- Error handling? ✓
- Docstrings? ✓
- Follows patterns? ✓
```

**Step 4: Add Rule (AI-Assisted)**
```
"Now add the Snakemake rule in workflow/rules/process_bowtie.smk
following process_strings.smk. Integrate with config system."
```

**Step 5: Add Tests (AI-Assisted)**
```
"Add tests in tests/test_bowtie_matcher.py following
test_string_matcher.py patterns."
```

**Step 6: Update Documentation (AI-Assisted)**
```
"Update docs/ARCHITECTURE.md to document the new method.
Add to System Components and Extension Points sections."
```

**Step 7: Complete Task (Human)**
```markdown
# tasks/tasks.md
- [x] TASK-001
  - **Completed:** 2024-01-20
```

```markdown
# docs/status.md
## Recent Completions
- ✅ Added Bowtie2 detection method (TASK-001)
```

---

## 🎯 Key Takeaways

1. **Always provide context** - Reference docs, tasks, and existing code
2. **Work incrementally** - Small, focused changes are better
3. **Review everything** - AI is a tool, not a replacement for understanding
4. **Use the documentation** - It's there to help both you and AI
5. **Update as you go** - Keep tasks and status current
6. **Test frequently** - Validate each change before moving on
7. **Follow patterns** - Consistency makes code maintainable

---

## 📚 Quick Reference

**Starting a new task:**
1. Add to `tasks/tasks.md`
2. Review relevant docs
3. Orient AI with context
4. Implement incrementally

**Continuing work:**
1. Check `docs/status.md`
2. Review task in `tasks/tasks.md`
3. Continue from last point

**Completing work:**
1. Verify acceptance criteria
2. Run tests
3. Update `tasks/tasks.md`
4. Update `docs/status.md`
5. Update docs if needed

**Getting help:**
1. Reference specific docs
2. Show similar code
3. Provide error messages
4. State what you've tried

---

This guide works in conjunction with the project management structure to create an effective AI-assisted development workflow. The key is providing context, working incrementally, and always validating AI output against project standards.
