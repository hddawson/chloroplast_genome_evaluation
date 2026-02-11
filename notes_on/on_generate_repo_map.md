    ## Run: 2026-02-10T09:22:53

    **Script**: `generate_repo_map.py`
    **Interpreter**: `python`
    **Hash**: `77b98c3968a1`

    ### Summary
    This script generates a comprehensive repository map by scanning a project directory and creating a markdown file (REPO_MAP.md) that helps LLMs navigate the codebase. It classifies files into three categories (readable text files, viewable binaries like images/PDFs, and opaque binaries) based on file extensions, then estimates token counts for text files using a character-to-token ratio. For readable files, it calls the Claude API to generate concise summaries describing each file's purpose, inputs/outputs, and role in the computational pipeline. The final output includes file statistics, context window usage estimates for various LLMs, a categorized file tree table, and individual file summaries to provide both structural overview and detailed content descriptions.

    ### Script
    ```
    #!/usr/bin/env python3
"""
generate_repo_map.py

Scans a project repo, classifies files by LLM-readability, estimates token counts,
calls Claude API to summarize each readable file, and writes REPO_MAP.md.

Usage:
    export ANTHROPIC_API_KEY="sk-..."
    python generate_repo_map.py /path/to/repo

Or:
    python generate_repo_map.py /path/to/repo --no-summarize   # skip API calls, structure only
"""

import os
import sys
import json
import argparse
import time
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

READABLE_EXTENSIONS = {
    # Code
    ".r", ".R", ".py", ".sh", ".jl", ".cpp", ".c", ".h", ".java", ".rs",
    # Notebooks (JSON-based, readable)
    ".ipynb",
    # Data (text-based)
    ".csv", ".tsv", ".txt", ".json", ".yaml", ".yml", ".toml", ".xml", ".bed",
    ".fasta", ".fa", ".fna", ".faa", ".fastq", ".fq", ".gff", ".gff3", ".gtf",
    ".vcf", ".sam", ".nwk", ".nex", ".phy",
    # Docs
    ".md", ".rst", ".tex", ".html", ".htm", ".bib",
    # Config
    ".cfg", ".ini", ".conf", ".env",
    # Log
    ".log",
}

BINARY_BUT_VIEWABLE = {".pdf", ".png", ".jpg", ".jpeg", ".gif", ".webp", ".svg"}

SKIP_DIRS = {".git", "__pycache__", "node_modules", ".pixi", ".snakemake", "env", "venv", ".venv"}

TOKENS_PER_CHAR = 0.25  # rough estimate

MODEL = "claude-sonnet-4-20250514"
MAX_SUMMARY_TOKENS = 300

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def estimate_tokens(text: str) -> int:
    return int(len(text) * TOKENS_PER_CHAR)


def classify_file(path: Path) -> str:
    """Returns 'readable', 'viewable', or 'opaque'."""
    ext = path.suffix.lower()
    if ext in READABLE_EXTENSIONS:
        return "readable"
    elif ext in BINARY_BUT_VIEWABLE:
        return "viewable"
    else:
        return "opaque"


def read_file_safe(path: Path, max_chars: int = 50_000) -> str | None:
    """Read a text file, return None if binary/too large."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
        if len(text) > max_chars:
            text = text[:max_chars] + f"\n\n... [truncated at {max_chars} chars]"
        return text
    except Exception:
        return None


def summarize_with_claude(filepath: str, content: str, api_key: str) -> str:
    """Call Claude API to get a one-paragraph summary of a file."""
    import urllib.request

    prompt = f"""You are summarizing a file in a computational genetics research repo for a REPO_MAP that helps LLMs navigate the project.

File: {filepath}

Respond with ONLY:
1. One sentence: what this file does / contains.
2. Key details (inputs, outputs, methods, dependencies) as a brief list.
3. If it's a script, note the rough pipeline stage (preprocessing, analysis, visualization, etc).

Keep it under 100 words total.

<file_content>
{content[:30000]}
</file_content>"""

    body = json.dumps({
        "model": MODEL,
        "max_tokens": MAX_SUMMARY_TOKENS,
        "messages": [{"role": "user", "content": prompt}],
    }).encode()

    req = urllib.request.Request(
        "https://api.anthropic.com/v1/messages",
        data=body,
        headers={
            "Content-Type": "application/json",
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
        },
    )

    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read())
            return data["content"][0]["text"].strip()
    except Exception as e:
        return f"[summary failed: {e}]"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def scan_repo(repo_root: Path) -> list[dict]:
    """Walk the repo and collect file metadata."""
    files = []
    for dirpath, dirnames, filenames in os.walk(repo_root):
        # Prune skipped dirs in-place
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fname in sorted(filenames):
            fpath = Path(dirpath) / fname
            rel = fpath.relative_to(repo_root)
            category = classify_file(fpath)
            size = fpath.stat().st_size

            entry = {
                "path": str(rel),
                "category": category,
                "size_bytes": size,
                "ext": fpath.suffix.lower(),
                "tokens_est": None,
                "summary": None,
            }

            if category == "readable":
                content = read_file_safe(fpath)
                if content is not None:
                    entry["tokens_est"] = estimate_tokens(content)
                    entry["_content"] = content  # stash for summarization

            files.append(entry)
    return files


def build_tree_string(files: list[dict]) -> str:
    """Build a markdown table of all files."""
    tag = {"readable": "✅", "viewable": "👁️", "opaque": "❌"}

    lines = [
        "| File | Directory | Type | LLM Access | Size (KB) | Est. Tokens |",
        "|------|-----------|------|------------|-----------|-------------|",
    ]

    for f in sorted(files, key=lambda x: x["path"]):
        parts = Path(f["path"]).parts
        name = parts[-1]
        d = str(Path(*parts[:-1])) if len(parts) > 1 else "."
        ext = f["ext"] or "—"
        access = tag[f["category"]]
        size_kb = f"{f['size_bytes'] / 1024:.1f}"
        tok = f"~{f['tokens_est']:,}" if f["tokens_est"] else "—"
        lines.append(f"| `{name}` | `{d}` | {ext} | {access} | {size_kb} | {tok} |")

    lines.append("")
    lines.append("*✅ = readable, 👁️ = viewable (image/PDF), ❌ = opaque binary*")

    return "\n".join(lines)


def build_repo_map(repo_root: Path, files: list[dict], summarize: bool) -> str:
    """Assemble the REPO_MAP.md content."""
    readable = [f for f in files if f["category"] == "readable"]
    viewable = [f for f in files if f["category"] == "viewable"]
    opaque = [f for f in files if f["category"] == "opaque"]

    total_tokens = sum(f["tokens_est"] or 0 for f in files)
    readable_tokens = sum(f["tokens_est"] or 0 for f in readable)

    md = []
    md.append(f"# REPO MAP — {repo_root.name}")
    md.append(f"\n*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}*\n")

    # Coverage summary
    md.append("## LLM Context Coverage\n")
    md.append(f"| Metric | Value |")
    md.append(f"|--------|-------|")
    md.append(f"| Total files | {len(files)} |")
    md.append(f"| LLM-readable (text) | {len(readable)} ({100*len(readable)/max(len(files),1):.0f}%) |")
    md.append(f"| Viewable (images/PDF) | {len(viewable)} |")
    md.append(f"| Opaque (binary) | {len(opaque)} |")
    md.append(f"| Est. readable tokens | ~{readable_tokens:,} |")
    md.append(f"| Est. total tokens | ~{total_tokens:,} |")

    # Context window estimates for common models
    context_windows = {
        "Claude Opus 4.6 / Sonnet 4.5": 200_000,
        "Claude Haiku 4.5": 200_000,
        "GPT-4o": 128_000,
        "Gemini 2.5 Pro": 1_000_000,
    }
    md.append(f"| | |")
    md.append(f"| **Model context usage** | |")
    for model_name, ctx in context_windows.items():
        pct = 100 * readable_tokens / max(ctx, 1)
        bar = "🟢" if pct < 50 else "🟡" if pct < 80 else "🔴"
        md.append(f"| {bar} {model_name} ({ctx//1000}k) | {pct:.1f}% used |")

    md.append("")

    if opaque:
        md.append("### ⚠️ Opaque files (LLM cannot read these)\n")
        for f in opaque:
            md.append(f"- `{f['path']}` ({f['ext']}, {f['size_bytes']/1024:.1f} KB)")
        md.append("")

    # Annotated tree
    md.append("## File Tree\n")
    md.append(build_tree_string(files))
    md.append("")

    # Summaries
    if summarize:
        md.append("## File Summaries\n")
        for f in sorted(readable, key=lambda x: x["path"]):
            md.append(f"### `{f['path']}`\n")
            md.append(f"{f['summary'] or '*No summary available.*'}\n")

    return "\n".join(md)


def main():
    parser = argparse.ArgumentParser(description="Generate REPO_MAP.md for LLM navigation")
    parser.add_argument("repo", help="Path to repo root")
    parser.add_argument("--no-summarize", action="store_true", help="Skip Claude API calls")
    parser.add_argument("--output", default=None, help="Output path (default: REPO_MAP.md in repo root)")
    args = parser.parse_args()

    repo_root = Path(args.repo).resolve()
    assert repo_root.is_dir(), f"Not a directory: {repo_root}"

    api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    summarize = not args.no_summarize

    if summarize:
        assert api_key, "Set ANTHROPIC_API_KEY env var or use --no-summarize"

    print(f"Scanning {repo_root} ...")
    files = scan_repo(repo_root)
    print(f"Found {len(files)} files.")

    if summarize:
        readable = [f for f in files if f["category"] == "readable" and f.get("_content")]
        print(f"Summarizing {len(readable)} readable files via Claude API...")
        for i, f in enumerate(readable):
            print(f"  [{i+1}/{len(readable)}] {f['path']}")
            f["summary"] = summarize_with_claude(f["path"], f["_content"], api_key)
            time.sleep(0.5)  # light rate limiting

    # Clean up stashed content
    for f in files:
        f.pop("_content", None)

    output_path = Path(args.output) if args.output else repo_root / "REPO_MAP.md"
    content = build_repo_map(repo_root, files, summarize)
    output_path.write_text(content)
    print(f"\nWrote {output_path}")


if __name__ == "__main__":
    main()
    ```

    ### Output
    ```

--- stderr ---
usage: generate_repo_map.py [-h] [--no-summarize] [--output OUTPUT] repo
generate_repo_map.py: error: the following arguments are required: repo

    ```

    ---

    ## Run: 2026-02-10T09:24:49

    **Script**: `generate_repo_map.py`
    **Interpreter**: `python`
    **Hash**: `77b98c3968a1`

    ### Summary
    This script generates a comprehensive repository map (REPO_MAP.md) for LLM navigation by scanning project files and classifying them as readable, viewable, or opaque based on file extensions. It estimates token counts for text files using a character-to-token ratio and optionally calls the Claude API to generate one-paragraph summaries of each readable file's purpose and key details. The output includes coverage statistics showing what percentage of the repository an LLM can access, a detailed file tree table with metadata, and AI-generated summaries to help LLMs understand the project structure and navigate effectively.

    ### Script
    ```
    #!/usr/bin/env python3
"""
generate_repo_map.py

Scans a project repo, classifies files by LLM-readability, estimates token counts,
calls Claude API to summarize each readable file, and writes REPO_MAP.md.

Usage:
    export ANTHROPIC_API_KEY="sk-..."
    python generate_repo_map.py /path/to/repo

Or:
    python generate_repo_map.py /path/to/repo --no-summarize   # skip API calls, structure only
"""

import os
import sys
import json
import argparse
import time
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

READABLE_EXTENSIONS = {
    # Code
    ".r", ".R", ".py", ".sh", ".jl", ".cpp", ".c", ".h", ".java", ".rs",
    # Notebooks (JSON-based, readable)
    ".ipynb",
    # Data (text-based)
    ".csv", ".tsv", ".txt", ".json", ".yaml", ".yml", ".toml", ".xml", ".bed",
    ".fasta", ".fa", ".fna", ".faa", ".fastq", ".fq", ".gff", ".gff3", ".gtf",
    ".vcf", ".sam", ".nwk", ".nex", ".phy",
    # Docs
    ".md", ".rst", ".tex", ".html", ".htm", ".bib",
    # Config
    ".cfg", ".ini", ".conf", ".env",
    # Log
    ".log",
}

BINARY_BUT_VIEWABLE = {".pdf", ".png", ".jpg", ".jpeg", ".gif", ".webp", ".svg"}

SKIP_DIRS = {".git", "__pycache__", "node_modules", ".pixi", ".snakemake", "env", "venv", ".venv"}

TOKENS_PER_CHAR = 0.25  # rough estimate

MODEL = "claude-sonnet-4-20250514"
MAX_SUMMARY_TOKENS = 300

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def estimate_tokens(text: str) -> int:
    return int(len(text) * TOKENS_PER_CHAR)


def classify_file(path: Path) -> str:
    """Returns 'readable', 'viewable', or 'opaque'."""
    ext = path.suffix.lower()
    if ext in READABLE_EXTENSIONS:
        return "readable"
    elif ext in BINARY_BUT_VIEWABLE:
        return "viewable"
    else:
        return "opaque"


def read_file_safe(path: Path, max_chars: int = 50_000) -> str | None:
    """Read a text file, return None if binary/too large."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
        if len(text) > max_chars:
            text = text[:max_chars] + f"\n\n... [truncated at {max_chars} chars]"
        return text
    except Exception:
        return None


def summarize_with_claude(filepath: str, content: str, api_key: str) -> str:
    """Call Claude API to get a one-paragraph summary of a file."""
    import urllib.request

    prompt = f"""You are summarizing a file in a computational genetics research repo for a REPO_MAP that helps LLMs navigate the project.

File: {filepath}

Respond with ONLY:
1. One sentence: what this file does / contains.
2. Key details (inputs, outputs, methods, dependencies) as a brief list.
3. If it's a script, note the rough pipeline stage (preprocessing, analysis, visualization, etc).

Keep it under 100 words total.

<file_content>
{content[:30000]}
</file_content>"""

    body = json.dumps({
        "model": MODEL,
        "max_tokens": MAX_SUMMARY_TOKENS,
        "messages": [{"role": "user", "content": prompt}],
    }).encode()

    req = urllib.request.Request(
        "https://api.anthropic.com/v1/messages",
        data=body,
        headers={
            "Content-Type": "application/json",
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
        },
    )

    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read())
            return data["content"][0]["text"].strip()
    except Exception as e:
        return f"[summary failed: {e}]"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def scan_repo(repo_root: Path) -> list[dict]:
    """Walk the repo and collect file metadata."""
    files = []
    for dirpath, dirnames, filenames in os.walk(repo_root):
        # Prune skipped dirs in-place
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fname in sorted(filenames):
            fpath = Path(dirpath) / fname
            rel = fpath.relative_to(repo_root)
            category = classify_file(fpath)
            size = fpath.stat().st_size

            entry = {
                "path": str(rel),
                "category": category,
                "size_bytes": size,
                "ext": fpath.suffix.lower(),
                "tokens_est": None,
                "summary": None,
            }

            if category == "readable":
                content = read_file_safe(fpath)
                if content is not None:
                    entry["tokens_est"] = estimate_tokens(content)
                    entry["_content"] = content  # stash for summarization

            files.append(entry)
    return files


def build_tree_string(files: list[dict]) -> str:
    """Build a markdown table of all files."""
    tag = {"readable": "✅", "viewable": "👁️", "opaque": "❌"}

    lines = [
        "| File | Directory | Type | LLM Access | Size (KB) | Est. Tokens |",
        "|------|-----------|------|------------|-----------|-------------|",
    ]

    for f in sorted(files, key=lambda x: x["path"]):
        parts = Path(f["path"]).parts
        name = parts[-1]
        d = str(Path(*parts[:-1])) if len(parts) > 1 else "."
        ext = f["ext"] or "—"
        access = tag[f["category"]]
        size_kb = f"{f['size_bytes'] / 1024:.1f}"
        tok = f"~{f['tokens_est']:,}" if f["tokens_est"] else "—"
        lines.append(f"| `{name}` | `{d}` | {ext} | {access} | {size_kb} | {tok} |")

    lines.append("")
    lines.append("*✅ = readable, 👁️ = viewable (image/PDF), ❌ = opaque binary*")

    return "\n".join(lines)


def build_repo_map(repo_root: Path, files: list[dict], summarize: bool) -> str:
    """Assemble the REPO_MAP.md content."""
    readable = [f for f in files if f["category"] == "readable"]
    viewable = [f for f in files if f["category"] == "viewable"]
    opaque = [f for f in files if f["category"] == "opaque"]

    total_tokens = sum(f["tokens_est"] or 0 for f in files)
    readable_tokens = sum(f["tokens_est"] or 0 for f in readable)

    md = []
    md.append(f"# REPO MAP — {repo_root.name}")
    md.append(f"\n*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}*\n")

    # Coverage summary
    md.append("## LLM Context Coverage\n")
    md.append(f"| Metric | Value |")
    md.append(f"|--------|-------|")
    md.append(f"| Total files | {len(files)} |")
    md.append(f"| LLM-readable (text) | {len(readable)} ({100*len(readable)/max(len(files),1):.0f}%) |")
    md.append(f"| Viewable (images/PDF) | {len(viewable)} |")
    md.append(f"| Opaque (binary) | {len(opaque)} |")
    md.append(f"| Est. readable tokens | ~{readable_tokens:,} |")
    md.append(f"| Est. total tokens | ~{total_tokens:,} |")

    # Context window estimates for common models
    context_windows = {
        "Claude Opus 4.6 / Sonnet 4.5": 200_000,
        "Claude Haiku 4.5": 200_000,
        "GPT-4o": 128_000,
        "Gemini 2.5 Pro": 1_000_000,
    }
    md.append(f"| | |")
    md.append(f"| **Model context usage** | |")
    for model_name, ctx in context_windows.items():
        pct = 100 * readable_tokens / max(ctx, 1)
        bar = "🟢" if pct < 50 else "🟡" if pct < 80 else "🔴"
        md.append(f"| {bar} {model_name} ({ctx//1000}k) | {pct:.1f}% used |")

    md.append("")

    if opaque:
        md.append("### ⚠️ Opaque files (LLM cannot read these)\n")
        for f in opaque:
            md.append(f"- `{f['path']}` ({f['ext']}, {f['size_bytes']/1024:.1f} KB)")
        md.append("")

    # Annotated tree
    md.append("## File Tree\n")
    md.append(build_tree_string(files))
    md.append("")

    # Summaries
    if summarize:
        md.append("## File Summaries\n")
        for f in sorted(readable, key=lambda x: x["path"]):
            md.append(f"### `{f['path']}`\n")
            md.append(f"{f['summary'] or '*No summary available.*'}\n")

    return "\n".join(md)


def main():
    parser = argparse.ArgumentParser(description="Generate REPO_MAP.md for LLM navigation")
    parser.add_argument("repo", help="Path to repo root")
    parser.add_argument("--no-summarize", action="store_true", help="Skip Claude API calls")
    parser.add_argument("--output", default=None, help="Output path (default: REPO_MAP.md in repo root)")
    args = parser.parse_args()

    repo_root = Path(args.repo).resolve()
    assert repo_root.is_dir(), f"Not a directory: {repo_root}"

    api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    summarize = not args.no_summarize

    if summarize:
        assert api_key, "Set ANTHROPIC_API_KEY env var or use --no-summarize"

    print(f"Scanning {repo_root} ...")
    files = scan_repo(repo_root)
    print(f"Found {len(files)} files.")

    if summarize:
        readable = [f for f in files if f["category"] == "readable" and f.get("_content")]
        print(f"Summarizing {len(readable)} readable files via Claude API...")
        for i, f in enumerate(readable):
            print(f"  [{i+1}/{len(readable)}] {f['path']}")
            f["summary"] = summarize_with_claude(f["path"], f["_content"], api_key)
            time.sleep(0.5)  # light rate limiting

    # Clean up stashed content
    for f in files:
        f.pop("_content", None)

    output_path = Path(args.output) if args.output else repo_root / "REPO_MAP.md"
    content = build_repo_map(repo_root, files, summarize)
    output_path.write_text(content)
    print(f"\nWrote {output_path}")


if __name__ == "__main__":
    main()
    ```

    ### Output
    ```

--- stderr ---
usage: generate_repo_map.py [-h] [--no-summarize] [--output OUTPUT] repo
generate_repo_map.py: error: the following arguments are required: repo

    ```

    ---

    ## Run: 2026-02-10T09:25:38

    **Script**: `generate_repo_map.py`
    **Interpreter**: `python`
    **Hash**: `77b98c3968a1`

    ### Summary
    This script scans a project repository to create an LLM-friendly navigation document called REPO_MAP.md. It classifies files into three categories (readable text files, viewable binaries like images/PDFs, and opaque binaries) based on file extensions, estimates token counts for readable files using a character-to-token ratio, and optionally calls the Claude API to generate one-paragraph summaries of each readable file's purpose and functionality. The output is a comprehensive markdown document containing file statistics, context window usage estimates for various LLM models, a structured file tree table, and AI-generated summaries to help LLMs understand the project structure and navigate the codebase effectively.

    ### Script
    ```
    #!/usr/bin/env python3
"""
generate_repo_map.py

Scans a project repo, classifies files by LLM-readability, estimates token counts,
calls Claude API to summarize each readable file, and writes REPO_MAP.md.

Usage:
    export ANTHROPIC_API_KEY="sk-..."
    python generate_repo_map.py /path/to/repo

Or:
    python generate_repo_map.py /path/to/repo --no-summarize   # skip API calls, structure only
"""

import os
import sys
import json
import argparse
import time
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

READABLE_EXTENSIONS = {
    # Code
    ".r", ".R", ".py", ".sh", ".jl", ".cpp", ".c", ".h", ".java", ".rs",
    # Notebooks (JSON-based, readable)
    ".ipynb",
    # Data (text-based)
    ".csv", ".tsv", ".txt", ".json", ".yaml", ".yml", ".toml", ".xml", ".bed",
    ".fasta", ".fa", ".fna", ".faa", ".fastq", ".fq", ".gff", ".gff3", ".gtf",
    ".vcf", ".sam", ".nwk", ".nex", ".phy",
    # Docs
    ".md", ".rst", ".tex", ".html", ".htm", ".bib",
    # Config
    ".cfg", ".ini", ".conf", ".env",
    # Log
    ".log",
}

BINARY_BUT_VIEWABLE = {".pdf", ".png", ".jpg", ".jpeg", ".gif", ".webp", ".svg"}

SKIP_DIRS = {".git", "__pycache__", "node_modules", ".pixi", ".snakemake", "env", "venv", ".venv"}

TOKENS_PER_CHAR = 0.25  # rough estimate

MODEL = "claude-sonnet-4-20250514"
MAX_SUMMARY_TOKENS = 300

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def estimate_tokens(text: str) -> int:
    return int(len(text) * TOKENS_PER_CHAR)


def classify_file(path: Path) -> str:
    """Returns 'readable', 'viewable', or 'opaque'."""
    ext = path.suffix.lower()
    if ext in READABLE_EXTENSIONS:
        return "readable"
    elif ext in BINARY_BUT_VIEWABLE:
        return "viewable"
    else:
        return "opaque"


def read_file_safe(path: Path, max_chars: int = 50_000) -> str | None:
    """Read a text file, return None if binary/too large."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
        if len(text) > max_chars:
            text = text[:max_chars] + f"\n\n... [truncated at {max_chars} chars]"
        return text
    except Exception:
        return None


def summarize_with_claude(filepath: str, content: str, api_key: str) -> str:
    """Call Claude API to get a one-paragraph summary of a file."""
    import urllib.request

    prompt = f"""You are summarizing a file in a computational genetics research repo for a REPO_MAP that helps LLMs navigate the project.

File: {filepath}

Respond with ONLY:
1. One sentence: what this file does / contains.
2. Key details (inputs, outputs, methods, dependencies) as a brief list.
3. If it's a script, note the rough pipeline stage (preprocessing, analysis, visualization, etc).

Keep it under 100 words total.

<file_content>
{content[:30000]}
</file_content>"""

    body = json.dumps({
        "model": MODEL,
        "max_tokens": MAX_SUMMARY_TOKENS,
        "messages": [{"role": "user", "content": prompt}],
    }).encode()

    req = urllib.request.Request(
        "https://api.anthropic.com/v1/messages",
        data=body,
        headers={
            "Content-Type": "application/json",
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
        },
    )

    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read())
            return data["content"][0]["text"].strip()
    except Exception as e:
        return f"[summary failed: {e}]"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def scan_repo(repo_root: Path) -> list[dict]:
    """Walk the repo and collect file metadata."""
    files = []
    for dirpath, dirnames, filenames in os.walk(repo_root):
        # Prune skipped dirs in-place
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fname in sorted(filenames):
            fpath = Path(dirpath) / fname
            rel = fpath.relative_to(repo_root)
            category = classify_file(fpath)
            size = fpath.stat().st_size

            entry = {
                "path": str(rel),
                "category": category,
                "size_bytes": size,
                "ext": fpath.suffix.lower(),
                "tokens_est": None,
                "summary": None,
            }

            if category == "readable":
                content = read_file_safe(fpath)
                if content is not None:
                    entry["tokens_est"] = estimate_tokens(content)
                    entry["_content"] = content  # stash for summarization

            files.append(entry)
    return files


def build_tree_string(files: list[dict]) -> str:
    """Build a markdown table of all files."""
    tag = {"readable": "✅", "viewable": "👁️", "opaque": "❌"}

    lines = [
        "| File | Directory | Type | LLM Access | Size (KB) | Est. Tokens |",
        "|------|-----------|------|------------|-----------|-------------|",
    ]

    for f in sorted(files, key=lambda x: x["path"]):
        parts = Path(f["path"]).parts
        name = parts[-1]
        d = str(Path(*parts[:-1])) if len(parts) > 1 else "."
        ext = f["ext"] or "—"
        access = tag[f["category"]]
        size_kb = f"{f['size_bytes'] / 1024:.1f}"
        tok = f"~{f['tokens_est']:,}" if f["tokens_est"] else "—"
        lines.append(f"| `{name}` | `{d}` | {ext} | {access} | {size_kb} | {tok} |")

    lines.append("")
    lines.append("*✅ = readable, 👁️ = viewable (image/PDF), ❌ = opaque binary*")

    return "\n".join(lines)


def build_repo_map(repo_root: Path, files: list[dict], summarize: bool) -> str:
    """Assemble the REPO_MAP.md content."""
    readable = [f for f in files if f["category"] == "readable"]
    viewable = [f for f in files if f["category"] == "viewable"]
    opaque = [f for f in files if f["category"] == "opaque"]

    total_tokens = sum(f["tokens_est"] or 0 for f in files)
    readable_tokens = sum(f["tokens_est"] or 0 for f in readable)

    md = []
    md.append(f"# REPO MAP — {repo_root.name}")
    md.append(f"\n*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}*\n")

    # Coverage summary
    md.append("## LLM Context Coverage\n")
    md.append(f"| Metric | Value |")
    md.append(f"|--------|-------|")
    md.append(f"| Total files | {len(files)} |")
    md.append(f"| LLM-readable (text) | {len(readable)} ({100*len(readable)/max(len(files),1):.0f}%) |")
    md.append(f"| Viewable (images/PDF) | {len(viewable)} |")
    md.append(f"| Opaque (binary) | {len(opaque)} |")
    md.append(f"| Est. readable tokens | ~{readable_tokens:,} |")
    md.append(f"| Est. total tokens | ~{total_tokens:,} |")

    # Context window estimates for common models
    context_windows = {
        "Claude Opus 4.6 / Sonnet 4.5": 200_000,
        "Claude Haiku 4.5": 200_000,
        "GPT-4o": 128_000,
        "Gemini 2.5 Pro": 1_000_000,
    }
    md.append(f"| | |")
    md.append(f"| **Model context usage** | |")
    for model_name, ctx in context_windows.items():
        pct = 100 * readable_tokens / max(ctx, 1)
        bar = "🟢" if pct < 50 else "🟡" if pct < 80 else "🔴"
        md.append(f"| {bar} {model_name} ({ctx//1000}k) | {pct:.1f}% used |")

    md.append("")

    if opaque:
        md.append("### ⚠️ Opaque files (LLM cannot read these)\n")
        for f in opaque:
            md.append(f"- `{f['path']}` ({f['ext']}, {f['size_bytes']/1024:.1f} KB)")
        md.append("")

    # Annotated tree
    md.append("## File Tree\n")
    md.append(build_tree_string(files))
    md.append("")

    # Summaries
    if summarize:
        md.append("## File Summaries\n")
        for f in sorted(readable, key=lambda x: x["path"]):
            md.append(f"### `{f['path']}`\n")
            md.append(f"{f['summary'] or '*No summary available.*'}\n")

    return "\n".join(md)


def main():
    parser = argparse.ArgumentParser(description="Generate REPO_MAP.md for LLM navigation")
    parser.add_argument("repo", help="Path to repo root")
    parser.add_argument("--no-summarize", action="store_true", help="Skip Claude API calls")
    parser.add_argument("--output", default=None, help="Output path (default: REPO_MAP.md in repo root)")
    args = parser.parse_args()

    repo_root = Path(args.repo).resolve()
    assert repo_root.is_dir(), f"Not a directory: {repo_root}"

    api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    summarize = not args.no_summarize

    if summarize:
        assert api_key, "Set ANTHROPIC_API_KEY env var or use --no-summarize"

    print(f"Scanning {repo_root} ...")
    files = scan_repo(repo_root)
    print(f"Found {len(files)} files.")

    if summarize:
        readable = [f for f in files if f["category"] == "readable" and f.get("_content")]
        print(f"Summarizing {len(readable)} readable files via Claude API...")
        for i, f in enumerate(readable):
            print(f"  [{i+1}/{len(readable)}] {f['path']}")
            f["summary"] = summarize_with_claude(f["path"], f["_content"], api_key)
            time.sleep(0.5)  # light rate limiting

    # Clean up stashed content
    for f in files:
        f.pop("_content", None)

    output_path = Path(args.output) if args.output else repo_root / "REPO_MAP.md"
    content = build_repo_map(repo_root, files, summarize)
    output_path.write_text(content)
    print(f"\nWrote {output_path}")


if __name__ == "__main__":
    main()
    ```

    ### Output
    ```

--- stderr ---
usage: generate_repo_map.py [-h] [--no-summarize] [--output OUTPUT] repo
generate_repo_map.py: error: the following arguments are required: repo

    ```

    ---

