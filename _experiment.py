"""
Experiment runner: wraps `pixi run <interpreter> <script>`, logs everything,
and calls the Claude API for a summary.

Usage:
    pixi run experiment <interpreter> <script> [--force]

Examples:
    pixi run experiment python analysis.py
    pixi run experiment Rscript model.R --force
"""

import sys
import os
import hashlib
import json
import subprocess
import pathlib
import datetime
import textwrap

# ── arg parsing ──────────────────────────────────────────────────────────────

args = sys.argv[1:]
force = "--force" in args
if force:
    args.remove("--force")

assert len(args) >= 2, (
    "Usage: pixi run experiment <interpreter> <script> [--force]\n"
    "  e.g. pixi run experiment python my_analysis.py"
)

interpreter = args[0]
script_path = pathlib.Path(args[1])
assert script_path.exists(), f"Script not found: {script_path}"

script_text = script_path.read_text()
script_hash = hashlib.sha256(script_text.encode()).hexdigest()

# ── duplicate detection ──────────────────────────────────────────────────────

HASH_FILE = pathlib.Path(".experiment_hashes")
hashes = {}
if HASH_FILE.exists():
    hashes = json.loads(HASH_FILE.read_text())

key = str(script_path)
if key in hashes and hashes[key] == script_hash and not force:
    print(
        f"ERROR: {script_path} has already been run with identical content.\n"
        f"Use --force to run anyway.",
        file=sys.stderr,
    )
    sys.exit(1)

# ── run the script ───────────────────────────────────────────────────────────

print(f"▶ Running: {interpreter} {script_path}")

# Force line-buffered/unbuffered output from child processes
env = os.environ.copy()
env["PYTHONUNBUFFERED"] = "1"  # for python scripts

# Use stdbuf to force line-buffering if available (helps R, etc.)
import shutil
cmd = [interpreter, str(script_path)]
if shutil.which("stdbuf") and interpreter != "python":
    cmd = ["stdbuf", "-oL", "-eL"] + cmd

proc = subprocess.Popen(
    cmd,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,  # merge stderr into stdout so we see everything live
    text=True,
    env=env,
)

output_lines = []
for line in proc.stdout:
    print(line, end="", flush=True)
    output_lines.append(line)
proc.wait()

output = "".join(output_lines)

if proc.returncode != 0:
    print(f"⚠ Script exited with code {proc.returncode}", file=sys.stderr)

# ── call Claude API for summary ─────────────────────────────────────────────

api_key = os.environ.get("ANTHROPIC_API_KEY")
assert api_key, (
    "ANTHROPIC_API_KEY not set. Set it in your environment or .env file."
)

import urllib.request

summary_prompt = (
    f"Briefly describe what this script does in 2-4 sentences. "
    f"Be specific about the methods/data involved. "
    f"Do not include code.\n\n```\n{script_text}\n```"
)

request_body = json.dumps({
    "model": "claude-sonnet-4-20250514",
    "max_tokens": 300,
    "messages": [{"role": "user", "content": summary_prompt}],
}).encode()

req = urllib.request.Request(
    "https://api.anthropic.com/v1/messages",
    data=request_body,
    headers={
        "Content-Type": "application/json",
        "x-api-key": api_key,
        "anthropic-version": "2023-06-01",
    },
    method="POST",
)

try:
    with urllib.request.urlopen(req) as resp:
        data = json.loads(resp.read())
    summary = data["content"][0]["text"]
except urllib.error.HTTPError as e:
    error_body = e.read().decode()
    summary = f"(Claude API call failed: {e} - {error_body})"
    print(f"⚠ {summary}", file=sys.stderr)
except Exception as e:
    summary = f"(Claude API call failed: {e})"
    print(f"⚠ {summary}", file=sys.stderr)

# ── write notes ──────────────────────────────────────────────────────────────

notes_dir = pathlib.Path("notes_on")
notes_dir.mkdir(exist_ok=True)

note_name = f"on_{script_path.stem}.md"
note_path = notes_dir / note_name

timestamp = datetime.datetime.now().isoformat(timespec="seconds")

entry = textwrap.dedent(f"""\
    ## Run: {timestamp}

    **Script**: `{script_path}`
    **Interpreter**: `{interpreter}`
    **Hash**: `{script_hash[:12]}`

    ### Summary
    {summary}

    ### Script
    ```
    {script_text}
    ```

    ### Output
    ```
    {output}
    ```

    ---

    """)

with open(note_path, "a") as f:
    f.write(entry)

print(f"✓ Notes appended to {note_path}")

# ── save hash ────────────────────────────────────────────────────────────────

hashes[key] = script_hash
HASH_FILE.write_text(json.dumps(hashes, indent=2))