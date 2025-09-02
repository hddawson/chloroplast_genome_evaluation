#!/usr/bin/env python3
"""
Tmux Pane GitHub Reporter

Reports tmux pane output to GitHub via:
1. GitHub Issues (one issue per job, comments for updates)
2. GitHub Gists (one gist per job, updated with new content)
3. Repository commits (commits new logs to a repo)

Usage:
    python tmux_watcher.py --add session_name:pane_id "Job Description" --method issues
    python tmux_watcher.py --list
    python tmux_watcher.py --daemon
"""

import subprocess
import json
import time
import requests
import argparse
import os
import hashlib
from datetime import datetime
from pathlib import Path
import base64

CONFIG_FILE = os.path.expanduser("~/.tmux_github_watcher.json")

class GitHubReporter:
    def __init__(self):
        self.config = self.load_config()
        self.running = True
        
    def load_config(self):
        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as f:
                return json.load(f)
        return {
            "watched_panes": {},
            "github_token": "",
            "github_repo": "",  # format: "username/repo-name"
            "default_method": "gists"  # issues, gists, or commits
        }
    
    def save_config(self):
        with open(CONFIG_FILE, 'w') as f:
            json.dump(self.config, f, indent=2)
    
    def get_headers(self):
        return {
            "Authorization": f"token {self.config['github_token']}",
            "Accept": "application/vnd.github.v3+json"
        }
    
    def get_pane_content(self, session, pane_id, lines=50):
        """Get recent content from a tmux pane"""
        try:
            result = subprocess.run([
                'tmux', 'capture-pane', '-t', f"{session}:{pane_id}", 
                '-p', '-S', f'-{lines}'
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                return result.stdout.strip()
            else:
                return None
        except Exception as e:
            print(f"Error capturing pane {session}:{pane_id}: {e}")
            return None
    
    def find_new_lines(self, old_content, new_content):
        """Find new lines added since last check"""
        if not old_content:
            return new_content
        
        old_lines = old_content.split('\n')
        new_lines = new_content.split('\n')
        
        if len(new_lines) <= len(old_lines):
            return ""
        
        # Find where old content ends in new content
        old_tail = old_lines[-3:] if len(old_lines) >= 3 else old_lines
        
        for i in range(len(new_lines) - len(old_tail) + 1):
            if new_lines[i:i+len(old_tail)] == old_tail:
                return '\n'.join(new_lines[i+len(old_tail):])
        
        return '\n'.join(new_lines[-20:])  # Last 20 lines as fallback
    
    def create_github_issue(self, job_name, session, pane_id, initial_content):
        """Create a new GitHub issue for job tracking"""
        url = f"https://api.github.com/repos/{self.config['github_repo']}/issues"
        
        body = f"""# Job Monitor: {job_name}

**Tmux Pane:** `{session}:{pane_id}`  
**Started:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  
**Host:** `{os.uname().nodename}`

## Initial Output
```
{initial_content[-1500:]}  
```

---
*This issue will be updated automatically with new output from the tmux pane.*
"""
        
        payload = {
            "title": f"🔬 {job_name}",
            "body": body,
            "labels": ["tmux-monitor", "automated"]
        }
        
        try:
            response = requests.post(url, json=payload, headers=self.get_headers())
            if response.status_code == 201:
                issue_data = response.json()
                print(f"✅ Created GitHub issue #{issue_data['number']} for {job_name}")
                return issue_data['number']
            else:
                print(f"❌ Failed to create issue: {response.status_code}")
                return None
        except Exception as e:
            print(f"❌ Error creating issue: {e}")
            return None
    
    def comment_on_issue(self, issue_number, new_content, timestamp):
        """Add a comment to the GitHub issue with new output"""
        url = f"https://api.github.com/repos/{self.config['github_repo']}/issues/{issue_number}/comments"
        
        body = f"""**Update:** {timestamp}

```
{new_content}
```"""
        
        payload = {"body": body}
        
        try:
            response = requests.post(url, json=payload, headers=self.get_headers())
            if response.status_code == 201:
                print(f"✅ Added comment to issue #{issue_number}")
                return True
            else:
                print(f"❌ Failed to comment: {response.status_code}")
                return False
        except Exception as e:
            print(f"❌ Error commenting: {e}")
            return False
    
    def create_gist(self, job_name, session, pane_id, content):
        """Create a GitHub gist for the job"""
        url = "https://api.github.com/gists"
        
        filename = f"{job_name.replace(' ', '_')}_{session}_{pane_id}.log"
        
        payload = {
            "description": f"Tmux Monitor: {job_name} ({session}:{pane_id})",
            "public": False,
            "files": {
                filename: {
                    "content": f"# {job_name}\n# Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n# Pane: {session}:{pane_id}\n\n{content}"
                }
            }
        }
        
        try:
            response = requests.post(url, json=payload, headers=self.get_headers())
            if response.status_code == 201:
                gist_data = response.json()
                print(f"✅ Created gist for {job_name}: {gist_data['html_url']}")
                return gist_data['id']
            else:
                print(f"❌ Failed to create gist: {response.status_code}")
                return None
        except Exception as e:
            print(f"❌ Error creating gist: {e}")
            return None
    
    def update_gist(self, gist_id, job_name, session, pane_id, full_content):
        """Update the GitHub gist with new content"""
        url = f"https://api.github.com/gists/{gist_id}"
        
        filename = f"{job_name.replace(' ', '_')}_{session}_{pane_id}.log"
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        
        payload = {
            "files": {
                filename: {
                    "content": f"# {job_name}\n# Last Updated: {timestamp}\n# Pane: {session}:{pane_id}\n\n{full_content}"
                }
            }
        }
        
        try:
            response = requests.patch(url, json=payload, headers=self.get_headers())
            if response.status_code == 200:
                print(f"✅ Updated gist for {job_name}")
                return True
            else:
                print(f"❌ Failed to update gist: {response.status_code}")
                return False
        except Exception as e:
            print(f"❌ Error updating gist: {e}")
            return False
    
    def commit_to_repo(self, job_name, session, pane_id, new_content):
        """Commit new output to a file in the repository"""
        filename = f"logs/{job_name.replace(' ', '_')}_{session}_{pane_id}.log"
        url = f"https://api.github.com/repos/{self.config['github_repo']}/contents/{filename}"
        
        # Get current file to append to it
        try:
            response = requests.get(url, headers=self.get_headers())
            if response.status_code == 200:
                file_data = response.json()
                current_content = base64.b64decode(file_data['content']).decode('utf-8')
                sha = file_data['sha']
            else:
                current_content = f"# {job_name}\n# Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n# Pane: {session}:{pane_id}\n\n"
                sha = None
        except:
            current_content = f"# {job_name}\n# Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n# Pane: {session}:{pane_id}\n\n"
            sha = None
        
        # Append new content with timestamp
        updated_content = current_content + f"\n\n--- {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ---\n{new_content}"
        
        payload = {
            "message": f"Update {job_name} logs",
            "content": base64.b64encode(updated_content.encode()).decode(),
        }
        
        if sha:
            payload["sha"] = sha
        
        try:
            response = requests.put(url, json=payload, headers=self.get_headers())
            if response.status_code in [200, 201]:
                print(f"✅ Committed update for {job_name}")
                return True
            else:
                print(f"❌ Failed to commit: {response.status_code}")
                return False
        except Exception as e:
            print(f"❌ Error committing: {e}")
            return False
    
    def send_update(self, pane_key, pane_info, new_content, full_content):
        """Send update using the configured method"""
        job_name = pane_info["name"]
        session, pane_id = pane_key.split(':', 1)
        method = pane_info.get("method", self.config["default_method"])
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        
        if method == "issues":
            issue_number = pane_info.get("issue_number")
            if issue_number:
                self.comment_on_issue(issue_number, new_content, timestamp)
        
        elif method == "gists":
            gist_id = pane_info.get("gist_id")
            if gist_id:
                self.update_gist(gist_id, job_name, session, pane_id, full_content)
        
        elif method == "commits":
            self.commit_to_repo(job_name, session, pane_id, new_content)
    
    def add_pane(self, session_pane, job_name, method=None, mode="periodic"):
        """Add a pane to watch list"""
        if ':' not in session_pane:
            print("❌ Format should be session_name:pane_id")
            return
        
        if not self.config["github_token"] or not self.config["github_repo"]:
            print("❌ Please configure GitHub token and repo first")
            return
        
        session, pane_id = session_pane.split(':', 1)
        method = method or self.config["default_method"]
        
        # Verify pane exists
        content = self.get_pane_content(session, pane_id, lines=50)
        if content is None:
            print(f"❌ Pane {session_pane} not found")
            return
        
        pane_info = {
            "name": job_name,
            "method": method,
            "added": datetime.now().isoformat(),
            "update_mode": mode,  # "changes" or "periodic"
            "last_content": content,
            "last_update": datetime.now().isoformat()
        }
        
        # Initialize the GitHub resource
        if method == "issues":
            issue_number = self.create_github_issue(job_name, session, pane_id, content)
            if issue_number:
                pane_info["issue_number"] = issue_number
            else:
                return
        
        elif method == "gists":
            gist_id = self.create_gist(job_name, session, pane_id, content)
            if gist_id:
                pane_info["gist_id"] = gist_id
            else:
                return
        
        elif method == "commits":
            # Initial commit
            if not self.commit_to_repo(job_name, session, pane_id, f"Initial content:\n{content}"):
                return
        
        self.config["watched_panes"][session_pane] = pane_info
        self.save_config()
        print(f"✅ Now watching {session_pane} ({job_name}) via {method}")
    
    def watch_panes(self):
        """Main monitoring loop"""
        print("🔍 Starting tmux pane watcher for GitHub...")
        
        while self.running:
            current_time = time.time()

            for pane_key, pane_info in self.config["watched_panes"].items():
                session, pane_id = pane_key.split(':', 1)
                job_name = pane_info["name"]
                last_content = pane_info.get("last_content", "")
                last_update_time = pane_info.get("last_update_time", 0)
                update_mode = pane_info.get("update_mode", "periodic")  # "changes" or "periodic"
                # Get current content
                current_content = self.get_pane_content(session, pane_id, lines=200)
                
                if current_content is None:
                    print(f"⚠️ Pane {pane_key} not found, skipping...")
                    continue

                should_update = False
                update_content = ""
                
                if update_mode == "periodic":
                    # Update every 10 minutes regardless of changes
                    if current_time - last_update_time > 900:  # 15 minutes
                        should_update = True
                        # Get last 25 lines for periodic updates
                        recent_lines = '\n'.join(current_content.split('\n')[-10:])
                        update_content = recent_lines
                        print(f"⏰ Periodic update for {job_name} (every 10 min)")
                else:
                    # Original change-detection mode
                    new_lines = self.find_new_lines(last_content, current_content)
                    if new_lines.strip():
                        should_update = True
                        update_content = new_lines
                        print(f"📝 New output detected in {job_name}")
                
                if should_update:
                    self.send_update(pane_key, pane_info, update_content, current_content)
                    
                    # Update stored content and timestamp
                    self.config["watched_panes"][pane_key]["last_content"] = current_content
                    self.config["watched_panes"][pane_key]["last_update"] = datetime.now().isoformat()
                    self.config["watched_panes"][pane_key]["last_update_time"] = current_time
            
            self.save_config()
            time.sleep(15)  # Check every 15 seconds
    
    def list_panes(self):
        """List all watched panes"""
        if not self.config["watched_panes"]:
            print("No panes being watched")
            return
        
        print("\n🔍 Watched Panes:")
        print("-" * 80)
        for pane_key, info in self.config["watched_panes"].items():
            method = info.get("method", "unknown")
            last_update = info.get("last_update", "Never")
            if last_update != "Never":
                last_update = datetime.fromisoformat(last_update).strftime("%m/%d %H:%M")
            
            # Show GitHub resource info
            resource_info = ""
            if method == "issues" and "issue_number" in info:
                resource_info = f"Issue #{info['issue_number']}"
            elif method == "gists" and "gist_id" in info:
                resource_info = f"Gist {info['gist_id'][:8]}..."
            elif method == "commits":
                resource_info = "Repository logs/"
            
            print(f"{pane_key:20} | {info['name']:20} | {method:8} | {resource_info:15} | {last_update}")

def main():
    parser = argparse.ArgumentParser(description='Watch tmux panes and report to GitHub')
    parser.add_argument('--add', nargs=2, metavar=('PANE', 'NAME'), 
                       help='Add pane to watch (format: session:pane_id "Job Name")')
    parser.add_argument('--method', choices=['issues', 'gists', 'commits'], 
                       help='GitHub reporting method')
    parser.add_argument('--mode', choices=['changes', 'periodic'], default='changes',
                       help='Update mode: "changes" (on new output) or "periodic" (every 10 min)')
    parser.add_argument('--list', action='store_true', help='List watched panes')
    parser.add_argument('--daemon', action='store_true', help='Start monitoring daemon')
    parser.add_argument('--token', metavar='TOKEN', help='Set GitHub personal access token')
    parser.add_argument('--repo', metavar='USER/REPO', help='Set GitHub repository')
    
    args = parser.parse_args()
    watcher = GitHubReporter()
    
    if args.token:
        watcher.config["github_token"] = args.token
        watcher.save_config()
        print(f"✅ GitHub token updated")
        return
        
    if args.repo:
        watcher.config["github_repo"] = args.repo
        watcher.save_config()
        print(f"✅ GitHub repo set to {args.repo}")
        return

    if args.add:
        watcher.add_pane(args.add[0], args.add[1], args.method, args.mode)
    elif args.list:
        watcher.list_panes()
    elif args.daemon:
        try:
            watcher.watch_panes()
        except KeyboardInterrupt:
            print("\n🛑 Stopping watcher...")
    else:
        parser.print_help()
        print("\n💡 Setup:")
        print("1. python tmux_watcher.py --token 'your_github_token'")
        print("2. python tmux_watcher.py --repo 'username/repo-name'")
        print("3. python tmux_watcher.py --add 'session:pane' 'Job Name' --method gists --mode periodic")
        print("   (--mode periodic = update every 10 min regardless of changes)")
        print("   (--mode changes = only update when new output detected)")

if __name__ == "__main__":
    main()