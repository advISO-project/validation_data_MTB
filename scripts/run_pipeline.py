#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path("scripts")

def run_step(name, script):
    print(f"\n🔹 Running step: {name}")
    try:
        subprocess.run(
            ["python", str(SCRIPTS / script)], 
            check=True
        )
        print(f"✅ {name} completed successfully")
    except subprocess.CalledProcessError:
        print(f"❌ {name} failed. Stopping pipeline.")
        sys.exit(1)

def main():
    print("\n🚀 Starting mini validation workflow...\n")

    run_step("Extract truth set", "extract_samples_from_tbtamr.py")
    run_step("Download / collect FASTQs", "extract_fastqs_for_tbtamr.py")
    run_step("Generate manifests (CSV + YAML)", "generate_manifest.py")

    print("\n🎉 All steps completed successfully!\n")
    print("Outputs are available in tbtamr/")

if __name__ == "__main__":
    main()
