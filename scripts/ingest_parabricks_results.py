"""
After running run_parabricks_benchmark.sh on the SageMaker host,
copy the timing file back and run this to update benchmark_summary.json.

Usage:
  # Copy results from SageMaker host (run locally):
  scp -P 17022 ec2-user@localhost:/home/ec2-user/SageMaker/parabricks_bench/logs/timing.txt \
      results/sagemaker/timing.txt
  scp -P 17022 ec2-user@localhost:/home/ec2-user/SageMaker/parabricks_bench/logs/benchmark.log \
      results/sagemaker/benchmark.log

  python scripts/ingest_parabricks_results.py
"""
import json
import sys
from pathlib import Path

TIMING_FILE = Path("results/sagemaker/timing.txt")
SUMMARY_FILE = Path("results/benchmark_summary.json")

if not TIMING_FILE.exists():
    print(f"ERROR: {TIMING_FILE} not found. Copy it from SageMaker first.")
    sys.exit(1)

timing = {}
for line in TIMING_FILE.read_text().splitlines():
    if "=" in line:
        k, v = line.strip().split("=", 1)
        try:
            timing[k.strip()] = float(v.strip())
        except ValueError:
            timing[k.strip()] = v.strip()

print("Timing data:", timing)

cpu_s = timing.get("cpu_mutect2_sec")
gpu_s = timing.get("gpu_mutect2_sec")
fq_s  = timing.get("gpu_fq2bam_sec")

speedup_mutect2 = round(cpu_s / gpu_s, 1) if cpu_s and gpu_s else None

summary = {
    "gpu_model": "NVIDIA A10G (23GB VRAM)",
    "parabricks_version": "4.7.1-1",
    "dataset": "HCC1143 tumor-normal (Broad GATK tutorial, chr17 subset)",
    "reference": "GRCh38 chr17",
    "cpu_mutect2_sec": cpu_s,
    "gpu_mutect2_sec": gpu_s,
    "speedup_mutect2": speedup_mutect2,
    "gpu_fq2bam_sec": fq_s,
    "cpu_variants": int(timing.get("cpu_variants", 0)),
    "gpu_variants": int(timing.get("gpu_variants", 0)),
    "notes": "Real tumor-normal matched BAMs. Chr17 subset for feasibility.",
}

SUMMARY_FILE.parent.mkdir(parents=True, exist_ok=True)
existing = {}
if SUMMARY_FILE.exists():
    existing = json.loads(SUMMARY_FILE.read_text())

existing.update({"parabricks_real_run": summary})
SUMMARY_FILE.write_text(json.dumps(existing, indent=2))
print(f"Updated {SUMMARY_FILE}")

if speedup_mutect2:
    print(f"\n=== Results ===")
    print(f"CPU Mutect2 (chr17):  {cpu_s}s")
    print(f"GPU Parabricks:       {gpu_s}s")
    print(f"Speedup (Mutect2):    {speedup_mutect2}x")
    if fq_s:
        print(f"GPU fq2bam:           {fq_s}s")
