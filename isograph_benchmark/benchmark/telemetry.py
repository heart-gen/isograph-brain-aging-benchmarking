from __future__ import annotations

import importlib.metadata
import json
import os
import platform
import resource
import shutil
import socket
import subprocess
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator


def _read_first_cpu_model() -> str | None:
    cpuinfo = Path("/proc/cpuinfo")
    if not cpuinfo.exists():
        return None
    for line in cpuinfo.read_text(errors="ignore").splitlines():
        if line.startswith("model name"):
            return line.split(":", 1)[1].strip()
    return None


def _mem_total_mb() -> int | None:
    meminfo = Path("/proc/meminfo")
    if not meminfo.exists():
        return None
    for line in meminfo.read_text(errors="ignore").splitlines():
        if line.startswith("MemTotal:"):
            return int(line.split()[1]) // 1024
    return None


def _version(package: str) -> str | None:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return None


def _run_short(cmd: list[str], timeout: int = 5) -> str | None:
    if shutil.which(cmd[0]) is None:
        return None
    try:
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=timeout)
    except Exception:
        return None
    text = (proc.stdout or proc.stderr).strip()
    return text or None


def software_versions(include_r: bool = False) -> dict[str, Any]:
    versions: dict[str, Any] = {
        "python": sys.version.split()[0],
        "isograph": _version("isograph"),
        "numpy": _version("numpy"),
        "pandas": _version("pandas"),
        "pyarrow": _version("pyarrow"),
        "torch": _version("torch"),
    }
    if include_r:
        versions["r"] = _run_short(["Rscript", "--version"])
        versions["wgcna"] = _run_short(
            ["Rscript", "-e", "cat(as.character(packageVersion('WGCNA')))"],
            timeout=10,
        )
    return versions


def hardware_info() -> dict[str, Any]:
    info: dict[str, Any] = {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "cpu_model": _read_first_cpu_model(),
        "logical_cpus": os.cpu_count(),
        "total_ram_mb": _mem_total_mb(),
    }
    if shutil.which("nvidia-smi") is not None:
        query = "name,memory.total,driver_version"
        gpu = _run_short(
            ["nvidia-smi", f"--query-gpu={query}", "--format=csv,noheader,nounits"],
            timeout=5,
        )
        info["nvidia_smi_gpu_query"] = gpu
    return info


def slurm_info() -> dict[str, Any]:
    keys = [
        "SLURM_JOB_ID",
        "SLURM_ARRAY_JOB_ID",
        "SLURM_ARRAY_TASK_ID",
        "SLURM_JOB_PARTITION",
        "SLURM_JOB_NODELIST",
        "SLURM_CPUS_PER_TASK",
        "SLURM_MEM_PER_NODE",
        "SLURM_MEM_PER_CPU",
        "SLURM_GPUS",
        "SLURM_JOB_GPUS",
        "SLURM_TIMELIMIT",
    ]
    return {key.lower(): os.environ.get(key) for key in keys}


def torch_info() -> dict[str, Any]:
    try:
        import torch
    except Exception:
        return {"torch_available": False}
    info: dict[str, Any] = {
        "torch_available": True,
        "cuda_available": bool(torch.cuda.is_available()),
        "torch_cuda_version": getattr(torch.version, "cuda", None),
    }
    if torch.cuda.is_available():
        info["cuda_device_count"] = torch.cuda.device_count()
        info["cuda_device_name"] = torch.cuda.get_device_name(0)
        info["cuda_total_memory_mb"] = int(torch.cuda.get_device_properties(0).total_memory // 1024**2)
    return info


def reset_torch_peak_memory() -> None:
    try:
        import torch
    except Exception:
        return
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()


def torch_peak_memory() -> dict[str, Any]:
    try:
        import torch
    except Exception:
        return {}
    if not torch.cuda.is_available():
        return {}
    return {
        "gpu_peak_allocated_mb": int(torch.cuda.max_memory_allocated() // 1024**2),
        "gpu_peak_reserved_mb": int(torch.cuda.max_memory_reserved() // 1024**2),
    }


def max_rss_mb() -> float:
    own = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    scale = 1024**2 if sys.platform == "darwin" else 1024
    return float(max(own, child) / scale)


@contextmanager
def measured_run() -> Iterator[dict[str, Any]]:
    start = time.time()
    start_perf = time.perf_counter()
    start_usage = resource.getrusage(resource.RUSAGE_SELF)
    payload: dict[str, Any] = {
        "start_time_unix": start,
        "start_time_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(start)),
    }
    try:
        yield payload
    finally:
        end = time.time()
        end_usage = resource.getrusage(resource.RUSAGE_SELF)
        payload.update(
            {
                "end_time_unix": end,
                "end_time_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(end)),
                "elapsed_sec": time.perf_counter() - start_perf,
                "user_sec": end_usage.ru_utime - start_usage.ru_utime,
                "sys_sec": end_usage.ru_stime - start_usage.ru_stime,
                "max_rss_mb": max_rss_mb(),
            }
        )


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
