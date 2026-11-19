"""Run the computational experiments of the paper (restartable).

    python scripts/run_experiments.py aggregation   # Sections 6.1 and 6.2
    python scripts/run_experiments.py scaling       # Section 6.3
    python scripts/run_experiments.py all

aggregation: dominance audits at 50 locations (6 cities) and the aggregation study with
             50/75/100/200/300/400 locations, 20 potential stations and a budget of 10.
             Results: <study>/raw/aggregation_cases, <study>/runtime/aggregation.
scaling:     300 locations; the budget grid (20/80/140 potential stations x 25/50/75 % budgets)
             and the station-count scaling study (25 % budget).
             Results: <study>/raw/scaling_cases, <study>/runtime/scaling.

Every case has a total time limit of one hour (enumeration, model construction and
optimization) and a memory guard of 16 GiB. Stop rules: in the aggregation study, larger
sizes of a city are skipped after its first failure; in the scaling study, larger station
counts are skipped for the same city and budget fraction. Completed cases (existing
result.json / summary.json) are reused when the script is started again.

A long-running Julia worker (scripts/worker.jl) is started with the JuBiC project
environment; it runs the station-path unit tests on start-up and refuses to start if one fails.
"""
from pathlib import Path
import json, os, subprocess, sys, time
from paths import B, CODE, REPO, JULIA, SETTINGS as C

if len(sys.argv) != 2 or sys.argv[1] not in ("aggregation", "scaling", "all"):
    sys.exit(__doc__)
GROUPS = ["aggregation", "scaling"] if sys.argv[1] == "all" else [sys.argv[1]]
THREADS = int(C["threads"])
MEMORY_LIMIT = C["memory_guard_gib"] * 2**30
TIME_LIMIT = C["total_case_seconds"]
GROUP_LABEL = {"aggregation": "aggregation", "scaling": "candidate_budget"}  # labels read by the analysis


def read(p):
    try:
        return json.loads(p.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError):
        return {}


def write(p, x):
    p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_suffix(".tmp")
    temp.write_text(json.dumps(x, indent=2), encoding="utf-8")
    # A reader can briefly deny replacement on Windows; retry the atomic handoff.
    for attempt in range(300):
        try:
            os.replace(temp, p)
            return
        except PermissionError:
            if attempt == 299:
                raise
            time.sleep(.1)


def memory(proc):
    """Private memory of the worker process in bytes (0 if it cannot be measured)."""
    if os.name == "nt":
        import ctypes, ctypes.wintypes
        class M(ctypes.Structure):
            _fields_ = [("cb", ctypes.wintypes.DWORD), ("faults", ctypes.wintypes.DWORD)] + [
                (n, ctypes.c_size_t) for n in ["peak_working_set", "working_set", "peak_paged_pool", "paged_pool",
                                               "peak_nonpaged_pool", "nonpaged_pool", "pagefile", "peak_pagefile",
                                               "private_bytes"]]
        m = M(); m.cb = ctypes.sizeof(m)
        fn = ctypes.windll.psapi.GetProcessMemoryInfo
        fn.argtypes = [ctypes.wintypes.HANDLE, ctypes.c_void_p, ctypes.wintypes.DWORD]
        return m.private_bytes if fn(int(proc._handle), ctypes.byref(m), m.cb) else 0
    try:
        import psutil
        return psutil.Process(proc.pid).memory_info().rss
    except Exception:
        return 0


def case(jobs, group, city, k, s, b):
    name = f"{city}_{k}_s{s}_b{b}"
    if any(x["id"] == name for x in jobs):
        return
    jobs.append(dict(id=name, kind="solve", group=group, groups=[GROUP_LABEL[group]], city=city, size=k, dominance=True,
                     candidate_count=s, station_budget=b,
                     directory=str(B / f"raw/{group}_cases" / name),
                     aggregation_root=str(B / "inputs/aggregation"), demand_root=str(B / "inputs/demand")))


def job_list(group):
    jobs = []
    if group == "aggregation":
        s, b = C["aggregation_candidates"], C["aggregation_budget"]
        for city in C["cities"]:
            case(jobs, group, city, C["targets"][0], s, b)
        for city in C["cities"]:
            k = C["targets"][0]
            jobs.append(dict(id=f"audit_{city}_{k}", kind="audit", group=group, groups=["dominance"], city=city, size=k,
                             directory=str(B / "runtime/aggregation" / f"audit_{city}_{k}")))
        for k in C["targets"][1:]:
            for city in C["cities"]:
                case(jobs, group, city, k, s, b)
    else:
        k = C["budget_study_locations"]
        pairs = {(s, f) for s in C["candidate_counts"] for f in C["budget_fractions"]}
        pairs |= {(s, C["scaling_budget_fraction"]) for s in C["scaling_candidate_counts"]}
        for s, f in sorted(pairs):
            for city in C["cities"]:
                case(jobs, group, city, k, s, round(s * f))
    return jobs


def run(group):
    root = B / "runtime" / group
    root.mkdir(parents=True, exist_ok=True)
    jobs = job_list(group)
    write(B / f"jobs_{group}.json", jobs)
    progress = B / f"progress_{group}.json"
    env = os.environ.copy()
    env.update(SLP_STUDY_DIR=str(B), SLP_RUNTIME=str(root), JUBIC_REPO=str(REPO))
    env["TEMP"] = env["TMP"] = str(root / "temp"); (root / "temp").mkdir(exist_ok=True)
    state = dict(proc=None, log=None, generation=0)

    def start_worker():
        state["generation"] += 1
        for name in ["command.json", "worker_ready.json", "done.json"]:
            (root / name).unlink(missing_ok=True)
        state["log"] = (root / f"worker_{state['generation']}.log").open("w", encoding="utf-8")
        cmd = [JULIA, f"--project={REPO}", "--compiled-modules=no", "--pkgimages=no", "-O2",
               f"--threads={THREADS}", str(CODE / "scripts" / "worker.jl")]
        flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        state["proc"] = proc = subprocess.Popen(cmd, cwd=root, env=env, stdout=state["log"],
                                                stderr=subprocess.STDOUT, creationflags=flags)
        write(root / "process.json", dict(pid=proc.pid, generation=state["generation"], command=cmd,
                                          started_unix=time.time()))
        start = time.time()
        while not read(root / "worker_ready.json").get("ready"):
            if proc.poll() is not None:
                raise RuntimeError("Worker startup failed; inspect " + state["log"].name)
            if time.time() - start > 2400:
                raise RuntimeError("Worker startup exceeded 40 minutes")
            write(root / "status.json", dict(phase="worker_startup", pid=proc.pid, seconds=time.time() - start))
            time.sleep(2)

    def stop_worker():
        proc = state["proc"]
        if proc and proc.poll() is None:
            proc.kill(); proc.wait(timeout=30)
        if state["log"]:
            state["log"].close()
        state["proc"] = None

    blocked_city, blocked_fraction, completed, skipped = set(), set(), [], []
    try:
        for job in jobs:
            directory = Path(job["directory"]); audit = job["kind"] == "audit"
            resultfile = directory / ("summary.json" if audit else "result.json")
            existing = read(resultfile)
            if existing:
                completed.append(dict(id=job["id"], result=existing)); status = existing.get("status", "complete")
            else:
                reason = None
                if not audit and group == "aggregation" and job["city"] in blocked_city:
                    reason = "larger aggregation skipped after earlier timeout/failure"
                if not audit and group == "scaling" and \
                        (job["city"], job["station_budget"] / job["candidate_count"]) in blocked_fraction:
                    reason = "larger candidate pool skipped after earlier timeout/failure at this budget fraction"
                if reason:
                    skipped.append(dict(id=job["id"], reason=reason))
                    write(progress, dict(phase="running", completed=len(completed), skipped=skipped))
                    continue
                if not (B / "inputs/aggregation" / job["city"] / "COMPLETE.json").exists():
                    raise RuntimeError(f"Missing prepared inputs for {job['city']}; run prepare_instances.jl first")
                if state["proc"] is None:
                    start_worker()
                proc = state["proc"]
                directory.mkdir(parents=True, exist_ok=True)
                write(root / "status.json", dict(phase="running", case=job["id"], completed=len(completed), skipped=len(skipped)))
                write(root / "command.json", job); dispatch = time.time(); peak = 0; failure = None
                while read(root / "done.json").get("id") != job["id"]:
                    if proc.poll() is not None:
                        failure = "worker_exit"; break
                    st = read(directory / "state.json"); phase = st.get("phase", "loading_instance")
                    begin = st.get("benchmark_start_unix"); mem = memory(proc); peak = max(peak, mem)
                    if mem > MEMORY_LIMIT:
                        failure = "memory_limit_during_" + phase; break
                    if begin and phase not in ("postprocessing", "done") and time.time() - begin > TIME_LIMIT + 15:
                        failure = "time_limit_during_" + phase; break
                    if not begin and time.time() - dispatch > TIME_LIMIT + 600:
                        failure = "loading_or_audit_guard"; break
                    if begin and time.time() - begin > TIME_LIMIT + 1200:
                        failure = "postprocessing_timeout"; break
                    time.sleep(1)
                result = read(resultfile)
                if failure:
                    stop_worker()
                    result.update(status=failure, city=job["city"], size=job["size"], time_budget_seconds=TIME_LIMIT)
                result["peak_private_bytes"] = peak
                write(resultfile, result)
                status = result.get("status", "complete"); completed.append(dict(id=job["id"], result=result))
                print(job["id"], status, result.get("model_seconds", result.get("wall_seconds")), flush=True)
            if not audit and status not in ["OPTIMAL", "ALMOST_OPTIMAL"]:
                if group == "aggregation":
                    blocked_city.add(job["city"])
                else:
                    blocked_fraction.add((job["city"], job["station_budget"] / job["candidate_count"]))
            write(progress, dict(phase="running", completed=len(completed), total_jobs=len(jobs),
                                 last_case=job["id"], last_status=status, skipped=skipped))
        write(progress, dict(phase="complete", completed=len(completed), total_jobs=len(jobs), skipped=skipped,
                             finished_unix=time.time()))
    except BaseException as e:
        write(progress, dict(phase="error", error=str(e), completed=len(completed), total_jobs=len(jobs), skipped=skipped))
        raise
    finally:
        stop_worker()


for g in GROUPS:
    run(g)
