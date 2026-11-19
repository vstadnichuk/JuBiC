"""Normalize the results of `run_experiments.py scaling` into data/scaling/ (Section 6.3).

Reads jobs_scaling.json, progress_scaling.json, raw/scaling_cases and runtime/scaling. Main outputs:
  solver_cases.csv              status, enumeration/assembly/solve times, model size and paths
  user_enumeration_summary.csv  per-user enumeration time statistics (median, 99th percentile, max)
  optimized_modal_shares.csv    modal split of the optimal designs
"""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import csv, json

from paths import B, CODE
OUT = B / "data/scaling"
OUT.mkdir(parents=True, exist_ok=True)

def load(path):
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError):
        return {}

def frac(value):
    return Fraction(str(value).replace("//", "/"))

def write(name, rows):
    if not rows:
        return
    with (OUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def true(x):
    return str(x).strip().lower() == "true"

def main():
    jobs = load(B / "jobs_scaling.json")
    progress = load(B / "progress_scaling.json")
    skipped = {x["id"]: x["reason"] for x in progress.get("skipped", [])}
    cases, modes, auditrows, warnings, eventstats = [], [], [], [], []
    for job in jobs:
        directory = B / ("runtime/aggregation" if job["kind"] == "audit" else "raw/scaling_cases") / job["id"]
        if job["kind"] == "audit":
            result = load(directory / "summary.json")
            if not result:
                continue
            users = []
            if (directory / "users.csv").exists():
                users = list(csv.DictReader((directory / "users.csv").open()))
            complete = [u for u in users if not true(u["raw_timeout"]) and not true(u["dominance_timeout"])]
            raw = sum(int(u["raw_paths"]) for u in complete)
            dom = sum(int(u["dominance_paths"]) for u in complete)
            legacy = sum(int(u["legacy_arc_dominance_paths"]) for u in complete)
            certificates = sum(int(u["dominance_not_in_raw"]) + int(u["certificate_failures"]) for u in complete)
            auditrows.append(dict(
                city=job["city"], users=2*job["size"]*(job["size"]-1), complete_users=len(complete),
                incomplete_users=2*job["size"]*(job["size"]-1)-len(complete), raw_paths=raw,
                dominance_paths=dom, legacy_arc_dominance_paths=legacy,
                reduction_pct=100*(1-dom/raw) if raw else None,
                users_station_below_legacy=sum(int(u["dominance_paths"]) < int(u["legacy_arc_dominance_paths"]) for u in complete),
                user_summed_raw_seconds=sum(float(u["raw_seconds"]) for u in complete),
                user_summed_dominance_seconds=sum(float(u["dominance_seconds"]) for u in complete),
                user_summed_legacy_seconds=sum(float(u["legacy_seconds"]) for u in complete),
                comparison_wall_seconds=result.get("wall_seconds"),
                one_path_users=sum(int(u["dominance_paths"]) == 1 for u in complete),
                one_path_pct=100*sum(int(u["dominance_paths"]) == 1 for u in complete)/len(complete) if complete else None,
                paths_per_user=dom/len(complete) if complete else None,
                certificate_failures=certificates,
                summary_certificate_failures=result.get("certificate_failures"),
                summary_incomplete_users=result.get("incomplete_users"),
                status=result.get("status", "complete")))
            continue
        result = load(directory / "result.json")
        instance = load(directory / "instance.json")
        state = load(directory / "state.json")
        total = result.get("model_seconds")
        time_source = "worker_measurement"
        if total is None and result and state.get("benchmark_start_unix"):
            total = max(0, (directory / "result.json").stat().st_mtime-state["benchmark_start_unix"])
            time_source = "guard_termination_checkpoint"
        row = dict(
            id=job["id"], city=job["city"], locations=job["size"],
            candidates=job["candidate_count"], budget=job["station_budget"],
            budget_fraction=job["station_budget"]/job["candidate_count"],
            groups="+".join(job["groups"]),
            status=result.get("status", "skipped" if job["id"] in skipped else "pending"),
            skipped_reason=skipped.get(job["id"], ""),
            users=result.get("users", instance.get("nusers")), layered_vertices=instance.get("nnodes"), arcs=instance.get("narcs"), path_count=result.get("path_count"),
            active_users=result.get("modeled_path_users"),
            fixed_singletons=result.get("eliminated_fixed_path_users"),
            variables=result.get("variables"), constraints=result.get("constraints"),
            enumeration_seconds=result.get("enumeration_seconds"),
            assembly_seconds=max(0, result["build_seconds"]-result["enumeration_seconds"])
            if "build_seconds" in result and "enumeration_seconds" in result else None,
            solve_seconds=result.get("solve_wall_seconds"),
            total_seconds=total, time_source=time_source,
            final_phase=state.get("phase"), users_finished=state.get("users_finished"),
            users_assembled=state.get("users_assembled"),
            gap=result.get("gap"), objective=result.get("objective"),
            opened_stations=result.get("opened_stations"),
            peak_private_gib=result.get("peak_private_bytes", 0)/2**30,
            intrazonal_pct=100*float(frac(instance.get("intrazonal_owner_mass", 0))
                                    + frac(instance.get("intrazonal_nonowner_mass", 0))))
        cases.append(row)
        path = directory / "solution_modes.csv"
        if path.exists():
            values = list(csv.DictReader(path.open()))
            mass = sum(float(x["demand"]) for x in values)
            expected = 1-row["intrazonal_pct"]/100
            assert abs(mass-expected) < 1e-5, (job["id"], mass, expected)
            families = defaultdict(float)
            for x in values:
                families[x["mode"]] += float(x["demand"])
            for mode in ["walking", "car", "transit", "bike_sharing", "multi_powered"]:
                modes.append(dict(id=job["id"], city=job["city"], locations=job["size"],
                                  candidates=job["candidate_count"], budget=job["station_budget"],
                                  mode=mode, original_demand_pct=100*families[mode],
                                  conditional_pct=100*families[mode]/mass,
                                  intrazonal_pct=row["intrazonal_pct"], status=row["status"]))
        path = directory / "gurobi.log"
        if path.exists():
            for line in path.read_text(errors="replace").splitlines():
                if any(word in line.lower() for word in
                       ["warning", "numerical", "unscaled primal", "small coefficients",
                        "large matrix coefficient", "large objective coefficient"]):
                    warnings.append(dict(case=job["id"], message=line.strip()))
        path = directory / "user_events.csv"
        if path.exists():
            byclass = defaultdict(list)
            finished_events = []
            with path.open() as stream:
                for event in csv.DictReader(stream):
                    if event["event"] == "finished":
                        byclass[event["car_available"]].append(event)
                        finished_events.append(event)
            if (row["enumeration_seconds"] is None and row["users"] is not None
                    and len(finished_events) == int(row["users"])):
                row["enumeration_seconds"] = max(float(x["case_elapsed_seconds"])
                                                 for x in finished_events)
                if row["total_seconds"] is not None:
                    row["assembly_seconds"] = max(0, row["total_seconds"]-row["enumeration_seconds"])
                row["time_source"] += "+enumeration_completion_trace"
            for car, entries in byclass.items():
                times = sorted(float(x["enumeration_seconds"]) for x in entries)
                eventstats.append(dict(
                    case=job["id"], car_available=car, completed_users=len(entries),
                    one_path_users=sum(int(x["path_count"]) == 1 for x in entries),
                    paths=sum(int(x["path_count"]) for x in entries),
                    median_enumeration_seconds=times[len(times)//2],
                    p99_enumeration_seconds=times[min(len(times)-1, int(.99*len(times)))],
                    maximum_enumeration_seconds=max(times),
                    timed_out_users=sum(x["timed_out"] == "true" for x in entries)))
    for log in (B / "runtime/scaling").glob("worker_*.log"):
        current = None
        for line in log.read_text(errors="replace").splitlines():
            if line.startswith("CASE_START "):
                current = line.split(" ",1)[1].strip()
            elif line.startswith("CASE_DONE "):
                current = None
            elif current and any(word in line.lower() for word in
                                 ["warning", "numerical trouble", "small coefficients",
                                  "large matrix coefficient", "large objective coefficient"]):
                warnings.append(dict(case=current, message=line.strip()))
    write("solver_cases.csv", cases)
    write("optimized_modal_shares.csv", modes)
    write("dominance50.csv", auditrows)
    warnings = list({(w["case"],w["message"]):w for w in warnings}.values())
    write("solver_warnings.csv", warnings)
    write("user_enumeration_summary.csv", eventstats)


    summary = dict(
        fresh_cost_eur_km=.22, configured_solver_cases=len(cases),
        results=sum(x["status"] not in ["pending", "skipped"] for x in cases),
        optimal=sum(x["status"] == "OPTIMAL" for x in cases),
        skipped=sum(x["status"] == "skipped" for x in cases),
        pending=sum(x["status"] == "pending" for x in cases),
        statuses={s: sum(x["status"] == s for x in cases) for s in sorted({x["status"] for x in cases})},
        audit_cities=len(auditrows),
        certificate_failures=sum(x["certificate_failures"] for x in auditrows),
        incomplete_audit_users=sum(x["incomplete_users"] for x in auditrows),
        legacy_equals_station_paths=all(x["legacy_arc_dominance_paths"] == x["dominance_paths"] for x in auditrows),
        warning_lines=len(warnings), controller_phase=progress.get("phase"),
        strict_checks_passed=all(x["certificate_failures"] == 0 for x in auditrows))
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
