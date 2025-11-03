from typing import List, Dict, Tuple

def _key_amr(r: Dict[str,str]) -> Tuple:
    return (r.get("determinant",""), r.get("contig",""), r.get("start",""), r.get("end",""))

def _key_vf(r: Dict[str,str]) -> Tuple:
    return (r.get("gene",""), r.get("contig",""), r.get("start",""), r.get("end",""))

def merge_for_sample(sample: str, amr_rows: List[Dict[str,str]], vf_rows: List[Dict[str,str]]):
    # Deduplicate AMR: prefer AMRFinderPlus when overlapping
    amr_seen = {}
    amr_merged: List[Dict[str,str]] = []
    for r in amr_rows:
        k = _key_amr(r)
        if k in amr_seen:
            i = amr_seen[k]
            if amr_merged[i]["tool"] != "AMRFinderPlus" and r["tool"] == "AMRFinderPlus":
                amr_merged[i] = r
        else:
            amr_seen[k] = len(amr_merged)
            amr_merged.append(r)

    # VF already preferred in traits_vf
    vf_merged = vf_rows

    return {"amr": amr_merged, "vf": vf_merged}
