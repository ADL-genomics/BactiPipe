# caps.py
# Central registry for organism capabilities (grows over time)

CAPS = {
    "salmonella": {
        "mlst": True,
        "cgmlst": True,
        "mlst_scheme": "senterica",     # what your mlst.py expects via -s
        "cgmlst_scheme": "salmonella",  # what cgMLST.py expects via -s
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": "cgmlstfinder_db",
        "serotyper": "seqsero2",
    },
    "ecoli": {
        "mlst": True, "cgmlst": True,
        "mlst_scheme": "ecoli",
        "cgmlst_scheme": "ecoli",
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": "cgmlstfinder_db",
        "serotyper": "serotypefinder",
    },
    "klebsiella": {
        "mlst": True, "cgmlst": False,
        "mlst_scheme": "kpneumoniae",
        "cgmlst_scheme": None,
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": None,
        "serotyper": "kleborate",
    },
    "clostridium": {
        "mlst": True, "cgmlst": True,
        "mlst_scheme": "cperfringens",
        "cgmlst_scheme": "clostridium",
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir":"cgmlstfinder_db",
        "serotyper": "kleborate",
    },
}


def supported_rows():
    def as_str(x: object) -> str:
        return x if isinstance(x, str) and x else "-"
    rows = []
    for org, spec in CAPS.items():
        rows.append([
            org,
            "yes" if spec.get("mlst") else "no",
            as_str(spec.get("mlst_scheme")),
            "yes" if spec.get("cgmlst") else "no",
            as_str(spec.get("cgmlst_scheme")),
            as_str(spec.get("serotyper")),
        ])
    return rows
