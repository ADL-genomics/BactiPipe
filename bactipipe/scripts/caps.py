# caps.py
# Central registry for organism capabilities (grows over time)

CAPS = {
    "abaumannii": {
        "mlst": True, "cgmlst": False,
        "mlst_scheme": "abaumannii",
        "cgmlst_scheme": None,
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": None,
        "serotyper":"",
    },
    "clostridium": {
        "mlst": True, "cgmlst": True,
        "mlst_scheme": "cperfringens",
        "cgmlst_scheme": "clostridium",
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": "cgmlstfinder_db",
        "serotyper": "",
    },
    "cfetus": {
        "mlst": True, "cgmlst": True,
        "mlst_scheme": "cfetus",
        "cgmlst_scheme": "campylobacter",
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": "cgmlstfinder_db",
        "serotyper": "",
    },
    "cjejuni": {
        "mlst": True, "cgmlst": True,
        "mlst_scheme": "cjejuni",
        "cgmlst_scheme": "campylobacter",
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": "cgmlstfinder_db",
        "serotyper": "",
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
    "listeria": {
        "mlst": True, "cgmlst": False,
        "mlst_scheme": "lmonocytogenes",
        "cgmlst_scheme": None,
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": None,
        "serotyper": "",
    },
    "pmultocida": {
        "mlst": True, "cgmlst": False,
        "mlst_scheme": "pmultocida",
        "cgmlst_scheme": None,
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": None,
        "serotyper": "",
    },
    "salmonella": {
        "mlst": True,
        "cgmlst": True,
        "mlst_scheme": "senterica",
        "cgmlst_scheme": "salmonella",
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": "cgmlstfinder_db",
        "serotyper": "seqsero2",
    },
    "sepidermidis": {
        "mlst": True, "cgmlst": False,
        "mlst_scheme": "sepidermidis",
        "cgmlst_scheme": None,
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": None,
        "serotyper": "",
    },
    "saureus": {
        "mlst": True, "cgmlst": False,
        "mlst_scheme": "saureus",
        "cgmlst_scheme": None,
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": None,
        "serotyper": "",
    },
    "spseudintermedius": {
        "mlst": True, "cgmlst": False,
        "mlst_scheme": "spseudintermedius",
        "cgmlst_scheme": None,
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": None,
        "serotyper": "",
    },
    "szooepidemicus": {
        "mlst": True, "cgmlst": False,
        "mlst_scheme": "szooepidemicus",
        "cgmlst_scheme": None,
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": None,
        "serotyper": "",
    },
    "orhinotracheale": {
        "mlst": True, "cgmlst": False,
        "mlst_scheme": "orhinotracheale",
        "cgmlst_scheme": None,
        "mlst_db_subdir": "mlst_db",
        "cgmlst_db_subdir": None,
        "serotyper": "",
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
