import csv
import logging
from pathlib import Path

import pytest

pytest.importorskip("boto3")
pytest.importorskip("pandas")
pytest.importorskip("reportlab")

from bactipipe.__version__ import __version__
from bactipipe.downstream_manifest import (
    DownstreamManifestError,
    read_downstream_manifest,
    write_command_sample_sheet,
)
from bactipipe.scripts import core, run_traits, traits_reconcile, traits_report, type_genomes, utils


def _logger():
    logger = logging.getLogger("bactipipe-test")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    return logger


def test_reported_version_uses_package_version():
    assert utils.bactipipe_version == f"v{__version__}"


def test_detect_parser_accepts_hyphenated_web_arguments():
    args = run_traits.parse_args(
        [
            "--analyst", "Analyst",
            "--genomes-dir", "/assemblies",
            "--sample-sheet", "/manifest.tsv",
            "--run-name", "run-1",
            "--organism", "Listeria",
            "--min-identity", "91",
            "--min-coverage", "61",
        ]
    )
    assert args.genomes_dir == "/assemblies"
    assert args.sample_sheet == "/manifest.tsv"
    assert args.min_identity == 91
    assert args.min_coverage == 61


def test_detect_sample_sheet_keeps_stable_id_display_name_and_specimen(tmp_path):
    genomes = tmp_path / "assemblies"
    genomes.mkdir()
    (genomes / "S01.fasta").write_text(">contig\nACGT\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("S01\tDisplay 01\tbrain\n", encoding="utf-8")

    assert run_traits._read_sample_sheet(str(sheet), str(genomes)) == [
        ("S01", "Display 01", "brain", str(genomes / "S01.fasta"))
    ]


def test_detect_rejects_duplicate_sample_ids(tmp_path):
    genomes = tmp_path / "assemblies"
    genomes.mkdir()
    (genomes / "S01.fasta").write_text(">contig\nACGT\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("S01\tFirst\tbrain\nS01\tSecond\tblood\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Duplicate isolate IDs"):
        run_traits._read_sample_sheet(str(sheet), str(genomes))


def test_relate_sample_sheet_retains_missing_assemblies_and_rejects_duplicates(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample\tisolate\tspecimen\tdate\nS01\tDisplay 01\tbrain\t2026-08-14\n",
        encoding="utf-8",
    )
    rows = type_genomes._read_sample_sheet(sheet, _logger())
    assert rows == [("S01", "Display 01", "brain")]

    manifest = type_genomes._build_manifest(rows, tmp_path / "assemblies", _logger())
    assert manifest[0].sample == "S01"
    assert manifest[0].isolate == "Display 01"
    assert manifest[0].exists is False

    sheet.write_text("S01\tFirst\tbrain\nS01\tSecond\tblood\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate sample IDs"):
        type_genomes._read_sample_sheet(sheet, _logger())


def test_amrfinder_is_preferred_during_duplicate_reconciliation():
    common = {"determinant": "blaTEM", "contig": "c1", "start": "1", "end": "100"}
    merged = traits_reconcile.merge_for_sample(
        "S01",
        [
            {**common, "tool": "abricate", "identity": "99"},
            {**common, "tool": "amrfinder", "identity": "100"},
        ],
        [],
    )
    assert merged["amr"] == [{**common, "tool": "amrfinder", "identity": "100"}]


def test_sample_and_consolidated_reports_keep_detection_fields(tmp_path):
    merged = {
        "amr": [
            {"sample": "S01", "determinant": "blaTEM", "type": "AMR", "phenotype": "beta-lactam"},
            {"sample": "S01", "determinant": "gyrA_S83F", "type": "mutation", "phenotype": "quinolone"},
        ],
        "vf": [
            {"sample": "S01", "gene": "stx2a", "base_gene": "stx2", "category": "Toxin", "note": "Subtype a"}
        ],
    }
    traits_report.write_sample_tsvs("S01", merged, str(tmp_path))
    with (tmp_path / "S01.vf.tsv").open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["base_gene"] == "stx2"
    assert row["note"] == "Subtype a"

    output = Path(
        traits_report.write_consolidated_tsv(
            {"Display 01": merged}, str(tmp_path), "cohort.tsv"
        )
    ).read_text(encoding="utf-8")
    assert "Acquired\tblaTEM\tbeta-lactam\t1" in output
    assert "Mutation\tgyrA_S83F\tquinolone\t1" in output


def test_conda_probe_uses_direct_argv(monkeypatch):
    observed = {}

    def fake_run(command, **kwargs):
        observed["command"] = command
        return object()

    monkeypatch.setattr(run_traits.subprocess, "run", fake_run)
    run_traits._conda_run("viramr", ["amrfinder", "--list_organisms"])
    assert observed["command"] == [
        "conda", "run", "-n", "viramr", "amrfinder", "--list_organisms"
    ]


def test_big_sdb_scheme_is_written_in_legacy_cge_layout(tmp_path, monkeypatch):
    scheme_url = "https://example.test/db/species/schemes/1"
    locus_urls = [
        "https://example.test/db/species/loci/geneA",
        "https://example.test/db/species/loci/geneB",
    ]
    responses = {
        scheme_url: (
            '{"description":"MLST","last_updated":"2026-08-14",'
            '"records":1,"loci":['
            '"https://example.test/db/species/loci/geneA",'
            '"https://example.test/db/species/loci/geneB"],'
            '"profiles_csv":"https://example.test/profiles"}'
        ).encode(),
        "https://example.test/profiles": b"ST\tgeneA\tgeneB\n1\t1\t2\n",
        f"{locus_urls[0]}/alleles_fasta": b">geneA_1\nACGT\n",
        f"{locus_urls[1]}/alleles_fasta": b">geneB_2\nTGCA\n",
    }
    monkeypatch.setattr(core, "_download_url", responses.__getitem__)

    source = core._download_legacy_mlst_scheme(
        str(tmp_path), "example", "Example species", scheme_url
    )

    assert source["loci"] == ["geneA", "geneB"]
    assert (tmp_path / "example" / "example.tsv").read_text() == (
        "ST\tgeneA\tgeneB\n1\t1\t2\n"
    )
    assert (tmp_path / "example" / "example.fsa").read_text() == (
        ">geneA_1\nACGT\n>geneB_2\nTGCA\n"
    )


def test_downstream_manifest_groups_independent_commands_and_reference(tmp_path):
    manifest = tmp_path / "run_downstream.tsv"
    manifest.write_text(
        "cohort_id\tcommand\torganism\tsample_id\tdisplay_name\tspecimen\trole\tmin_identity\n"
        "outbreak_a\trelate\tsalmonella\tS01\tCase 1\tfeces\treference\t\n"
        "outbreak_a\trelate\tsalmonella\tS02\tCase 2\tfeces\tsample\t\n"
        "traits_a\tdetect\tEscherichia\tS03\tCase 3\tbrain\tsample\t91\n",
        encoding="utf-8",
    )

    cohorts = read_downstream_manifest(manifest)
    assert [cohort.cohort_id for cohort in cohorts] == ["outbreak_a", "traits_a"]
    assert cohorts[0].reference.sample_id == "S01"
    assert cohorts[1].options["min_identity"] == 91.0
    assert cohorts[0].signature != cohorts[1].signature

    relate_sheet = write_command_sample_sheet(cohorts[0], tmp_path / "relate.tsv")
    assert relate_sheet.read_text(encoding="utf-8").splitlines()[0] == (
        "sample\tisolate\tspecimen\tdate"
    )
    detect_sheet = write_command_sample_sheet(cohorts[1], tmp_path / "detect.tsv")
    assert detect_sheet.read_text(encoding="utf-8") == "S03\tCase 3\tbrain\n"


def test_downstream_manifest_rejects_path_like_cohort_and_conflicts(tmp_path):
    manifest = tmp_path / "run_downstream.tsv"
    manifest.write_text(
        "cohort_id\tcommand\torganism\tsample_id\n"
        "../unsafe\trelate\tsalmonella\tS01\n"
        "../unsafe\tdetect\tsalmonella\tS02\n",
        encoding="utf-8",
    )
    with pytest.raises(DownstreamManifestError, match="Cohort ID"):
        read_downstream_manifest(manifest)
