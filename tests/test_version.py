import sys

import pytest

from bactipipe.__version__ import __version__
from bactipipe.cli import main


@pytest.mark.parametrize("flag", ["-v", "--version"])
def test_cli_reports_authoritative_package_version(monkeypatch, capsys, flag):
    monkeypatch.setattr(sys, "argv", ["bactipipe", flag])

    with pytest.raises(SystemExit) as exc_info:
        main()

    assert exc_info.value.code == 0
    assert capsys.readouterr().out.strip() == f"bactipipe {__version__}"
