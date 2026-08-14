from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from typing import Iterable, List, Optional, Sequence, Tuple, Union

import pandas as pd

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


# ──────────────────────────────────────────────────────────────────────────────
# Layout configuration
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class LayoutConfig:
    # Page & margins
    page_size: Tuple[float, float] = landscape(letter)
    left_margin: float = 36
    right_margin: float = 36
    top_margin: float = 48
    bottom_margin: float = 36

    # Dual header band split (portion of usable width for left vs right panels)
    left_panel_ratio: float = 0.65
    right_panel_ratio: float = 0.35

    # Column width strategy for main QC table
    reserved_cols: Tuple[int, int] = (2, 3)     # Expected & Identified organism (0-based)
    reserved_width: float = 2.0 * inch
    last3_narrow_width: float = 0.9 * inch

    # Typography (separate sizes for header blocks vs the table)
    header_block_font_size: int = 10            # used for Run Metadata + Tools blocks
    table_header_font_size: int = 8             # table header row
    table_body_font_size: int = 8               # table data rows

    # Header wrapping threshold for table columns
    wrap_threshold: int = 12

    # Body alignment: columns to LEFT-align (0-based). Others are centered.
    left_align_cols: Tuple[int, ...] = (0, 2, 3)

    # Conditional coloring
    pass_color: colors.Color = colors.HexColor("#1b5e20")
    fail_color: colors.Color = colors.HexColor("#b71c1c")


# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

def render_pdf_from_rows(
    rows: Sequence[Sequence[object]],
    header_labels: Sequence[str],
    pdf_path: str,
    run_metadata: Optional[Iterable[str]] = None,
    legend: Optional[Iterable[str]] = None,
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
    tools_header: Optional[Sequence[Tuple[str, str]]] = None,  # list of (tool, version/cutoff)
    layout: LayoutConfig = LayoutConfig(),
    run_name: Optional[str] = None,  # NEW: explicit run name override
) -> str:
    """Render a QC report PDF from in-memory rows."""
    df = pd.DataFrame(rows, columns=list(header_labels))
    # If no explicit run_name, try to infer from run_metadata
    if not run_name and run_metadata:
        run_name = _infer_run_name_from_metadata(list(run_metadata))
    return _render_pdf_common(
        df=df,
        pdf_path=pdf_path,
        run_metadata=list(run_metadata) if run_metadata else None,
        legend=list(legend) if legend else None,
        title=title,
        subtitle=subtitle,
        tools_header=list(tools_header) if tools_header else None,
        layout=layout,
        run_name=run_name,
    )


def render_pdf_from_tsv(
    tsv_path: str,
    pdf_path: str,
    run_metadata: Optional[Iterable[str]] = None,
    legend: Optional[Iterable[str]] = None,
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
    layout: LayoutConfig = LayoutConfig(),
    run_name: Optional[str] = None,  # NEW: explicit run name override
) -> str:
    """Render a QC report PDF directly from an already written TSV file."""
    lines = _read_lines(tsv_path)
    header_idx = _find_table_header_idx(lines)
    if header_idx is None:
        raise RuntimeError("Could not locate the table header line containing 'Sample ID'.")

    # Parse region above the table: left metadata lines (with any '====' preserved)
    # and a tools span (tool >> value) surrounded by separators if present.
    top_block, tools_pairs = _parse_header_blocks(lines, header_idx)

    table_text = "\n".join(lines[header_idx:])
    df = pd.read_csv(StringIO(table_text), sep="\t")

    run_meta = list(run_metadata) if run_metadata else (top_block or None)

    # Infer run name if not provided
    if not run_name:
        run_name = _infer_run_name_from_metadata(run_meta)

    return _render_pdf_common(
        df=df,
        pdf_path=pdf_path,
        run_metadata=run_meta,
        legend=list(legend) if legend else None,
        title=title,
        subtitle=subtitle,
        tools_header=tools_pairs or None,
        layout=layout,
        run_name=run_name,
    )


# ──────────────────────────────────────────────────────────────────────────────
# Parsing helpers
# ──────────────────────────────────────────────────────────────────────────────

def _read_lines(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as f:
        return f.read().splitlines()


def _find_table_header_idx(lines: List[str]) -> Optional[int]:
    for i, line in enumerate(lines):
        if "Sample ID" in line and "\t" in line:
            return i
    return None


def _is_separator(line: str) -> bool:
    # Be permissive: allow lines that are mostly '=' after stripping whitespace.
    s = line.strip()
    if len(s) < 3:
        return False
    eq_count = sum(1 for ch in s if ch == '=')
    return eq_count >= max(3, int(0.8 * len(s)))


def _parse_header_blocks(lines: List[str], header_idx: int) -> Tuple[List[str], List[Tuple[str, str]]]:
    """
    Parse pre-table region (lines[:header_idx]) into:
      - top_block: metadata lines (preserve as-is, including any '====' lines)
      - tools_pairs: list of (tool, value) parsed from lines containing '>>'
                     that are inside the last [==== ... ====] span in the region.
    If no such tools span is found, everything is metadata and tools_pairs is empty.
    """
    region = lines[:header_idx]

    # Identify all separator positions
    sep_positions = [i for i, ln in enumerate(region) if _is_separator(ln)]

    # Collect candidate spans [start, end] bounded by separators
    spans: List[Tuple[int, int]] = []
    for i in range(len(sep_positions) - 1):
        start = sep_positions[i]
        end = sep_positions[i + 1]
        if end > start + 1:
            spans.append((start, end))

    # Choose the last span that contains at least one '>>' line as the tools block
    tools_span: Optional[Tuple[int, int]] = None
    for start, end in spans:
        if any(">>" in ln for ln in region[start + 1:end]):
            tools_span = (start, end)

    if tools_span:
        start, end = tools_span
        # Metadata: everything BEFORE the tools span (preserve including any '====')
        top_block = [ln for ln in region[:start] if ln.strip()]
        # Tools: parse '>>' lines only
        tools_pairs = []
        for ln in region[start + 1:end]:
            if ">>" in ln:
                k, v = ln.split(">>", 1)
                tools_pairs.append((k.strip(), v.strip()))
        return top_block, tools_pairs

    # No tools span found → treat all as metadata
    top_block = [ln for ln in region if ln.strip()]
    return top_block, []


def _infer_run_name_from_metadata(run_meta: Optional[List[str]]) -> Optional[str]:
    if not run_meta:
        return None
    for ln in run_meta:
        if ln.lower().startswith("run:"):
            return ln.split(":", 1)[1].strip()
    return None


# ──────────────────────────────────────────────────────────────────────────────
# Styles & layout helpers
# ──────────────────────────────────────────────────────────────────────────────

def _styles(layout: LayoutConfig):
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="TitleCenter", parent=styles["Title"], alignment=TA_CENTER))
    styles.add(ParagraphStyle(
        name="HeaderBlock",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=layout.header_block_font_size,
        leading=layout.header_block_font_size + 2,
    ))
    styles.add(ParagraphStyle(name="Header", parent=styles["Heading2"], spaceAfter=4))
    styles.add(ParagraphStyle(name="HeaderRight", parent=styles["Heading2"], spaceAfter=4))
    styles.add(ParagraphStyle(
        name="Separator",
        parent=styles["Normal"],
        fontName="Courier",
        fontSize=layout.header_block_font_size,
        leading=layout.header_block_font_size + 2,
    ))
    # NEW: style for wrapped organism names in table body
    styles.add(ParagraphStyle(
        name="TableBodyPara",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=layout.table_body_font_size,
        leading=layout.table_body_font_size + 2,
    ))
    return styles


def _page_footer_factory(run_name: Optional[str]):
    """Returns a ReportLab onPage callback that draws the desired footer."""
    prefix = "BactiPipe QC Report"
    if run_name:
        label = f"{prefix} - Run {run_name}"
    else:
        label = prefix

    def _footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        text = f"{label}    |    Page {doc.page}"
        # Right-align to the page's right margin
        canvas.drawRightString(doc.pagesize[0] - 36, 20, text)
        canvas.restoreState()

    return _footer


def _wrap_header(text: str, threshold: int = 12) -> str:
    if len(text) <= threshold:
        return text
    parts = text.split()
    if len(parts) > 1:
        mid = len(parts) // 2
        return "\n".join([" ".join(parts[:mid]), " ".join(parts[mid:])])
    mid = len(text) // 2
    return text[:mid] + "\n" + text[mid:]


def _compute_col_widths(df: pd.DataFrame, layout: LayoutConfig) -> List[float]:
    page_w, _ = layout.page_size
    usable = page_w - (layout.left_margin + layout.right_margin)

    n_cols = df.shape[1]
    last3_total = 3 * layout.last3_narrow_width
    reserved_total = 2 * layout.reserved_width
    remaining = max(usable - (last3_total + reserved_total), 1.0)

    other_cols = n_cols - 5  # all except last 3 and reserved 2
    other_width = remaining / other_cols if other_cols > 0 else remaining

    widths: List[float] = []
    for i in range(n_cols):
        if i in layout.reserved_cols:
            widths.append(layout.reserved_width)
        elif i >= n_cols - 3:
            widths.append(layout.last3_narrow_width)
        else:
            widths.append(other_width)
    return widths


# ──────────────────────────────────────────────────────────────────────────────
# Building sections
# ──────────────────────────────────────────────────────────────────────────────

def _build_dual_header_band(
    run_metadata_lines: Optional[List[str]],
    tools_header: Optional[Sequence[Tuple[str, str]]],
    layout: LayoutConfig,
    styles,
) -> Optional[Table]:
    """
    Two-column band:
      Left panel:  "Run Metadata" + lines (preserves any '====' lines from TSV)
      Right panel: "Tools & Cutoffs" + '====' under title + (tool,value) rows + '====' after list
    Right panel sits near the right edge via layout ratios.
    """
    if not run_metadata_lines and not tools_header:
        return None

    # LEFT: Run Metadata (preserve separators from TSV)
    left_flow: List[Union[Paragraph, Table]] = []
    if run_metadata_lines:
        left_flow.append(Paragraph("Run Metadata", styles["Header"]))
        for ln in run_metadata_lines:
            if _is_separator(ln):
                left_flow.append(Paragraph(ln.strip(), styles["Separator"]))
            else:
                left_flow.append(Paragraph(str(ln), styles["HeaderBlock"]))

    # RIGHT: Tools & Cutoffs (explicit separators under title and after the list)
    right_flow: List[Union[Paragraph, Table]] = []
    if tools_header:
        right_flow.append(Paragraph("Tool version Numbers & QC Cutoffs", styles["HeaderRight"]))
        right_flow.append(Paragraph("=" * 35, styles["Separator"]))  # under title

        # Build tool table with Paragraph cells (same leading as metadata)
        tool_data = [
            [Paragraph(str(k), styles["HeaderBlock"]),
             Paragraph(str(v), styles["HeaderBlock"])]
            for (k, v) in tools_header
        ]
        tool_table = Table(tool_data, hAlign="LEFT")
        tool_table.setStyle(TableStyle([
            ("ALIGN", (0,0), (-1,-1), "LEFT"),
            ("VALIGN", (0,0), (-1,-1), "TOP"),
            ("LEFTPADDING", (0,0), (-1,-1), 0),
            ("RIGHTPADDING", (0,0), (-1,-1), 0),
            ("TOPPADDING", (0,0), (-1,-1), 0),
            ("BOTTOMPADDING", (0,0), (-1,-1), 0),  # compact vertically
        ]))
        right_flow.append(tool_table)
        right_flow.append(Paragraph("=" * 35, styles["Separator"]))  # after list

    # Parent 2-col table pushed to the right edge for tools
    page_w, _ = layout.page_size
    usable = page_w - (layout.left_margin + layout.right_margin)

    left_ratio = max(0.05, min(0.95, layout.left_panel_ratio))
    right_ratio = max(0.05, min(0.95, layout.right_panel_ratio))
    total = left_ratio + right_ratio or 1.0
    left_ratio /= total
    right_ratio /= total

    col_widths = [usable * left_ratio, usable * right_ratio]

    parent = Table([[left_flow, right_flow]], colWidths=col_widths)
    parent.setStyle(TableStyle([
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("ALIGN", (0,0), (-1,-1), "LEFT"),
        ("LEFTPADDING", (0,0), (-1,-1), 0),
        ("RIGHTPADDING", (0,0), (-1,-1), 0),
        ("TOPPADDING", (0,0), (-1,-1), 0),
        ("BOTTOMPADDING", (0,0), (-1,-1), 0),
    ]))
    return parent


def _build_table(df: pd.DataFrame, layout: LayoutConfig) -> Table:
    styles = _styles(layout)

    # Headers: wrap except reserved organism columns
    wrapped_headers = []
    for i, c in enumerate(df.columns):
        txt = str(c)
        if i in layout.reserved_cols:
            wrapped_headers.append(txt)  # no wrapping for these headers
        else:
            wrapped_headers.append(_wrap_header(txt, layout.wrap_threshold))

    # Body rows: wrap only organism columns (reserved_cols)
    data = [wrapped_headers]
    org_cols = set(layout.reserved_cols)  # usually {2, 3}
    for row in df.astype(str).values.tolist():
        r = []
        for i, cell in enumerate(row):
            if i in org_cols:
                r.append(Paragraph(cell, styles["TableBodyPara"]))
            else:
                r.append(cell)
        data.append(r)

    # Column widths
    col_widths = _compute_col_widths(df, layout)
    t = Table(data, repeatRows=1, colWidths=col_widths)

    style = TableStyle([
        # Header row
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#f0f0f0")),
        ("TEXTCOLOR", (0,0), (-1,0), colors.black),
        ("VALIGN", (0,0), (-1,0), "MIDDLE"),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTSIZE", (0,0), (-1,0), layout.table_header_font_size),
        ("BOTTOMPADDING", (0,0), (-1,0), 6),

        # Body defaults
        ("ALIGN", (0,1), (-1,-1), "CENTER"),
        ("VALIGN", (0,1), (-1,-1), "MIDDLE"),
        ("GRID", (0,0), (-1,-1), 0.25, colors.grey),
        ("FONTSIZE", (0,1), (-1,-1), layout.table_body_font_size),
    ])

    # Alternating row backgrounds
    for i in range(1, len(data)):
        style.add("BACKGROUND", (0,i), (-1,i),
                  colors.HexColor("#fcfcfc") if i % 2 == 1 else colors.white)

    # Match header alignment to body
    n_cols = df.shape[1]
    style.add("ALIGN", (0,0), (-1,0), "CENTER")  # default header center
    for col in layout.left_align_cols:
        if 0 <= col < n_cols:
            style.add("ALIGN", (col, 0), (col, 0), "LEFT")   # header col
            style.add("ALIGN", (col, 1), (col, -1), "LEFT")  # body col

    # Conditional coloring for Overall Quality
    if "Overall Quality" in df.columns:
        oq_idx = df.columns.get_loc("Overall Quality")
        for r_idx, val in enumerate(df["Overall Quality"].astype(str), start=1):
            v = val.strip().lower()
            if v == "pass":
                style.add("TEXTCOLOR", (oq_idx, r_idx), (oq_idx, r_idx), layout.pass_color)
            elif v == "fail":
                style.add("TEXTCOLOR", (oq_idx, r_idx), (oq_idx, r_idx), layout.fail_color)

    t.setStyle(style)
    return t


# ──────────────────────────────────────────────────────────────────────────────
# Document assembly
# ──────────────────────────────────────────────────────────────────────────────

def _render_pdf_common(
    df: pd.DataFrame,
    pdf_path: str,
    run_metadata: Optional[List[str]],
    legend: Optional[List[str]],
    title: Optional[str],
    subtitle: Optional[str],
    tools_header: Optional[Sequence[Tuple[str, str]]],
    layout: LayoutConfig,
    run_name: Optional[str],  # used in footer
) -> str:
    styles = _styles(layout)

    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=layout.page_size,
        rightMargin=layout.right_margin,
        leftMargin=layout.left_margin,
        topMargin=layout.top_margin,
        bottomMargin=layout.bottom_margin,
    )

    story: List[object] = []

    # Title / subtitle
    if title:
        story.append(Paragraph(title, styles["TitleCenter"]))
    if subtitle:
        story.append(Paragraph(subtitle, styles["HeaderBlock"]))
    if title or subtitle:
        story.append(Spacer(1, 0.2 * inch))

    # Dual header band
    band = _build_dual_header_band(run_metadata, tools_header, layout, styles)
    if band:
        story.append(band)
        story.append(Spacer(1, 0.15 * inch))

    # Legend (optional)
    if legend:
        story.append(Paragraph("Legend", styles["Header"]))
        for l in legend:
            story.append(Paragraph("• " + str(l), styles["HeaderBlock"]))
        story.append(Spacer(1, 0.2 * inch))

    # Main QC table
    story.append(Paragraph("Quality Summary Table", styles["Header"]))
    story.append(_build_table(df, layout))

    # Footer with run name + page number
    footer_func = _page_footer_factory(run_name)

    doc.build(story, onFirstPage=footer_func, onLaterPages=footer_func)
    return pdf_path