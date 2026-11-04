#!/usr/bin/env python3
"""Write TSV files as an .xlsx with one sheet per file, formatted as a supplementary table. No dependencies.

    python3 scripts/tsv_to_xlsx.py table.tsv "Table S7.xlsx"
    python3 scripts/tsv_to_xlsx.py VDR.tsv GC.tsv "Table S6.xlsx"   # sheets named VDR and GC
    python3 scripts/tsv_to_xlsx.py --group "Genome search=Assembly,Species" ... "Table S9.xlsx"

A single file gives one sheet named Sheet1. With several, each sheet is named after its file (without .tsv).

Arial 10 throughout, a bold header row, and thin borders around every cell of the table. The table ends at the first
empty row or the first row with only its first cell filled (a table note), and the notes have no borders.
--group SHEET=COL,... (repeatable) names leading columns of a sheet that group the rows: a value is written once
for each run of rows that share it (and the values of the grouping columns to its left), with one border around the
run. For a single sheet, SHEET is Sheet1.

Numbers are written as numbers, everything else as text. Column widths follow the widest value (notes excluded),
with capitals counted wider than other characters.
A row of empty cells is kept as a blank row.
"""
import argparse
import html
import os
import re
import zipfile

NUM = re.compile(r'-?\d+(\.\d+)?([eE][-+]?\d+)?$')
# cellXfs: 0 plain (notes), 1 header, 2 body, 3 to 6 grouping cell (first, middle, last row of its run, single row)
STYLES = (
    '<fonts count="2"><font><sz val="10"/><name val="Arial"/><family val="2"/></font>'
    '<font><b/><sz val="10"/><name val="Arial"/><family val="2"/></font></fonts>'
    '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'
    '<borders count="5">' + ''.join(
        '<border>' + ''.join(f'<{s} style="thin"><color auto="1"/></{s}>' if s in sides else f'<{s}/>'
                             for s in ('left', 'right', 'top', 'bottom')) + '<diagonal/></border>'
        for sides in ('', 'left right top bottom', 'left right top', 'left right', 'left right bottom')) +
    '</borders>'
    '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
    '<cellXfs count="7"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
    '<xf numFmtId="0" fontId="1" fillId="0" borderId="1" xfId="0" applyFont="1" applyBorder="1"/>'
    + ''.join(f'<xf numFmtId="0" fontId="0" fillId="0" borderId="{b}" xfId="0" applyBorder="1"/>' for b in (1, 2, 3, 4, 1)) +
    '</cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>')


def text_width(s):  # in Arial digit widths: capitals count wider (accessions), and never less than 0.9 per character
    return max(0.9 * len(s), sum(1.2 if c.isupper() else 1 if c.isdigit() else 0.9 if c.islower() else 0.5 for c in s))


def col_name(i):
    s = ''
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def cell(ref, value, style):
    s = f' s="{style}"' if style else ''
    if value == '':
        return f'<c r="{ref}"{s}/>' if style else ''
    if NUM.match(value):
        return f'<c r="{ref}"{s}><v>{value}</v></c>'
    return f'<c r="{ref}" t="inlineStr"{s}><is><t xml:space="preserve">{html.escape(value)}</t></is></c>'


def sheet_xml(src, group):
    rows = [line.rstrip('\r\n').split('\t') for line in open(src, encoding='utf-8') if line.rstrip('\r\n')]
    ncol = len(rows[0])
    end = next((r for r in range(1, len(rows)) if not any(rows[r][1:])), len(rows))  # first blank or note row
    body = [row + [''] * (ncol - len(row)) for row in rows[1:end]]
    n_group = 0
    while n_group < ncol and rows[0][n_group] in group:
        n_group += 1
    if missing := set(group) - set(rows[0][:n_group]):
        raise SystemExit(f'{src}: grouping columns must be the leading columns, not found there: {sorted(missing)}')
    sized = [r for r in rows if any(r[1:])]
    widths = [min(130, round(max(text_width(r[i]) for r in sized if i < len(r))) + 3) for i in range(ncol)]
    cols = ''.join(f'<col min="{i + 1}" max="{i + 1}" width="{w}" customWidth="1"/>' for i, w in enumerate(widths))

    def key(r, c):  # a grouping cell continues the run above when it and the grouping cells to its left match
        return tuple(body[r][:c + 1]) if 0 <= r < len(body) else None

    out = []
    for r, row in enumerate(rows):
        if r == 0:
            cells = [cell(f'{col_name(c)}1', v, 1) for c, v in enumerate(row)]
        elif r <= len(body):
            b, cells = r - 1, []
            for c, v in enumerate(body[b]):
                style = 2
                if c < n_group:
                    first, last = key(b - 1, c) != key(b, c), key(b + 1, c) != key(b, c)
                    style = 6 if first and last else 3 if first else 5 if last else 4
                    v = v if first else ''
                cells.append(cell(f'{col_name(c)}{r + 1}', v, style))
        else:
            cells = [cell(f'{col_name(c)}{r + 1}', v, 0) for c, v in enumerate(row)]
        out.append(f'<row r="{r + 1}">' + ''.join(cells) + '</row>')
    return cols, ''.join(out)


def main(srcs, dst, groups):
    names = ['Sheet1'] if len(srcs) == 1 else [os.path.basename(f).removesuffix('.tsv') for f in srcs]
    if unknown := set(groups) - set(names):
        raise SystemExit(f'--group names unknown sheets: {sorted(unknown)}')
    sheets = [sheet_xml(f, groups.get(n, [])) for f, n in zip(srcs, names)]
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    rel = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
    parts = {
        '[Content_Types].xml':
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            + ''.join(f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" ContentType="application/vnd.openxmlformats-'
                      'officedocument.spreadsheetml.worksheet+xml"/>' for i in range(len(sheets))) +
            '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            '</Types>',
        '_rels/.rels':
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rId1" Type="{rel}/officeDocument" Target="xl/workbook.xml"/></Relationships>',
        'xl/workbook.xml':
            f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="{ns}" xmlns:r="{rel}">'
            '<sheets>' + ''.join(f'<sheet name="{html.escape(n)}" sheetId="{i + 1}" r:id="rId{i + 1}"/>'
                                 for i, n in enumerate(names)) + '</sheets></workbook>',
        'xl/_rels/workbook.xml.rels':
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            + ''.join(f'<Relationship Id="rId{i + 1}" Type="{rel}/worksheet" Target="worksheets/sheet{i + 1}.xml"/>'
                      for i in range(len(sheets))) +
            f'<Relationship Id="rId{len(sheets) + 1}" Type="{rel}/styles" Target="styles.xml"/></Relationships>',
        'xl/styles.xml': f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><styleSheet xmlns="{ns}">{STYLES}</styleSheet>',
    }
    for i, (cols, data) in enumerate(sheets):
        parts[f'xl/worksheets/sheet{i + 1}.xml'] = (
            f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="{ns}">'
            f'<cols>{cols}</cols><sheetData>{data}</sheetData></worksheet>')
    with zipfile.ZipFile(dst, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, xml in parts.items():
            z.writestr(name, xml)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--group', action='append', default=[], metavar='SHEET=COL,...')
    ap.add_argument('files', nargs='+', help='TSV files, then the .xlsx to write')
    a = ap.parse_args()
    if len(a.files) < 2:
        ap.error('give at least one TSV file and the .xlsx')
    groups = {}
    for g in a.group:
        sheet, _, columns = g.partition('=')
        groups[sheet] = columns.split(',')
    main(a.files[:-1], a.files[-1], groups)
