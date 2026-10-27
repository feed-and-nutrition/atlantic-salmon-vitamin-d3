#!/usr/bin/env python3
"""Write TSV files as an .xlsx with one sheet per file and a bold header row. No dependencies.

    python3 scripts/tsv_to_xlsx.py table.tsv "Table S8 ....xlsx"
    python3 scripts/tsv_to_xlsx.py VDR.tsv GC.tsv "Table S7 ....xlsx"   # sheets named VDR and GC

A single file gives one sheet named Sheet1. With several, each sheet is named after its file (without .tsv).

Numbers are written as numbers, everything else as text. Column widths follow the longest value. A row of empty
cells is kept as a blank row, and a row with only its first cell filled (a table note) does not widen the column.
"""
import html
import os
import re
import sys
import zipfile

NUM = re.compile(r'-?\d+(\.\d+)?([eE][-+]?\d+)?$')


def col_name(i):
    s = ''
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def cell(ref, value, style):
    if NUM.match(value):
        return f'<c r="{ref}"{style}><v>{value}</v></c>'
    return f'<c r="{ref}" t="inlineStr"{style}><is><t xml:space="preserve">{html.escape(value)}</t></is></c>'


def sheet_xml(src):
    rows = [line.rstrip('\r\n').split('\t') for line in open(src, encoding='utf-8') if line.rstrip('\r\n')]
    sized = [r for r in rows if any(r[1:])]
    widths = [min(60, max(len(r[i]) for r in sized if i < len(r)) + 2) for i in range(len(rows[0]))]
    cols = ''.join(f'<col min="{i + 1}" max="{i + 1}" width="{w}" customWidth="1"/>' for i, w in enumerate(widths))
    data = ''.join(
        f'<row r="{r + 1}">' + ''.join(cell(f'{col_name(c)}{r + 1}', v, ' s="1"' if r == 0 else '')
                                       for c, v in enumerate(row)) + '</row>'
        for r, row in enumerate(rows))
    return cols, data


def main(srcs, dst):
    names = ['Sheet1'] if len(srcs) == 1 else [os.path.basename(f).removesuffix('.tsv') for f in srcs]
    sheets = [sheet_xml(f) for f in srcs]
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
        'xl/styles.xml':
            f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><styleSheet xmlns="{ns}">'
            '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
            '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
            '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'
            '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
            '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
            '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>'
            '</styleSheet>',
    }
    for i, (cols, data) in enumerate(sheets):
        parts[f'xl/worksheets/sheet{i + 1}.xml'] = (
            f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="{ns}">'
            f'<cols>{cols}</cols><sheetData>{data}</sheetData></worksheet>')
    with zipfile.ZipFile(dst, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, xml in parts.items():
            z.writestr(name, xml)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1:-1], sys.argv[-1])
