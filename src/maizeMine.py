"""
MaizeMine client for retrieving gene sequences and Pfam domains.
Translated from Kotlin, with sequence/Pfam additions.
"""

import requests
from typing import Optional
from dataclasses import dataclass
import html


BASE_URL = "https://maizemine.rnet.missouri.edu/maizemine/service"


@dataclass
class GeneSummary:
    query: str
    primary_identifier: str
    symbol: Optional[str] = None
    description: Optional[str] = None


@dataclass
class PfamDomain:
    identifier: str
    name: Optional[str] = None
    description: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None


@dataclass
class GeneSequence:
    primary_identifier: str
    residues: str
    length: int


def _escape_xml(s: str) -> str:
    return html.escape(s, quote=True)


def _post_tab_query(xml_query: str, max_rows: Optional[int] = None) -> str:
    url = f"{BASE_URL.rstrip('/')}/query/results"
    data = {
        "format": "tab",
        "query": xml_query,
        "columnheaders": "1",
    }
    if max_rows is not None:
        data["size"] = str(max_rows)
    
    resp = requests.post(url, data=data)
    resp.raise_for_status()
    return resp.text


def _parse_tab(tab: str, expected_headers: list[str]) -> list[dict[str, str]]:
    lines = [l for l in tab.replace("\r\n", "\n").replace("\r", "\n").split("\n") if l.strip()]
    if not lines:
        return []
    
    first_cols = lines[0].split("\t")
    looks_like_header = any(">" in c for c in first_cols)
    
    if looks_like_header:
        headers = first_cols
        data_start = 1
    else:
        headers = expected_headers
        data_start = 0
    
    out = []
    for line in lines[data_start:]:
        cols = line.split("\t")
        if all(c.strip() == "" for c in cols):
            continue
        row = {}
        for j, h in enumerate(headers):
            if j < len(cols):
                val = cols[j].strip()
                if val.startswith('"') and val.endswith('"'):
                    val = val[1:-1]
                row[h] = val
        out.append(row)
    return out


def resolve_to_primary_identifier(gene_input: str) -> Optional[str]:
    """Resolve gene ID/symbol/alias to canonical primaryIdentifier."""
    q = gene_input.strip()
    if not q:
        return None
    
    # Try exact primaryIdentifier match
    summary = gene_summary_by_primary_identifier(q)
    if summary:
        return summary.primary_identifier
    
    # Fallback to LOOKUP
    summary = gene_summary_by_lookup(q)
    if summary:
        return summary.primary_identifier
    
    return None


def gene_summary_by_primary_identifier(primary_id: str) -> Optional[GeneSummary]:
    pid = primary_id.strip()
    if not pid:
        return None
    
    view = "Gene.primaryIdentifier Gene.symbol Gene.description"
    expected = ["Gene > Gene ID", "Gene > Symbol", "Gene > Description"]
    
    xml = f'''<query model="genomic" view="{view}">
      <constraint path="Gene.primaryIdentifier" op="=" value="{_escape_xml(pid)}"/>
    </query>'''
    
    tab = _post_tab_query(xml)
    rows = _parse_tab(tab, expected)
    if not rows:
        return None
    
    row = rows[0]
    gid = row.get("Gene > Gene ID", "").strip()
    if not gid:
        return None
    
    return GeneSummary(
        query=pid,
        primary_identifier=gid,
        symbol=row.get("Gene > Symbol") or None,
        description=row.get("Gene > Description") or None,
    )


def gene_summary_by_lookup(gene_input: str) -> Optional[GeneSummary]:
    q = gene_input.strip()
    if not q:
        return None
    
    view = "Gene.primaryIdentifier Gene.symbol Gene.description"
    expected = ["Gene > Gene ID", "Gene > Symbol", "Gene > Description"]
    
    xml = f'''<query model="genomic" view="{view}">
      <constraint path="Gene" op="LOOKUP" value="{_escape_xml(q)}"/>
    </query>'''
    
    tab = _post_tab_query(xml)
    rows = _parse_tab(tab, expected)
    if not rows:
        return None
    
    row = rows[0]
    gid = row.get("Gene > Gene ID", "").strip()
    if not gid:
        return None
    
    return GeneSummary(
        query=q,
        primary_identifier=gid,
        symbol=row.get("Gene > Symbol") or None,
        description=row.get("Gene > Description") or None,
    )


def get_gene_sequence(gene_input: str) -> Optional[GeneSequence]:
    """Get protein or transcript sequence for a gene."""
    resolved = resolve_to_primary_identifier(gene_input)
    if not resolved:
        return None
    
    # Try protein sequence first
    view = "Gene.primaryIdentifier Gene.proteins.sequence.residues Gene.proteins.sequence.length"
    expected = ["Gene > Gene ID", "Gene > Proteins > Sequence > Residues", "Gene > Proteins > Sequence > Length"]
    
    xml = f'''<query model="genomic" view="{view}">
      <constraint path="Gene.primaryIdentifier" op="=" value="{_escape_xml(resolved)}"/>
    </query>'''
    
    tab = _post_tab_query(xml)
    rows = _parse_tab(tab, expected)
    
    if rows:
        row = rows[0]
        residues = row.get("Gene > Proteins > Sequence > Residues", "").strip()
        length_str = row.get("Gene > Proteins > Sequence > Length", "").strip()
        if residues:
            length = int(length_str) if length_str.isdigit() else len(residues)
            return GeneSequence(primary_identifier=resolved, residues=residues, length=length)
    
    return None


def get_pfam_domains(gene_input: str, limit_rows: int = 500) -> list[PfamDomain]:
    """Get Pfam protein domains for a gene."""
    resolved = resolve_to_primary_identifier(gene_input)
    if not resolved:
        return []
    
    view = (
        "Gene.primaryIdentifier "
        "Gene.proteins.proteinDomainRegions.proteinDomain.primaryIdentifier "
        "Gene.proteins.proteinDomainRegions.proteinDomain.name "
        "Gene.proteins.proteinDomainRegions.proteinDomain.description "
        "Gene.proteins.proteinDomainRegions.start "
        "Gene.proteins.proteinDomainRegions.end"
    )
    expected = [
        "Gene > Gene ID",
        "Gene > Proteins > Protein Domain Regions > Protein Domain > DB Identifier",
        "Gene > Proteins > Protein Domain Regions > Protein Domain > Name",
        "Gene > Proteins > Protein Domain Regions > Protein Domain > Description",
        "Gene > Proteins > Protein Domain Regions > Start",
        "Gene > Proteins > Protein Domain Regions > End",
    ]
    
    xml = f'''<query model="genomic" view="{view}">
      <constraint path="Gene.primaryIdentifier" op="=" value="{_escape_xml(resolved)}"/>
    </query>'''
    
    tab = _post_tab_query(xml, max_rows=limit_rows)
    rows = _parse_tab(tab, expected)
    
    domains = []
    seen = set()
    for row in rows:
        pfam_id = row.get("Gene > Proteins > Protein Domain Regions > Protein Domain > DB Identifier", "").strip()
        if not pfam_id:
            continue
        
        start_str = row.get("Gene > Proteins > Protein Domain Regions > Start", "").strip()
        end_str = row.get("Gene > Proteins > Protein Domain Regions > End", "").strip()
        
        key = (pfam_id, start_str, end_str)
        if key in seen:
            continue
        seen.add(key)
        
        domains.append(PfamDomain(
            identifier=pfam_id,
            name=row.get("Gene > Proteins > Protein Domain Regions > Protein Domain > Name") or None,
            description=row.get("Gene > Proteins > Protein Domain Regions > Protein Domain > Description") or None,
            start=int(start_str) if start_str.isdigit() else None,
            end=int(end_str) if end_str.isdigit() else None,
        ))
    
    return domains


# --- Quick test / example usage ---
if __name__ == "__main__":
    test_gene = "Zm00001eb000010"
    
    print(f"Testing with gene: {test_gene}")
    
    # Test resolution
    resolved = resolve_to_primary_identifier(test_gene)
    assert resolved is not None, f"Failed to resolve {test_gene}"
    print(f"Resolved: {resolved}")
    
    # Test sequence
    seq = get_gene_sequence(test_gene)
    if seq:
        print(f"Sequence length: {seq.length}")
        print(f"Sequence preview: {seq.residues[:50]}...")
    else:
        print("No sequence found")
    
    # Test Pfam domains
    domains = get_pfam_domains(test_gene)
    print(f"Pfam domains found: {len(domains)}")
    for d in domains[:5]:
        print(f"  {d.identifier}: {d.name} ({d.start}-{d.end})")