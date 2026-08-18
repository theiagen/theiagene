# theiagene.lib.feature

The feature data model shared by every theiagene command. It turns the flat,
line-oriented rows of a GFF3 annotation into a navigable tree of genes,
their RNAs, and the CDS/exon segments beneath them — then lets you reach any
class of feature by name, look one up by ID, sort the whole thing into a
canonical order, and serialize it back out to GFF3.

## Introduction

A GFF file is a flat list of rows. Each row has a `Parent` attribute pointing
at the `ID` of another row, so a gene, its mRNA, and that mRNA's several CDS
segments arrive as unrelated lines that only *imply* a hierarchy. This module
reifies that hierarchy.

There are two classes:

- **`Feature`** — one annotated interval (a gene, an mRNA, a CDS segment, …).
  It holds the parsed GFF columns, a link up to its `parent`, and a list of its
  `descendants`. It knows how to serialize itself (and its subtree) back to
  GFF3.
- **`FeatureCol`** — an ordered *collection* of `Feature`s. On construction it
  wires up the `parent`/`descendants` links from `Parent`↔`ID` relationships,
  then exposes the features both by **class** (`.genes`, `.rnas`, `.cds`,
  `.exons`, or `col["mRNA"]`) and by **ID** (`col.get("rna-1")`).

A typical flow: a parser builds a `list[Feature]`, hands it to `FeatureCol`,
and downstream code walks `gene → rna → cds` through `.descendants` or pulls a
named class straight off the collection.

```python
from theiagene.lib.feature import Feature, FeatureCol

# usually built by a parser (see theiagene.lib.parsers.assimilate_gff);
# shown here constructed by hand. ingest=True is what derives each feature's
# fid/pid from its ID/Parent attributes -- without it fid stays None and the
# Parent<->ID wiring below has nothing to match against
gene = Feature(seqid="chr1", source="example", type="gene", start=0, end=900,
               attributes={"ID": "gene-FKS1", "Name": "FKS1"}, ingest=True)
rna  = Feature(seqid="chr1", source="example", type="mRNA", start=0, end=900,
               attributes={"ID": "rna-1", "Parent": "gene-FKS1"}, ingest=True)
cds  = Feature(seqid="chr1", source="example", type="CDS", start=100, end=800,
               attributes={"ID": "cds-1", "Parent": "rna-1"}, ingest=True)

col = FeatureCol([gene, rna, cds])   # groups by Parent<->ID on construction

repr(col)                            # 'FeatureCol(1 genes, 1 RNAs, 1 CDS, 0 exons)'
col["mRNA"]                          # [rna-1]  (raw GFF type resolves to the 'rna' class)
col.rnas[0].descendants              # [cds-1]  (hierarchy wired up)
col.get("cds-1").parent.fid          # 'rna-1'
print(gene.to_gff())                 # serialize the gene + its subtree to GFF3
```

`Feature` defines no `__repr__`, so the lists above print as the default
`<theiagene.lib.feature.Feature object at 0x...>`; they are shown by `fid` for
readability. The `to_gff` call emits the subtree, converted back to GFF3's
1-based, both-inclusive coordinates:

```
chr1	example	gene	1	900	.	.	.	ID=gene-FKS1;Name=FKS1
chr1	example	mRNA	1	900	.	.	.	ID=rna-1;Parent=gene-FKS1
chr1	example	CDS	101	800	.	.	.	ID=cds-1;Parent=rna-1
```

Coordinates are stored **0-based, half-open** internally (BED-style), regardless
of the 1-based, both-inclusive convention GFF3 uses on disk; conversion happens
only at serialization time.

---

## API reference

### `class Feature`

A single cross-annotation feature — one interval on one contig.

```python
Feature(fid=None, pid=None, seqid=None, source=None, type=None,
        start=None, end=None, score=None, strand=None, phase=None,
        attributes=None, descendants=None, parent=None,
        sequence="", origin_sequence="", ingest=False)
```

**Parameters**

| name | type | description |
| --- | --- | --- |
| `fid` | `str` | Feature ID. If falsy and `ingest=True`, derived from the `ID`/`Id`/`id` attribute. A feature with no resolvable ID is left with `fid=None` rather than raising — GFF3 only requires `ID` on features that have children or span multiple lines, so callers that need an ID for every record (as `FeatureCol` grouping does) call `synthesize_id` to mint one. |
| `pid` | `str` | Parent ID. If falsy and `ingest=True`, derived from the `Parent`/`parent` attribute. |
| `seqid` | `str` | Contig/sequence name (GFF column 1). |
| `source` | `str` | Annotation source (GFF column 2). |
| `type` | `str` | Feature type such as `gene`, `mRNA`, `CDS`, `exon` (GFF column 3). |
| `start`, `end` | `int` | Interval bounds, **0-based half-open**. Required: both are coerced with `_as_int` and stored sorted so `start <= end` always holds, but an undefined placeholder (`.`/`?`/`""`/`None`) raises `ValueError` — a feature without coordinates cannot be built. `0` is a valid coordinate and is not treated as undefined. |
| `score` | `int` \| placeholder | GFF column 6; undefined → `None`. |
| `strand` | `+`/`-`/`1`/`-1`/placeholder | Coerced to a BioPython-style strand integer (`1`, `-1`) or `None`. |
| `phase` | `int` \| placeholder | GFF column 8; undefined → `None`. |
| `attributes` | `dict` \| `str` | Column-9 attributes. A dict is used as-is; a raw GFF3 column-9 string is parsed only when `ingest=True`. Defaults to a fresh `{}`. |
| `descendants` | `list[Feature]` | Child features. Defaults to a fresh `[]` (never a shared mutable default). |
| `parent` | `Feature` | Link to the parent feature. |
| `sequence` | `str` | Explicit nucleotide sequence; must be at least `end - start` long or `ValueError` is raised. |
| `origin_sequence` | `str` | Full contig sequence to slice `[start:end]` from, populating `sequence`. |
| `ingest` | `bool` | When `True`, run `_ingest()` to parse a string `attributes` and backfill `fid`/`pid` from them. |

**Attributes** — the constructor stores normalized `fid`, `pid`, `seqid`,
`source`, `type`, `start`, `end`, `score`, `strand`, `phase`, `attributes`
(dict), `descendants`, `parent`, and `sequence`.

**Methods**

- **`to_gff() -> str`**
  Serialize this feature *and all of its descendants* to a GFF3 string, emitted
  depth-first (each parent before its children), joined by newlines with no
  trailing newline. Coordinates are converted back to GFF3's 1-based,
  both-inclusive columns; undefined numeric/strand fields (`score`, `phase`,
  `strand`) render as `.`. The string columns (`seqid`, `source`, `type`) are
  written as-is and have no placeholder — serializing a feature that is missing
  one raises `TypeError`, so set them on any feature you intend to write back
  out.

- **`synthesize_id(count) -> str`**
  Assign and return a synthetic `fid` for a record that carried no `ID`, built
  as `{type}{count}` (lowercased type plus a per-type occurrence count) and
  suffixed with `_{pid}` when a parent is known, which keeps multi-segment rows
  sharing one parent distinct. The value is written back into
  `attributes["ID"]` so it survives a round-trip through `to_gff`.
  `theiagene.lib.parsers.iter_gff_features` calls this for every ID-less record
  it reads.

  > `_ingest()`, `_to_gff_line()` (single-line serialization, no descendants)
  > are internal helpers.

**Raises**

- `ValueError` — `start` or `end` is undefined; a defined
  `start`/`end`/`score`/`phase` is not an integer; an invalid strand token; or a
  `sequence` shorter than the coordinate span.

---

### `class FeatureCol`

An ordered collection of `Feature`s linked by `Parent`↔`ID` relationships.

```python
FeatureCol(features=None, group=True)
```

**Parameters**

| name | type | description |
| --- | --- | --- |
| `features` | `list[Feature]` | Features to collect. Copied into `self.features`. Defaults to empty. |
| `group` | `bool` | When `True` (default), run `group_features` to wire up `parent`/`descendants` in place. Pass `False` to preserve links already established (see `roots`). |

On construction the features are (optionally) grouped and then bucketed by
canonical class, and a `fid → Feature` index is built.

**Attributes**

| name | type | description |
| --- | --- | --- |
| `features` | `list[Feature]` | All features, in collection order. |
| `genes`, `rnas`, `cds`, `exons` | `list[Feature]` | Features of each canonical class. Any GFF type ending in `rna`/`transcript` (mRNA, tRNA, ncRNA, primary_transcript, …) buckets into `rnas`. |

**Methods**

- **`by_id(fid) -> Feature`** — return the feature whose `fid` matches;
  raises `KeyError` if absent. Only features in `self.features` are indexed
  (descendants reachable only via `.descendants` are not).
- **`get(fid, default=None) -> Feature | default`** — non-raising counterpart
  of `by_id`.
- **`sort() -> FeatureCol`** — order features by contig (`seqid`), then
  parent-before-descendant, then `start`; sorts every `descendants` list in
  place so this walk and `Feature.to_gff` share one sibling order. Undefined
  `seqid`/`start` sort to the end of their group. Returns `self`.
- **`roots() -> FeatureCol`** — a new `FeatureCol` of only the root features
  (no parent within this collection), built with `group=False` so the existing
  hierarchy stays reachable through `.descendants`.

**Item access** (`__getitem__`)

| key | result |
| --- | --- |
| `col["gene"]`, `col["rna"]`, `col["cds"]`, `col["exon"]` | that class's list |
| `col["genes"]` (plural) | same as the singular |
| `col["mRNA"]`, `col["tRNA"]`, … (raw GFF type) | resolves to the matching class list |
| `col[0]`, `col[2:5]` (int / slice) | index or slice into `self.features` |

Keys are case-insensitive; an unrecognized string raises `KeyError`, a
non-str/int/slice key raises `TypeError`.

**Other protocols** — `iter(col)` yields `self.features`; `len(col)` is the
feature count; `repr(col)` summarizes class counts, e.g.
`FeatureCol(3 genes, 3 RNAs, 8 CDS, 8 exons)`.

**Raises**

- `KeyError` — two features share an `fid` after grouping (surfaced rather than
  silently collapsed), or an unrecognized string key.

---

### `group_features(features, parent_ids=["Parent", "parent"], ids=["ID", "Id", "id"]) -> dict`

Hierarchically group a list of features by `Parent`↔`ID`, linking each feature
to its `parent` and appending it to its parent's `descendants` **in place**.
Duplicate IDs are de-replicated by appending a count (`id`, `id_1`, `id_2`, …)
so multi-segment CDS features index distinctly; a `Parent` pointing at a
duplicated (ambiguous) ID raises `KeyError`. Returns a `{seqid: [Feature, …]}`
dict.

`FeatureCol` calls this for you on construction — call it directly only when
working with raw feature lists outside a collection.

---

### Coordinate & attribute conventions

- **Coordinates** are 0-based, half-open in memory; `_to_gff_line` converts back
  to GFF3's 1-based, both-inclusive columns on output.
- **Strand** is stored as `1` / `-1` / `None`; `+` / `-` on the way in, `.` for
  `None` on the way out.
- **Undefined** GFF fields (`.`, `?`, `""`, `None`) normalize to `None` (or `{}`
  for attributes).
- **Attribute (column 9)** parsing/serialization percent-encodes the GFF3
  reserved characters (`% ; = & ,` tab/newline/CR); `%` is handled first so
  escapes are not double-encoded.
