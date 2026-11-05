"""Unit tests for theiagene.report_variants (SnpEff ANN reporting)."""

import pytest

from theiagene import report_variants as rv
from theiagene.lib.parsers import assimilate_gff


_HEADER = "\n".join(
    [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1>",
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
        "'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | "
        "Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | "
        "CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' \">",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
        '##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Ref obs">',
        '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alt obs">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
    ]
)


def _ann(allele, annotation, feature_id="rna-x", hgvsc="", hgvsp="", messages="",
         feature_type="transcript"):
    """One SnpEff ANN entry, its unused subfields left empty"""
    fields = dict.fromkeys(rv._ANN_FIELDS, "")
    fields.update({
        "Allele": allele, "Annotation": annotation, "Feature_Type": feature_type,
        "Feature_ID": feature_id, "HGVS.c": hgvsc, "HGVS.p": hgvsp,
        "Messages": messages,
    })
    return "|".join(fields[name] for name in rv._ANN_FIELDS)


def _record(pos, ref, alt, anns, fmt="GT:AD", sample="1:0,562"):
    info = "ANN=" + ",".join(anns) if anns else "."
    return f"chr1\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t{info}\t{fmt}\t{sample}"


def _write_vcf(tmp_path, records, header=_HEADER, name="annotated.vcf"):
    path = tmp_path / name
    path.write_text("\n".join([header] + records) + "\n")
    return str(path)


# --------------------------------------------------------------------------- #
# iter_annotations
# --------------------------------------------------------------------------- #

def test_iter_annotations_yields_one_entry_per_ann_item(tmp_path):
    vcf = _write_vcf(tmp_path, [
        _record(100, "T", "C", [
            _ann("C", "missense_variant", "rna-x", "c.428A>G", "p.Lys143Arg"),
            _ann("C", "synonymous_variant", "rna-y", "c.12A>G", "p.Lys4Lys"),
        ]),
        # a record SnpEff left unannotated yields nothing
        _record(200, "G", "A", []),
    ])
    entries = [entry for _, entry in rv.iter_annotations(vcf)]
    assert [entry["Feature_ID"] for entry in entries] == ["rna-x", "rna-y"]
    assert entries[0]["HGVS.p"] == "p.Lys143Arg"


def test_iter_annotations_rejects_vcf_without_ann(tmp_path):
    header = "\n".join(
        line for line in _HEADER.splitlines() if "ID=ANN" not in line
    )
    vcf = _write_vcf(tmp_path, [_record(100, "T", "C", [])], header=header)
    with pytest.raises(ValueError, match="ANN"):
        list(rv.iter_annotations(vcf))


# --------------------------------------------------------------------------- #
# _depth_suffix
# --------------------------------------------------------------------------- #

def _first_record(tmp_path, record):
    vcf = _write_vcf(tmp_path, [record])
    return next(r for r, _ in rv.iter_annotations(vcf))


def test_depth_suffix_renders_ref_then_alt(tmp_path):
    record = _first_record(tmp_path, _record(100, "T", "C", [_ann("C", "x")]))
    assert rv._depth_suffix(record, "C") == "T:0 C:562"


def test_depth_suffix_selects_matching_alt_of_multiallelic(tmp_path):
    record = _first_record(
        tmp_path, _record(100, "A", "G,T", [_ann("T", "x")], sample="1:5,15,10")
    )
    assert rv._depth_suffix(record, "T") == "A:5 T:10"


def test_depth_suffix_falls_back_to_ro_ao_and_scrubs_scientific_gq(tmp_path):
    # no AD, so RO/AO supply the depths; the float GQ is scrubbed by import_vcf
    record = _first_record(
        tmp_path,
        _record(100, "GA", "GAA", [_ann("GAA", "x")], fmt="GT:GQ:RO:AO",
                sample="1:3.5e+01:2:23"),
    )
    assert rv._depth_suffix(record, "GAA") == "GA:2 GAA:23"


def test_depth_suffix_none_without_allele_depths(tmp_path):
    record = _first_record(
        tmp_path, _record(100, "T", "C", [_ann("C", "x")], fmt="GT", sample="1")
    )
    assert rv._depth_suffix(record, "C") is None


def test_depth_suffix_none_for_allele_absent_from_alts(tmp_path):
    record = _first_record(tmp_path, _record(100, "T", "C", [_ann("G", "x")]))
    assert rv._depth_suffix(record, "G") is None


# --------------------------------------------------------------------------- #
# report_line / _abbreviate
# --------------------------------------------------------------------------- #

def test_report_line_appends_depths_after_semicolon():
    line = rv.report_line(
        "ERG11", "lanosterol 14-alpha demethylase", "missense_variant",
        "c.428A>G", "p.Lys143Arg", "T:0 C:562",
    )
    assert line == (
        'ERG11: "lanosterol 14-alpha demethylase" '
        "(missense_variant c.428A>G p.Lys143Arg; T:0 C:562)"
    )


def test_report_line_omits_absent_pieces():
    line = rv.report_line("FKS1", "product", "splice_region_variant", "c.100+5G>A")
    assert line == 'FKS1: "product" (splice_region_variant c.100+5G>A)'


@pytest.mark.parametrize(
    "hgvs,expected",
    [
        ("p.Lys143Arg", "K143R"),
        ("p.Arg13545Ser", "R13545S"),
        ("p.Asp164Asp", "D164D"),
        ("p.Arg458*", "R458*"),
        ("p.Lys128fs", "K128fs"),
        ("p.Ter117Trpext*?", "*117Wext*?"),
        ("p.Phe432_Gly436del", "F432_G436del"),
        ("p.Ala292_Ser293insSerGly", "A292_S293insSG"),
        ("c.428A>G", "428A>G"),
        ("c.100+5G>T", "100+5G>T"),
        ("c.-14G>C", "-14G>C"),
        ("c.*32G>A", "*32G>A"),
        ("c.383delA", "383delA"),
        ("c.1222dupT", "1222dupT"),
        ("n.42A>T", "42A>T"),
    ],
)
def test_abbreviate(hgvs, expected):
    assert rv._abbreviate(hgvs) == expected


# --------------------------------------------------------------------------- #
# report_variants (end to end, through the GFF feature hierarchy)
# --------------------------------------------------------------------------- #

# gene -> mRNA -> CDS; the product's ',' is percent-encoded as GFF3 column 9
# requires
_GFF = "\n".join(
    [
        "##gff-version 3",
        "chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=gene-FKS1;Name=FKS1",
        "chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=rna-x;Parent=gene-FKS1",
        "chr1\t.\tCDS\t1\t1000\t.\t+\t0\tID=cds-x;Parent=rna-x;protein_id=prot-x;"
        "product=1%2C3-beta-glucan synthase component FKS1",
    ]
) + "\n"

_FKS1 = '"1,3-beta-glucan synthase component FKS1"'


def _report(tmp_path, records, gff=_GFF, suppress=(), qualifiers=("product",), **kwargs):
    gff_path = tmp_path / "reference.gff"
    gff_path.write_text(gff)
    vcf = _write_vcf(tmp_path, records)
    return rv.report_variants(
        vcf, assimilate_gff(str(gff_path)), set(suppress), "CDS", list(qualifiers), **kwargs
    )


_SYNONYMOUS = _record(100, "T", "C", [
    _ann("C", "synonymous_variant", "rna-x", "c.492C>T", "p.Asp164Asp"),
])


def test_report_variants_labels_lines_by_query_gene(tmp_path):
    records = _report(tmp_path, [_SYNONYMOUS], query_list=["FKS1"])
    assert records == [[
        "FKS1",
        "rna-x:c.492C>T",
        "prot-x:p.Asp164Asp",
        "492C>T",
        "D164D",
        f"FKS1: {_FKS1} (synonymous_variant c.492C>T p.Asp164Asp; T:0 C:562)",
    ]]


def test_report_variants_falls_back_to_product_label(tmp_path):
    # no query matches this feature, so the product names the line instead
    records = _report(tmp_path, [_SYNONYMOUS], query_list=["ERG11"])
    assert records[0][0] == "1.3-beta-glucan.synthase.component.FKS1"


def test_report_variants_exact_match_rejects_substring_query(tmp_path):
    # 'FKS1' is only a substring of the product/gene id, so --exact_match drops it
    records = _report(tmp_path, [_SYNONYMOUS], query_list=["FKS1"], exact_match=True)
    assert records[0][-1].startswith("1.3-beta-glucan.synthase.component.FKS1: ")


def test_report_variants_marks_missing_protein_change_na(tmp_path):
    (record,) = _report(tmp_path, [_record(200, "G", "A", [
        _ann("A", "splice_region_variant&intron_variant", "rna-x", "c.100+5G>A"),
    ])])
    assert record[1:5] == ["rna-x:c.100+5G>A", "NA", "100+5G>A", "NA"]
    assert record[5].endswith("(splice_region_variant&intron_variant c.100+5G>A; G:0 A:562)")


def test_report_variants_drops_non_transcript_annotations(tmp_path):
    # SnpEff gives intergenic regions an HGVS.c, so the feature type has to drop them
    assert _report(tmp_path, [_record(100, "T", "C", [
        _ann("C", "intergenic_region", "gene-a-gene-b", "n.100T>C",
             feature_type="intergenic_region"),
    ])]) == []


def test_report_variants_drops_annotation_without_hgvs(tmp_path):
    # SnpEff run with -noHgvs leaves both HGVS subfields empty
    assert _report(tmp_path, [_record(100, "T", "C", [
        _ann("C", "missense_variant", "rna-x"),
    ])]) == []


def test_report_variants_suppresses_on_any_consequence_term(tmp_path):
    records = [
        _record(100, "T", "C", [
            _ann("C", "missense_variant&splice_region_variant", "rna-x",
                 "c.428A>G", "p.Lys143Arg"),
        ]),
        _SYNONYMOUS.replace("\t100\t", "\t200\t"),
    ]
    reported = _report(tmp_path, records, suppress=["splice_region_variant"])
    assert [record[1] for record in reported] == ["rna-x:c.492C>T"]


@pytest.mark.parametrize(
    "warning",
    ["WARNING_TRANSCRIPT_NO_START_CODON", "WARNING_TRANSCRIPT_INCOMPLETE",
     "INFO_REALIGN_3_PRIME&WARNING_TRANSCRIPT_NO_START_CODON"],
)
def test_report_variants_treats_partial_transcript_as_noncoding(tmp_path, warning):
    (record,) = _report(tmp_path, [_record(100, "A", "C", [
        _ann("C", "missense_variant", "rna-x", "c.42A>C", "p.Lys14Asn", warning),
    ])])
    assert record[1:5] == ["rna-x:c.42A>C", "NA", "42A>C", "NA"]
    assert "p.Lys14Asn" not in record[5]


def test_report_variants_keeps_protein_change_under_other_messages(tmp_path):
    (record,) = _report(tmp_path, [_record(100, "AT", "A", [
        _ann("A", "frameshift_variant", "rna-x", "c.383delT", "p.Lys128fs",
             "INFO_REALIGN_3_PRIME"),
    ], sample="1:0,40")])
    assert record[2] == "prot-x:p.Lys128fs"


def test_run_cli_writes_headed_tsv(tmp_path):
    gff = tmp_path / "reference.gff"
    gff.write_text(_GFF)
    vcf = _write_vcf(tmp_path, [_SYNONYMOUS])
    output = tmp_path / "report.tsv"
    rv.main([
        "--vcf", vcf, "--reference_gff", str(gff),
        "--query_genes", "FKS1", "--output", str(output),
    ])
    header, record = output.read_text().splitlines()
    assert header == "#GENE\tHGVSc\tHGVSp\tNT\tAA\tREPORT"
    assert record.split("\t")[:5] == [
        "FKS1", "rna-x:c.492C>T", "prot-x:p.Asp164Asp", "492C>T", "D164D"
    ]


# --------------------------------------------------------------------------- #
# unit selection (SnpEff's 'Feature_ID' as a tiebreak, not as an ID lookup)
# --------------------------------------------------------------------------- #

# two transcripts overlapping one locus, each naming itself with a `Name` that
# differs from its `ID`, and neither CDS carrying a protein_id
_OVERLAP_GFF = "\n".join(
    [
        "##gff-version 3",
        "chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=gene-a;Name=FKS1",
        "chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=rna-a;Name=FKS1-T1;Parent=gene-a",
        "chr1\t.\tCDS\t1\t1000\t.\t+\t0\tID=cds-a;Parent=rna-a;"
        "product=1%2C3-beta-glucan synthase component FKS1",
        "chr1\t.\tgene\t1\t1000\t.\t-\t.\tID=gene-b;Name=ERG11",
        "chr1\t.\tmRNA\t1\t1000\t.\t-\t.\tID=rna-b;Name=ERG11-T1;Parent=gene-b",
        "chr1\t.\tCDS\t1\t1000\t.\t-\t0\tID=cds-b;Parent=rna-b;"
        "product=lanosterol 14-alpha demethylase",
    ]
) + "\n"


@pytest.fixture
def overlapping(tmp_path):
    """Report the one variant over the overlapping locus, under a given
    ``Feature_ID``, returning each record's HGVSp and REPORT columns."""

    def _run(feature_id, qualifiers=("product",), **kwargs):
        records = _report(tmp_path, [_record(100, "T", "C", [
            _ann("C", "missense_variant", feature_id, "c.428A>G", "p.Lys143Arg"),
        ])], gff=_OVERLAP_GFF, qualifiers=qualifiers, **kwargs)
        return [(record[2], record[-1]) for record in records]

    return _run


def test_report_variants_matches_feature_id_against_id(overlapping):
    # with no protein_id on the CDS, HGVS.p is prefixed with the transcript
    assert overlapping("rna-b") == [(
        "rna-b:p.Lys143Arg",
        'lanosterol.14-alpha.demethylase: "lanosterol 14-alpha demethylase" '
        "(missense_variant c.428A>G p.Lys143Arg; T:0 C:562)",
    )]


def test_report_variants_matches_feature_id_against_name_attribute(overlapping):
    # 'FKS1-T1' is the mRNA's Name and never its ID; the tiebreak still
    # discriminates rather than returning whichever unit happened to sort first
    assert overlapping("FKS1-T1") == [(
        "FKS1-T1:p.Lys143Arg",
        f"1.3-beta-glucan.synthase.component.FKS1: {_FKS1} "
        "(missense_variant c.428A>G p.Lys143Arg; T:0 C:562)",
    )]


def test_report_variants_drops_ambiguous_overlap_no_feature_answers_to(overlapping):
    # two units overlap and neither answers to the Feature_ID, so the entry
    # cannot be attributed to either and is dropped rather than guessed at
    assert overlapping("rna-c") == []


def test_report_variants_selected_unit_drives_the_query_label(overlapping):
    # the unit the Feature_ID picked is also the one the --query_genes term is
    # matched against
    (report,) = overlapping("ERG11-T1", qualifiers=("product", "Name"),
                            query_list=["ERG11"])
    assert report[1].startswith('ERG11: "lanosterol 14-alpha demethylase" ')


# Ensembl-style identifiers carry their own ':'; the CDS lacks a protein_id, so
# HGVS.p falls back to the (also colon-bearing) Feature_ID
_ENSEMBL_GFF = "\n".join(
    [
        "##gff-version 3",
        "chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=gene:GA;Name=GA",
        "chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=transcript:TA;Parent=gene:GA",
        "chr1\t.\tCDS\t1\t1000\t.\t+\t0\tID=CDS:TA;Parent=transcript:TA;product=protein A",
    ]
) + "\n"


def test_report_variants_quotes_colon_bearing_identifiers(tmp_path):
    (record,) = _report(tmp_path, [_record(100, "G", "A", [
        _ann("A", "stop_gained", "transcript:TA", "c.9G>A", "p.Trp3*"),
    ])], gff=_ENSEMBL_GFF)
    assert record[1:5] == ['"transcript:TA":c.9G>A', '"transcript:TA":p.Trp3*', "9G>A", "W3*"]
    # the REPORT line carries no identifier, so it is untouched
    assert "(stop_gained c.9G>A p.Trp3*;" in record[5]


@pytest.mark.parametrize(
    "identifier, expected",
    [("rna-x", "rna-x:c.9G>A"), ("transcript:TA", '"transcript:TA":c.9G>A')],
)
def test_prefix_hgvs_quotes_only_colon_bearing_identifiers(identifier, expected):
    assert rv._prefix_hgvs(identifier, "c.9G>A") == expected


# --------------------------------------------------------------------------- #
# unit selection: an ambiguous overlap without a Feature_ID
# --------------------------------------------------------------------------- #

def test_report_variants_drops_ambiguous_overlap_without_feature_id(overlapping):
    # two units overlap and the entry carries no Feature_ID to arbitrate, so it
    # is dropped rather than attributed to whichever unit sorted first
    assert overlapping("") == []


def test_report_variants_takes_sole_unit_without_feature_id(tmp_path):
    # a sole overlapping unit needs no tiebreak, Feature_ID or not
    (record,) = _report(tmp_path, [_record(100, "T", "C", [
        _ann("C", "synonymous_variant", "", "c.492C>T", "p.Asp164Asp"),
    ])])
    assert record[0] == "1.3-beta-glucan.synthase.component.FKS1"


# --------------------------------------------------------------------------- #
# unit fallbacks: gene -> CDS (no transcript level) and bare CDS
# --------------------------------------------------------------------------- #

# gene -> CDS, as e.g. many prokaryotic GFFs are laid out; SnpEff invents a
# transcript for each gene, so its Feature_ID matches no GFF identifier
_GENE_CDS_GFF = "\n".join(
    [
        "##gff-version 3",
        "chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=gene-a;Name=geneA",
        "chr1\t.\tCDS\t1\t1000\t.\t+\t0\tID=cds-a;Parent=gene-a;protein_id=prot-a;"
        "product=protein A",
        "chr1\t.\tgene\t2001\t3000\t.\t+\t.\tID=gene-b;Name=geneB",
        "chr1\t.\tCDS\t2001\t3000\t.\t+\t0\tID=cds-b;Parent=gene-b;protein_id=prot-b;"
        "product=protein B",
    ]
) + "\n"

# CDS records owning no gene or transcript at all
_BARE_CDS_GFF = "\n".join(
    [
        "##gff-version 3",
        "chr1\t.\tCDS\t1\t1000\t.\t+\t0\tID=cds-a;protein_id=prot-a;product=protein A",
        "chr1\t.\tCDS\t2001\t3000\t.\t+\t0\tID=cds-b;protein_id=prot-b;product=protein B",
    ]
) + "\n"


@pytest.mark.parametrize("gff", [_GENE_CDS_GFF, _BARE_CDS_GFF], ids=["gene-cds", "bare-cds"])
@pytest.mark.parametrize("pos, product, protein_id", [
    (100, "protein A", "prot-a"),
    (2100, "protein B", "prot-b"),
])
def test_report_variants_falls_back_below_the_transcript_level(
    tmp_path, gff, pos, product, protein_id
):
    # with no mRNA level the variant resolves through the gene (else the bare
    # CDS) it overlaps, even under a Feature_ID SnpEff made up itself
    (record,) = _report(tmp_path, [_record(pos, "T", "C", [
        _ann("C", "missense_variant", "TRANSCRIPT_made_up", "c.100T>C", "p.Phe34Leu"),
    ])], gff=gff)
    assert record[0] == product.replace(" ", ".")
    assert record[1:3] == ["TRANSCRIPT_made_up:c.100T>C", f"{protein_id}:p.Phe34Leu"]
    assert record[5].startswith(f'{product.replace(" ", ".")}: "{product}" ')


# --------------------------------------------------------------------------- #
# up-/downstream entries resolve against the units within the SnpEff buffer
# --------------------------------------------------------------------------- #

# geneA spans the variant; geneB sits ~1.7 kb downstream of it and does not
_NEIGHBOUR_GFF = "\n".join(
    [
        "##gff-version 3",
        "chr1\t.\tgene\t1\t300\t.\t+\t.\tID=gene-A;Name=GENEA",
        "chr1\t.\tmRNA\t1\t300\t.\t+\t.\tID=rna-A;Parent=gene-A",
        "chr1\t.\tCDS\t1\t300\t.\t+\t0\tID=cds-A;Parent=rna-A;protein_id=prot-A;"
        "product=protein A",
        "chr1\t.\tgene\t2001\t3000\t.\t+\t.\tID=gene-B;Name=GENEB",
        "chr1\t.\tmRNA\t2001\t3000\t.\t+\t.\tID=rna-B;Parent=gene-B",
        "chr1\t.\tCDS\t2001\t3000\t.\t+\t0\tID=cds-B;Parent=rna-B;protein_id=prot-B;"
        "product=protein B",
    ]
) + "\n"

_NEIGHBOUR_VARIANT = _record(9, "G", "A", [
    _ann("A", "stop_gained", "rna-A", "c.9G>A", "p.Trp3*"),
    _ann("A", "upstream_gene_variant", "rna-B", "c.-1992G>A"),
], sample="1:2,40")


def test_report_variants_attributes_upstream_entry_to_the_neighbouring_gene(tmp_path):
    # the upstream entry belongs to geneB, which the variant does not overlap;
    # it must not inherit geneA's label and product
    records = _report(tmp_path, [_NEIGHBOUR_VARIANT], gff=_NEIGHBOUR_GFF)
    assert [(record[0], record[1]) for record in records] == [
        ("protein.A", "rna-A:c.9G>A"),
        ("protein.B", "rna-B:c.-1992G>A"),
    ]


# --------------------------------------------------------------------------- #
# run_cli: BED-derived query labels and stdout output
# --------------------------------------------------------------------------- #

def _cli_inputs(tmp_path, records=(_SYNONYMOUS,), gff=_GFF):
    gff_path = tmp_path / "reference.gff"
    gff_path.write_text(gff)
    return ["--vcf", _write_vcf(tmp_path, list(records)), "--reference_gff", str(gff_path)]


def test_run_cli_labels_by_bedfile_names(tmp_path, capsys):
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t0\t1000\tFKS1\n")
    rv.main(_cli_inputs(tmp_path) + ["--bedfile", str(bed)])
    _, record = capsys.readouterr().out.splitlines()
    assert record.split("\t")[0] == "FKS1"


def test_run_cli_query_genes_take_precedence_over_bedfile(tmp_path, capsys):
    # the BED names a term that matches nothing; --query_genes wins regardless
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t0\t1000\tERG11\n")
    rv.main(_cli_inputs(tmp_path) + ["--bedfile", str(bed), "--query_genes", "FKS1"])
    _, record = capsys.readouterr().out.splitlines()
    assert record.split("\t")[0] == "FKS1"


def test_run_cli_prints_tsv_to_stdout_and_suppresses_up_downstream_by_default(
    tmp_path, capsys
):
    # no --output: the TSV goes to stdout, and with --suppress left at its
    # default the upstream entry for geneB is dropped
    rv.main(_cli_inputs(tmp_path, [_NEIGHBOUR_VARIANT], gff=_NEIGHBOUR_GFF))
    assert capsys.readouterr().out.splitlines() == [
        "#GENE\tHGVSc\tHGVSp\tNT\tAA\tREPORT",
        "\t".join([
            "protein.A", "rna-A:c.9G>A", "prot-A:p.Trp3*", "9G>A", "W3*",
            'protein.A: "protein A" (stop_gained c.9G>A p.Trp3*; G:2 A:40)',
        ]),
    ]


# --------------------------------------------------------------------------- #
# ANN subfield count
# --------------------------------------------------------------------------- #

def test_ann_header_fields_parses_snpeff_description():
    fields = rv.ann_header_fields(
        "Functional annotations: 'Allele | Annotation | Annotation_Impact' "
    )
    assert fields == ["Allele", "Annotation", "Annotation_Impact"]
    assert rv.ann_header_fields("Functional annotations") == []


def test_iter_annotations_rejects_entry_shifted_by_pipe_in_identifier(tmp_path):
    # SnpEff copies an NCBI gnl| transcript ID unescaped into Feature_ID, which
    # would otherwise read Transcript_BioType/Rank as HGVS.c/HGVS.p
    shifted = _ann("GC", "frameshift_variant", "rna-gnl|WGS:PGLS|mrna_CJI97_000402T0",
                   "c.10dupC", "p.Leu4fs")
    vcf = _write_vcf(tmp_path, [_record(100, "G", "GC", [shifted])])
    with pytest.raises(ValueError, match="18 subfields; the header declares 16"):
        list(rv.iter_annotations(vcf))


def test_iter_annotations_rejects_header_without_subfields(tmp_path):
    header = _HEADER.replace(
        _HEADER.splitlines()[2],
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">',
    )
    vcf = _write_vcf(tmp_path, [_record(100, "T", "C", [])], header=header)
    with pytest.raises(ValueError, match="declares no subfields"):
        list(rv.iter_annotations(vcf))


def test_iter_annotations_rejects_header_with_unexpected_subfield_count(tmp_path):
    header = _HEADER.replace(" | ERRORS / WARNINGS / INFO", "")
    vcf = _write_vcf(tmp_path, [_record(100, "T", "C", [])], header=header)
    with pytest.raises(ValueError, match="declares 15 subfields; 16 expected"):
        list(rv.iter_annotations(vcf))
