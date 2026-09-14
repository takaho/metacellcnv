#!/usr/bin/env python3
"""Generate markers/*.csv (dog, human, mouse) from one shared panel definition,
so the three species tables stay structurally consistent.

Provenance and validation status
---------------------------------
- The panel structure was built against Case1 of Dog10K_Boxer_Tasha in this
  project. Macrophage / Endothelial / Fibroblast / Epithelial were confirmed
  specific by depth-matched testing (validation=depth_matched_case1_dog).
- All other panels, and the human/mouse tables, are produced by nomenclature
  conversion only and have not been validated on this project's data
  (validation=orthology_conversion); the loader warns if these are used for
  labeling.
- Genes with no 1:1 mouse ortholog are dropped: GNLY, MNDA. For LYZ we use
  Lyz2 (the myeloid-predominant paralog of Lyz1/Lyz2); for FCGR3A we use
  Fcgr3.
"""
from __future__ import annotations
import csv
from pathlib import Path

HERE = Path(__file__).parent

# cell_type -> (role, use_as_normal_reference, validated_in_dog, [genes])
#   role: primary (normal calling) / fallback (used only if all primaries are absent)
#   use_as_normal_reference: whether this type may serve as the CNV normal reference
#     Fibroblast / Epithelial can themselves be tumor, so they are excluded by default
PANELS: dict[str, tuple[str, bool, bool, list[str]]] = {
    "T/NK": ("primary", True, False, [
        "CD3D", "CD3E", "CD3G", "CD2", "LCK", "ITK", "THEMIS", "SKAP1", "CD247",
        "IL7R", "CD28", "GZMA", "GZMB", "NKG7", "GNLY", "KLRD1", "PRF1", "EOMES"]),
    "B": ("primary", True, False, [
        "CD79A", "CD79B", "MS4A1", "CD19", "PAX5", "BANK1", "EBF1", "BLK",
        "CD22", "FCRL1"]),
    "Plasma": ("primary", True, False, [
        "JCHAIN", "MZB1", "XBP1", "DERL3", "PRDM1", "TNFRSF17", "SDC1"]),
    "Myeloid": ("primary", True, False, [
        "CD68", "CSF1R", "MRC1", "CD163", "C1QA", "C1QB", "C1QC", "MSR1",
        "AIF1", "TYROBP", "LYZ", "CD14", "ITGAM", "FCGR3A", "MARCO", "VSIG4",
        "FCER1G", "MNDA", "IRF8", "ZBTB46", "BATF3", "FLT3", "S100A8", "S100A9",
        "MMP9", "CSF3R"]),
    "Macrophage": ("primary", True, True, [
        "CD68", "AIF1", "TYROBP", "CD14", "FCER1G", "LYZ",
        "C1QA", "C1QB", "C1QC", "CTSS"]),
    "Mast": ("primary", True, False, [
        "KIT", "CPA3", "MS4A2", "TPSAB1", "CMA1", "GATA2", "HDC"]),
    "Endothelial": ("primary", True, True, [
        "PECAM1", "CDH5", "VWF", "KDR", "CLDN5", "EGFL7", "TEK", "ERG",
        "FLT1", "ESAM", "RAMP2", "PLVAP"]),
    "Fibroblast": ("primary", False, True, [
        "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "DCN", "LUM",
        "FBN1", "POSTN", "SPARC", "THY1", "FAP", "PDGFRB"]),
    "SmoothMuscle": ("primary", True, False, [
        "ACTA2", "TAGLN", "MYH11", "RGS5", "NOTCH3", "CNN1", "DES"]),
    "Epithelial": ("primary", False, True, [
        "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "SFN", "KRT5",
        "KRT14", "CLDN4"]),
    "Leukocyte": ("fallback", True, False, [
        "PTPRC", "LAPTM5", "CORO1A", "CD52", "SRGN"]),
}

# Genes with no 1:1 mouse ortholog (dropped)
MOUSE_ABSENT = {"GNLY", "MNDA"}
# Mouse symbols that aren't a simple Title case of the human symbol
MOUSE_SPECIAL = {"LYZ": "Lyz2", "FCGR3A": "Fcgr3"}


def to_mouse(sym: str) -> str | None:
    if sym in MOUSE_ABSENT:
        return None
    if sym in MOUSE_SPECIAL:
        return MOUSE_SPECIAL[sym]
    return sym[0].upper() + sym[1:].lower()


FIELDS = ["cell_type", "gene", "role", "use_as_normal_reference",
          "use_for_labels", "validation"]


def write(species: str, conv) -> None:
    rows = []
    for ct, (role, as_ref, validated, genes) in PANELS.items():
        for g in genes:
            sym = conv(g)
            if sym is None:
                continue
            rows.append({
                "cell_type": ct, "gene": sym, "role": role,
                "use_as_normal_reference": "TRUE" if as_ref else "FALSE",
                "use_for_labels": "TRUE" if validated else "FALSE",
                "validation": ("depth_matched_case1_dog"
                               if (validated and species == "dog")
                               else "orthology_conversion"),
            })
    path = HERE / f"markers_{species}.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"{path.name}: {len(rows)} rows / {len({r['cell_type'] for r in rows})} types")


if __name__ == "__main__":
    write("dog", lambda s: s)     # CellRanger's dog reference uses HGNC-like symbols
    write("human", lambda s: s)
    write("mouse", to_mouse)
