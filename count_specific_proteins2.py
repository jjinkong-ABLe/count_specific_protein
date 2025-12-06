#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OrthoFinder 결과에서 종별 Shared / Specific 단백질 수 계산

입력:
  1) Orthogroups.tsv (또는 Orthogroups_excel.xlsx)
  2) Orthogroups_UnassignedGenes.tsv (또는 엑셀 변환본)

출력:
  - shared_specific_protein_counts_per_species.tsv

정의:
  - Proteins_in_OGs:
      Orthogroups.tsv에 배정된 해당 종 유전자(단백질) 총합
  - Shared_OGs:
      해당 종 유전자가 포함되어 있고, 2종 이상이 함께 포함된 orthogroup 수
  - Shared_proteins:
      위 Shared_OGs에 속한 해당 종 유전자 수 합
  - Specific_OGs:
      해당 종 유전자만 포함된 orthogroup 수 (다른 종 컬럼은 모두 비어 있음)
  - Specific_proteins_in_OGs:
      종 특이 orthogroup에 속한 해당 종 유전자 수 합
  - Unassigned:
      Orthogroups_UnassignedGenes.tsv에서 해당 종에 속한 유전자 수
  - Total_proteins:
      Proteins_in_OGs + Unassigned
  - Specific_proteins_total:
      Specific_proteins_in_OGs + Unassigned
"""

import os
import re
import argparse
import pandas as pd

def read_table(path: str) -> pd.DataFrame:
    """tsv/xlsx 자동 판별 로더"""
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path)
    # 기본 tsv
    return pd.read_csv(path, sep="\t")

def count_genes(cell) -> int:
    """Orthogroups 셀의 유전자 개수 계산"""
    if pd.isna(cell):
        return 0
    s = str(cell).strip()
    if not s:
        return 0
    # OrthoFinder는 보통 "gene1, gene2, gene3" 형식
    parts = [p.strip() for p in s.split(",")]
    parts = [p for p in parts if p]
    return len(parts)

def nonempty_count(series: pd.Series) -> int:
    """Unassigned 파일에서 비어있지 않은 셀 개수"""
    # NaN 제외 + 공백 문자열 제외
    s = series.dropna().astype(str).str.strip()
    return (s != "").sum()

def infer_species_columns(df: pd.DataFrame) -> list:
    """
    Orthogroups 파일에서 종 컬럼 자동 추정
    - 첫 컬럼이 'Orthogroup'이거나 유사 이름이면 제외
    """
    cols = list(df.columns)
    if not cols:
        return []

    first = cols[0].lower()
    if "orthogroup" in first:
        return cols[1:]
    # 혹시 첫 컬럼명이 없거나 이상한 경우를 대비해
    # 'Orthogroup'라는 컬럼이 존재하면 그걸 기준으로
    if "Orthogroup" in cols:
        idx = cols.index("Orthogroup")
        return [c for i, c in enumerate(cols) if i != idx]
    # fallback: 전체 컬럼 중 문자열 리스트가 들어가는 컬럼을 종으로 가정
    # (대부분은 첫 컬럼만 ID이므로 이렇게 써도 안전)
    return cols[1:]

def main():
    parser = argparse.ArgumentParser(description="OrthoFinder Shared/Specific protein counter")
    parser.add_argument(
        "--orthogroups",
        default="Orthogroups.tsv",
        help="Orthogroups.tsv 또는 Orthogroups_excel.xlsx 경로"
    )
    parser.add_argument(
        "--unassigned",
        default="Orthogroups_UnassignedGenes.tsv",
        help="Orthogroups_UnassignedGenes.tsv 또는 엑셀 변환본 경로"
    )
    parser.add_argument(
        "--out",
        default="shared_specific_protein_counts_per_species.tsv",
        help="출력 TSV 경로"
    )
    args = parser.parse_args()

    # --- Load ---
    og = read_table(args.orthogroups)
    un = read_table(args.unassigned)

    # --- Species columns ---
    species_cols = infer_species_columns(og)
    if not species_cols:
        raise ValueError("Orthogroups 파일에서 종 컬럼을 찾지 못했습니다.")

    # Unassigned에도 동일한 종 컬럼이 있는지 교집합으로 정리
    un_species_cols = [c for c in species_cols if c in un.columns]
    if not un_species_cols:
        # unassigned 파일 컬럼명이 다를 수 있으니,
        # 공통 매칭이 안 되면 unassigned의 전체 컬럼을 종으로 사용
        # (단, 첫 컬럼이 의미 없는 경우를 대비해 Orthogroup 같은 이름 제외)
        un_species_cols = [c for c in un.columns if "orthogroup" not in str(c).lower()]

    # --- Count matrix for Orthogroups ---
    # 각 셀에 포함된 유전자 수로 변환
    counts = og[species_cols].applymap(count_genes)

    # 해당 종이 그 OG에 존재하는지
    present = counts > 0
    num_species_present = present.sum(axis=1)

    # --- Unassigned counts per species ---
    unassigned_counts = {}
    for sp in species_cols:
        if sp in un.columns:
            unassigned_counts[sp] = nonempty_count(un[sp])
        else:
            # 컬럼명이 없으면 0 처리
            unassigned_counts[sp] = 0

    # --- Aggregate per species ---
    rows = []
    for sp in species_cols:
        sp_counts = counts[sp]

        proteins_in_ogs = int(sp_counts.sum())

        # shared: 2종 이상 포함된 OG에서의 해당 종 단백질
        shared_mask = (num_species_present >= 2) & (sp_counts > 0)
        shared_ogs = int(shared_mask.sum())
        shared_proteins = int(sp_counts[shared_mask].sum())

        # specific OG: 오직 1종만 존재하는 OG에서 해당 종이 가진 경우
        specific_mask = (num_species_present == 1) & (sp_counts > 0)
        specific_ogs = int(specific_mask.sum())
        specific_proteins_in_ogs = int(sp_counts[specific_mask].sum())

        unassigned_n = int(unassigned_counts.get(sp, 0))

        total_proteins = proteins_in_ogs + unassigned_n
        specific_total = specific_proteins_in_ogs + unassigned_n

        rows.append({
            "Species": sp,
            "Proteins_in_OGs": proteins_in_ogs,
            "Unassigned": unassigned_n,
            "Total_proteins": total_proteins,

            "Shared_OGs": shared_ogs,
            "Shared_proteins": shared_proteins,

            "Specific_OGs": specific_ogs,
            "Specific_proteins_in_OGs": specific_proteins_in_ogs,
            "Specific_proteins_total": specific_total
        })

    out_df = pd.DataFrame(rows)

    # 보기 좋게 정렬(원하면 대문자/알파벳 순)
    out_df = out_df.sort_values("Species").reset_index(drop=True)

    # --- Save TSV ---
    out_df.to_csv(args.out, sep="\t", index=False, encoding="utf-8")

    print(f"[OK] Saved: {args.out}")
    print(out_df)

if __name__ == "__main__":
    main()

