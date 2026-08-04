#!/usr/bin/env python3
"""One-time documentation migration for the completed TamI/TxtC pilot."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def update_todo() -> None:
    path = ROOT / "todo.md"
    text = path.read_text(encoding="utf-8")
    replacements = {
        "- [ ] 最初の解析対象を決める。TamI のみで始めるか、TamI と TxtC を並行して進めるかを選ぶ。":
            "- [x] 最初の解析対象を TamI と TxtC の並行 pilot に決定する。",
        "- [ ] homolog search strategy を定義する: database、E-value threshold、identity range、coverage cutoff、maximum sequence count。":
            "- [x] homolog search strategy を定義する: local Swiss-Prot/Bacteria、E-value 1e-20、identity 25–95%、coverage 70%以上、最大500件。",
        "- [ ] full-length MSA と catalytic-domain MSA を作成する。":
            "- [x] full-length MSA と tree-building 用 smart-gap trimmed MSA を作成する。\n- [ ] catalytic-domain MSA を追加する。",
        "- [ ] MSA quality metrics を記録する: sequence count、Neff、gap fraction、coverage、conserved P450 motif alignment。":
            "- [x] pilot MSA quality metrics を記録する: sequence count、gap fraction、coverage、conserved P450 motif alignment、reference residue mapping。\n- [ ] 最終共進化用セットを拡張した後に Neff を算出する。",
        "- [ ] すべての解析コマンドと parameter を、再実行可能な script または notebook に保存する。":
            "- [x] pilot 解析コマンドと parameter を、再実行可能な script、config、IQ-TREE log に保存する。",
    }
    for old, new in replacements.items():
        if old not in text and new not in text:
            raise RuntimeError(f"Expected todo item not found: {old}")
        text = text.replace(old, new)
    text = text.replace("最終更新日: 2026-07-09", "最終更新日: 2026-08-04", 1)
    path.write_text(text, encoding="utf-8", newline="\n")


def update_readme() -> None:
    path = ROOT / "README.md"
    text = path.read_text(encoding="utf-8")
    marker = "## TamI / TxtC pilot MSA・系統樹 (2026-08-04)"
    if marker not in text:
        text += """

## TamI / TxtC pilot MSA・系統樹 (2026-08-04)

TamI と TxtC について、reference 確認、local Swiss-Prot/Bacteria homolog 検索、filter、MAFFT full-length MSA、残基・P450 motif QC、smart-gap tree MSA、IQ-TREE ModelFinder + SH-aLRT 1000 + UFBoot 1000 を実行した。結果と制約は `reports/TamI_TxtC_MSA_tree_workflow_v0.1.md`、共通 parameter は `config/p450_phylo_pilot_v0.1.json`、集計値は各 system の `metadata/msa/<System>/msa_tree_summary.csv` に記録している。

TamI は reference を含む 69 配列、TxtC は 44 配列である。TxtC は手順書の pilot 推奨 50 配列を下回るため、最終共進化解析へ進む前により広いデータベースで homolog set を拡張する。
"""
    text = text.replace("最終更新日: 2026-07-09", "最終更新日: 2026-08-04", 1)
    path.write_text(text, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    update_todo()
    update_readme()
