# BEEP

BEEP (BGC-Heme Enzyme Evolution Project) は、生合成遺伝子クラスター (BGC) に含まれる heme enzyme、特に細菌および真菌由来 cytochrome P450 を対象とした、標準化データ整理および共進化解析プロジェクトである。

本プロジェクトの目的は、BGC レベルの代謝経路情報、P450 配列の共進化シグナル、三次元構造情報、文献証拠を統合し、機能決定残基 Functional Determinant Residues (FDetR) を推定するための再利用可能な解析基盤を構築することである。

最終更新日: 2026-08-04

## 現在の範囲

現段階では、benchmark データセットの標準化を優先する。まず、出典、ファイル、ハッシュ、配列対応、構造対応を追跡できるローカルデータパッケージを作成し、その後に MSA、MIp/DCA/EVcouplings、構造マッピング、候補残基ランキングへ進む。

現在の core dataset v0.1 統合版に入っている benchmark 系は以下の通り。

| Dataset | BGC / Product | Main P450 | Status | Notes |
| --- | --- | --- | --- | --- |
| `tami_project` | tirandamycin BGC / tirandamycin | TamI | v0.1 complete | 第一 benchmark。NCBI と UniProt の配列は一致し、Leu101, Leu244, Leu295 を既知参照残基として記録済み。 |
| `txtc_project` | thaxtomin BGC / thaxtomin A | TxtC | v0.1 complete | 第二 benchmark。BGC/NCBI 配列と UniProt/PDB 構造配列の由来が異なるため、配列差分表を作成済み。 |
| `metadata/per_system` | rosamicin BGC / rosamicin | RosC | integrated in master v0.1 | FG-loop と変異体情報を含む loop dynamics benchmark。P107S/L176Q、Leu176 FG loop、PDB 9VGM/9VGN/9VGQ/9VGP/9VGO を記録済み。 |
| `metadata/per_system` | mycophenolic acid BGC / mycophenolic acid | MpaDE / MpaD | integrated in master v0.1 | 真菌由来 P450-hydrolase fusion benchmark。MpaDE fusion と MpaD/MpaE domain cross-reference を一行で記録する。 |

今後の候補系として、さらに thalianol、avenacin、ginkgolide などの比較進化系を想定している。

統合版 master CSV は `metadata/core_bgc_p450_dataset_v0.1_master.csv` に置く。体系別の normalized row は `metadata/per_system/` に保存し、既存の `tami_project/` と `txtc_project/` の原始データパッケージは移動せず保持する。

## リポジトリ構成

```text
BEEP/
  README.md
  todo.md
  concept.md
  metadata/
    core_bgc_p450_dataset_v0.1_master.csv
    data_collection_log_master_v0.1.csv
    validation_checklist_master_v0.1.md
    per_system/
  tami_project/
    data/
      raw/
      processed/
    metadata/
    literature_notes/
    manuals/
    reports/
  txtc_project/
    data/
      raw/
      processed/
    metadata/
    literature_notes/
    manuals/
    reports/
  backup/
```

各 `*_project` ディレクトリは、できるだけ同じ構造で管理する。

- `data/raw/`: MIBiG、NCBI、UniProt、PDB/mmCIF などの未加工ダウンロードファイル。
- `data/processed/`: 後続処理結果、または legacy reference file。
- `metadata/`: コアデータ表、ダウンロードログ、検証チェックリスト、配列対応表。
- `literature_notes/`: 文献証拠、反応機構、benchmark residue の整理。
- `manuals/`: データ整理の操作手順書。
- `reports/`: 進捗報告、会議資料、解析レポート。

## データセットの状態

### TamI

`tami_project` は v0.1 データパッケージとして整備済み。

- MIBiG: BGC0001052.5
- Product: tirandamycin
- Organism: Streptomyces sp. 307-9
- Protein IDs: NCBI ADC79647.1; UniProt D3Y1J3
- Structures: PDB 6XA2, 6XA3
- Benchmark residues: Leu101, Leu244, Leu295

TamI は、residue-level coevolution が既知の機能部位や基質ポケット周辺を回収できるかを確認するための第一 benchmark として用いる。

### TxtC

`txtc_project` は v0.1 データパッケージとして整備済み。

- MIBiG: BGC0000444.5
- Product: thaxtomin A
- BGC organism: Streptomyces turgidiscabies Car8
- BGC protein: NCBI AAW49320.1
- Structure reference: UniProt Q8VS75, Streptomyces acidiscabies 84-104
- Structures: PDB 6F0B, 6F0C
- Candidate mechanism residues: Q8VS75/PDB 番号系における Ser280, Thr385

重要な注意点: AAW49320.1 と Q8VS75 はどちらも 395 aa だが、40 箇所でアミノ酸が異なる。Ser280 は Q8VS75/PDB 番号系では Ser だが、AAW49320.1 の同じ位置では Ala である。残基レベルで解釈する前に、`txtc_project/metadata/TxtC_AAW49320.1_vs_Q8VS75_sequence_mismatches_v0.1.csv` を参照する。

### RosC

RosC は `metadata/core_bgc_p450_dataset_v0.1_master.csv` の v0.1 統合行として整備済み。

- MIBiG: BGC0002086.3
- Product: rosamicin
- Main literature source: Micromonospora rosaria IFO 13697
- Protein IDs: RosC BAM35932.1; RosD BAM35933.1
- Structures: PDB 9VGM, 9VGN, 9VGQ, 9VGP, 9VGO
- Loop/mutant evidence: Pro107, Leu176 in FG loop, P107S/L176Q mutant

注意点: MIBiG 側の taxonomy は Micromonospora carbonacea で、RosC/RosD の機能・変異体文献は Micromonospora rosaria に基づく。BGC source と RosC sequence source を同一視する前に、この source split を確認する。

### MpaDE_MpaD

MpaDE_MpaD は `metadata/core_bgc_p450_dataset_v0.1_master.csv` の v0.1 統合行として整備済み。

- MIBiG: BGC0000104.4
- Product: mycophenolic acid
- Organism: Penicillium brevicompactum IBT 23078
- Fusion reference: BK008023.1 / DAA35181.1
- Legacy domain references: mpaD ADY00131.1; mpaE ADY00132.1 from HQ731031.1
- Domain architecture: N-terminal MpaD/CYP631B5 P450 + C-terminal MpaE Zn-dependent hydrolase/lactone synthase

注意点: MpaDE_MpaD は一つの比較体系として記録する。後続解析では full-fusion numbering と P450-domain-only numbering を分けて扱う必要がある。

## 解析方針

想定する解析パイプラインは以下の通り。

1. BGC/P450 benchmark データセットを整理する。
2. 配列出典、残基番号、構造ファイル、文献証拠を確認する。
3. 各 P450 について homolog set と MSA を構築する。
4. MIp/APC、DCA、EVcouplings を実行する。
5. 共進化スコアを PDB または AlphaFold 構造へマッピングする。
6. BGC context、構造位置、文献証拠、系統特異的置換を統合する。
7. FDetR 候補表および mutation-priority report を作成する。

## データ管理ルール

- ダウンロードした raw file は `data/raw/` に保存し、手動編集しない。
- 各 raw file について、source URL、collection date、file size、SHA256 hash、validation status を `metadata/data_collection_log_*` に記録する。
- 残基番号を benchmark 証拠として使う前に、配列出典の違いや番号対応を明記する。
- 候補残基は保守的に分類する。direct mutant evidence、structure-supported residue、mechanism candidate を区別する。
- 新しい benchmark 系を追加する場合、実際の schema gap が見つかるまでは TamI/TxtC の構造を再利用する。


## TamI / TxtC pilot MSA・系統樹 (2026-08-04)

TamI と TxtC について、reference 確認、local Swiss-Prot/Bacteria homolog 検索、filter、MAFFT full-length MSA、残基・P450 motif QC、smart-gap tree MSA、IQ-TREE ModelFinder + SH-aLRT 1000 + UFBoot 1000 を実行した。結果と制約は `reports/TamI_TxtC_MSA_tree_workflow_v0.1.md`、共通 parameter は `config/p450_phylo_pilot_v0.1.json`、集計値は各 system の `metadata/msa/<System>/msa_tree_summary.csv` に記録している。

TamI は reference を含む 69 配列、TxtC は 44 配列である。TxtC は手順書の pilot 推奨 50 配列を下回るため、最終共進化解析へ進む前により広いデータベースで homolog set を拡張する。
