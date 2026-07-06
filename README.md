# BEEP

BEEP (BGC-Heme Enzyme Evolution Project) は、生合成遺伝子クラスター (BGC) に含まれる heme enzyme、特に細菌由来 cytochrome P450 を対象とした、標準化データ整理および共進化解析プロジェクトである。

本プロジェクトの目的は、BGC レベルの代謝経路情報、P450 配列の共進化シグナル、三次元構造情報、文献証拠を統合し、機能決定残基 Functional Determinant Residues (FDetR) を推定するための再利用可能な解析基盤を構築することである。

最終更新日: 2026-07-06

## 現在の範囲

現段階では、benchmark データセットの標準化を優先する。まず、出典、ファイル、ハッシュ、配列対応、構造対応を追跡できるローカルデータパッケージを作成し、その後に MSA、MIp/DCA/EVcouplings、構造マッピング、候補残基ランキングへ進む。

現在整備済みの benchmark 系は以下の通り。

| Dataset | BGC / Product | Main P450 | Status | Notes |
| --- | --- | --- | --- | --- |
| `tami_project` | tirandamycin BGC / tirandamycin | TamI | v0.1 complete | 第一 benchmark。NCBI と UniProt の配列は一致し、Leu101, Leu244, Leu295 を既知参照残基として記録済み。 |
| `txtc_project` | thaxtomin BGC / thaxtomin A | TxtC | v0.1 complete | 第二 benchmark。BGC/NCBI 配列と UniProt/PDB 構造配列の由来が異なるため、配列差分表を作成済み。 |

今後の候補系として、RosC/RosD、MpaD/MpaE/MpaDE、さらに thalianol、avenacin、ginkgolide などの比較進化系を想定している。

## リポジトリ構成

```text
BEEP/
  README.md
  todo.md
  concept.md
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