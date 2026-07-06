# BEEP To-Do

最終更新日: 2026-07-06

## Phase 0: リポジトリ整理

- [x] 文字化けしていた root `README.md` を、読めるプロジェクト概要に置き換える。
- [x] 文字化けしていた root `todo.md` を、現在の進捗に合わせたタスクリストに置き換える。
- [x] root の主要 Markdown を日本語で統一する。
- [ ] `concept.md` を今後どの粒度で維持するか決める。短い研究概要にするか、詳細な研究計画書にするかを選ぶ。
- [ ] TamI/TxtC を横断比較する summary table が必要になった時点で、root に `metadata/` または `docs/` を追加する。

## Phase 1: Benchmark データセット基盤

- [x] `tami_project` を標準化データセット構造で作成する。
- [x] TamI の raw data を収集する: MIBiG、NCBI FASTA、UniProt FASTA、PDB/mmCIF。
- [x] TamI の metadata CSV、data collection log、literature note、validation checklist を作成する。
- [x] `txtc_project` を `tami_project` と同じ基本構造で作成する。
- [x] TxtC の raw data を収集する: MIBiG、NCBI FASTA、UniProt FASTA、PDB/mmCIF。
- [x] TxtC の metadata CSV、data collection log、literature note、validation checklist を作成する。
- [x] TxtC の AAW49320.1 と Q8VS75 の配列差分を記録する。
- [ ] TamI と TxtC の主要項目を比較する root-level benchmark summary table を作成する。
- [ ] RosC/RosD、MpaD/MpaE/MpaDE を追加する前に、最低限必要な metadata fields を確定する。

## Phase 2: データセット検証

- [x] TamI: NCBI と UniProt の配列一致を確認する。
- [x] TamI: Leu101、Leu244、Leu295 が 6XA2/6XA3 構造中で確認できることを検証する。
- [x] TxtC: MIBiG BGC0000444.5 の parse result と source status を確認する。
- [x] TxtC: AAW49320.1 と Q8VS75 の長さ、および mismatch positions を確認する。
- [x] TxtC: Ser280 と Thr385 が 6F0B/6F0C 構造中で確認できることを検証する。
- [ ] TxtC の下流解析で使う残基番号体系を決める。候補は Q8VS75/PDB numbering、AAW49320.1 numbering、または二列 mapping。
- [ ] candidate residue を benchmark residue に昇格できる条件を短く定義する。

## Phase 3: 次の Benchmark 系

- [ ] RosC/RosD の対象 BGC、protein IDs、PDB/AlphaFold availability、文献証拠を確認する。
- [ ] `rosc_project` の初期構造を作成する。
- [ ] MpaD、MpaE、MpaDE を一つの系として扱うか、分けて扱うか決める。
- [ ] `mpa_project` または `mpade_project` の初期構造を作成する。
- [ ] 新規系ごとに、MIBiG、NCBI、UniProt、PDB の source URLs を先に記録してからファイルを取得する。
- [ ] 新規系ごとに、metadata CSV、data collection log、literature note、validation checklist を共進化解析より前に作成する。

## Phase 4: MSA と共進化 MVP

- [ ] 最初の解析対象を決める。TamI のみで始めるか、TamI と TxtC を並行して進めるかを選ぶ。
- [ ] homolog search strategy を定義する: database、E-value threshold、identity range、coverage cutoff、maximum sequence count。
- [ ] full-length MSA と catalytic-domain MSA を作成する。
- [ ] MSA quality metrics を記録する: sequence count、Neff、gap fraction、coverage、conserved P450 motif alignment。
- [ ] 最初の軽量 baseline として MIp/APC を実行する。
- [ ] MSA quality が許容範囲に入った後、DCA または EVcouplings を実行する。
- [ ] すべての解析コマンドと parameter を、再実行可能な script または notebook に保存する。

## Phase 5: 構造マッピング

- [ ] TamI の共進化スコアを 6XA2 と 6XA3 にマッピングする。
- [ ] TxtC の共進化スコアを Q8VS75/PDB numbering で 6F0B と 6F0C にマッピングする。
- [ ] region、PDB availability、ligand/heme proximity、literature support を含む residue annotation table を作成する。
- [ ] SRS、BC loop、FG loop、substrate channel、heme-proximal residues をどのように score 化するか決める。

## Phase 6: 候補残基ランキング

- [ ] DCA/MIp rank、structure priority、literature support、uncertainty penalty を含む初期 FDetR scoring table を設計する。
- [ ] TamI の既知残基を用いて、ranking が既知の機能残基または pocket residue を回収できるか確認する。
- [ ] TxtC を用いて、dynamic/conformational-switching residue を過剰解釈せずに扱えるか確認する。
- [ ] 最初の candidate-residue report を Markdown 形式で作成する。

## 現在の優先事項

1. TamI と TxtC の cross-dataset summary table を完成させる。
2. TxtC の下流解析における residue numbering convention を決定する。
3. 新しい系を追加する前に、TamI の MSA/coevolution MVP を開始する。