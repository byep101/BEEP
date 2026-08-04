# BEEP To-Do

最終更新日: 2026-08-04

## Phase 0: リポジトリ整理

- [x] 文字化けしていた root `README.md` を、読めるプロジェクト概要に置き換える。
- [x] 文字化けしていた root `todo.md` を、現在の進捗に合わせたタスクリストに置き換える。
- [x] root の主要 Markdown を日本語で統一する。
- [ ] `concept.md` を今後どの粒度で維持するか決める。短い研究概要にするか、詳細な研究計画書にするかを選ぶ。
- [x] root `metadata/` に core dataset v0.1 の master CSV と体系別 normalized row を追加する。

## Phase 1: Benchmark データセット基盤

- [x] `tami_project` を標準化データセット構造で作成する。
- [x] TamI の raw data を収集する: MIBiG、NCBI FASTA、UniProt FASTA、PDB/mmCIF。
- [x] TamI の metadata CSV、data collection log、literature note、validation checklist を作成する。
- [x] `txtc_project` を `tami_project` と同じ基本構造で作成する。
- [x] TxtC の raw data を収集する: MIBiG、NCBI FASTA、UniProt FASTA、PDB/mmCIF。
- [x] TxtC の metadata CSV、data collection log、literature note、validation checklist を作成する。
- [x] TxtC の AAW49320.1 と Q8VS75 の配列差分を記録する。
- [x] TamI、TxtC、RosC、MpaDE_MpaD を横断する root-level master CSV を作成する。
- [x] RosC/RosD、MpaD/MpaE/MpaDE を追加するための最低限の metadata fields を確定する。

## Phase 2: データセット検証

- [x] TamI: NCBI と UniProt の配列一致を確認する。
- [x] TamI: Leu101、Leu244、Leu295 が 6XA2/6XA3 構造中で確認できることを検証する。
- [x] TxtC: MIBiG BGC0000444.5 の parse result と source status を確認する。
- [x] TxtC: AAW49320.1 と Q8VS75 の長さ、および mismatch positions を確認する。
- [x] TxtC: Ser280 と Thr385 が 6F0B/6F0C 構造中で確認できることを検証する。
- [x] TxtC の下流解析で使う残基番号体系を決める。Q8VS75/PDB numbering を構造マッピング基準とし、AAW49320.1 は BGC/NCBI anchor として二列 mapping で保持する。
- [ ] candidate residue を benchmark residue に昇格できる条件を短く定義する。

## Phase 3: 次の Benchmark 系

- [x] RosC/RosD の対象 BGC、protein IDs、PDB availability、文献証拠を確認し、master CSV に統合する。
- [ ] `rosc_project` の初期構造を作成する。
- [x] MpaD、MpaE、MpaDE を一つの fusion/domain 比較体系として扱う方針を決め、master CSV に統合する。
- [ ] `mpa_project` または `mpade_project` の初期構造を作成する。
- [x] 新規系ごとに、MIBiG、NCBI、PDB または literature の source URLs を master log に記録する。
- [x] 新規系ごとに、metadata CSV、data collection log、validation checklist を共進化解析より前に作成する。

## Phase 4: MSA と共進化 MVP

- [x] 最初の解析対象を TamI と TxtC の並行 pilot に決定する。
- [x] homolog search strategy を定義する: local Swiss-Prot/Bacteria、E-value 1e-20、identity 25–95%、coverage 70%以上、最大500件。
- [x] full-length MSA と tree-building 用 smart-gap trimmed MSA を作成する。
- [ ] catalytic-domain MSA を追加する。
- [x] pilot MSA quality metrics を記録する: sequence count、gap fraction、coverage、conserved P450 motif alignment、reference residue mapping。
- [ ] 最終共進化用セットを拡張した後に Neff を算出する。
- [ ] 最初の軽量 baseline として MIp/APC を実行する。
- [ ] MSA quality が許容範囲に入った後、DCA または EVcouplings を実行する。
- [x] pilot 解析コマンドと parameter を、再実行可能な script、config、IQ-TREE log に保存する。

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

1. TamI の MSA/coevolution MVP を開始する。
2. TxtC では Q8VS75/PDB numbering と AAW49320.1 BGC anchor の二列 mapping を維持して解析する。
3. RosC と MpaDE_MpaD は、Issues_to_check に残した source split、chain ID、domain boundary を解消してから構造マッピングへ進める。
