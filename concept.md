# BEEP 研究コンセプト

## 概要

BEEP (BGC-Heme Enzyme Evolution Project) は、生合成遺伝子クラスター (BGC) に含まれる heme enzyme、特に cytochrome P450 の機能分化を、BGC context、配列共進化、三次元構造、文献証拠を統合して解析する研究プロジェクトである。

中心となる問いは次の通りである。

> なぜ、ある P450 は特定の BGC の中で特定の基質選択性、反応順序、酸化位置、または遠隔的な構造制御を獲得したのか。

この問いに対して、単なる配列類似性ではなく、BGC 内での役割、隣接遺伝子との関係、残基間共進化、構造上の位置、既知の実験証拠を組み合わせて答えることを目指す。

## 背景

天然物生合成に関わる P450 は、保存された heme 結合モチーフを持つ一方で、基質選択性や酸化順序を決める領域が多様である。重要な情報は、活性中心近傍だけでなく、substrate recognition sites、BC loop、FG loop、基質チャネル、水分子ネットワーク、さらには fusion enzyme の interface にも現れる可能性がある。

従来の配列類似性検索は「どの酵素に似ているか」を示すには有効だが、「どの残基が機能分化を支えたか」を説明するには不十分である。BEEP では、BGC という代謝経路単位の文脈を加え、P450 内部の残基レベル共進化と構造情報を統合する。

## 作業仮説

1. BGC 内に取り込まれた P450 は、周辺遺伝子、基質供給、反応順序による制約を受けている。
2. 機能分化に関わる残基は、単純な全体保存度ではなく、系統特異的置換や残基間共進化として現れやすい。
3. 重要残基は heme 近傍だけでなく、基質チャネル、loop、SRS、水分子ネットワーク、動的構造制御領域にも分布しうる。
4. MIp/DCA/EVcouplings の結果は、構造マッピングと文献証拠を組み合わせることで、実験候補として解釈しやすくなる。
5. ProteinMPNN などの構造条件付き配列設計は、機能証明ではなく、構造的妥当性を補助的に評価する手段として扱うべきである。

## 初期 Benchmark 系

| System | BGC / Product | P450 | Benchmark value |
| --- | --- | --- | --- |
| TamI | tirandamycin BGC / tirandamycin | TamI | 既知の構造、反応機構、文献証拠、候補残基が揃っている第一 benchmark。 |
| TxtC | thaxtomin BGC / thaxtomin A | TxtC | substrate conformational switching を扱う第二 benchmark。構造配列と BGC 配列の差異も、番号対応の検証例として有用。 |
| RosC/RosD | rosamicin 関連 BGC | RosC/RosD | 多段階酸化や loop dynamics の検証候補。 |
| MpaD/MpaE/MpaDE | mycophenolic acid BGC | MpaD/MpaE/MpaDE | P450-hydrolase fusion や domain interface の解析候補。 |

## 解析パイプライン

```text
BGC / P450 benchmark selection
  -> raw data collection and metadata logging
  -> sequence-source validation and residue-numbering check
  -> homolog search and MSA construction
  -> MIp/APC, DCA, EVcouplings
  -> structural mapping to PDB or AlphaFold models
  -> integration with BGC context and literature evidence
  -> FDetR candidate ranking
  -> mutation-priority report
```

## FDetR の考え方

FDetR は Functional Determinant Residue の略称として使う。これは false discovery rate の FDR と混同しないようにするための作業上の表記である。

FDetR candidate は、以下の要素を統合して優先順位をつける。

- DCA / EVcouplings rank
- MIp/APC rank
- BGC context support
- structural priority
- lineage-specific substitution
- literature support
- sequence-numbering uncertainty
- MSA quality uncertainty

初期スコアの概念式は以下のように置く。

```text
FDetR_score =
  w1 * z(DCA_rank_score)
+ w2 * z(MIp_rank_score)
+ w3 * gene_context_score
+ w4 * structural_priority_score
+ w5 * lineage_specific_substitution_score
+ w6 * literature_support_score
- w7 * uncertainty_penalty
```

このスコアは最初から絶対的な答えとして扱わず、benchmark 系で検証しながら重みと解釈を調整する。

## 直近の研究方針

短期的には、TamI と TxtC を用いて、データ整理、配列検証、構造マッピング、候補残基整理の標準手順を固める。次に TamI の MSA/coevolution MVP を実行し、既知残基 Leu101、Leu244、Leu295 がどの程度回収されるかを確認する。

TxtC では、Q8VS75/PDB 番号系と AAW49320.1/BGC 番号系の違いを明示的に扱い、Ser280/Thr385 のような candidate mechanism residues を過剰に benchmark residue として扱わないようにする。

中期的には、RosC/RosD と MpaD/MpaE/MpaDE を追加し、loop dynamics、多段階酸化、fusion enzyme interface など、TamI/TxtC だけでは見えにくい機能分化のパターンを比較する。

最終的には、未知または未注釈 BGC 中の P450 に対して、根拠付きの候補機能残基リストと mutation-priority report を出力できる解析基盤を構築する。