#!/usr/bin/env python3
"""Rewrite the pilot search-deviation note with the final adopted workflow."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEXT = """# Pilot parameter deviation (2026-08-04)

手順書既定の `nr + Bacteria[Organism] + blastp + max_target_seqs 500` を TamI/TxtC に実行したが、NCBI remote service が両 query とも CPU usage limit で終了した。

その後、次の remote search を段階的に試行した。

- `nr + Actinomycetota + blastp-fast`
- `refseq_protein + Bacteria/Actinomycetota + blastp-fast`
- `Swiss-Prot + Bacteria + blastp-fast`

これらは接続失敗、CPU 制限、または実用的な待機時間を超える queue/stall となったため、pilot dataset には採用しなかった。失敗時の stderr log と空の raw output は provenance として各 system の `logs/` と `data/raw/homologs/` に保持する。

最終的に採用した検索は、NCBI preformatted BLAST database の local Swiss-Prot snapshot（更新日 2026-07-31）に対する `blastp-fast`、taxonomy filter `taxid:2` (Bacteria) である。E-value 1e-20、最大 500 hits、identity 25–95%、query coverage 70%以上、長さ比 0.75–1.25、X 1%以下という主要 filter は維持した。配列は同 database から `blastdbcmd` で抽出した。

制約: Swiss-Prot は高品質だが nr より収録範囲が狭い。TamI は reference 込み 69 配列を確保した一方、TxtC は 44 配列で手順書の pilot 推奨 50 配列を下回る。したがって本成果物は MSA/tree 手順と mapping の pilot 検証用であり、最終 coevolution 用 homolog set は local nr、分割 remote search、または MMseqs2 などで拡張する。
"""

(ROOT / "metadata" / "msa" / "pilot_parameter_deviation_2026-08-04.md").write_text(
    TEXT, encoding="utf-8", newline="\n"
)
