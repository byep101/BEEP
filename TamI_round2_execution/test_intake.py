import unittest, tempfile
from pathlib import Path
from collect_and_filter import select_hits, parse_genbank

class IntakeTests(unittest.TestCase):
    def select(self,lines):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'hits.tsv'; p.write_text('\n'.join(lines)+'\n')
            return select_hits(p)
    def row(self,acc='WP_1.1',identity=25,qend=290,slen=310,bitscore=100):
        return f'ADC79647.1\t{acc}\t{identity}\t{qend}\t1\t{qend}\t1\t{qend}\t1e-30\t{bitscore}\t413\t{slen}'
    def test_query_span_boundary(self):
        _,audit,selected=self.select([self.row(qend=289)])
        self.assertEqual(selected,[])
        self.assertIn('query_span_coverage_lt_0.7',audit[0]['reason'])
        self.assertEqual(self.select([self.row(qend=290)])[2],['WP_1.1'])
    def test_identity_and_length_boundaries(self):
        self.assertEqual(self.select([self.row(identity=95,slen=516)])[2],['WP_1.1'])
        for kw in ({'identity':95.01},{'identity':24.99},{'slen':309},{'slen':517}):
            self.assertEqual(self.select([self.row(**kw)])[2],[])
    def test_one_best_HSP_not_summed_coverage(self):
        n,a,selected=self.select([self.row(qend=290,bitscore=90),self.row(qend=200,bitscore=110)])
        self.assertEqual(n,2); self.assertEqual(len(a),1); self.assertEqual(selected,[])
    def test_xml_rejects_error_and_preserves_metadata(self):
        with self.assertRaises(ValueError): parse_genbank(b'<ERROR>rate limit</ERROR>')
        xml=b'<GBSet><GBSeq><GBSeq_accession-version>WP_1.1</GBSeq_accession-version><GBSeq_sequence>acde</GBSeq_sequence><GBSeq_definition>P450</GBSeq_definition><GBSeq_taxonomy>Bacteria; Actinomycetota</GBSeq_taxonomy></GBSeq></GBSet>'
        result=parse_genbank(xml)['WP_1.1']
        self.assertEqual(result['sequence'],'ACDE')
        self.assertEqual(result['annotation'],'P450')
        self.assertEqual(result['taxonomy'],'Bacteria; Actinomycetota')

if __name__=='__main__': unittest.main()
