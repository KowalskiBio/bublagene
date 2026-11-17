// A small built-in library of short, well-known elements. Long elements
// (resistance genes, origins, viral promoters) are deliberately not guessed at:
// they come from the reader's own saved elements and from the features already
// annotated in other open files (see detect.js).
//
// kind 'dna'     — matched on either strand, allowing a few mismatches.
// kind 'protein' — matched in all six reading frames, so any codon choice finds it.

const dna = (name, type, seq, note = '') => ({ name, type, kind: 'dna', seq, note });
const aa = (name, seq, note = '') => ({ name, type: 'misc_feature', kind: 'protein', seq, note });

export const BUILTIN_LIBRARY = [
  dna('T7 promoter', 'promoter', 'TAATACGACTCACTATAGG', 'Bacteriophage T7 RNA polymerase promoter'),
  dna('T3 promoter', 'promoter', 'AATTAACCCTCACTAAAGG', 'Bacteriophage T3 RNA polymerase promoter'),
  dna('SP6 promoter', 'promoter', 'ATTTAGGTGACACTATAG', 'Bacteriophage SP6 RNA polymerase promoter'),
  dna('lac operator', 'protein_bind', 'AATTGTGAGCGGATAACAATT', 'LacI binding site'),
  dna('M13 fwd', 'primer_bind', 'GTAAAACGACGGCCAGT', 'Common sequencing primer site'),
  dna('M13 rev', 'primer_bind', 'CAGGAAACAGCTATGAC', 'Common sequencing primer site'),
  dna('T7 terminator', 'terminator', 'CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG', 'Bacteriophage T7 transcription terminator'),
  dna('Kozak', 'regulatory', 'GCCACCATGG', 'Translation initiation context in eukaryotes'),
  dna('loxP', 'misc_feature', 'ATAACTTCGTATAGCATACATTATACGAAGTTAT', 'Cre recombinase site'),
  dna('FRT', 'misc_feature', 'GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC', 'Flp recombinase site'),
  dna('attB1', 'misc_feature', 'ACAAGTTTGTACAAAAAAGCAGGCT', 'Gateway recombination site'),
  dna('attB2', 'misc_feature', 'ACCACTTTGTACAAGAAAGCTGGGT', 'Gateway recombination site'),
  dna('I-SceI site', 'misc_feature', 'TAGGGATAACAGGGTAAT', 'I-SceI homing endonuclease site'),

  aa('His tag', 'HHHHHH', 'Polyhistidine tag'),
  aa('FLAG tag', 'DYKDDDDK'),
  aa('HA tag', 'YPYDVPDYA'),
  aa('Myc tag', 'EQKLISEEDL'),
  aa('V5 tag', 'GKPIPNPLLGLDST'),
  aa('Strep-tag II', 'WSHPQFEK'),
  aa('T7 tag', 'MASMTGGQQMG'),
  aa('S-tag', 'KETAAAKFERQHMDS'),
  aa('E-tag', 'GAPVPYPDPLEPR'),
  aa('Avi tag', 'GLNDIFEAQKIEWHE', 'Biotinylation site'),
  aa('Spy tag', 'AHIVMVDAYKPTK'),
  aa('HiBiT', 'VSGWRLFKKIS'),
  aa('SV40 NLS', 'PKKKRKV', 'Nuclear localisation signal'),
  aa('TEV site', 'ENLYFQG', 'TEV protease cleavage site'),
  aa('HRV 3C site', 'LEVLFQG', 'PreScission protease cleavage site'),
  aa('Thrombin site', 'LVPRGS', 'Thrombin cleavage site'),
  aa('Factor Xa site', 'IEGR', 'Factor Xa cleavage site'),
  aa('Enterokinase site', 'DDDDK', 'Enterokinase cleavage site'),
  aa('P2A', 'ATNFSLLKQAGDVEENPGP', 'Ribosome-skipping peptide'),
  aa('T2A', 'EGRGSLLTCGDVEENPGP', 'Ribosome-skipping peptide'),
];
