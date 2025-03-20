# 3주차 수업
## 4.1 생물학적 서열 표현과 기본 알고리즘

### DNA 유효성 검사 함수  
유효성 검사는 항상 습관적으로 해야한다.

```
def validate_dna (dna_seq) :
    seqm = dna_seq.upper()
    valid = seqm.count('A') + seqm.count('C') + seqm.count('T') + seqm.count('G')
    if valid == len(seqm) : return True
    else : return False
```

### 문자 빈도 분석 함수

```
def frequency (seq) :
    dic = {}
    for s in seq.upper() :
        if s in dic : dic[s] += 1
        else : dic[s] = 1
    return dic
```

### GC 비율 계산 함수
 
```
def gc_content (dna_seq) :
    gc_count = 0
    for s in dna_seq :
        if s in "GCgc" : gc_count += 1
    return gc_count / len(dna_seq)
```


```
# 구간별 gc content
def gc_content_subseq (dna_seq, k=10):
    res = []
    for i in range(0,  len(dna_seq)-k+1, k) :
        subseq = dna_seq[i:i+k]
        gc = gc_content(subseq)
        res.append(gc)
    return res
```

### 부분서열 (subseq) 찾기

```
def find_subseq(dna_seq, pattern) :
    positions = []
    start = 0
    while True :
        idx = dna_seq.find(pattern, start)
        if idx == -1 :
            break
        positions.append(idx)
        start = idx + 1
    return positions
```

## 4.2 전사와 역상보

### 전사

```
def transcription (dna_seq) :
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    return dna_seq.upper().replace("T","U") #T를 U로 바꾸기
```

### 역상보
```
def reverse_complement (dna_seq) :
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    comp = ""
    for c in dna_seq.upper() :
        if c == 'A' :
            comp = "T" + comp
        elif c == "T" :
            comp = "A" + comp
        elif c == "C" :
            comp = "G" + comp
        elif c == "G" :
            comp = "C" + comp
    return comp
```

## 4.3 번역

### 코돈 -> 아미노산 매핑
```
def translate_codon (cod) :
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A", 
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", 
      "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    if cod in tc : return tc[cod]
    else : return None
```

### DNA 서열을 아미노산 서열로 변환

```
def translate_seq (dna_seq, ini_pos = 0) :
    assert validate_dna (dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    seq_aa = ""
    for pos in range(ini_pos, len(seqm)-2,3) :
        cod = seqm[pos:pos+3]
        seq_aa += translate_codon(cod)
    return seq_aa
```

### 주어진 DNA 서열에서 특정 아미노산을 암호화하는 각 코돈의 사용 빈도를 계산하는 함수

```
#변형되어 시험 출제될 수 있음
def codon_usage(dna_seq, aa) :
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    dic = {}
    total = 0
    for i in range(0, len(seqm)-2, 3) :
        cod = seqm[i:i+3]
        if translate_codon(cod) == aa:
            if cod in dic :
                dic[cod] += 1
            else: dic[cod] = 1
            total += 1
    if total >0 :
        for k in dic :
            dic[k] /= total
    return dic
```
