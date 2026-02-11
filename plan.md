# Alzheimer's Disease QSP Model Implementation Plan

## 1. Project Overview

**목표**: Madrasi et al. (2021) 논문의 알츠하이머병 Aβ plaque QSP 모델을 R로 재현

**논문**: "Systematic in silico analysis of clinically tested drugs for reducing amyloid-beta plaque accumulation in Alzheimer's disease"

**모델링 대상 약물 (7종)**:
- Anti-Aβ mAbs: Aducanumab, Crenezumab, Solanezumab, Bapineuzumab
- BACE inhibitors: Elenbecestat, Verubecestat
- γ-secretase inhibitor: Semagacestat

---

## 2. Model Structure

### 2.1 Compartments (4개)

| Compartment | Volume | Reference |
|-------------|--------|-----------|
| Plasma | 3 L | Pearson et al., 2008 |
| Peripheral | 3 L | Assumed same as plasma |
| CSF | 0.139 L | Nau et al., 2010 |
| Brain ISF | 0.261 L | Shah & Betts, 2012 |

### 2.2 Species (State Variables)

**Aβ Biology Species:**
- APP (Amyloid Precursor Protein)
- BACE (β-secretase, membrane-bound)
- BACEs (soluble BACE)
- γ-secretase (Gamma)
- CTFβ (C-terminal fragment β)
- Aβ monomer
- Aβ oligomer (Aolig)
- Aβ plaque (Aplaq) - Brain ISF only
- FcR (Fc receptor) - Brain ISF only

**Drug Species:**
- mAb (free antibody)
- mAb:Aβ complex
- mAb:Aolig complex
- mAb:Aplaq complex
- mAb:Aolig:FcR complex
- mAb:Aplaq:FcR complex
- BI (BACE inhibitor)
- BI:BACE complex
- GI (γ-secretase inhibitor)
- GI:γ complex

---

## 3. Key Reactions & ODEs

### 3.1 Aβ Production Pathway
```
APP + BACE ⟷ APP:BACE → CTFβ + BACE (kcatBACE)
CTFβ + γ ⟷ CTFβ:γ → Aβ + γ (kcatGamma)
```

### 3.2 Aβ Aggregation
```
Aβ monomer ⟷ Aβ oligomer (kM2G / kG2M)
Aβ oligomer ⟷ Aβ plaque (kG2P / kP2G) - Brain ISF only
```

### 3.3 Drug Mechanisms

**mAb Binding:**
```
mAb + Aβ ⟷ mAb:Aβ (konAb / koffma0)
mAb + Aolig ⟷ mAb:Aolig (konAb / koffma1)
mAb + Aplaq ⟷ mAb:Aplaq (konPD / koffma2)
```

**ADCP (Antibody-Dependent Cellular Phagocytosis):**
```
mAb:Aolig + FcR ⟷ mAb:Aolig:FcR → clearance (kcatADCP)
mAb:Aplaq + FcR ⟷ mAb:Aplaq:FcR → clearance (kcatADCP)
```

**BACE Inhibition:**
```
BI + BACE ⟷ BI:BACE (competitive inhibition)
```

**γ-secretase Inhibition:**
```
GI + γ ⟷ GI:γ (competitive inhibition)
```

---

## 4. Parameters

### 4.1 Core Biology Parameters (Table S3)

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| APP-BACE association | konPP | 1E-3 | 1/(nM·s) |
| BACE turnover | kcatBACE | 7.2E-3 | 1/s |
| γ-secretase turnover | kcatGamma | 1.2E-3 | 1/s |
| Monomer→Oligomer | kM2G | 1.4E-5 | 1/s |
| Oligomer→Monomer (plasma) | kG2M_plasma | 1.4E5 | 1/s |
| Oligomer→Monomer (brain) | kG2M_bisf | 1.4E-8 | 1/s |
| Oligomer→Plaque | kG2P | 7.0E-8 | 1/s |
| Plaque→Oligomer | kP2G_bisf | 7E-11 | 1/s |
| Aβ clearance (plasma) | kclearAbeta_plasma | 1.9E-4 | 1/s |
| Aβ clearance (brain) | kclearAbeta_bisf | 5.5E-5 | 1/s |
| Plaque clearance | kclearAplaq | 8.0E-9 | 1/s |
| APP clearance | kclearAPP | 6.93E-5 | 1/s |
| BACE clearance | kclearBACE | 1.2E-5 | 1/s |
| γ-secretase clearance | kclearGamma | 8.0E-6 | 1/s |
| FcR clearance | kclearFcR | 9.6E-5 | 1/s |

### 4.2 Transport Parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Plasma→Peripheral (Aβ) | k12Abeta | 3E-4 | 1/s |
| Plasma→CSF (Aβ) | k13Abeta | 1.72E-9 | 1/s |
| Plasma→Brain ISF (Aβ) | k14Abeta | 1.48E-7 | 1/s |
| Peripheral→Plasma (Aβ) | k21Abeta | 2.0E-5 | 1/s |
| CSF→Plasma (Aβ) | k31Abeta | 4.5E-5 | 1/s |
| Brain ISF→Plasma (Aβ) | k41Abeta | 1.48E-8 | 1/s |
| Brain ISF→CSF (Aβ) | k43Abeta | 7.5E-5 | 1/s |

### 4.3 mAb-Specific Parameters (Table S4)

| Parameter | Aducanumab | Bapineuzumab | Crenezumab | Solanezumab |
|-----------|------------|--------------|------------|-------------|
| MW (g/mol) | 150000 | 148800 | 148800 | 148800 |
| kclearmAb (1/s) | 1.5E-7 | 8.4E-7 | 1.0E-6 | 8.9E-7 |
| Kd_monomer (nM) | 8000 | 4 | 24 | 1.8 |
| Kd_oligomer (nM) | 4 | 1.5 | 6 | 4 |
| Kd_plaque (nM) | 4 | 1.5 | 6 | NA |
| Kd_FcR (nM) | 15.6 | 15.6 | 62.4 | 15.6 |
| kcatADCP (1/s) | 3.7E-5 | 3.7E-5 | 1.2E-5 | 3.7E-5 |

### 4.4 Small Molecule Parameters (Table S5)

| Parameter | Elenbecestat | Verubecestat | Semagacestat |
|-----------|--------------|--------------|--------------|
| MW (g/mol) | 437.4 | 409.4 | 361.4 |
| Oral scaling | 0.01 | 0.02 | 0.04 |
| kabs (1/s) | 1.50E-4 | 1.50E-4 | 2.00E-4 |
| kclear (1/s) | 2.5E-4 | 8.3E-5 | 1.0E-3 |
| Ki (nM) | 2 | 7.8 | 15 |

---

## 5. Implementation Steps

### Phase 1: Core Infrastructure
- [ ] R project 설정 및 패키지 설치 (deSolve, ggplot2, tidyverse)
- [ ] Parameter 파일 생성 (parameters.R)
- [ ] Compartment volumes 정의

### Phase 2: Aβ Biology Model
- [ ] State variables 정의
- [ ] ODE system 구현 (Aβ production, aggregation, transport)
- [ ] Steady state 검증 (Table S2 값과 비교)
- [ ] SILK experiment 시뮬레이션으로 검증 (Figure S1)

### Phase 3: Drug Models
- [ ] mAb PK model 구현
- [ ] mAb-Aβ binding 및 ADCP 구현
- [ ] BACE inhibitor model 구현
- [ ] γ-secretase inhibitor model 구현

### Phase 4: Simulation & Validation
- [ ] Aducanumab 시뮬레이션 (Figure 2A 재현)
- [ ] Bapineuzumab 시뮬레이션 (Figure 2B 재현)
- [ ] Crenezumab 시뮬레이션 (Figure 2C 재현)
- [ ] Solanezumab 시뮬레이션 (Figure 2D 재현)
- [ ] Elenbecestat 시뮬레이션 (Figure 3A 재현)
- [ ] Verubecestat 시뮬레이션 (Figure 3B 재현)
- [ ] Semagacestat 시뮬레이션 (Figure 3C 재현)

### Phase 5: Analysis
- [ ] Plaque reduction 비교 분석 (Figure 4)
- [ ] Binding profile별 분석 (Figure 5)
- [ ] Plaque turnover 민감도 분석 (Figure 6A)
- [ ] Peripheral sink hypothesis 검증 (Figure 6B)

---

## 6. File Structure

```
AD/
├── plan.md                    # This file
├── R/
│   ├── 01_parameters.R        # All model parameters
│   ├── 02_model_odes.R        # ODE system definition
│   ├── 03_simulation.R        # Simulation functions
│   ├── 04_drug_models.R       # Drug-specific models
│   ├── 05_validation.R        # Model validation
│   └── 06_analysis.R          # Analysis and figures
├── data/
│   └── clinical_data.csv      # Digitized clinical data
├── output/
│   └── figures/               # Generated figures
└── main.R                     # Main execution script
```

---

## 7. Key Equations (ODE Form)

### Aβ Monomer in Brain ISF:
```
d[Aβ_bisf]/dt = kcatGamma * [CTFβ:γ_bisf]
              - kM2G * [Aβ_bisf]
              + kG2M_bisf * [Aolig_bisf]
              - kclearAbeta_bisf * [Aβ_bisf]
              - k41Abeta * [Aβ_bisf] * (V_bisf/V_plasma)
              + k14Abeta * [Aβ_plasma] * (V_plasma/V_bisf)
              - k43Abeta * [Aβ_bisf] * (V_bisf/V_csf)
              - konAb * [mAb_bisf] * [Aβ_bisf]
              + koffma0 * [mAb:Aβ_bisf]
```

### Plaque in Brain ISF:
```
d[Aplaq_bisf]/dt = kG2P * [Aolig_bisf]
                 - kP2G_bisf * [Aplaq_bisf]
                 - kclearAplaq * [Aplaq_bisf]
                 - konPD * [mAb_bisf] * [Aplaq_bisf]
                 + koffma2 * [mAb:Aplaq_bisf]
```

### ADCP Clearance:
```
d[mAb:Aplaq:FcR]/dt = konPF * [mAb:Aplaq] * [FcR]
                    - koffPF * [mAb:Aplaq:FcR]
                    - kcatADCP * [mAb:Aplaq:FcR]
```

---

## 8. Expected Steady State (Table S2)

| Species | Plasma | CSF | Brain ISF |
|---------|--------|-----|-----------|
| Aβ Monomer (nM) | 0.066 | 2.8 | 0.75 |
| Aβ Oligomer (nM) | 6.9E-12 | 0.02 | 4.4 |
| Aβ Plaque (nM) | 0 | 0 | 39 |

---

## 9. Validation Criteria

1. **SILK Data Fit**: CSF Aβ dynamics (Figure S1)
2. **Steady State**: Match Table S2 values
3. **Drug PK**: Match single-dose PK profiles (Figures 2, 3)
4. **Drug PD**: Match plaque reduction curves
5. **Plaque Half-life**: ~2.75 years (calibrated value)

---

## 10. Notes

- 모든 속도 상수는 초(second) 단위
- 시뮬레이션 시 체중 70kg 가정
- AD 환자 vs 건강인: k31Abeta가 건강인에서 1.4배 빠름
- Plaque 정상상태 도달: 67년 시뮬레이션 후 약물 투여 시작
- SUVR 데이터 정규화: cutoff 값 차감 후 baseline 대비 % 변화 계산
