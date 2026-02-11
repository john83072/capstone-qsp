# Alzheimer's Disease QSP Model (R Implementation)

Madrasi et al. (2021)의 알츠하이머병 QSP 모델을 R로 재현한 프로젝트입니다.

## 1. Model Validation (Steady State)
[cite_start]67년 시뮬레이션을 통해 약물 투여 전 환자의 Steady State 농도를 검증한 결과입니다.

| Species | Compartment | Expected (nM) | Simulated (nM) | Status |
| :--- | :--- | :--- | :--- | :--- |
| Aβ Monomer | Plasma | 0.066 | 0.0663 | ✅ OK |
| Aβ Monomer | Brain ISF | 0.75 | 0.731 | ✅ OK |
| Aβ Plaque | Brain ISF | 39.0 | 37.4 | ✅ OK |

* [cite_start]**Plaque Half-life:** 2.75 years (Paper value: ~2.75 years) [cite: 87]
* [cite_start]**Aβ plasma t1/2:** 1.01 hours (Paper value: ~1 hour) [cite: 87]

## 2. Methodology
- **R Version:** 4.x
- **Packages:** `deSolve`, `ggplot2`
- **Simulation Time:** 67 years (for baseline calibration)
