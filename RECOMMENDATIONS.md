# Chemistry Model Analysis & Recommendations

## Executive Summary

**Current Status**: Te-dependent rates implemented for electron-impact and recombination reactions (57 reactions). All rates correctly calibrated at Te=1 eV.

**Key Findings**:
1. Group 7 (Neutral-Neutral) reactions have NO Tg dependence → CRITICAL issue for Tg studies
2. H3+ species exists but has NO formation pathway → orphaned species
3. Group 8 (Termolecular) reactions have NO Tg dependence → moderate issue

---

## Priority Recommendations

### 🔴 CRITICAL PRIORITY

#### 1. H3+ Formation Chemistry (IMMEDIATE)

**Issue**: H3+ appears in species list with mobility/drift/sticking, but has **ZERO formation pathways**!

**Missing Reactions**:
```
H2+ + H2 → H3+ + H        k ~ 2.0e-9 cm³/s  (very fast!)
H3+ + CH4 → CH5+ + H2      k ~ 1.5e-9 cm³/s  (proton transfer)
```

**Impact**: If your simulations show any H3+ density, it's numerically incorrect since there's no way to form it.

**Effort**: LOW (2 reactions, 5 minutes)

**Action Required**: Check if H2+ ionization exists. If not, add:
```
e + H2 → H2+ + 2e          k ~ scale_ionization(Te)
```

---

#### 2. Group 7 Tg-Dependence (HIGH PRIORITY if varying Tg)

**Issue**: ALL 65 neutral-neutral reactions are constants calibrated at Tg=400K

**Critical Reactions Needing Tg-dependence**:
```
H + CH4 → CH3 + H2         Ea ~ 10 kcal/mol  (rate × 10 from 300→500K)
CH3 + H → CH2 + H2         Ea ~ 15 kcal/mol  (rate × 20 from 300→500K)
CH + H2 → CH2 + H          Ea ~ 1.3 kcal/mol (rate × 2 from 300→500K)
CH3 + CH3 → C2H6           k ~ Tg^(-0.5)     (rate × 0.8 from 300→500K)
```

**Impact by Tg Range**:
- **Tg = 350-450K**: Moderate (20-30% error)
- **Tg = 300-500K**: LARGE (factor of 2-10× error)
- **Fixed Tg = 400K**: No issue (current calibration)

**Arrhenius Form**: k(Tg) = A × (Tg/400)^n × exp(-Ea/R × (1/Tg - 1/400))

**Effort**: MEDIUM (~20 key reactions need Ea and n parameters from literature)

**Priority**:
- If you vary Tg → **IMMEDIATE**
- If Tg fixed at 400K → **SKIP**

---

### 🟡 MEDIUM PRIORITY

#### 3. Group 8 Termolecular Tg-Dependence

**Current**: 3 reactions constant
```
H + H + M → H2 + M
CH3 + CH3 + M → C2H6 + M
CH3 + H + M → CH4 + M
```

**Should be**: k ~ Tg^(-1.5 to -2.5)

**Impact**: Factor of 2-3× for Tg = 300-500K

**Effort**: LOW (3 reactions)

---

#### 4. CH2+ Formation Chemistry

**Current**: Only one formation pathway:
```
CH2 + CH3+ → CH3 + CH2+    (charge transfer)
```

**Missing**:
```
e + CH2 → CH2+ + 2e        (electron-impact ionization)
Ar+ + CH2 → CH2+ + Ar      (charge transfer)
```

**Impact**: Low (CH2+ is minor ion in most regimes)

**Effort**: LOW (2-3 reactions)

---

### 🟢 LOW PRIORITY

#### 5. Group 9 Wall Sticking Tg^0.5

**Current**: All constant

**Should scale**: k_stick ~ sqrt(8kTg/πm) (thermal velocity)

**Impact**: Factor of 1.3× for 300-500K (small)

**Effort**: LOW

---

#### 6. Groups 3-4 (Ar* & Penning) Tg^0.5

**Impact**: Factor of 1.3× for 300-500K (small)

**Effort**: LOW

**Note**: These are collision-limited, already near max rate

---

#### 7. Group 5 (Ion-Neutral) Tg^(-0.16)

**Impact**: Factor of 1.1× for 300-500K (very small)

**Effort**: LOW

**Note**: Langevin-type reactions, already at collision rate

---

## Missing Chemistry Pathways

### A. Excited State Quenching (Optional)

**Present species**: H2*, CH4*, C2H2*, C3H4* (appear in reactions)

**Missing explicit quenching**:
```
H2* + M → H2 + M           k ~ 1e-11 cm³/s
CH4* + M → CH4 + M         k ~ 1e-11 cm³/s
```

**Priority**: Only if you track excited state densities explicitly

---

### B. Collisional Detachment (Optional)

**Present**: Attachment reactions
```
e + CH4 → CH3- + H         (already included)
e + H2 → H- + H            (already included)
```

**Missing**: Thermal detachment
```
CH3- + M → CH3 + e + M     k ~ 1e-13 cm³/s × exp(-Ea/kTe)
H- + M → H + e + M         k ~ 1e-13 cm³/s × exp(-Ea/kTe)
```

**Priority**: Only important at high Tg (>600K) or if negative ion balance is critical

---

### C. Larger Hydrocarbon Ions (Optional)

**Missing**: C3Hx+, C4Hx+ formation/destruction

**Priority**: Only for long-residence-time plasmas (>10 ms) or if you observe significant C3+/C4+ signals

---

## Questions to Guide Implementation

### 1. Temperature Range
**Q**: What is your Tg range?
- **Fixed at 400K** → No Group 7 changes needed
- **300-500K range** → Must implement Group 7 Tg-dependence
- **>500K** → Also need detachment reactions

### 2. Parameter Sweeps
**Q**: Do you vary Tg in parameter sweeps?
- **YES** → Implement Groups 7 & 8 Tg-dependence (HIGH PRIORITY)
- **NO** → Can skip

### 3. H3+ Observations
**Q**: Do you observe H3+ density in your results?
- **YES** → MUST add H3+ formation immediately
- **NO** → Still recommend adding (species is defined)

### 4. Experimental Validation
**Q**: Are you matching experimental data at fixed conditions?
- **YES at Te~1eV, Tg~400K** → Current model is correctly calibrated, no changes needed
- **YES at other conditions** → Need temperature dependencies
- **NO (predictive)** → Should add all dependencies for robustness

### 5. Simulation Type
**Q**: What kind of simulations?
- **Single-point steady-state** → Minimal changes needed
- **Parameter sweeps (Te, Tg)** → Need all temperature dependencies
- **Transient dynamics** → Need temperature dependencies + excited state chemistry

### 6. Residence Time
**Q**: What is your residence time/pressure?
- **<1 ms, <1 Torr** → Current chemistry sufficient
- **>10 ms, >1 Torr** → May need larger hydrocarbon ion chemistry

---

## Recommended Implementation Order

### Phase 1: Critical Fixes (30 minutes)
1. ✅ Te-dependent rates (DONE)
2. ⬜ Add H3+ formation chemistry (2 reactions)
3. ⬜ Validate H3+ density in results

### Phase 2: If Varying Tg (4-6 hours)
1. ⬜ Literature search for Ea & n for ~20 key Group 7 reactions
2. ⬜ Implement Arrhenius forms for Group 7
3. ⬜ Add Tg^(-2) scaling for Group 8 termolecular
4. ⬜ Validate against experimental Tg trends

### Phase 3: If Needed (1-2 hours)
1. ⬜ Add CH2+ formation pathways
2. ⬜ Add wall sticking Tg^0.5 scaling
3. ⬜ Add excited state quenching (if tracking excited states)

---

## Quick Decision Matrix

| Your Situation | Action Required |
|----------------|-----------------|
| Fixed Tg=400K, matching experiments | ✅ Done! Current model is good |
| Fixed Tg≠400K | Add H3+, consider Group 7 at new Tg |
| Varying Tg in sweeps | Add H3+ + Group 7 Tg-dependence (CRITICAL) |
| Long residence time | Add H3+ + larger hydrocarbon ions |
| Tracking excited states | Add excited state quenching |
| High Tg (>600K) | Add H3+ + Group 7 + detachment |

---

## Implementation Support

I can implement any of these priorities for you. Just let me know:
1. Your Tg range
2. Whether you vary Tg in simulations
3. Which priorities you'd like to tackle

The Te-dependent rates are already perfectly implemented and validated! ✓
