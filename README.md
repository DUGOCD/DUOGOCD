# EIEAPT in MATLAB

## File Structure Overview

### Core Execution
- `main.   m` - Main entry point (run this file)

### Initialization
- `chaotic_initialization.   m` - Chaotic population initialization (Logistic map)
- `initialize.   m` - Basic population initialization
- `initialize1.   m` -Special case initialization
- `initialize2.   m` - Special case initialization
- `UniformPoint.   m` - Reference vector generation

### Genetic Operators
- `cross_mutation1.   m`
- `cross_mutation6.   m`

### Selection Mechanisms
- `F_mating.   m` - Mating pool selection
- `F_select11.   m` - Elite selection
- `elitsm.   m` - Elite preservation
- `non_domination_sort.   m` - Fast non-dominated sorting

### Evaluation & Visualization
- `object_fun.   m` - Objective function evaluation
- `TURE_PF.   m` - True Pareto front generation
- `Draw.   m` - Pareto front visualization

2.    Run in MATLAB:

## Key Parameters (configurable in main.m)

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| `pop` | Population size | 100 |
| `maxFE` | Max function evaluations | 100,000 |
| `f_num` | Number of objectives | 2-3 |
| `x_num` | Decision variables | 500-2000 |
| `pc` | Crossover probability | 1.0 |
| `pm` | Mutation probability | 1/D |

## Customization Guide

### Change Test Problem
Modify in `main.m`:
```matlab
fun = 'LSMOP3';     % Switch to LSMOP3 problem
```

### Adjust Chaotic Mapping
Edit `chaotic_initialization.m`:
```matlab
map_type = 'tent';     % Change to Tent map
```

### Modify Crossover
Tune in `cross_mutation6.m`:
```matlab
yita1 = 15;     % SBX distribution index
```

## Visualization
View results with:
```matlab
Draw(chromo(:,end-f_num+1:end));     % Show final Pareto front
```

## Troubleshooting
- **"Undefined function" errors**: Ensure all files are in MATLAB path
- **Dimension mismatches**: Verify `x_num` matches problem dimensions
