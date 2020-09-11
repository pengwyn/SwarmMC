# SwarmMC.jl

Documentation for SwarmMC.jl

Simple API list in order of use.

## Setup of simulation

```@docs
SwarmMC.PARAMS
SwarmMC.Finalise
SwarmMC.PARTTYPE
SwarmMC.GAS
SwarmMC.CreateCollFreq
SwarmMC.CF_STYLE_ENUM
SwarmMC.ISS_ENUM
SwarmMC.GNS_ENUM
SwarmMC.MEAS_BIN
SwarmMC.MEAS_QUANT
SwarmMC.MakeSaveName
```

## Running of simulation

```@docs
SwarmMC.Save
SwarmMC.LoopMaxTime
SwarmMC.BunchedPropagate
SwarmMC.PROPS_OUT
```

## Loading/Combining results

```@docs
SwarmMC.ReadAll
SwarmMC.GetAll
SwarmMC.CollectUp
SwarmMC.NameElements
SwarmMC.PrefixSets

```

## Analysing results

```@docs
SwarmMC.Quants
SwarmMC.EpsFromVel

SwarmMC.ΔR
SwarmMC.RGrid

SwarmMC.TGrid
SwarmMC.TInstGrid
SwarmMC.ΔT

SwarmMC.EpsGrid
SwarmMC.ΔEps

SwarmMC.ZGrid
SwarmMC.ΔZ

SwarmMC.CosthGrid
SwarmMC.ΔCosth
```

## Advanced

### Definining new measurement types

```@docs
SwarmMC.Log2Fac
SwarmMC.InstValue
```

### User callbacks

```@docs
SwarmMC.INTEGRATOR
```
