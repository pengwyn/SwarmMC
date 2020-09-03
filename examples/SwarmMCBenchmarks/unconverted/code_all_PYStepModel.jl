
using SwarmMC
using Generic

using PYStepModel

params = SetupParams()
props = LoopMaxTime(params, 100)
Save(params, props)
