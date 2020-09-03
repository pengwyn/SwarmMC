
using SwarmMC

p = DefaultParams()

props = LoopBunchedPropagate(p, 1, 100, 1e8, 1e10)

@profile LoopBunchedPropagate(p, 10, 1000, 1e8, 1e10)
