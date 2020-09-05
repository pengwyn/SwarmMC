
union!(LOAD_PATH, ["."])
using MagneticField

#Blist =  (1:9)' .* 10.^(0:4)
#Blist =  [1,5]' .* 10.^(0:4)
#Blist = [0 ; Blist[:]]
Blist = [1000.]*Hx
@CheckTurns for BHx in Blist,
    Bθ in [0,30,45,60,90]

    p = SetupParams(;BHx, Bθ)

    props = LoopMaxTime(p, 1)

    Save(p, props)
end
