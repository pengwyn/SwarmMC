
using SwarmMC, StaticArrays
using PyPlot
using NamedArrays

quant_list = [:eps,:flux_W,:flux_DT]

model_list = ['A','B','C','D']
method_list = [:func, :grid, :legendre]

#quant_comp = Dict(key => Matrix(undef, 4,3) for key in quant_list)
#quant_comp = Dict(key => NamedArray(Matrix(undef, 4,3), (model_list, method_list), ("Model", "MC Method")) for key in quant_list)
quant_comp = NamedArray(zeros(4,3,3), (model_list, method_list, quant_list), ("Model", "MC Method", "Quant"))

for model in model_list
    #figure()
    #axeps = PyPlot.axes()
    #figure()
    #axW = PyPlot.axes()

    for method in method_list
        params,props = ReadAll("ReidAniso:model=$model:method=$(method)")

        #axeps.plot(TGrid(params), props[:eps])
        #axW.plot(TGrid(params), props[(:vel,:z)] * params.len_unit / params.time_unit)

        quants = Quants(params,props,include_errors=false)

        for key in quant_list
            #quant_comp[key][ind_model,ind_method] = quants[key]
            #quant_comp[key][model,method] = quants[key]
            quant_comp[model,method,key] = quants[key]
        end
    end
end

