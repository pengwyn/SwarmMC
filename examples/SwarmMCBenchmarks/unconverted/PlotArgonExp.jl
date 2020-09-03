

ParseLine(line) = float.(split(strip(line), ','))

function ReadData(filename)
    set = []
    
    open(filename) do file
        while !eof(file)
            name = readline(file) |> strip
            E = ParseLine(readline(file))
            val = ParseLine(readline(file))

            push!(set, (name, E, val))
        end
    end

    set
end


using Plots
function PlotDrifts()
    set = ReadData("ArExp_W.txt")

    p = plot(reuse=false)
    for item in set
        name,E,val = item
        plot!(E,val,label=name, marker=:circle)
    end
    xaxis!(scale=:log10)
    title!("Drifts")
    display(p)
end

function PlotTownsend()
    set = ReadData("ArExp_alpha.txt")

    p = plot(reuse=false)
    for item in set
        name,E,val = item
        plot!(E,val,label=name, marker=:circle)
    end
    xaxis!(scale=:log10)
    title!("Alpha")
    display(p)
end
function PlotChE()
    set = ReadData("ArExp_ChEnDL.txt")

    p = plot(reuse=false)
    for item in set
        name,E,val = item
        plot!(E,val,label=name, marker=:circle)
    end
    xaxis!(scale=:log10)
    title!("Drifts")
    display(p)

    set = ReadData("ArExp_ChEnDT.txt")

    for item in set
        name,E,val = item
        plot!(E,val,label=name, marker=:circle)
    end
    xaxis!(scale=:log10)
    title!("Drifts")
    display(p)
end
