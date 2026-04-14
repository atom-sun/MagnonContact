function process_resonance_parameters(filename)
    s = read(filename, String)
    s_clean = replace(s, r"\s" => "")
    params = eval(Meta.parse(s_clean))
    params_near = map(params) do (alpha, Jz)
        Jz_near = round(Jz + 0.2, RoundUp, digits=1)
        [alpha, Jz_near]
    end
    return params_near
end
