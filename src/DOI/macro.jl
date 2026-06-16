macro variable(model, var::Symbol, phase)
    m = esc(model); v = esc(var); p = esc(phase)
    name = String(var)
    return quote
        local _m = $m
        _m isa Interesso.Optimizer ||
            throw(ArgumentError("@variable: model must be Interesso.Optimizer, got $(typeof(_m))"))
        local _p = $p
        _p isa DOI.PhaseIndex ||
            throw(ArgumentError("@variable: phase must be DOI.PhaseIndex, got $(typeof(_p))"))
        MOI.is_valid(_m, _p) || throw(ArgumentError("@variable: invalid phase index for this model"))

        local _dv = DOI.add_dynamic_variable(_m, _p)
        MOI.set(_m, DOI.DynamicVariableName(), _dv, $name)
        $v = _dv
        _dv
    end
end

macro control(model, var::Symbol, phase)
    m = esc(model); v = esc(var); p = esc(phase)
    name = String(var)
    return quote
        local _m = $m
        _m isa Interesso.Optimizer ||
            throw(ArgumentError("@control: model must be Interesso.Optimizer, got $(typeof(_m))"))
        local _p = $p
        _p isa DOI.PhaseIndex ||
            throw(ArgumentError("@control: phase must be DOI.PhaseIndex, got $(typeof(_p))"))
        MOI.is_valid(_m, _p) || throw(ArgumentError("@control: invalid phase index for this model"))

        local _dv = DOI.add_dynamic_variable(_m, _p)
        MOI.set(_m, DOI.DynamicVariableName(), _dv, $name)
        Interesso.mark_control!(_m, _dv)
        $v = _dv
        _dv
    end
end
