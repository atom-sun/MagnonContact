
function time_step_and_nsteps(t, time_step, nsteps)
    """t is the total unitary evolution time.
    One can also specify `time_step` or `nsteps`.
    If they are both specified, they must satisfy `time_step * nsteps == t`.
    If neither are specified, the default is `nsteps=1`, which means that
    `time_step == t`.

    This function should work as ITensorTDVP.time_step_and_nsteps.

    """
    if !isnothing(time_step) && !isnothing(nsteps)
        if t !== time_step * nsteps
            throw(ArgumentError("t not match: $t, ($time_step, $nsteps)"))
        end
    elseif !isnothing(time_step)
        nsteps = t / time_step
    elseif !isnothing(nsteps)
        time_step = t / nsteps
    else
        time_step = t
        nsteps = 1
    end
    return (time_step, nsteps)
end
