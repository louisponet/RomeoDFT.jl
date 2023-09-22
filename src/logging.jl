function process_stacktrace()
    st = stacktrace()
    uid = findfirst(x -> occursin("Systems", string(x.file)), st)
    
    if uid === nothing
        uid = findfirst(x -> occursin("threadingconstructs.jl",string(x.file)), st)
        if uid === nothing
            return "At $(st[end-1])"
        else
            uid -= 1
        end
    end
    stackline = st[uid]
    curline = 1
    curupdate = ""
    open(string(stackline.file), "r") do f
        while curline < stackline.line
            line = readline(f)
            if occursin("update", line) && occursin("function", line)
                curupdate = line
            end
            curline += 1
        end
    end
    m = match(r"update\([\w]*::(\w+).", curupdate)
    if m !== nothing
        sys = m.captures[1]
        return "In $(split(string(stackline.file), "RomeoDFT/")[end]):$(stackline.line) [$sys]"
    else
        return "In $(split(string(stackline.file), "RomeoDFT/")[end]):$(stackline.line) "
    end
end

struct SafeLoggingComponent{T, CT<:Overseer.AbstractComponent{T}} <: Overseer.AbstractComponent{T}
    c::CT
    lc::Component{Log}
    lock::ReentrantLock
end

function Base.lock(f, c::SafeLoggingComponent)
    lock(c.lock)
    try
        f()
    finally
        unlock(c.lock)
    end
end

function Base.setindex!(c::SafeLoggingComponent{T}, v::T, e::Union{Overseer.AbstractEntity, Int}) where {T}
    @debugv 3 if e in c.lc
        stmsg = process_stacktrace()
        push!(c.lc[e].logs, "($(now())) -- $stmsg: Setting $T")
        "$stmsg: $(Entity(e)) - Setting $T"
    end
    if e ∉ c.c
        lock(c) do 
            return c.c[e] = v
        end
    else
        return c.c[e]=v
    end
end

function Base.pop!(c::SafeLoggingComponent{T}, e::AbstractEntity) where {T}
    @debugv 3 if e in c.lc
        stmsg = process_stacktrace()
        push!(c.lc[e].logs, "($(now())) -- $stmsg: Popping $T")
        "$stmsg: $(Entity(e)) - Popping $T"
    end
    lock(c) do 
        return pop!(c.c, e)
    end
end
function Base.pop!(c::SafeLoggingComponent{T}) where {T}
    e = lock(c) do
        pop!(c.c)
    end
    
    @debugv 3 if e in c.lc
        stmsg = process_stacktrace()
        push!(c.lc[e].logs, "($(now())) -- $stmsg: Popping $T")
        "$stmsg: $(Entity(e)) - Popping $T"
    end
    return e
end

Overseer.component(c::SafeLoggingComponent) = c.c

Overseer.parent(c::SafeLoggingComponent, args...) = Overseer.parent(c.c, args...)
Overseer.npools(c::SafeLoggingComponent, args...) = Overseer.npools(c.c, args...)
Overseer.pools(c::SafeLoggingComponent)           = Overseer.pools(c.c)
Overseer.pool(c::SafeLoggingComponent, args...)   = Overseer.pool(c.c, args...)

function log(e::Overseer.EntityState, msg::AbstractString)
    @debugv 3 begin
        stmsg = process_stacktrace()
        debmsg = "$stmsg: $(e.e) -- $msg"
        msg = "$stmsg: $msg"
        debmsg
    end
    @debugv 1 begin
        "$(e.e) -- $msg"
    end
    logc = e.components[1].lc
    push!(logc[e.e].logs, "$(now()) -- $msg")
end
log(m::Searcher, e::AbstractEntity, msg::AbstractString) = log(m[e], msg)
