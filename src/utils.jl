using MacroTools

macro error_capturing(expr)
    @capture(expr, for e_ in f__ body_ end)
    
    loop_expr = f[1] 
    if loop_expr.head == :call
        m = loop_expr.args[2]
    elseif loop_expr.head == :macrocall
        m = loop_expr.args[3]
    else
        error("Did not recognize $expr")
    end

    error_comp_name = gensym("err_comp")
    if !inexpr(loop_expr, :Error)
        error_expr = quote
            if $e in $error_comp_name
                continue
            end
        end
    else
        error_expr = ()
    end
    err_name = gensym("err")
    return esc(quote
        $error_comp_name = $m[Error]
        for $e in $loop_expr
            $error_expr
            try
                $body
            catch $err_name
                $error_comp_name[$e] = Error($e, $err_name, stacktrace(catch_backtrace()))
            end
        end
    end)
end

macro error_capturing_threaded(expr)
    @capture(expr, for e_ in f__ body_ end)
    
    loop_expr = f[1] 
    if loop_expr.head == :call
        m = loop_expr.args[2]
    elseif loop_expr.head == :macrocall
        m = loop_expr.args[3]
    else
        error("Did not recognize $expr")
    end

    error_comp_name = gensym("err_comp")
    if !inexpr(loop_expr, :Error)
        error_expr = quote
            if $e in $error_comp_name
                continue
            end
        end
    else
        error_expr = ()
    end
    err_name = gensym("err")
    body = MacroTools.postwalk(x-> x isa Expr && x.head in (:continue, :break) ? :(return) : x, body)
    return esc(quote
        $error_comp_name = $m[Error]
        @sync for $e in $loop_expr
            $error_expr
            Threads.@spawn try
                $body
            catch $err_name
                $error_comp_name[$e] = Error($e, $err_name, stacktrace(catch_backtrace()))
            end
        end
    end)
end

"Calculates the order of eigenvectors of `occ1` which corresponds most to the order of eigenvectors in `occ2`."
function eigvec_order(occ1, occ2)
    eigvals1, eigvecs1 = eigen(occ1)
    eigvals2, eigvecs2 = eigen(occ2)
    dim = size(eigvecs1, 1)
    order = zeros(Int, dim)
    for c1 in 1:dim
        best = 0.0
        for c2 in 1:dim
            t = abs(dot(eigvecs1[:, c1], eigvecs2[:, c2]))
            if t > best && !(c2 in order)
                best = t
                order[c1] = c2
            end
        end
    end
    return (; eigvals1, eigvals2, eigvecs1, eigvecs2, order)
end

maximum_generation(v)                 = maximum(x -> x.generation, v, init=0)
maximum_generation(m::AbstractLedger) = maximum_generation(m[Generation])

ismagnetic(at::Atom) = sum(at.magnetization) != 0 && at.dftu.U != 0
