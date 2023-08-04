# Helper Macros for checking that a function doesn't allocate memory
# https://forem.julialang.org/dpinol/detecting-test-allocated-gotchas-34op

macro testnoallocations(expressions...)
    # Uncomment the following line if non-const globals should be prohibited
    _failIfNonConstGlobalsInExpressions(__module__, expressions...)
    return esc(
        quote
            !@isCalledFromFunction() &&
                @warn "Since not called from a function @allocated could be imprecise"
            
            # Executes the code twice to exclude the allocation cost of the first call
            $(expressions...)
            @test (@allocated $(expressions...)) === 0
        end
    )
end

macro testnoallocations_noprecompile(expressions...)
    # Uncomment the following line if non-const globals should be prohibited
    _failIfNonConstGlobalsInExpressions(__module__, expressions...)
    return esc(
        quote
            !@isCalledFromFunction() &&
                @warn "Since not called from a function @allocated could be imprecise"
            @test (@allocated $(expressions...)) === 0
        end
    )
end

macro isCalledFromFunction()
    expr = esc(:(
        try
            currentFunctionName = nameof(var"#self#")
            true
        catch
            false
        end
    ))
    return expr
end

function _failIfNonConstGlobalsInExpressions(mod::Module, expressions...)
    for e in expressions
        nonConstGlobals = (
            arg for arg in expressionsymbols(e) if isdefined(mod, arg) && !isconst(mod, arg)
        )
        if !isempty(nonConstGlobals)
            error(
                "testnoallocations called with expression containing non const global symbols $(collect(
            nonConstGlobals
        ))",
            )
        end
    end
end

" Return all the symbols that make up an expression (or itself if a symbol is passed)"
function expressionsymbols(e::Union{Expr, Symbol, Number})
    !isa(e, Expr) && return (ex for ex in (e,) if isa(e, Symbol))
    topSymbols = (arg for arg in e.args if isa(arg, Symbol))
    subExpressions = (expressionsymbols(arg) for arg in e.args if isa(arg, Expr))
    return (topSymbols..., Iterators.flatten(subExpressions)...)
end
